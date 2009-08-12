#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

/* definition of rest mass: Eq. (11) in arXiv:0812.2245v2 M = \int \rho W \sqrt\gamma d^3x
 * this meshes with the EOM for hydro which gives: D_{,t} + (D [\alpha v^i - \beta^i])_{,\i} = 0 
 * (eg. Font lrr-2008-7 Eq. (28) ff) where D=D_{whisky} = \sqrt\gamma*D_{lrr}
 *
 * with this we can easyly get M_{,t} = - \oint_{\partial V} D [\alpha v^i - \beta^i] d\Sigma_i
 * accessing only HydroBase and ADMBase variables this means we need: 
 *   D = \sqrt\gamma W \rho, \gamma = \sqrt{\det g_{ij}}, W^2 = 1/(1-v^2)
 * this calculation is performed in get_gab_ja_onto_detector below, although I think one should 
 * really consider using a generalized get_vars_onto_detector coupled to an expression evaluator 
 * later on (roland)
 */

// scheduled routines
void outflow_init(CCTK_ARGUMENTS);
void outflow_postrecovery(CCTK_ARGUMENTS);
void outflow (CCTK_ARGUMENTS);

#define DIM 3
#define NUM_INPUT_ARRAYS 6+4+4
#define NUM_OUTPUT_ARRAYS 6+4+4

#define MAX_NUMBER_DETECTORS 100
static CCTK_INT file_created[MAX_NUMBER_DETECTORS];

static inline CCTK_REAL pow2(CCTK_REAL x) {return x*x;}
static int upstairs_metric(CCTK_REAL gg[3][3],CCTK_REAL gu[3][3]);
static int drdth_drdph(int i, int j, int ntheta, int nphi,
                int sn,
                CCTK_REAL dth, CCTK_REAL dph,
                CCTK_INT verbose,
                CCTK_INT maxntheta, CCTK_INT maxnphi,
                CCTK_REAL *sf_radius,
                CCTK_REAL *ht, CCTK_REAL *hp);
static int get_gab_ja_onto_detector(CCTK_ARGUMENTS,
                             CCTK_INT sn,
                             CCTK_REAL *g11, CCTK_REAL *g12, CCTK_REAL *g13,
                             CCTK_REAL *g22, CCTK_REAL *g23, CCTK_REAL *g33,
                             CCTK_REAL *jx, CCTK_REAL *jy, CCTK_REAL *jz);
static int get_g_j_local(int i, int j, int ntheta,
                  CCTK_REAL *g11_det,CCTK_REAL *g12_det,CCTK_REAL *g13_det,
                  CCTK_REAL *g22_det,CCTK_REAL *g23_det,CCTK_REAL *g33_det,
                  CCTK_REAL *j1_det,CCTK_REAL *j2_det,CCTK_REAL *j3_det,
                  CCTK_REAL gloc[3][3],CCTK_REAL jloc[3]);
static CCTK_INT outflow_get_local_memory(CCTK_INT npoints);
static int Outflow_write_output(CCTK_ARGUMENTS, CCTK_INT det, CCTK_REAL sum);
static CCTK_REAL *outflow_allocate_array(CCTK_INT npoints, const char *name);


/* IO */
static int Outflow_write_output(CCTK_ARGUMENTS, CCTK_INT det, CCTK_REAL sum)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const char *fmode = (file_created[det]>0) ? "a" : "w";  // file mode: append if already written
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[2048]; // XXX fixed size
  FILE *file;

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "writing output");
  }

  // filename
  sprintf(varname, "outflow_det_%d",det);

  filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
                              strlen (file_extension) +2);
  assert(filename);
  sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);

  // open file
  file = fopen (filename, fmode);
  if (!file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "write_outflow: Could not open scalar output file '%s'",
                filename);
    return -1;
  }

  // write header on startup
  if (file_created[det]<=0) {
    fprintf(file,"# Outflow\n");
    fprintf(file,"# detector no.=%d\n",det);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:it 2:t 3:flux\n");
  }

  // write data
  sprintf (format_str_real,
           "%%d\t%%%s\t%%%s\n", 
           out_format,out_format);
  fprintf(file, format_str_real, cctk_iteration, cctk_time, sum);

  fclose(file); 
  free(filename);

  if (det>=MAX_NUMBER_DETECTORS) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "warn: det=%d, but MAX_NUMBER_DETECTORS=%d, increase",
               det,MAX_NUMBER_DETECTORS);
  }

  if (file_created[det]==0) {
    file_created[det]=1;
  }
  return 1;
}

/* fills g11...g33, j1...j3 with the interpolated numbers */
static int get_gab_ja_onto_detector(CCTK_ARGUMENTS,
                             CCTK_INT sn,
                             CCTK_REAL *g11, CCTK_REAL *g12, CCTK_REAL *g13,
                             CCTK_REAL *g22, CCTK_REAL *g23, CCTK_REAL *g33,
                             CCTK_REAL *jx, CCTK_REAL *jy, CCTK_REAL *jz)
{
  DECLARE_CCTK_ARGUMENTS; 
  DECLARE_CCTK_PARAMETERS; 
  int ierr, retval = 1;
  CCTK_INT ind,ind2d;
  CCTK_REAL th,ph,ct,st,cp,sp;
  CCTK_INT ntheta,nphi;
  // auxilliary variables used in constructing j
  static CCTK_REAL *rho = NULL, *velx = NULL, *vely = NULL, *velz = NULL;
  static CCTK_REAL *beta1 = NULL, *beta2 = NULL, *beta3 = NULL, *alpha = NULL;

  assert(sn>=0);
  assert(g11); assert(g12); assert(g13);
  assert(g22); assert(g23); assert(g33);
  assert(jx); assert(jy); assert(jz);

  ntheta=sf_ntheta[sn];
  nphi=sf_nphi[sn];

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"surface %d (%g,%g,%g) nth,nph (%d,%d)",
	       sn,sf_centroid_x[sn],sf_centroid_y[sn],sf_centroid_z[sn],
	       ntheta,nphi);
  }

  const CCTK_REAL pi=acos(-1.0);
  const CCTK_REAL dth=pi/(ntheta-1);
  const CCTK_REAL dph=2.0*pi/(nphi-1);
  const CCTK_INT npoints=ntheta*nphi;

  CCTK_INT interp_npoints=npoints;
  // uni-processor code - only work on CPU 0
  const CCTK_INT myproc= CCTK_MyProc(cctkGH);
  if ( myproc != 0 ) {
    interp_npoints=0;
  }

  // allocate memory for auxilliary arrays
# define ALLOCATE_TEMP(name) \
  if(name == NULL) \
    name = outflow_allocate_array(npoints, #name); \
  assert(name)
  
  ALLOCATE_TEMP(rho); 
  ALLOCATE_TEMP(velx);
  ALLOCATE_TEMP(vely);
  ALLOCATE_TEMP(velz);
  ALLOCATE_TEMP(beta1);
  ALLOCATE_TEMP(beta2);
  ALLOCATE_TEMP(beta3);
  ALLOCATE_TEMP(alpha);
# undef ALLOCATE_TEMP

  // coordinates setup
  CCTK_REAL det_x[npoints], det_y[npoints], det_z[npoints];
  for (int i=0;i<ntheta;i++) {
    th=i*dth;
    ct=cos(th);
    st=sin(th);
    for (int j=0;j<nphi;j++) {
      ind=i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
      ph=j*dph;
      cp=cos(ph);
      sp=sin(ph);
      ind2d=i + ntheta*j;
      det_x[ind2d]=sf_centroid_x[sn]+sf_radius[ind]*cp*st;
      det_y[ind2d]=sf_centroid_y[sn]+sf_radius[ind]*sp*st;
      det_z[ind2d]=sf_centroid_z[sn]+sf_radius[ind]*ct;
    }
  }

  const void* interp_coords[3] 
    = { (const void *) det_x,
        (const void *) det_y,
        (const void *) det_z };

  // 3d input arrays
  const CCTK_INT input_array_indices[NUM_INPUT_ARRAYS]
    = { CCTK_VarIndex("ADMBase::gxx"),
        CCTK_VarIndex("ADMBase::gxy"),
        CCTK_VarIndex("ADMBase::gxz"),
        CCTK_VarIndex("ADMBase::gyy"),
        CCTK_VarIndex("ADMBase::gyz"),
        CCTK_VarIndex("ADMBase::gzz"),

        CCTK_VarIndex("HydroBase::velx"),
        CCTK_VarIndex("HydroBase::vely"),
        CCTK_VarIndex("HydroBase::velz"),
        CCTK_VarIndex("HydroBase::rho"),

        CCTK_VarIndex("ADMBase::betax"),
        CCTK_VarIndex("ADMBase::betay"),
        CCTK_VarIndex("ADMBase::betaz"),
        CCTK_VarIndex("ADMBase::alp"),
      };
  for(unsigned int i = 0 ; i < sizeof(input_array_indices)/sizeof(input_array_indices[0]) ; i++) {
    if(input_array_indices[i] < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "couldn't find variable '%s'",
        CCTK_VarName(input_array_indices[i]));
        return -1; /*NOTREACHED*/
    }
  }

  const CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS]
    = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
      };

  // 2d output arrays
  void * output_arrays[NUM_OUTPUT_ARRAYS]
    = { (void *) g11, 
        (void *) g12,
        (void *) g13,
        (void *) g22,
        (void *) g23,
        (void *) g33,

        (void *) velx, 
        (void *) vely,
        (void *) velz,
        (void *) rho, 

        (void *) beta1, 
        (void *) beta2,
        (void *) beta3,
        (void *) alpha, 
      };

  const CCTK_INT operand_indices[NUM_OUTPUT_ARRAYS]
    = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };

  const CCTK_INT opcodes[NUM_OUTPUT_ARRAYS]
    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0 };

  // handles setup
  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "couldn't find interpolator \"%s\"!",
               interpolator_name);

  int param_table_handle = Util_TableCreateFromString(interpolator_pars);

  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                        operand_indices, "operand_indices");
  
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS, 
                        opcodes, "opcodes");
  
  if (param_table_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "bad interpolator parameter(s) \"%s\"!",
               interpolator_pars);
  
  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0)
    CCTK_VWarn(-1, __LINE__, __FILE__, CCTK_THORNSTRING,
        "can't get coordinate system handle for coordinate system \"%s\"!",
               coord_system);

  // actual interpolation call
  ierr = CCTK_InterpGridArrays(cctkGH,
                               DIM, // number of dimensions 
                               operator_handle,
                               param_table_handle,
                               coord_system_handle,
                               interp_npoints,
                               CCTK_VARIABLE_REAL,
                               interp_coords,
                               NUM_INPUT_ARRAYS, // Number of input arrays
                               input_array_indices,
                               NUM_OUTPUT_ARRAYS, // Number of output arrays
                               output_array_types,
                               output_arrays);

  // compute current from primitive values
  for(int i = 0 ; i < npoints ; i++) {
    CCTK_REAL detg, dens, v2, w_lorentz;

    detg = 2*g12[i]*g13[i]*g23[i] + g33[i]*(g11[i]*g22[i] - pow2(g12[i])) - g22[i]*pow2(g13[i]) - g11[i]*pow2(g23[i]);

    v2 = g11[i]*pow2(velx[i]) + g22[i]*pow2(vely[i]) + g33[i]*pow2(velz[i]) + 
         2*g12[i]*velx[i]*vely[i] + 2*g13[i]*velx[i]*velz[i] + 2*g23[i]*vely[i]*velz[i];
    w_lorentz = sqrt(1. / (1. - v2));
    assert(w_lorentz >= 1.);
    dens = sqrt(detg)*rho[i]*w_lorentz;

    jx[i] = - dens * (alpha[i]*velx[i] - beta1[i]);
    jy[i] = - dens * (alpha[i]*vely[i] - beta2[i]);
    jz[i] = - dens * (alpha[i]*velz[i] - beta3[i]);
  }

  if (ierr<0) {
    CCTK_WARN(1,"interpolation screwed up");
    retval = -1;
  }

  ierr = Util_TableDestroy(param_table_handle);

  if (ierr != 0) {
    CCTK_WARN(1,"Could not destroy table");
    retval = -1;
  }

  return retval;

}

/* invert 3-metric to get g^ij - see TAT/TGRtensor/src/tensor.F90*/
static int upstairs_metric(CCTK_REAL gg[3][3],CCTK_REAL gu[3][3])
{
  const CCTK_REAL dtg=(gg[0][0]*gg[1][1]-gg[0][1]*gg[0][1])*gg[2][2]
    -gg[0][0]*gg[1][2]*gg[1][2]
    +2.0*gg[0][1]*gg[0][2]*gg[1][2]
    -gg[0][2]*gg[0][2]*gg[1][1];

  if (dtg==0) {
      return -1;
  }
  if (dtg<0) {
      return -2;
  }

  gu[0][0] = (gg[1][1] * gg[2][2] - gg[1][2] *gg[1][2] ) / dtg;
  gu[0][1] = (gg[0][2] * gg[1][2] - gg[0][1] *gg[2][2] ) / dtg;
  gu[0][2] = (gg[0][1] * gg[1][2] - gg[0][2] *gg[1][1] ) / dtg;
  gu[1][1] = (gg[0][0] * gg[2][2] - gg[0][2] *gg[0][2] ) / dtg;
  gu[1][2] = (gg[0][2] * gg[0][1] - gg[1][2] *gg[0][0] ) / dtg;
  gu[2][2] = (gg[0][0] * gg[1][1] - gg[0][1] *gg[0][1] ) / dtg;

  gu[1][0] = gu[0][1];
  gu[2][0] = gu[0][2];
  gu[2][1] = gu[1][2];

#if 0
  CCTK_REAL tmp;
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      tmp=0;
      for (int a=0;a<3;a++) {
        tmp+=gu[i][a]*gg[a][j];
      }
      fprintf(stderr,"test: %d %d: %g\n",i,j,tmp);
    }
  }
#endif

  return 1;
}



static int get_g_j_local(int i, int j, int ntheta,
                  CCTK_REAL *g11_det,CCTK_REAL *g12_det,CCTK_REAL *g13_det,
                  CCTK_REAL *g22_det,CCTK_REAL *g23_det,CCTK_REAL *g33_det,
                  CCTK_REAL *j1_det,CCTK_REAL *j2_det,CCTK_REAL *j3_det,
                  CCTK_REAL gloc[3][3],CCTK_REAL jloc[3])
{
  /* gloc_ij = g_ij, kloc_ij=K_ij - downstairs index */
  CCTK_INT ind2d=i + ntheta*j;
  gloc[0][0]=g11_det[ind2d];
  gloc[0][1]=g12_det[ind2d];
  gloc[0][2]=g13_det[ind2d];
  gloc[1][0]=g12_det[ind2d];
  gloc[1][1]=g22_det[ind2d];
  gloc[1][2]=g23_det[ind2d];
  gloc[2][0]=g13_det[ind2d];
  gloc[2][1]=g23_det[ind2d];
  gloc[2][2]=g33_det[ind2d];
  
  /* jloc_i - upstairs index */
  jloc[0]=j1_det[ind2d];
  jloc[1]=j2_det[ind2d];
  jloc[2]=j3_det[ind2d];

  return 1;
}

static int drdth_drdph(int i, int j, int ntheta, int nphi,
                int sn,
                CCTK_REAL dth, CCTK_REAL dph,
                CCTK_INT verbose,
                CCTK_INT maxntheta, CCTK_INT maxnphi,
                CCTK_REAL *sf_radius,
                CCTK_REAL *ht, CCTK_REAL *hp)
{
  CCTK_INT ind=0;
  CCTK_REAL htp1,htm1,hpp1,hpm1;

  /* dr/dth dr/dph */
  if (i != 0 && i != ntheta-1) {
    ind=(i+1) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
    htp1=sf_radius[ind];
    ind=(i-1) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
    htm1=sf_radius[ind];
    *ht=(htp1-htm1)/(2.*dth);
    if (verbose>5) {
      fprintf(stderr,"  normal : i=%d j=%d ht=%g\n",i,j,*ht);
    }
  }
  else {
    if (i==0) {
      ind=maxntheta *(j+maxnphi*sn);
      CCTK_REAL ht1=sf_radius[ind];
      ind=1 + maxntheta *(j+maxnphi*sn);
      CCTK_REAL ht2=sf_radius[ind];
      ind=2 + maxntheta *(j+maxnphi*sn);
      CCTK_REAL ht3=sf_radius[ind];
      *ht=-1./(2.*dth)*(3.*ht1-4.*ht2+ht3);
      if (verbose>5) { 
        fprintf(stderr,"  one-sided for i=0: %d,%d ht=%g\n",i,j,*ht);            
      }
    }
    else if (i==ntheta-1) {
      ind= ntheta-1 + maxntheta *(j+maxnphi*sn);
      CCTK_REAL htn=sf_radius[ind];
      ind= ntheta-2 + maxntheta *(j+maxnphi*sn);
      CCTK_REAL htn1=sf_radius[ind];
      ind= ntheta-3 + maxntheta *(j+maxnphi*sn);
      CCTK_REAL htn2=sf_radius[ind];
      *ht=-1./(2.*dth)*(3.*htn-4.*htn1+htn2);
      if (verbose>5) { 
        fprintf(stderr,"  one-sided for i=ntheta-1: %d,%d ht=%g\n",i,j,*ht);            
      }
    }
    else {
      if (verbose>2) CCTK_WARN(1,"unknown theta position");
      return -1;
    }
  }
  if (j !=0 && j != nphi-1) {
    ind=i + maxntheta *((j+1)+maxnphi*sn); // XXX not sf_ntheta!
    hpp1=sf_radius[ind];
    ind=i + maxntheta *((j-1)+maxnphi*sn); // XXX not sf_ntheta!
    hpm1=sf_radius[ind];
    *hp=(hpp1-hpm1)/(2.*dph);
    if (verbose>5) { 
      fprintf(stderr,"  normal %d %d hp=%g\n",i,j,*hp);            
    }
  }
  else {
    if (j==0) {
      ind=i+maxntheta *(maxnphi*sn);
      CCTK_REAL hp1=sf_radius[ind];
      ind=i + maxntheta *(1+maxnphi*sn);
      CCTK_REAL hp2=sf_radius[ind];
      ind=i + maxntheta *(2+maxnphi*sn);
      CCTK_REAL hp3=sf_radius[ind];
      *hp=-1./(2.*dph)*(3.*hp1-4.*hp2+hp3);
      if (verbose>5) { 
        fprintf(stderr,"  one-sided for j=0: %d %d hp=%g\n",i,j,*hp);            
      }
    }
    else if (j==nphi-1) {
      ind= i + maxntheta *(nphi-1+maxnphi*sn);
      CCTK_REAL hpn=sf_radius[ind];
      ind= i + maxntheta *(nphi-2+maxnphi*sn);
      CCTK_REAL hpn1=sf_radius[ind];
      ind= i + maxntheta *(nphi-3+maxnphi*sn);
      CCTK_REAL hpn2=sf_radius[ind];
      *hp=-1./(2.*dph)*(3.*hpn-4.*hpn1+hpn2);
      if (verbose>5) { 
        fprintf(stderr,"  one-sided for j=nphi-1: %d %d hp=%g\n",i,j,*hp);            
      }
    }
    else {
      if (verbose>2) CCTK_WARN(1,"unknown phi position");
      return -2;
    }
  }
  return 1;
}


void outflow_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
 
  if (verbose>0) {
    CCTK_INFO("initialize outflow stuff");
  }
  for (int i=0;i<num_detectors;i++) {
    outflow_flux[i]=0;
  }

  for (int i=0;i<MAX_NUMBER_DETECTORS;i++) {
    file_created[i]=0;
  }
}

void outflow_postrecovery(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose>0) CCTK_INFO("postrecovery: set flag to recreate header");

  for (int i=0;i<MAX_NUMBER_DETECTORS;i++) {
    file_created[i]=0;
  }
}


CCTK_REAL *outflow_allocate_array(CCTK_INT npoints, const char *name)
{
  DECLARE_CCTK_PARAMETERS;

  if (npoints<=0) {
    CCTK_WARN(0,"can't allocate array with npoints <=0");
  }
  if (name==NULL) {
    CCTK_WARN(0,"give a name");
  }

  CCTK_REAL *res;
  res=(CCTK_REAL *) malloc(sizeof(CCTK_REAL)*npoints);
  if (res==NULL) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "%s allocation for npoints=%d failed",name,npoints);
  }
  for (int i=0;i<npoints;i++) res[i]=0;
  const CCTK_REAL mbyt=npoints*sizeof(CCTK_REAL)/(1024.0*1024.0);
  if (verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"allocated array %s with %d elements -> %g MB",name,npoints,mbyt);
  }
  return res;
}


static CCTK_INT have_integrand_memory=0;
static CCTK_REAL *g11_det, *g12_det, *g13_det, *g22_det, *g23_det, *g33_det;
static CCTK_REAL *j1_det, *j2_det, *j3_det;
static CCTK_INT outflow_get_local_memory(CCTK_INT npoints)
{
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose>1) CCTK_INFO("in allocate_memory");

  if (have_integrand_memory==0) {
    if (verbose>0) CCTK_INFO("allocating new memory");
    // metric on detector (covariant rank 2 tensor)
    g11_det=outflow_allocate_array(npoints,"g11_det");
    g12_det=outflow_allocate_array(npoints,"g12_det");
    g13_det=outflow_allocate_array(npoints,"g13_det");
    g22_det=outflow_allocate_array(npoints,"g22_det");
    g23_det=outflow_allocate_array(npoints,"g23_det");
    g33_det=outflow_allocate_array(npoints,"g33_det");
    // current density on detector (vector)
    j1_det=outflow_allocate_array(npoints,"j1_det");
    j2_det=outflow_allocate_array(npoints,"j2_det");
    j3_det=outflow_allocate_array(npoints,"j3_det");
    // update memory allocation flag
    have_integrand_memory=1;
  }
  else {
    if (verbose>1) CCTK_INFO("already allocated memory");
    return 2;
  }

  return 1;
}



void outflow (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT sn, ind, ierr;
  CCTK_REAL ht,hp;
  CCTK_REAL sint,sinp,cost,cosp,rp;

  if (cctk_iteration % compute_every != 0) {
    return;
  }

  /* local memory allocation */
  CCTK_INT npoints=maxntheta*maxnphi;
  ierr=outflow_get_local_memory(npoints);
  if (ierr<0) {
    CCTK_WARN(1,"failed to allocate memory");
    return;
  }

  /* loop over detectors */
  for (int det=0;det<num_detectors;det++)
  {
    sn=surface_index[det];
    if (sn>=sphericalsurfaces_nsurfaces) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                 sn,sphericalsurfaces_nsurfaces);
      continue;
    }
    if (sf_valid[sn]<=0) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "didn't find valid detector surface for sn=%d, det=%d",sn,det);
      continue;
    }
    const CCTK_INT ntheta=sf_ntheta[sn];
    const CCTK_INT nphi=sf_nphi[sn];
    const CCTK_REAL pi=acos(-1.0);
    const CCTK_REAL dth=pi/(ntheta-1);
    const CCTK_REAL dph=2.0*pi/(nphi-1);
    const CCTK_REAL dtp=dth*dph;

    if (verbose>2) {
      CCTK_VInfo(CCTK_THORNSTRING,"ntheta=%d nphi=%d dth=%g dph=%g",
                 ntheta,nphi,dth,dph);
    }

    ierr=get_gab_ja_onto_detector(CCTK_PASS_CTOC, sn,
                              g11_det, g12_det, g13_det,
                              g22_det, g23_det, g33_det,
                              j1_det, j2_det, j3_det);
    if (ierr<0) {
      CCTK_WARN(1,"unable to get g_ab and j^a on the detector. not doing anything.");
      return;
    }

    CCTK_REAL /*phi[3],rvec[3],*/rdn[3]/*,rup[3]*/;
/*    CCTK_REAL phi_x[3], phi_y[3], phi_z[3];*/
    CCTK_REAL gloc[3][3], jloc[3], gup[3][3];
    CCTK_REAL th,ph,dA;
    CCTK_REAL sum=0; // the value of the flux integral

    const CCTK_INT myproc= CCTK_MyProc(cctkGH);

    CCTK_REAL iwtheta=0,iwphi=0,intweight=0;
    /* loop over detector surface */
    for (int i=0;i<ntheta;i++)
    {
      th=i*dth;
      cost=cos(th);
      sint=sin(th);
      if      (i==0 || i==ntheta-1) iwtheta=3.0/8.0;
      else if (i==1 || i==ntheta-2) iwtheta=7.0/6.0;
      else if (i==2 || i==ntheta-3) iwtheta=23.0/24.0;
      else iwtheta=1;

      for (int j=0;j<nphi;j++)
      {
        ph=j*dph;
        cosp=cos(ph);
        sinp=sin(ph);
        // ok I am fairly confused as to why this is necessary in the phi direction which is periodic (roland)
        // in particular see this article: http://www.math.lsa.umich.edu/~rauch/Paper4.pdf
        // frankly a simple trapezoid rule which sums up all terms with the same weight should be perfectly enough for this
	if      (j==0 || j==nphi-1) iwphi=3.0/8.0;
	else if (j==1 || j==nphi-2) iwphi=7.0/6.0;
	else if (j==2 || j==nphi-3) iwphi=23.0/24.0;
	else iwphi=1;
	intweight=iwphi*iwtheta;

        ind=i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
        rp=sf_radius[ind];

        get_g_j_local(i,j,ntheta,
                      g11_det,g12_det,g13_det,
                      g22_det,g23_det,g33_det,
                      j1_det,j2_det,j3_det,
                      gloc,jloc);

        /* invert 3-metric to get g^ij */
        ierr=upstairs_metric(gloc,gup);
        if (ierr<0) {
          if (verbose>2) {
            CCTK_VInfo(CCTK_THORNSTRING,"unable to compute inverse metric at th=%g ph=%g",
                       th,ph);
            if (ierr==-1) {
                CCTK_INFO("detg=0!");
            }
            if (ierr==-2) {
                CCTK_INFO("metric negative");
            }
          }
          continue; // skip this point
        }

        /* get derivatives of r in theta and phi direction */
        ierr=drdth_drdph(i, j, ntheta, nphi, 
                         sn,
                         dth,dph,
                         verbose,
                         maxntheta, maxnphi,
                         sf_radius,
                         &ht, &hp);
        if (ierr<0) {
          CCTK_WARN(1,"derivative computation failed");
          continue;
        }

        /* STEP 1 : COMPUTATION OF APPROXIMATE KILLING VECTOR */
        /* dx^a/dr^b */
#if 0
        CCTK_REAL dxdr=sint*cosp;
        CCTK_REAL dydr=sint*sinp;
        CCTK_REAL dzdr=cost;
#endif
        CCTK_REAL dxdp=rp*sint*(-sinp);
        CCTK_REAL dydp=rp*sint*cosp;
        CCTK_REAL dzdp=0;
        CCTK_REAL dxdt=rp*cosp*cost;
        CCTK_REAL dydt=rp*sinp*cost;
        CCTK_REAL dzdt=rp*(-sint);

#if 0
        /* phi^x - not killing, simply coordinate, but this vector has to be in the surface */
        phi[0]=dxdr*hp+dxdp;
        phi[1]=dydr*hp+dydp;
        phi[2]=dzdr*hp+dzdp;
#endif

        /* STEP 2 : COMPUTATION OF SURFACE NORMAL VECTOR */
        /* d_i theta, d_i phi */
        CCTK_REAL X=rp*cosp*sint;
        CCTK_REAL Y=rp*sinp*sint;
        CCTK_REAL Z=rp*cost;
        /* stay away from Z-axis */
        if (X*X+Y*Y == 0) {
          X=0;
          Y=1e-8;
        }

#if 0
	/* flat space rotational Killing vectors */
	phi_x[0]=0;
	phi_x[1]=-Z;
	phi_x[2]=Y;

	phi_y[0]=Z;
	phi_y[1]=0;
	phi_y[2]=-X;

	phi_z[0]=-Y;
	phi_z[1]=X;
	phi_z[2]=0;
#endif

        CCTK_REAL dtdx=X*1.0*Z*1.0/(sqrt(Y*Y+X*X)*(Z*Z+Y*Y+X*X));
        CCTK_REAL dtdy=Y*1.0*Z*1.0/(sqrt(Y*Y+X*X)*(Z*Z+Y*Y+X*X));
        CCTK_REAL dtdz=-sqrt(Y*Y+X*X)*1.0/(Z*Z+Y*Y+X*X);
        CCTK_REAL dpdx=-Y*1.0/(Y*Y+X*X);
        CCTK_REAL dpdy=X*1.0/(Y*Y+X*X);
        CCTK_REAL dpdz=0;
        CCTK_REAL dhdi[3];
        dhdi[0]=ht*dtdx+hp*dpdx;
        dhdi[1]=ht*dtdy+hp*dpdy;
        dhdi[2]=ht*dtdz+hp*dpdz;
        CCTK_REAL ni[3],fi[3];
        /* ni: unit normal vector */
        ni[0]=cosp*sint;
        ni[1]=sinp*sint;
        ni[2]=     cost;
        fi[0]=ni[0];
        fi[1]=ni[1];
        fi[2]=ni[2];
        for (int idir=0;idir<3;idir++) {
          fi[idir]=fi[idir]-dhdi[idir];
        }

        CCTK_REAL metnorm=0.;
        for (int idir=0;idir<3;idir++) {
          for (int jdir=0;jdir<3;jdir++) {
            metnorm = metnorm+fi[idir]*fi[jdir]*gup[idir][jdir];
          }
        }
        metnorm=sqrt(metnorm);
        for (int idir=0;idir<3;idir++) {
          rdn[idir]=fi[idir]/metnorm;
        }
#if 0
        for (int idir=0;idir<3;idir++) {
          rup[idir]=0;
          for (int jdir=0;jdir<3;jdir++) {
            rup[idir]=rup[idir]+gup[idir][jdir]*rdn[jdir];
          }
        }
        rvec[0]=rup[0];
        rvec[1]=rup[1];
        rvec[2]=rup[2];

        /* check rvec[] and phi[] are orthogonal */
        if (verbose>5) {
          CCTK_REAL tmp=0;
          for (int a=0;a<3;a++) {
            for (int b=0;b<3;b++) {
              tmp+=gloc[a][b]*rvec[a]*phi[b];
            }
          }
          fprintf(stderr,"i=%d j=%d normalization test=%g ==0?\n",i,j,tmp);
        }
#endif

        /* STEP 3 : COMPUTATION OF AREA ELEMENT dA */
        /* two metric tdn_ij: g_ab e_i^a e_j^b
           e^a_b=dx^a/dth^b */
        CCTK_REAL ee[3][2];
        ee[0][0] = dxdt;
        ee[0][1] = dxdp;
        ee[1][0] = dydt;
        ee[1][1] = dydp;
        ee[2][0] = dzdt;
        ee[2][1] = dzdp;

        CCTK_REAL tdn[2][2];
        for (int a=0;a<2;a++) {
          for (int b=0;b<2;b++) {
            tdn[a][b]=0;
            for (int c=0;c<3;c++) {
              for (int d=0;d<3;d++) {
                tdn[a][b]+=gloc[c][d]*ee[c][a]*ee[d][b];
              }
            }
          }
        }

        /* determinant */
        CCTK_REAL tdetg;
        tdetg=tdn[0][0]*tdn[1][1]-tdn[1][0]*tdn[0][1];

        if (tdetg>0) {
          dA=sqrt(tdetg)*dtp*intweight; // area element  (weighted for 4th order)
        }
        else {
          dA=0;
        }

#if 0
        CCTK_REAL xvec[3],yvec[3],zvec[3];
        xvec[0]=1;xvec[1]=0;xvec[2]=0;
        yvec[0]=0;yvec[1]=1;yvec[2]=0;
        zvec[0]=0;zvec[1]=0;zvec[2]=1;
#endif

        /* STEP 4 : FLUX FORMULA */
        for (int a=0;a<3;a++) {
          sum   += jloc[a] * rdn[a] * dA;
        }

        if (verbose>4) {
          fprintf(stderr,"sum=%g\n",sum);
        }
      } // j : phi
    } // i : theta

    if (verbose>0) {
      CCTK_VInfo(CCTK_THORNSTRING,"flux value=%g on detector %d", sum,det);
    }

    outflow_flux[det]=sum;

    /* IO on CPU 0 */
    if (myproc == 0) {
      ierr=Outflow_write_output(CCTK_PASS_CTOC,det, sum);
      if (ierr<0) {
	CCTK_WARN(1,"writing of information to files failed");
      }
    }
  } // det loop over detector number
}
