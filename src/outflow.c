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
#define NGHOSTS 2

#define MAX_NUMBER_DETECTORS 100
static CCTK_INT file_created[MAX_NUMBER_DETECTORS];

static inline CCTK_REAL pow2(CCTK_REAL x) {return x*x;}
static int drdth_drdph(int i, int j,
                int sn,
                CCTK_REAL dth, CCTK_REAL dph,
                CCTK_INT verbose,
                CCTK_INT maxntheta, CCTK_INT maxnphi,
                CCTK_REAL *sf_radius,
                CCTK_REAL *ht, CCTK_REAL *hp);
static int get_ja_onto_detector(CCTK_ARGUMENTS,
                             CCTK_INT sn,
                             CCTK_REAL *jx, CCTK_REAL *jy, CCTK_REAL *jz);
static int get_j_local(int i, int j, int ntheta,
                  CCTK_REAL *j1_det,CCTK_REAL *j2_det,CCTK_REAL *j3_det,
                  CCTK_REAL jloc[3]);
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

/* fills j1...j3 with the interpolated numbers */
static int get_ja_onto_detector(CCTK_ARGUMENTS,
                             CCTK_INT sn,
                             CCTK_REAL *jx, CCTK_REAL *jy, CCTK_REAL *jz)
{
  DECLARE_CCTK_ARGUMENTS; 
  DECLARE_CCTK_PARAMETERS; 
  int ierr, retval = 1;
  CCTK_INT ind,ind2d;
  CCTK_REAL th,ph,ct,st,cp,sp;
  CCTK_INT ntheta,nphi,npoints;
  // auxilliary variables used in constructing j
  static CCTK_REAL *rho = NULL, *velx = NULL, *vely = NULL, *velz = NULL;
  static CCTK_REAL *beta1 = NULL, *beta2 = NULL, *beta3 = NULL, *alpha = NULL;
  static CCTK_REAL *g11 = NULL, *g12 = NULL, *g13 = NULL, *g22 = NULL;
  static CCTK_REAL *g23 =  NULL, *g33 = NULL;

  assert(sn>=0);
  assert(jx); assert(jy); assert(jz);

  ntheta=sf_ntheta[sn]-2*nghoststheta[sn];
  nphi=sf_nphi[sn]-2*nghostsphi[sn];
  npoints=ntheta*nphi;

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"surface %d (%g,%g,%g) nth,nph (%d,%d)",
	       sn,sf_centroid_x[sn],sf_centroid_y[sn],sf_centroid_z[sn],
	       ntheta,nphi);
  }

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
  
  ALLOCATE_TEMP(g11);
  ALLOCATE_TEMP(g12);
  ALLOCATE_TEMP(g13);
  ALLOCATE_TEMP(g22);
  ALLOCATE_TEMP(g23);
  ALLOCATE_TEMP(g33);
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
  const CCTK_INT imin=nghoststheta[sn], imax=sf_ntheta[sn]-nghoststheta[sn]-1;
  const CCTK_INT jmin=nghostsphi[sn], jmax=sf_nphi[sn]-nghostsphi[sn]-1;
  const CCTK_REAL oth=sf_origin_theta[sn];
  const CCTK_REAL oph=sf_origin_phi[sn];
  const CCTK_REAL dth=sf_delta_theta[sn];
  const CCTK_REAL dph=sf_delta_phi[sn];

  CCTK_REAL det_x[npoints], det_y[npoints], det_z[npoints];
  for (int i=imin,n=0;i<=imax;i++,n++) { // theta in [0.5 delta_th, pi-0.5 delta_th]
    th=oth + i * dth;
    ct=cos(th);
    st=sin(th);
    for (int j=jmin,m=0;j<=jmax;j++,m++) { // phi in [0,2pi-delta_phi]
      ind=i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
      ph=oph + j * dph;
      cp=cos(ph);
      sp=sin(ph);
      ind2d=n + ntheta*m;
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

    jx[i] = dens * (alpha[i]*velx[i] - beta1[i]);
    jy[i] = dens * (alpha[i]*vely[i] - beta2[i]);
    jz[i] = dens * (alpha[i]*velz[i] - beta3[i]);
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

static int get_j_local(int i, int j, int ntheta,
                  CCTK_REAL *j1_det,CCTK_REAL *j2_det,CCTK_REAL *j3_det,
                  CCTK_REAL jloc[3])
{
  CCTK_INT ind2d=i + ntheta*j;
  /* jloc_i - upstairs index */
  jloc[0]=j1_det[ind2d];
  jloc[1]=j2_det[ind2d];
  jloc[2]=j3_det[ind2d];

  return 1;
}

static int drdth_drdph(int i, int j, 
                int sn,
                CCTK_REAL dth, CCTK_REAL dph,
                CCTK_INT verbose,
                CCTK_INT maxntheta, CCTK_INT maxnphi,
                CCTK_REAL *sf_radius,
                CCTK_REAL *ht, CCTK_REAL *hp)
{
  CCTK_INT ind;
  CCTK_REAL htp1,htm1,hpp1,hpm1;
  CCTK_REAL htp2,htm2,hpp2,hpm2;

  /* dr/dth dr/dph */
  ind=(i+1) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htp1=sf_radius[ind];
  ind=(i-1) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htm1=sf_radius[ind];
  ind=(i+2) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htp2=sf_radius[ind];
  ind=(i-2) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htm2=sf_radius[ind];
  *ht = (1./12.*htm2-2./3.*htm1+2./3.*htp1-1./12.*htp2)/dth;
  if (verbose>5) {
    fprintf(stderr,"  normal : i=%d j=%d ht=%g\n",i,j,*ht);
  }

  ind=i + maxntheta *((j+1)+maxnphi*sn); // XXX not sf_ntheta!
  hpp1=sf_radius[ind];
  ind=i + maxntheta *((j-1)+maxnphi*sn); // XXX not sf_ntheta!
  hpm1=sf_radius[ind];
  ind=i + maxntheta *((j+2)+maxnphi*sn); // XXX not sf_ntheta!
  hpp2=sf_radius[ind];
  ind=i + maxntheta *((j-2)+maxnphi*sn); // XXX not sf_ntheta!
  hpm2=sf_radius[ind];
  *hp = (1./12.*hpm2-2./3.*hpm1+2./3.*hpp1-1./12.*hpp2)/dph;
  if (verbose>5) { 
    fprintf(stderr,"  normal : i=%d j=%d hp=%g\n",i,j,*hp);            
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
static CCTK_REAL *j1_det, *j2_det, *j3_det;
static CCTK_INT outflow_get_local_memory(CCTK_INT npoints)
{
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose>1) CCTK_INFO("in allocate_memory");

  if (have_integrand_memory==0) {
    if (verbose>0) CCTK_INFO("allocating new memory");
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
    if ( cctk_iteration % compute_every_det[det] != 0 ) {
      continue;
    }
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
    if(nghoststheta[sn] < NGHOSTS || nghostsphi[sn] < NGHOSTS) { // we need at least NGHOSTS ghost zones 
      CCTK_VInfo(CCTK_THORNSTRING,
                 "number of ghost zones for spherical surface %d must be at least %d.",sn,NGHOSTS);
      continue;
    }

    const CCTK_INT imin=nghoststheta[sn], imax=sf_ntheta[sn]-nghoststheta[sn]-1;
    const CCTK_INT jmin=nghostsphi[sn], jmax=sf_nphi[sn]-nghostsphi[sn]-1;
    const CCTK_INT ntheta = imax-imin+1, nphi = jmax-jmin+1;
    const CCTK_REAL oth=sf_origin_theta[sn];
    const CCTK_REAL oph=sf_origin_phi[sn];
    const CCTK_REAL dth=sf_delta_theta[sn];
    const CCTK_REAL dph=sf_delta_phi[sn];
    const CCTK_REAL dtp=dth*dph;

    if (verbose>2) {
      CCTK_VInfo(CCTK_THORNSTRING,"ntheta=%d nphi=%d dth=%g dph=%g",
                 ntheta,nphi,dth,dph);
    }

    ierr=get_ja_onto_detector(CCTK_PASS_CTOC, sn,
                              j1_det, j2_det, j3_det);
    if (ierr<0) {
      CCTK_WARN(1,"unable to get g_ab and j^a on the detector. not doing anything.");
      return;
    }

    CCTK_REAL rdn[3], rhat[3], phihat[3], thetahat[3];
    CCTK_REAL jloc[3];
    CCTK_REAL th,ph;
    CCTK_REAL sum; // the value of the flux integral

    const CCTK_INT myproc= CCTK_MyProc(cctkGH);

    CCTK_REAL iwtheta,iwphi,intweight;
    /* loop over detector surface */
    sum = 0.;
    for (int i=imin,n=0;i<=imax;i++,n++) // theta in [0.5 delta_th, pi-0.5 delta_th]
    {
      th=oth + i * dth;
      cost=cos(th);
      sint=sin(th);
      // weigths from NR(C++,2) 4.1.14 plus extrapolated 1/2 integral (see Maple worksheet)
      if      (i==imin+0 || i==imax-0) iwtheta=13.0/12.0;
      else if (i==imin+1 || i==imax-1) iwtheta= 7.0/ 8.0;
      else if (i==imin+2 || i==imax-2) iwtheta=25.0/24.0;
      else iwtheta=1;

      for (int j=jmin,m=0;j<=jmax;j++,m++) // phi in [0,2pi-delta_phi]
      {
        ph=oph + j * dph;
        cosp=cos(ph);
        sinp=sin(ph);
	iwphi=1; // trapezoid rule
	intweight=iwphi*iwtheta;

        ind=i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
        rp=sf_radius[ind];

        if (verbose>5) {
          fprintf(stderr,"r=%g theta=%g phi=%g\n",rp,th,ph);
        }

        // operates on interpolated values
        get_j_local(n,m,ntheta,
                      j1_det,j2_det,j3_det,
                      jloc);

        // the flat space-like 3d unit vectors
        rhat    [0] =  cosp*sint;rhat    [1] =  sinp*sint;rhat    [2] =  cost;
        thetahat[0] =  cosp*cost;thetahat[1] =  sinp*cost;thetahat[2] = -sint;
        phihat  [0] = -sinp     ;phihat  [1] =  cosp     ;phihat  [2] =     0;

        /* get derivatives of r in theta and phi direction */
        ierr=drdth_drdph(i, j,
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

        // the vector surface element
        for(int idir = 0 ; idir < 3 ; idir++)
        {
          rdn[idir] = pow2(rp)*sint*rhat[idir] - ht*rp*thetahat[idir]
                      - hp*rp*phihat[idir];
        }

        // sum the integral
        for (int a=0;a<3;a++) {
          sum   += jloc[a] * rdn[a] * intweight * dtp;
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
