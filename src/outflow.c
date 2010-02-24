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
#include "util_String.h"

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
#define MAX_NUMBER_DETECTORS 100
#define MAX_NUMBER_EXTRAS 20
#define MAX_NUMBER_TRESHOLDS 20
#define NUM_INPUT_ARRAYS 6+4+4
#define NUM_OUTPUT_ARRAYS 6+4+4
#define NGHOSTS 2
#define ARRAY_INIT_VALUE 0.

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279
#endif

static CCTK_INT file_created[MAX_NUMBER_DETECTORS];
static CCTK_INT fluxdens_file_created[MAX_NUMBER_DETECTORS];
static CCTK_INT extras_file_created[MAX_NUMBER_DETECTORS*MAX_NUMBER_EXTRAS];

static inline CCTK_REAL pow2(CCTK_REAL x) {return x*x;}

/* copied from Multipole */
static void fill_variable(int idx, const char *optstring, void *callback_arg);

static int drdth_drdph(int i, int j,
                int sn,
                CCTK_REAL dth, CCTK_REAL dph,
                CCTK_INT verbose,
                CCTK_INT maxntheta, CCTK_INT maxnphi,
                CCTK_REAL *sf_radius,
                CCTK_REAL *ht, CCTK_REAL *hp);
static int get_ja_w_and_extras_onto_detector(CCTK_ARGUMENTS, CCTK_INT det,
        CCTK_INT num_extras, CCTK_INT extras_ind[MAX_NUMBER_EXTRAS], CCTK_REAL
        *jx, CCTK_REAL *jy, CCTK_REAL *jz, CCTK_REAL *w, CCTK_REAL **extras);
static int get_j_and_w_local(int i, int j, int ntheta, CCTK_REAL
        *j1_det,CCTK_REAL *j2_det,CCTK_REAL *j3_det, CCTK_REAL *w_det, CCTK_REAL
        jloc[3], CCTK_REAL *wloc);
static CCTK_INT outflow_get_local_memory(CCTK_INT npoints);
static int Outflow_write_output(CCTK_ARGUMENTS, CCTK_INT det, CCTK_REAL flux,
        CCTK_REAL w_lorentz, CCTK_REAL *threshold_fluxes);
static int Outflow_write_2d_output(CCTK_ARGUMENTS, const char *varname, CCTK_INT
        det, CCTK_INT *file_created_2d, CCTK_REAL *data_det, CCTK_REAL *w_det,
        CCTK_REAL *surfaceelement_det);
static CCTK_REAL *outflow_allocate_array(CCTK_INT npoints, const char *name);
static CCTK_REAL *get_surface_projection(CCTK_ARGUMENTS, int extra_num);

/* IO */
static int Outflow_write_2d_output(CCTK_ARGUMENTS, const char *varname, CCTK_INT
        det, CCTK_INT *file_created_2d, CCTK_REAL *data_det, CCTK_REAL *w_det,
        CCTK_REAL *surfaceelement_det)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char const *fmode;
  char *filename;
  char format_str_real[2048]; // XXX fixed size
  FILE *file;

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "writing output");
  }

  // check input data
  if (det>=MAX_NUMBER_DETECTORS) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "warn: det=%d, but MAX_NUMBER_DETECTORS=%d, increase",
               det,MAX_NUMBER_DETECTORS);
  }
  assert(surface_index[det] >= 0);
  
  // file mode: append if already written
  fmode = (*file_created_2d>0) ? "a" : "w";

  // filename
  Util_asprintf (&filename, "%s/outflow_surface_det_%d_%s.asc", out_dir, det, varname);
  assert(filename);

  // open file
  file = fopen (filename, fmode);
  if (!file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "write_outflow: Could not open scalar output file '%s'",
                filename);
    free(filename);
    return -1;
  }

  // write header on startup
  if (*file_created_2d<=0) {
    const CCTK_INT sn = surface_index[det];
    const CCTK_INT ntheta=sf_ntheta[sn]-2*nghoststheta[sn];
    const CCTK_INT nphi=sf_nphi[sn]-2*nghostsphi[sn];
    fprintf(file,"# 2d Outflow\n");
    fprintf(file,"# detector no.=%d ntheta=%d nphi=%d\n",det,ntheta,nphi);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:it 2:t 3:x 4:y 5:z 6:%s 7:w_lorentz 8:surface_element\n", varname);
  }

  // write data
  sprintf (format_str_real,
           "%%d\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n", 
           out_format, out_format, out_format, out_format, out_format,
           out_format, out_format);

  const CCTK_INT sn = surface_index[det];
  const CCTK_INT ntheta=sf_ntheta[sn]-2*nghoststheta[sn];
  const CCTK_INT imin=nghoststheta[sn], imax=sf_ntheta[sn]-nghoststheta[sn]-1;
  const CCTK_INT jmin=nghostsphi[sn], jmax=sf_nphi[sn]-nghostsphi[sn]-1;
  const CCTK_REAL oth=sf_origin_theta[sn];
  const CCTK_REAL oph=sf_origin_phi[sn];
  const CCTK_REAL dth=sf_delta_theta[sn];
  const CCTK_REAL dph=sf_delta_phi[sn];

  CCTK_REAL th, ph, ct,st, cp,sp,rp;
  CCTK_REAL det_x, det_y, det_z;
  for (int i=imin,n=0;i<=imax;i++,n++) { // theta in [0.5 delta_th, pi-0.5 delta_th]
    th=oth + i * dth;
    ct=cos(th);
    st=sin(th);

    for (int j=jmin,m=0;j<=jmax;j++,m++) { // phi in [0,2pi-delta_phi]
      int ind = i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
      int ind2d = n + ntheta*m;

      ph=oph + j * dph;
      cp=cos(ph);
      sp=sin(ph);
      rp=rad_rescale[det]*sf_radius[ind];

      if(output_relative_coordinates) {
        det_x=rp*cp*st;
        det_y=rp*sp*st;
        det_z=rp*ct;
      } else {
        det_x=sf_centroid_x[sn]+rp*cp*st;
        det_y=sf_centroid_y[sn]+rp*sp*st;
        det_z=sf_centroid_z[sn]+rp*ct;
      }

      fprintf(file, format_str_real, cctk_iteration, cctk_time,
              det_x,det_y,det_z, data_det[ind2d], w_det[ind2d],
              surfaceelement_det[ind2d]);
    }
    /* repeat first angle to ge a closed surface in gnuplot */
    int ind = i + maxntheta *(jmin+maxnphi*sn); // XXX not sf_ntheta!
    int ind2d = n + ntheta*0;
    ph=oph + jmin * dph;
    cp=cos(ph);
    sp=sin(ph);
    rp=rad_rescale[det]*sf_radius[ind];
    if(output_relative_coordinates) {
      det_x=rp*cp*st;
      det_y=rp*sp*st;
      det_z=rp*ct;
    } else {
      det_x=sf_centroid_x[sn]+rp*cp*st;
      det_y=sf_centroid_y[sn]+rp*sp*st;
      det_z=sf_centroid_z[sn]+rp*ct;
    }

    fprintf(file, format_str_real, cctk_iteration, cctk_time, det_x,det_y,det_z,
            data_det[ind2d], w_det[ind2d], surfaceelement_det[ind2d]);

    fprintf(file, "\n"); /* create a grid edge for gnuplot */
  }
  fprintf(file, "\n"); /* create a block for gnuplot */

  fclose(file); 
  free(filename);
  *file_created_2d = 1;

  return 1;
}

static int Outflow_write_output(CCTK_ARGUMENTS, CCTK_INT det, CCTK_REAL flux,
        CCTK_REAL w_lorentz, CCTK_REAL *threshold_fluxes)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char const *fmode;
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[2048]; // XXX fixed size
  int thresh, col;
  FILE *file;

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "writing output");
  }

  // file mode: append if already written
  if (det>=MAX_NUMBER_DETECTORS) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "warn: det=%d, but MAX_NUMBER_DETECTORS=%d, increase",
               det,MAX_NUMBER_DETECTORS);
  }
  fmode = (file_created[det]>0) ? "a" : "w";

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
    free(filename);
    return -1;
  }

  // write header on startup
  if (file_created[det]<=0) {
    fprintf(file,"# Outflow\n");
    fprintf(file,"# detector no.=%d\n",det);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:it 2:t 3:flux 4:avg(w_lorentz)");
    if(num_thresholds > 0) {
      col = 5;
      for(thresh = 0 ; thresh < num_thresholds ; thresh++) {
        fprintf(file," %d:w>=%g",col++,threshold[thresh]);
      }
    }
    fprintf(file,"\n");
  }

  // write data
  sprintf (format_str_real,
           "%%d\t%%%s\t%%%s\t%%%s", 
           out_format,out_format,out_format);
  fprintf(file, format_str_real, cctk_iteration, cctk_time, flux, w_lorentz);
  sprintf (format_str_real, "\t%%%s", out_format);
  for(thresh = 0 ; thresh < num_thresholds ; thresh++) {
    fprintf(file,format_str_real,threshold_fluxes[thresh]);
  }
  fprintf(file,"\n");

  fclose(file); 
  free(filename);

  if (file_created[det]==0) {
    file_created[det]=1;
  }
  return 1;
}

/* return pointer to surface_projection_<extra_num> */
static CCTK_REAL *get_surface_projection(CCTK_ARGUMENTS, int extra_num)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_REAL *retval = NULL;
  switch(extra_num)
  {
    case  0: retval = surface_projection_0; break;
    case  1: retval = surface_projection_1; break;
    case  2: retval = surface_projection_2; break;
    case  3: retval = surface_projection_3; break;
    case  4: retval = surface_projection_4; break;
    case  5: retval = surface_projection_5; break;
    case  6: retval = surface_projection_6; break;
    case  7: retval = surface_projection_7; break;
    case  8: retval = surface_projection_8; break;
    case  9: retval = surface_projection_9; break;
    case 10: retval = surface_projection_10; break;
    case 11: retval = surface_projection_11; break;
    case 12: retval = surface_projection_12; break;
    case 13: retval = surface_projection_13; break;
    case 14: retval = surface_projection_14; break;
    case 15: retval = surface_projection_15; break;
    case 16: retval = surface_projection_16; break;
    case 17: retval = surface_projection_17; break;
    case 18: retval = surface_projection_18; break;
    case 19: retval = surface_projection_19; break;
    default: CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "invalid extra variable number %d passed. Must be less than 20",
                extra_num);
             break;
  }

  return retval;
}

/* fills j1...j3,w and the extras with the interpolated numbers */
static int get_ja_w_and_extras_onto_detector(CCTK_ARGUMENTS, CCTK_INT det,
        CCTK_INT num_extras, CCTK_INT extras_ind[MAX_NUMBER_EXTRAS], CCTK_REAL
        *jx, CCTK_REAL *jy, CCTK_REAL *jz, CCTK_REAL *w, CCTK_REAL **extras)
{
  DECLARE_CCTK_ARGUMENTS; 
  DECLARE_CCTK_PARAMETERS; 
  int ierr;
  CCTK_INT sn,ind,ind2d;
  CCTK_REAL th,ph,ct,st,cp,sp,rp;
  CCTK_INT ntheta,nphi,npoints;
  // auxilliary variables used in constructing j
  static CCTK_REAL *rho = NULL, *velx = NULL, *vely = NULL, *velz = NULL;
  static CCTK_REAL *beta1 = NULL, *beta2 = NULL, *beta3 = NULL, *alpha = NULL;
  static CCTK_REAL *g11 = NULL, *g12 = NULL, *g13 = NULL, *g22 = NULL;
  static CCTK_REAL *g23 =  NULL, *g33 = NULL;

  assert(det>=0);
  assert(det<num_detectors);
  assert(jx); assert(jy); assert(jz);
  assert(w); assert(extras);
  for (int i=0; i<num_extras; i++)
      assert(extras[i]);

  sn = surface_index[det];
  assert(sn>=0);

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

  // allocate memory for auxilliary arrays (of maximum possible size)
# define ALLOCATE_TEMP(name) \
  if(name == NULL) \
    name = outflow_allocate_array(maxntheta*maxnphi, #name); \
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
      if (override_radius[det]) {
        rp = radius[det];
        assert(rp > 0.);
      } else {
        rp = rad_rescale[det]*sf_radius[ind];
      }
      det_x[ind2d]=sf_centroid_x[sn]+rp*cp*st;
      det_y[ind2d]=sf_centroid_y[sn]+rp*sp*st;
      det_z[ind2d]=sf_centroid_z[sn]+rp*ct;
    }
  }

  const void* interp_coords[3] 
    = { (const void *) det_x,
        (const void *) det_y,
        (const void *) det_z };

  // 3d input arrays
  CCTK_INT input_array_indices[NUM_INPUT_ARRAYS + MAX_NUMBER_EXTRAS]
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
  for(int i = 0 ; i < num_extras ; i++) {
     input_array_indices[NUM_INPUT_ARRAYS + i] = extras_ind[i];
  }
  for(int i = 0 ; i < NUM_INPUT_ARRAYS + num_extras ; i++) {
    if(input_array_indices[i] < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "couldn't find variable '%s'",
        CCTK_VarName(input_array_indices[i]));
        return -1; /*NOTREACHED*/
    }
  }

  CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS + num_extras ; i++) {
    output_array_types[i] = CCTK_VARIABLE_REAL;
  }

  // 2d output arrays
  void * output_arrays[NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS]
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
  for(int i = 0 ; i < num_extras ; i++) {
     output_arrays[NUM_OUTPUT_ARRAYS + i] = extras[i];
  }

  CCTK_INT operand_indices[NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS + num_extras  ; i++) {
    operand_indices[i] = i;
  }

  CCTK_INT opcodes[NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS + num_extras  ; i++) {
    opcodes[i] = 0;
  }

  // handles setup
  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "couldn't find interpolator \"%s\"!",
               interpolator_name);

  int param_table_handle = Util_TableCreateFromString(interpolator_pars);
  if (param_table_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "bad interpolator parameter(s) \"%s\"!",
               interpolator_pars);
  }
  
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS + num_extras,
                        operand_indices, "operand_indices");
  
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS + num_extras, 
                        opcodes, "opcodes");
  
  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "can't get coordinate system handle for coordinate system \"%s\"!",
               coord_system);
  }

  // actual interpolation call
  ierr = CCTK_InterpGridArrays(cctkGH,
                               DIM, // number of dimensions 
                               operator_handle,
                               param_table_handle,
                               coord_system_handle,
                               interp_npoints,
                               CCTK_VARIABLE_REAL,
                               interp_coords,
                               NUM_INPUT_ARRAYS + num_extras, // Number of input arrays
                               input_array_indices,
                               NUM_OUTPUT_ARRAYS + num_extras, // Number of output arrays
                               output_array_types,
                               output_arrays);

  if (ierr<0) {
    CCTK_WARN(1,"interpolation screwed up");
    Util_TableDestroy(param_table_handle);
    return -1;
  }

  ierr = Util_TableDestroy(param_table_handle);
  if (ierr != 0) {
    CCTK_WARN(1,"Could not destroy table");
    return -1;
  }

  // compute current from primitive values
  for(int i = 0 ; i < interp_npoints ; i++) {
    CCTK_REAL detg, dens, v2, w_lorentz;

    detg = 2*g12[i]*g13[i]*g23[i] + g33[i]*(g11[i]*g22[i] - pow2(g12[i])) -
        g22[i]*pow2(g13[i]) - g11[i]*pow2(g23[i]);
    if( detg < 0. ) 
    {
        static CCTK_INT last_warned = -1;

        if(verbose > 1 || (verbose > 0 && last_warned != cctk_iteration))
        {
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "%s: Metric determinant in iteration %6d at %15.6g,%15.6g,%15.6g is :" 
                "%15.6g from data g = [%15.6g,%15.6g,%15.6g,%15.6g,%15.6g,%15.6g]",
                __func__, cctk_iteration, det_x[i],det_y[i],det_z[i], detg, 
                g11[i],g12[i],g13[i],g22[i],g23[i],g33[i]);
          last_warned = cctk_iteration;
        }

        detg = 1.;
    }

    v2 = g11[i]*pow2(velx[i]) + g22[i]*pow2(vely[i]) + g33[i]*pow2(velz[i]) +
        2*g12[i]*velx[i]*vely[i] + 2*g13[i]*velx[i]*velz[i] +
        2*g23[i]*vely[i]*velz[i];

    w_lorentz = sqrt(1. / (1. - v2));
    if( w_lorentz < 1. || v2 > 1 ) 
    {
        static CCTK_INT last_warned = -1;

        if(verbose > 1 || (verbose > 0 && last_warned != cctk_iteration))
        {
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "%s: Unphysical Lorentz factor %15.6g, v2 = %15.6g for data "
                "g = [%15.6g,%15.6g,%15.6g,%15.6g,%15.6g,%15.6g] "
                "vel = [%15.6g,%15.6g,%15.6g] occured in iteration %d at location [%15.6g,%15.6g,%15.6g]",
                __func__, w_lorentz,v2, g11[i],g12[i],g13[i],g22[i],g23[i],g33[i],
                velx[i],vely[i],velz[i], cctk_iteration,
                det_x[i],det_y[i],det_z[i]);
          last_warned = cctk_iteration;
        }

        w_lorentz = 1.;
    }
    dens = sqrt(detg)*rho[i]*w_lorentz;

    jx[i] = dens * (alpha[i]*velx[i] - beta1[i]);
    jy[i] = dens * (alpha[i]*vely[i] - beta2[i]);
    jz[i] = dens * (alpha[i]*velz[i] - beta3[i]);
    w[i]  = w_lorentz;
  }

  return interp_npoints;

}

static int get_j_and_w_local(int i, int j, int ntheta,
                  CCTK_REAL *j1_det,CCTK_REAL *j2_det,CCTK_REAL *j3_det,
                  CCTK_REAL *w_det, CCTK_REAL jloc[3], CCTK_REAL *wloc)
{
  CCTK_INT ind2d=i + ntheta*j;
  /* jloc_i - upstairs index */
  jloc[0]=j1_det[ind2d];
  jloc[1]=j2_det[ind2d];
  jloc[2]=j3_det[ind2d];
  *wloc  =w_det[ind2d];

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
    fluxdens_file_created[i]=0;
    for (int j=0;j<MAX_NUMBER_EXTRAS;j++) {
      extras_file_created[i*MAX_NUMBER_EXTRAS+j]=0;
    }
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


static CCTK_REAL *j1_det, *j2_det, *j3_det, *w_det, *fluxdens_det, *surfaceelement_det;
static CCTK_INT outflow_get_local_memory(CCTK_INT npoints)
{
  DECLARE_CCTK_PARAMETERS;
  
  static CCTK_INT have_integrand_memory=0;

  if (verbose>1) CCTK_INFO("in allocate_memory");

  if (have_integrand_memory==0) {
    if (verbose>0) CCTK_INFO("allocating new memory");
    // current density on detector (vector)
    j1_det=outflow_allocate_array(npoints,"j1_det");
    j2_det=outflow_allocate_array(npoints,"j2_det");
    j3_det=outflow_allocate_array(npoints,"j3_det");
    w_det =outflow_allocate_array(npoints,"w_det");
    fluxdens_det =outflow_allocate_array(npoints,"fluxdens_det");
    surfaceelement_det =outflow_allocate_array(npoints,"surfaceelement_det");
    // update memory allocation flag
    have_integrand_memory=1;
  }
  else {
    if (verbose>1) CCTK_INFO("already allocated memory");
    return 2;
  }

  return 1;
}

/* callback routine to add one extra variable to be output */
static void fill_variable(int idx, const char *optstring, void *callback_arg)
{
  assert(idx >= 0);
  assert(callback_arg);

  CCTK_INT *extras_ind = (CCTK_INT * ) callback_arg;

  if(extras_ind[MAX_NUMBER_EXTRAS-1] != -1) { /* no more free slots */
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "too many extra variables, ignoring variable '%s'.",
               CCTK_VarName(idx));
    return;
  }

  /* find the first free slot in extras_ind and store the new index in it */
  for(int i = 0 ; i < MAX_NUMBER_EXTRAS ; i++)
  {
    if(extras_ind[i] == -1) {
      extras_ind[i] = idx;
      break;
    }
  }
}

void outflow (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT sn, ind, ind2d, ierr;
  CCTK_REAL ht,hp;
  CCTK_REAL sint,sinp,cost,cosp,rp;
  CCTK_INT interp_npoints;
  // variables related to the extra projected grid functions
  CCTK_INT num_extras;
  CCTK_INT extras_ind[MAX_NUMBER_EXTRAS];
  CCTK_REAL *extras[MAX_NUMBER_EXTRAS];

  /* local memory allocation */
  CCTK_INT maxnpoints=maxntheta*maxnphi;
  ierr=outflow_get_local_memory(maxnpoints);
  if (ierr<0) {
    CCTK_WARN(1,"failed to allocate memory");
    return;
  }

  /* clear the grid arrays */
  for(int i = 0 ; i < maxntheta*maxnphi*num_detectors ; i++)
    fluxdens_projected[i] = w_lorentz_projected[i] = ARRAY_INIT_VALUE;
  for(int e = 0 ; e < MAX_NUMBER_EXTRAS ; e++)
  {
    CCTK_REAL *surface_projection = get_surface_projection(CCTK_PASS_CTOC, e);
    for(int i = 0 ; i < maxntheta*maxnphi*num_detectors ; i++)
    {
      surface_projection[i] = ARRAY_INIT_VALUE;
    }
  }

  /* parse variables string and allocate temporary memory for them */
  if(!CCTK_Equals(extra_variables, "")) {
    for(int i = 0 ; i < MAX_NUMBER_EXTRAS ; i++) /* initialize so that we can count later */
    {
      extras_ind[i] = -1;
    }

    ierr = CCTK_TraverseString(extra_variables, fill_variable, extras_ind,
            CCTK_GROUP_OR_VAR);
    assert(ierr > 0);
  
    for(num_extras = 0 ; num_extras < MAX_NUMBER_EXTRAS ; num_extras++) /* count valid indices */
    {
      if(extras_ind[num_extras] == -1)
        break;
      extras[num_extras] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*maxnpoints);
      assert(extras[num_extras]);
    }
  } else {
    num_extras = 0;
  }

  /* loop over detectors */
  for (int det=0;det<num_detectors;det++)
  {
    CCTK_INT my_compute_every;

    /* check parameters and decide if we have to do anythin */
    if ( compute_every_det[det] >= 0 ) {
        my_compute_every = compute_every_det[det];
    } else {
        my_compute_every = compute_every;
    }
    if ( my_compute_every == 0 || cctk_iteration % my_compute_every != 0 ) {
      continue;
    }

    sn=surface_index[det];
    if (sn < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d is invalid for detector %d", sn,det);
      continue;
    } else if (sn>=sphericalsurfaces_nsurfaces) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                 sn,sphericalsurfaces_nsurfaces);
      continue;
    }
    if (sf_valid[sn]<=0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "didn't find valid detector surface for sn=%d, det=%d",sn,det);
      continue;
    }

    if(nghoststheta[sn] < NGHOSTS || nghostsphi[sn] < NGHOSTS) { // we need at least NGHOSTS ghost zones 
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
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

    interp_npoints=get_ja_w_and_extras_onto_detector(CCTK_PASS_CTOC, det, num_extras,
                    extras_ind, j1_det, j2_det, j3_det, w_det, extras);
    if (interp_npoints<0) {
      CCTK_WARN(1,"unable to get g_ab, j^a and the extra variables onto the detector. not doing anything.");
      continue;
    }
    if (interp_npoints==0) {
      /* nothing to do (we are not cpu 0) */
      if (verbose > 1) {
        CCTK_VInfo(CCTK_THORNSTRING, "I have nothing to do for detector %d", det);
      }
      continue;
    }
    if (verbose > 1) {
      CCTK_VInfo(CCTK_THORNSTRING, "integrating detector %d", det);
    }

    CCTK_REAL rdn[3], rhat[3], phihat[3], thetahat[3];
    CCTK_REAL jloc[3], wloc;
    CCTK_REAL th,ph;
    CCTK_REAL sum, sum_thresh[MAX_NUMBER_TRESHOLDS], sum_w_lorentz; // the value of the flux integral

    CCTK_REAL iwtheta,iwphi,intweight;
    /* init integration vars */
    sum = sum_w_lorentz = 0.;
    for (int t = 0 ; t < MAX_NUMBER_TRESHOLDS ; t++) {
      sum_thresh[t] = 0.;
    }

    /* loop over detector surface */
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
        if (override_radius[det]) {
            rp = radius[det];
            assert(rp > 0.);
        } else {
          rp=rad_rescale[det]*sf_radius[ind];
        }

        if (verbose>5) {
          fprintf(stderr,"r=%g theta=%g phi=%g\n",rp,th,ph);
        }

        // operates on interpolated values
        get_j_and_w_local(n,m,ntheta,
                      j1_det,j2_det,j3_det,
                      w_det,jloc,&wloc);

        // the flat space-like 3d unit vectors
        rhat    [0] =  cosp*sint;rhat    [1] =  sinp*sint;rhat    [2] =  cost;
        thetahat[0] =  cosp*cost;thetahat[1] =  sinp*cost;thetahat[2] = -sint;
        phihat  [0] = -sinp     ;phihat  [1] =  cosp     ;phihat  [2] =     0;

        /* get derivatives of r in theta and phi direction */
        if (override_radius[det]) {
          ht = hp = 0.; /* spherical */
        } else {
          ierr=drdth_drdph(i, j, sn, dth,dph, verbose, maxntheta, maxnphi,
                  sf_radius, &ht, &hp);
	  ht=rad_rescale[det]*ht;
	  hp=rad_rescale[det]*hp;
          if (ierr<0) {
            CCTK_WARN(1,"derivative computation failed");
            continue;
          }
        }

        // the vector surface element
        CCTK_REAL mag_rdn = 0;
        for(int idir = 0 ; idir < 3 ; idir++)
        {
          rdn[idir] = pow2(rp)*sint*rhat[idir] - ht*rp*thetahat[idir] -
                      hp*rp*phihat[idir];
          mag_rdn += pow2(rdn[idir]);
        }
        mag_rdn = sqrt(mag_rdn);

        // sum the integral
        CCTK_REAL df = 0, fluxdens_temp = 0;
        for (int a=0;a<3;a++) {
          df += jloc[a] * rdn[a] * intweight * dtp;
          fluxdens_temp += jloc[a] * rdn[a]/mag_rdn;
        }
        fluxdens_det[n + ntheta*m] = fluxdens_temp; /* store flux density for output */
        surfaceelement_det[n + ntheta*m] = mag_rdn * dtp; /* area element */

        // flux
        sum += df;
        if (verbose>4) {
          fprintf(stderr,"sum=%g\n",sum);
        }

        // Lorentz factor
        sum_w_lorentz += wloc * intweight * dtp;

        for(int t = 0 ; t < num_thresholds ; t++)
        {
          if(threshold[t] == -1) {
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "threshold %d is set to -1, ignoring.", t);
            continue;
          }
          if(wloc >= threshold[t]) {
            sum_thresh[t] += df;
          }
          if (verbose>4) {
            fprintf(stderr,"sum_thresh[%d]=%g\n",t,sum_thresh[t]);
          }
        }

      } // j : phi
    } // i : theta
    sum_w_lorentz /= 4*M_PI; // average w_lorentz

    if (verbose>0) {
      CCTK_VInfo(CCTK_THORNSTRING,"flux value=%g on detector %d", sum,det);
    }

    outflow_flux[det]=sum;

    /* store results in grid arrays, translating indices on the way */
    /* we fill in only the upper left corner of the grid array */
    for (int i=imin,n=0;i<=imax;i++,n++) // theta in [0.5 delta_th, pi-0.5 delta_th]
    {
      for (int j=jmin,m=0;j<=jmax;j++,m++) // phi in [0,2pi-delta_phi]
      {
        ind = i + maxntheta * (j + maxnphi*det);
        ind2d = n + ntheta * m;
        fluxdens_projected[ind] = fluxdens_det[ind2d];
        w_lorentz_projected[ind] = w_det[ind2d];
      }
    }
    for(int e = 0 ; e < num_extras ; e++)
    {
      CCTK_REAL *surface_projection = get_surface_projection(CCTK_PASS_CTOC, e);
      for (int i=imin,n=0;i<=imax;i++,n++) // theta in [0.5 delta_th, pi-0.5 delta_th]
      {
        for (int j=jmin,m=0;j<=jmax;j++,m++) // phi in [0,2pi-delta_phi]
        {
          ind = i + maxntheta * (j + maxnphi*det);
          ind2d = n + ntheta * m;
          surface_projection[ind] = extras[e][ind2d];
        }
      }
    }

    /* IO (we only get here if we are CPU #0) */
    ierr=Outflow_write_output(CCTK_PASS_CTOC,det, sum, sum_w_lorentz, sum_thresh);
    if (ierr<0) {
      CCTK_WARN(1,"writing of information to files failed");
    }
    if (output_2d_data) {
      ierr=Outflow_write_2d_output(CCTK_PASS_CTOC, "fluxdens", det,
              fluxdens_file_created, fluxdens_det, w_det, surfaceelement_det);
      if (ierr<0) {
        CCTK_WARN(1,"writing of fluxdens information to files failed");
      }
      for(int i = 0 ; i < num_extras ; i++)
      {
        ierr=Outflow_write_2d_output(CCTK_PASS_CTOC,
                CCTK_VarName(extras_ind[i]), det,
                &extras_file_created[det*MAX_NUMBER_EXTRAS+i], extras[i],
                w_det, surfaceelement_det);
        if (ierr<0) {
          CCTK_WARN(1,"writing of extras information to files failed");
        }
      }
    }
  } // det loop over detector number

  /* free temporary memory */
  for(int i = 0 ; i < num_extras ; i++)
  {
    free(extras[i]);
  }
}
