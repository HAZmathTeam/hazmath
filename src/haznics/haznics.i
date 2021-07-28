%module haznics 
%{
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION 
#define SWIG_FILE_WITH_INIT
//#define SWIG_PYTHON_SILENT_MEMLEAK
#include <numpy/arrayobject.h>
#include "hazmath.h"
%}

%include "numpy.ii" 

%init%{
import_array();
%}

%numpy_typemaps(short, NPY_SHORT, SHORT)
%numpy_typemaps(double, NPY_DOUBLE, REAL)

/*%typemap(in) void (*func_ptr)(double*, double*, void*) {
  $1 = $input;
}*/

%include "macro.h"
%include "mesh.h"
%include "amr.h"
%include "sparse.h"
%include "dense.h"
%include "vec.h"
%include "fem.h"
%include "solver.h"
%include "nonlinear.h"
%include "timestep.h"
%include "param.h"
//%include "eigen.h" 
%include "graphs.h"
%include "hazmath.h"


%nodefaultdtor dvector;
%extend dvector{
   ~dvector() {
      dvec_free(self);
   }
}

%nodefaultdtor dCSRmat;
%extend dCSRmat{
   ~dCSRmat() {
        dcsr_free(self);
   }
}

%nodefaultdtor AMG_data;
%extend AMG_data{
   ~AMG_data() {
        amg_data_free(self, NULL);
   }
}

%nodefaultdtor HX_curl_data;
%extend HX_curl_data{
   ~HX_curl_data() {
        HX_curl_data_free(self, 1);
   }
}

%nodefaultdtor HX_div_data;
%extend HX_div_data{
   ~HX_div_data() {
        HX_div_data_free(self, 1);
   }
}

%nodefaultdtor precond_data;
%extend precond_data{
   ~precond_data() {
        precond_data_free(self);
   }
}

%nodefaultdtor precond_ra_data;
%extend precond_ra_data{
   ~precond_ra_data() {
        precond_ra_data_free(self);
   }
}

/*
// not sure about this one yet
%nodefaultdtor AMG_param;
%extend AMG_param{
   ~AMG_param() {
        amli_coef_free(self);
   }
}
*/

%nodefaultctor input_param;
%extend input_param{
   input_param() {
        input_param *inparam = (input_param *) malloc(sizeof(input_param));
        param_input_init(inparam);
        return inparam;
   }
   input_param(const char *filenm) {
        input_param *inparam = (input_param *) malloc(sizeof(input_param));
        param_input(filenm, inparam);
        return inparam;
   }
}

/* need to add destructors for other structs!! */


%extend precond{
    void apply(dvector *lhs, dvector* rhs) { 
      printf("lhs-> row %d\n", lhs->row); 
      apply_precond(lhs->val, rhs->val, $self); 
    } 
    precond_data* precond_data(){
      return (precond_data*)$self->data;
    }
}

%extend block_dCSRmat{
  void init(int n, int m) {
    bdcsr_alloc(n,m,$self);
  } 
  void debugPrint(){
    printf("n, m = %d, %d \n", $self->bcol, $self->brow); 
  }
  dCSRmat* get(int i, int j) { 
      return ($self)->blocks[i + $self->brow*j]; 
    }
  void set(int i, int j, dCSRmat* mat) { 
      ($self)->blocks[i + $self->brow*j] = mat; 
    }

    
};  


%apply (double* IN_ARRAY1, int DIM1) {(double* A, int nnz)};
%apply (int* IN_ARRAY1, int DIM1) {(int* ja, int nnz2)};
%apply (int* IN_ARRAY1, int DIM1) {(int* ia, int n)};
dCSRmat* create_matrix(double *A, int nnz, int *ja, int nnz2, int *ia, int n);   
/* should clear the typemaps */ 


%apply (double* IN_ARRAY1, int DIM1) {(double* x, int n)};
dvector* create_dvector(double *x, int n); 
%clear (double* a, int n);


/* these three below should probably be removed */ 
// input_param* create_input_param();
// AMG_param* create_AMG_param(input_param *in_param);
// linear_itsolver_param* create_linear_itsolver_param(input_param *in_param);
/* this is here because helper functions seems to not be available in the library */
/* NB: it will produce warnings in haznicswrap - this should be fixed later */
void apply_precond(REAL *r, REAL *z, precond *pc);
precond* create_precond_amg(dCSRmat *A, AMG_param *amgparam);
precond* set_precond(void *data, void (*foo)(REAL*, REAL*, void*));
precond* create_precond(dCSRmat *A, AMG_param *amgparam);
precond* create_precond_famg(dCSRmat *A, dCSRmat *M, AMG_param *amgparam);
precond* create_precond_ra(dCSRmat *A, dCSRmat *M, REAL s_frac_power, REAL t_frac_power, REAL alpha, REAL beta, REAL scaling_a, REAL scaling_m, AMG_param *amgparam);
precond* create_precond_hxcurl(dCSRmat *Acurl, dCSRmat *Pcurl, dCSRmat *Grad, SHORT prectype, AMG_param *amgparam);
precond* create_precond_hxdiv_3D(dCSRmat *Adiv, dCSRmat *P_div, dCSRmat *Curl, dCSRmat *P_curl, SHORT prectype, AMG_param *amgparam);
precond* create_precond_hxdiv_2D(dCSRmat *Adiv,dCSRmat *P_div, dCSRmat *Curl, SHORT prectype, AMG_param *amgparam);
INT get_poles_no(precond *pc);

%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* vec1),
                                      (int len2, double* vec2)}
//%apply (precond* pc_) {(precond* pc)}
%rename (apply_precond) my_apply_precond;
%exception my_apply_precond {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
    void my_apply_precond(int len1, double* vec1, int len2, double* vec2, precond* pc) {
        if (len1 != len2) {
            PyErr_Format(PyExc_ValueError,
                         "Arrays of lengths (%d,%d) given",
                         len1, len2);
        }
        if (!pc) {
            PyErr_Format(PyExc_ValueError,
                         "Preconditioner class not set up");
        }
        if (!pc->data) {
            PyErr_Format(PyExc_ValueError,
                         "Preconditioner data not set up");
        }
        if (!pc->fct) {
            PyErr_Format(PyExc_ValueError,
                         "Preconditioner matvec not set up");
        }
        apply_precond(vec1, vec2, pc);
    }
%}

/* callback function as constant?
%constant void precond_amg(double*, double*, void*);
%constant void precond_amli(double*, double*, void*);
%constant void precond_nl_amli(double*, double*, void*);
%constant void precond_amg_add(double*, double*, void*);
 */

/* callback function directive
%callback("%s_cb");
void precond_amg(double*, double*, void*);
void precond_amli(double*, double*, void*);
void precond_nl_amli(double*, double*, void*);
void precond_amg_add(double*, double*, void*);
%nocallback;
 */




