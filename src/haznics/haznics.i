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

   // return numpy array
   PyObject* to_ndarray() {
        npy_intp dims[1]; dims[0] = self->row;
        PyObject *array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void*)(self->val));
        if(!array) {
            PyErr_SetString(PyExc_MemoryError, "Array not allocated!");
            // SWIG_fail; // for some reason swig says this label is not defined
        }
        return array;
   }
}

%nodefaultdtor ivector;
%extend ivector{
   ~ivector() {
        ivec_free(self);
   }

   // return numpy array
   PyObject* to_ndarray() {
        npy_intp dims[1]; dims[0] = self->row;
        PyObject *array = PyArray_SimpleNewFromData(1, dims, NPY_INT, (void*)(self->val));
        if(!array) {
            PyErr_SetString(PyExc_MemoryError, "Array not allocated!");
            // SWIG_fail; // for some reason swig says this label is not defined
        }
        return array;
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

%nodefaultctor AMG_param;
%extend AMG_param{
   AMG_param() {
        AMG_param *amgparam = (AMG_param *) malloc(sizeof(AMG_param));
        param_amg_init(amgparam);
        return amgparam;
   }
}

%nodefaultctor linear_itsolver_param;
%extend linear_itsolver_param{
   linear_itsolver_param() {
        linear_itsolver_param *itsparam = (linear_itsolver_param *) malloc(sizeof(linear_itsolver_param));
        param_linear_solver_init(itsparam);
        return itsparam;
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

/* need to add destructors for other structs!! */


%extend precond{
    void apply(dvector *lhs, dvector* rhs) {
      printf("lhs-> row %d\n", lhs->row);
      apply_precond(lhs->val, rhs->val, $self);
    }
    /*precond_ra_data* precond_data(){
      return (precond_ra_data*)$self->data;
    }*/
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
%apply (int ncol) {(int ncol)};
dCSRmat* create_matrix(double *A, int nnz, int *ja, int nnz2, int *ia, int n, int ncol);
%clear (double *A, int nnz);
%clear (int* ja, int nnz2);
%clear (int *ia, int n);
%clear (int ncol);

%apply (double* IN_ARRAY1, int DIM1) {(double* x, int n)};
dvector* create_dvector(double *x, int n);
%clear (double* x, int n);

%apply (int* IN_ARRAY1, int DIM1) {(int* x, int n)};
ivector* create_ivector(int *x, int n);
%clear (int* x, int n);

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
INT fenics_bsr_solver(INT block_size, dCSRmat *A, dvector *b, dvector *sol);
// dvector* compute_ra_aaa(REAL s_frac_power, REAL t_frac_power, REAL alpha, REAL beta, REAL scaling_a, REAL scaling_m);

%apply (int DIM1, double* IN_ARRAY1) {(int numval, double* z),
                                      (int numval2, double* f)};
//%apply (double AAA_tol) {(double AAA_tol)};
%rename (ra_aaa) my_ra_aaa;
%exception my_ra_aaa {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
dvector* my_ra_aaa(int numval, double* z, int numval2, double* f, double AAA_tol) {
    if (numval != numval2) {
        PyErr_Format(PyExc_ValueError,
                     "Arrays of lengths (%d,%d) given",
                     numval, numval2);
        return NULL;
    }
    return ra_aaa(numval, z, f, AAA_tol);
}
%}

//void print_precond_ra_amgparam(precond *pc);

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

%pythoncode %{
    def param_input_set_dict(d, inparam):
        assert isinstance(inparam, input_param), " Second argument should be of type haznics.input_param! "
        if isinstance(d, dict):
            for key in d:
                if isinstance(d[key], str):
                    exec("inparam.%s = \"%s\"" % (key, d[key]))
                else:
                    exec("inparam.%s = %s" % (key, d[key]))
        else:
            from warnings import warn
            from inspect import currentframe
            warn(" In function %s : "%(currentframe().f_code.co_name) + \
             "First argument is a not of type dict, input_param will not be set up! ", RuntimeWarning)
%}

%pythoncode %{
    def param_amg_set_dict(d, amgparam):
        assert isinstance(amgparam, AMG_param), " Second argument should be of type haznics.AMG_param! "
        if isinstance(d, dict):
            for key in d:
                if isinstance(d[key], str):
                    exec("amgparam.%s = \"%s\"" % (key, d[key]))
                else:
                    exec("amgparam.%s = %s" % (key, d[key]))
        else:
            from warnings import warn
            from inspect import currentframe
            warn(" In function %s : "%(currentframe().f_code.co_name) + \
             "First argument is a not of type dict, AMG_param will not be set up! ", RuntimeWarning)
%}

%pythoncode %{
    def param_linear_solver_set_dict(d, itsparam):
        assert isinstance(itsparam, linear_itsolver_param), " Second argument should be of type haznics.linear_itsolver_param! "
        if isinstance(d, dict):
            for key in d:
                if isinstance(d[key], str):
                    exec("itsparam.%s = \"%s\"" % (key, d[key]))
                else:
                    exec("itsparam.%s = %s" % (key, d[key]))
        else:
            from warnings import warn
            from inspect import currentframe
            warn(" In function %s : "%(currentframe().f_code.co_name) + \
             "First argument is a not of type dict, linear_itsolver_param will not be set up! ", RuntimeWarning)
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

INT wrapper_krylov_amg(dCSRmat *mat, dvector *rhs, dvector *sol, REAL tol);
INT fenics_metric_amg_solver(block_dCSRmat *A, dvector *b, dvector *x, block_dCSRmat *AD, block_dCSRmat *M, dCSRmat *interface_dof);
