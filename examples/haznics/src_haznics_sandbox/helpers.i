%module helpers
%{
#include "helpers.h"
#include "pyhelpers.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION 
#include "ndarrayobject.h"



%}

%init%{
import_array();
%}



%typemap(in) double* {
PyArrayObject* array = (PyArrayObject*)$input; 
const long int i=0; 
/* this should be improved, no checking at all here, 
   but numpy has a rich api, not sure if it is supposed to be i=1 up there */ 
$1 =  (double*) PyArray_GetPtr(array, &i);
}



%typemap(in) void (*func_ptr)(double*, double*) {
  $1 = $input; 
}

%typemap(in) void *data {
  $1 = $input; 
}; 


%include "helpers.h"
%include "pyhelpers.h"


