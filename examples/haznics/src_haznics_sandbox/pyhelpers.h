#ifndef __pyhelpers_H
#define __pyhelpers_H


#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION 
#include "ndarrayobject.h"
#include "helpers.h"


PyObject* callback(double* a, double *b, PyObject* pyfunc);
void dabla5(dabla_struct1 *d); 

#endif 
