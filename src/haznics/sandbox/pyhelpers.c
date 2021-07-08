#include "pyhelpers.h"

PyObject* callback(double* a, double *b, PyObject* pyfunc){


  printf("here I am inside the callback \n"); 

  import_array();


  printf("in C: a[2]=%e, b[2]=%e \n", a[2], b[2]); 
  npy_intp D[1]; D[0]=4; 
  // create new Python arrays 
  PyObject* aa = PyArray_SimpleNewFromData(1, D, NPY_DOUBLE, (void*)a);
  PyObject* bb = PyArray_SimpleNewFromData(1, D, NPY_DOUBLE, (void*)b);

  // memory management
  Py_INCREF(aa);  
  Py_INCREF(bb);

  // Build up the argument list...
  PyObject *arglist = Py_BuildValue("(OO)", aa, bb);

  // ...for calling the Python compare function
  PyObject *result = PyEval_CallObject(pyfunc, arglist);

  PyObject *return_obj = PyUnicode_FromString("everything is gone be ok");
  return return_obj; 

}


void dabla5(dabla_struct1 *d){
  if (d->type == 1) { // assume that 1 means Python object 
    PyObject* pyobj = (PyObject*) d->data;   
    callback(d->a, d->b, pyobj); 
  } else {
    (*d->func1)(d->a, d->b, d->data); 
  }
}
