#include "stdio.h"
#include "helpers.h"

void dabla(double *a, double* b){
   /* I know a, b are of size 4, but try to run further */ 
   for (int i=0; i< 7; i++) {
	   printf("i=%d a=%e b=%e\n", i, a[i], b[i]); 
   }
}

void dabla2(void* func_ptr, double* a, double* b){ 

	void (*func_ptr2)(double*, double*)= func_ptr;  

	(*func_ptr2)(a, b); 
}

func_ptr make_dabla() {
  return dabla;
}

void dabla4(dabla_struct *d){
  (*d->func)(d->a, d->b); 
}


