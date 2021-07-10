#include "helpers.h"
#include "stdio.h"

int main() {
	double a[4]; 
	double b[4]; 
	for (int i=0; i<4; i++) {
		a[i] = i; 
		b[i] = i*i; 
	}
	dabla(a, b); 

        printf("\n\ndealing with func ptr\n\n"); 
	void (*func_ptr0)(double*, double*)=&dabla;  
	(*func_ptr0)(a, b); 


        printf("\n\nsending in func ptr to function\n\n"); 
	dabla2(func_ptr0, a, b); 


        printf("\n\n more swig friendly version \n\n"); 
	void (*func_ptr1)(double*, double*)=make_dabla();  
	dabla2(func_ptr1, a, b); 

}


