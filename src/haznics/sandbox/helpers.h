#ifndef __helpers_H
#define __helpers_H



void dabla(double *a, double* b);

void dabla2(void* func_ptr0, double* a, double* b); 

typedef void (*func_ptr)(double*, double*);
func_ptr make_dabla();

typedef struct {
  double *a; 
  double *b; 
  func_ptr func;
} dabla_struct; 


typedef void (*func_ptr1)(double*, double*, void*);

typedef struct {
  int type; 
  void *data; 
  double *a; 
  double *b; 
  func_ptr1 func1;
} dabla_struct1; 



void dabla4(dabla_struct *d); 


#endif 


