#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

/* ===================================================== */
#include "macro.h"
#include "sparse.h"
#include "vec.h"
/*=============================================================*/

/*** Auxillary Files (from Ludmil) *******************************************************/

/****************************************************************************************/
void rveci_(FILE *fp, INT *vec, INT *nn)       
/* reads a vector of integers of size nn from a file fp*/
{
	
  INT n;
  INT *vec_end;
  n = *nn;
  vec_end  =  vec + n;
  for ( ; vec < vec_end; ++vec)
    fscanf(fp,"%i",vec);
  //fprintf(stdout,"Read %d INTEGERS", n);
  return;
}
/****************************************************************************************/

/****************************************************************************************/
void rvecd_(FILE *fp,  REAL *vec, INT *nn)
/* reads a vector of REALS of size nn from a file fp*/
{
  INT n;
  REAL *vec_end;  
  n= *nn;
  vec_end =  vec + n;
  for ( ; vec < vec_end; ++vec)
    fscanf(fp,"%lg",vec);
  //fprintf(stdout,"Read %d REALS", n);
  return;
}
/****************************************************************************************/
