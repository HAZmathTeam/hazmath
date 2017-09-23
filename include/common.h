/*this is included by all files in HAZMATH */
#ifndef COMMON_BLOCK_SIZE
#define COMMON_BLOCK_SIZE 8
/* 
   this has support for 8-dimensional simplices...
   flag_for_simplex[0] holds a flag for the current vertex (0-d
   simplex) flag_for_simplex[1] holds a flag for the current edge (1-d
   simplex) flag_for_simplex[0] holds a flag for the current vertex
   (2-d simplex) and so on.  these can be changed with in any loop
   (w.r.t. elements, vertices, edges so on and should be accessible
   from all functions).
*/
extern int flag_for_simplex[COMMON_BLOCK_SIZE];
#else
extern int flag_for_simplex[COMMON_BLOCK_SIZE];
#endif

/* 
   the following line should be in the c-file containing the main
   program #include "global.h" before including hazmath.h 
*/
