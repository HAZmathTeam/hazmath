/*

  This is included ONLY ONCE by the C-file containing the main program
  in a particular application. See examples/ The array defined below
  is global and accessible by all hazmath functions defined in files
  that include hazmath.h

THIS FILE MUST NOT BE INCLUDED BY hazmath.h

*/
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
int flag_for_simplex[COMMON_BLOCK_SIZE]={0,0,0,0,0,0,0};

