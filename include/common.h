/*this is included by all files in HAZMATH */
#ifndef COMMON_BLOCK_SIZE
#define COMMON_BLOCK_SIZE 16
/* 
   flag_for_simplex[] holds global variables needed in evaluating
   functions. for the current vertex (0-d

   flag_for_simplex[0] is the boundary code of the dof or the face
   under consideration. 

   flag_for_simplex[1] >= 0 indicates that this is a dof code (for
   spaces with more than component, e.g. mixed method and in such case
   it denotes the FE space number (0 or 1 or 2 etc, i.e. the block
   number).  

   flag_for_simplex[2] >= 0 indicates that this is a a face and the
   numer is the dimension of the face.

*/
extern int flag_for_simplex[COMMON_BLOCK_SIZE];
#else
extern int flag_for_simplex[COMMON_BLOCK_SIZE];
#endif

/* 
   the following line should be in the c-file containing the main
   program 
   
   #include "global.h" 

   before including "hazmath.h" 
*/
