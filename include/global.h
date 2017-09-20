/*

  This is included ONLY ONCE by the C-file containing the main program
  in a particular application. See examples/ The array defined below
  is global and accessible by all hazmath functions defined in files
  that include hazmath.h

THIS FILE MUST NOT BE INCLUDED BY hazmath.h

*/
#define COMMON_BLOCK_SIZE 16
/* 
   flag_for_simplex[] holds global variables needed in evaluating
   functions. for the current vertex (0-d

   flag_for_simplex[0] is the boundary code of the dof or the face
   under consideration. 

   flag_for_simplex[2] >= 0 indicates that this is a dof code (for
   spaces with more than component, e.g. mixed method and in such case
   it denotes the FE space number (0 or 1 or 2 etc, i.e. the block
   number).

   flag_for_simplex[1] > 0 indicates that this is a face.

*/
int flag_for_simplex[COMMON_BLOCK_SIZE]={-1,-1,-1,-1,	\
					 -1,-1,-1,-1,	\
					 -1,-1,-1,-1,	\
					 -1,-1,-1,-1};

