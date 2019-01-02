#ifndef INT
#define INT int
#endif

#ifndef REAL
#define REAL double
#endif

#ifndef FILENAMESIZE
#define FILENAMESIZE  1024
#endif

/* MAX LEVELS of refinement */
#define MAXREFLEVELS 5

/* 
   preceeds all file names, so it can be a directory name and a
   prefix, etc
*/

#define OPREFIX "./mesh_3d" 
#define IPREFIX "./macroel" 

#define VTKDO  1

/* spatial dimension */
#define DIM 3

#define USE_FEATURES 0
#define FEATURES_DIR  "./"

//#define FEATURES_FILE_IN "spiral2d.inp"
//#define FEATURES_FILE_IN "spiral.inp"
//#define FEATURES_FILE_IN "circle.inp"
#define FEATURES_FILE_IN "sphere.inp"

#define FEATURES_FILE_OUT "opts123.out"


/* 
   boundary conditions: boundary codes for essential BC and Robin type BC.
   the last two are ignored in 2D. If they are between 
1--16 Dirichlet; 17--32 Neumann; 33->64 Robin.
*/
// example with analytical solution, if set to 1;
#ifndef LEFTBC
#define LEFTBC 1
#endif
#ifndef RIGHTBC
#define RIGHTBC 2
#endif
#ifndef FRONTBC
#define FRONTBC 30
#endif
#ifndef BACKBC
#define BACKBC 31
#endif
#ifndef BOTTOMBC
#define BOTTOMBC 19
#endif
#ifndef TOPBC
#define TOPBC 22
#endif
/*********************************************************************/

