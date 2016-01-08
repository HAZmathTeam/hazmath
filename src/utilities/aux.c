#include "hazmat.h"

/*** Auxillary Files (some from Ludmil) *******************************************************/

/****************************************************************************************/
void constcoeff(REAL *val,REAL* x,REAL constval) 
{
    *val = constval;
}
/****************************************************************************************/


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

/****************************************************************************************/
void det3D(REAL *mydet,REAL* vec1,REAL* vec2,REAL* vec3)          
{
  /* gets determinant of 3 3D vectors */
  REAL dettmp;
  
  dettmp = vec1[0]*(vec2[1]*vec3[2]-vec2[2]*vec3[1]) - vec1[1]*(vec2[0]*vec3[2]-vec2[2]*vec3[0]) + vec1[2]*(vec2[0]*vec3[1]-vec2[1]*vec3[0]);

  *mydet = dettmp;

  return;
}
/****************************************************************************************/

/****************************************************************************************/
void baddimension()          
{
  // Print Error if dimension is not 2 or 3
  printf("\n!!!  You have now entered the Twilight Zone.  Your dimesnion is not 2 or 3!  !!!\n\n");
  exit(2);
  return;
}
/****************************************************************************************/

/**************** Read input file that sets parameters ***********************************/
void getinput(char* gridfile,REAL* fpar,INT* ipar)
{
  /* This routine opens file input.dat and grabs all parameters.  It puts INTeger parameters 
   * in ipar and REAL parameters in fpar
   */
	
  FILE* inputid;
  char infile[50] = "input.dat";
  char dummy[900];
	
				 
  // Open file			 
  inputid = fopen(infile,"r");
  if (inputid == NULL) {
    fprintf(stderr, "Can't open input file input.dat!\n");
    exit(1);
  }
  fscanf(inputid,"%s",dummy);                                // Header
  fscanf(inputid,"%d %s",&ipar[0],dummy);	             // Dimension 
  fscanf(inputid,"%d %s",&ipar[1],dummy);	             // nq1d quad nodes in each direction
  fscanf(inputid,"%d %d %s",&ipar[2],&ipar[3],dummy);        // order of velocity Stokes elements and order of pressure elements
    
  fscanf(inputid,"%lg %s",&fpar[0],dummy);                   // coefficient in front of B
  fscanf(inputid,"%lg %s",&fpar[1],dummy);                   // coefficient in front of curl curl B
  fscanf(inputid,"%lg %s",&fpar[2],dummy);                   // coefficient in front of u
  fscanf(inputid,"%lg %s",&fpar[3],dummy);                   // coefficient in front of Laplace u
    
  fscanf(inputid,"%s %s",gridfile,dummy);                    // Name of grid file
    
  fscanf(inputid,"%d %d %d %d %d %d %d %d %s",&ipar[4],\
	 &ipar[5],&ipar[6],&ipar[7],&ipar[8],&ipar[9],\
	 &ipar[10],&ipar[11],dummy);                         // Indicates whether true solution is known and given function (u1,u2,u3,p,B1,B2,B3)
  fscanf(inputid,"%d %d %d %d %d %d %d %d %s",&ipar[12],\
	 &ipar[13],&ipar[14],&ipar[15],&ipar[16],&ipar[17],\
	 &ipar[18],&ipar[19],dummy);                          // Indicates whether boundary conditions are timedependent and functions for u1 u2 u3 p B1 B2 B3
  fscanf(inputid,"%d %d %d %d %d %d %s",&ipar[20],&ipar[21],\
	 &ipar[22],&ipar[23],&ipar[24],&ipar[25],dummy);    // Indicates rhs functions for equations of u1 u2 u3 B1 B2 B3
  fscanf(inputid,"%d %d %d %d %d %d %d %s",&ipar[26],\
	 &ipar[27],&ipar[28],&ipar[29],&ipar[30],&ipar[31],\
	 &ipar[32],dummy);                                  // Indicates initial condition (u1,u2,u3,p,B1,B2,B3)
    
  fscanf(inputid,"%d %d %lg %s",&ipar[33],&ipar[34],\
	 &fpar[4],dummy);                                   // Timestepping: Type (0-CN,1-BDF1,2-BDF2) Time Steps, Time step size
    
  fscanf(inputid,"%d %lg %s",&ipar[35],&fpar[5],dummy);     // Newton:  Max Steps, Tolerance
    
  fscanf(inputid,"%d %d %d %d %s",&ipar[36],&ipar[37],\
	 &ipar[38],&ipar[39],dummy);                        // Dumping parameters: Sol(0 none 1 screen 2 file 3 both) Mesh (1 or 0) Matrix Iteration Output (0 or 1)
    
  fscanf(inputid,"%d %d %lg %d %s",&ipar[40],&ipar[41],\
	 &fpar[6],&ipar[42],dummy);                         // Solver: Type(0 CG 1 MINRES 2 GMRES 3 Precond MINRES 4 Precond GMRES) Max Iter. Tol Restart (-1=full)
    
  fscanf(inputid,"%lg %s",&fpar[7],dummy);                   // Maxwell disappearing solution parameter
    
  fscanf(inputid,"%d %d %d %d %s",&ipar[43],&ipar[44],&ipar[45],&ipar[46],dummy);     
// Multilevel solver stuff fine level, coarse level (level 0->h=1 1->h=0.5 etc), presmoothing steps, size of subdomains (# subdomains = (2^(param)-1)^2)
  
  fclose(inputid);

  return;
}
/****************************************************************************************/

