/*! \file examples/mg_geometric/gmg.c
 *
 *  Copyright 2007_HAZmath_. (Hu, Adler, Zikatanov).
 *
 * \brief This program solves the Poisson equation
 *
 *        -Delta(u) = f,
 *
 *        on the unit square in 2D or the unit cube in 3D  using V-cycle
 * 	  multigrid method with Dirichlet boundary conditions
 *
 * \note Created by Ludmil Zikatanov 20070109 (yyyymmdd)
 * \note With contributions from Johannes K. Kraus and Yunrong Zhu.
 *
 */
/***********************************************************************/
/*HHAZmath includes*/
/***********************************************************************/
//#include "hazmath.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
// look below for int and long int. 
#define __longint 1
#define __longdouble 0
///////////////////////////////////////////////
#ifdef INT
#undef INT
#endif
#ifdef REAL
#undef REAL
#endif
///////////////////////////////////////////////
#if( __longint == 1 )
#define INT long int
#else
#define INT int
#endif


#if( __longdouble == 1 )
#define REAL long double
#else 
#define REAL double
#endif
/* all routines */
#include "gmg_routines.h"
/* right hand side and boundary conditions */
#include "rhs_bc.h"

/* for the MG */
#define MAX_LEVELS 20
#define MIN_LEVELS 2
#define MIN_LEVELS 2
#define MAX_ITER 500
#define N_SMOOTH 1
/****************************************************/
void main()
{
  INT ke[MAX_LEVELS];
  REAL *u=NULL;
  REAL *f=NULL;
  /* WORKING */
  REAL *uk=NULL;
  REAL *rk=NULL;
  REAL *tk=NULL;
  REAL tol;
  REAL h;
  INT nspdim=3,n,levels,nx,ny,nz,ntotal;
  INT malls;
  INT mxitr, nsweeps, nvcycles;
  INT k;
  INT iterv,isfmg;
  REAL reduce, errrel;

  /* timing? */
  clock_t beg ; double time0;
  levels = 9; 
  fprintf(stdout, "\nINPUT NOTE: Number of points in one direction is = 2^(levels)+1.\n");
  fprintf(stdout, "           *Enter the desired number of levels (2-%2i): ",MAX_LEVELS);
  //
  FILE *fp;
  //  fp=stdin;
  fp=fopen("sample.input","r");
  //
#if( __longint == 1 )
  k=fscanf(fp,"%ld", &levels);
#else
  k=fscanf(fp,"%d", &levels);
#endif
  fprintf(stdout,"\n");
  if(levels > MAX_LEVELS) 
    levels=MAX_LEVELS;
  else if(levels < MIN_LEVELS) 
    levels=MIN_LEVELS;
  /* Number of pre- and post- smoothing steps */
  fprintf(stdout, "           *Enter spatial dimension (2 or 3): ");
#if( __longint == 1 )
  k=fscanf(fp,"%ld", &nspdim);
#else
  k=fscanf(fp,"%d", &nspdim);
#endif
  fprintf(stdout,"\n");
  if(nspdim > 3) 
    nspdim=3;
  else if(nspdim < 2) 
    nspdim=2;
  /* */
  fprintf(stdout, "\n           *Enter number of smoothing steps: ");
#if( __longint == 1 )
  k=fscanf(fp,"%ld", &nsweeps);
#else
  k=fscanf(fp,"%d", &nsweeps);
#endif
  fprintf(stdout,"\n");
  if(nsweeps < 1) 
    nsweeps=N_SMOOTH;
  /* */
  fprintf(stdout, "\n           *(V-cycles or FMG): [0 or 1]: ");
#if( __longint == 1 )
  k=fscanf(fp,"%ld", &isfmg);
#else
  k=fscanf(fp,"%d", &isfmg);
#endif
  fclose(fp);
  fprintf(stdout,"\n");
  if(isfmg) { 
    isfmg=1;
    nvcycles = 1;
  } else {
    isfmg=0;
    nvcycles=0;
  }
  /* */
  n=(INT )pow(2,levels)+1;
  nx = n;
  ny = n;
  if(nspdim == 2)   
    nz = 1;
  else {
    nspdim=3;
    nz = n;
  }
  tol = (REAL )pow(2.,-20);
  /* */

  memgmg(isfmg,levels,nspdim,nsweeps,&malls,ke,nx,ny,nz,&ntotal,&h,tol);
  u = (REAL *)malloc(ntotal*sizeof(REAL));
  f = (REAL *)malloc(ntotal*sizeof(REAL));
  tk = (REAL *)malloc(ntotal*sizeof(REAL));
  uk = (REAL *)malloc(malls*sizeof(REAL));
  rk = (REAL *)malloc(malls*sizeof(REAL));
  /*  u = new REAL [ntotal];
      f = new REAL [ntotal];
      tk = new REAL  [ntotal];
      uk = new REAL [malls];
      rk = new REAL [malls];
  */
  if(u && f && rk && tk && uk) {
    /* MAX num V-cycle iterations */
    mxitr=MAX_ITER; 
    /* 
       Set up the problem, that is get rsh and boundary conditions
    */
    for (k=0; k < ntotal; k++)
      {
	f[k]=0e+00;
	u[k]=0e+00;
      }
    if(nspdim == 3) {
      rhs3d(f,nx,ny,nz);
      bnd3d(u,f,nx,ny,nz);
      //
      beg=clock();
      gmg3D(isfmg,levels,nx,ny,nz,ntotal,f,u,tk,uk,rk,ke,tol,	\
	       nsweeps,nvcycles,mxitr,&iterv,&errrel,&reduce, h);
      time0 = (double)(clock() - beg) / CLOCKS_PER_SEC;
    } else {
      /*       call rhs2d(f,nx,ny)*/
      /* Set initial guess */
      rhs2d(f,nx,ny);
      bnd2d(u,f,nx,ny);
      beg=clock();      
      gmg2D(isfmg,levels,nx,ny,ntotal,f,u,tk,uk,rk,ke,tol,		\
	       nsweeps,nvcycles,mxitr,&iterv,&errrel,&reduce, h);
      time0 = (double)(clock() - beg) / CLOCKS_PER_SEC;
    }
    if(isfmg) {
      printf("\n FMG (%ld V-cycles); ||r||/||f||=%10.3e\n",(long int)nvcycles,errrel);
    } else {
      printf("\n %ld V-cycles; rel.residual=%10.3e ; damp.factor=%10.3e\n",\
	     (long int)iterv,errrel,reduce);
    }
    printf("\n ***CPU (calculated by clock())=%10.2f seconds\n\n",time0);    
    if(u) free(u);
    if(f) free(f);
    if(tk) free(tk);
    if(uk) free(uk);
    if(rk) free(rk);
  } else {
    fprintf(stdout, "\n****ERROR: (levels=%ld) is too big.\n****Could not allocate memory for",(long int)levels);
    if(!u) 
      fprintf(stdout, " the solution U()");
    else if(!f)
      fprintf(stdout, " the right hand side F()");
    else if(!rk)
      fprintf(stdout, " a working vector RK()");
    else
      fprintf(stdout, " a working vector TK()");
    fprintf(stdout, " ****\n\nExiting with code 255...\n\n");
    exit(255);
  }
  return;
} 
