/*! \file examples/mg_geometric/rhs_bc.h
 *
 *  Copyright 2007_HAZmath_. (Hu, Adler, Zikatanov).
 *
 * \brief Header with right hand side and the boundary conditions for
 *        the MG solver for Poisson equation on a lattice grid
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
/*==========================================================================
/* nomxyz should be replaced by num_lattice from amr/unigrid.c
 *
 * NOMXYZ gives the global number of the node (i,j,k). should be
 * replaced by num_lattice from hazmatth/src/amr/unigrid.c
==========================================================================*/
static INT nomxyz(INT i, INT j, INT k, INT nx, INT ny, INT nz)
{
/*
...  NOMXYZ gives the global number of the node (i,j,k)
*/
INT nxyz;

  if(i < 1 || j < 1 || k < 1 || i > nx || j > ny || k > nz)
    nxyz = 0;
  else
    nxyz = (k-1)*nx*ny + (j-1)*nx + i;

  return nxyz;
}
/*==========================================================================
/* nomxy should be replaced by num_lattice from amr/unigrid.c
 *
 * NOMXY gives the global number of the node (i,j). should be
 * replaced by num_lattice from hazmatth/src/amr/unigrid.c
==========================================================================*/
static INT nomxy(INT i, INT j, INT nx, INT ny)
{
INT nxy;

  if(i < 1 || j < 1 || i > nx || j > ny)
    nxy = 0;
  else
    nxy = (j-1)*nx + i;

  return nxy;
}
/*==========================================================================*/
/* 3D CASE                                                                  */
/*==========================================================================*/
REAL fxyz(REAL x, REAL y, REAL z)
{
  return 1.;
}
/*==========================================================================*/
void rhs3d(REAL *f, INT nx, INT ny, INT nz)
{
/* Forms right hand side. External function used: fxyz() */
INT i,j,k,ii,j0,k0,nxy,nxyz;
REAL h,hh,x,y,z;

  h=1./((double)(nx-1));
  nxy = nx*ny;
  z=0.;

  for (k=1; k<=nz; k++)
  {
    y=0.;
    k0 = (k-1)*nxy;
    for (j=1; j<=ny; j++)
    { 
      x=0.;
      j0 = (j-1)*nx;
      for (i=0; i < nx; i++)					/*!*/
      { 
        ii = i + j0 + k0;
        f[ii] = fxyz(x,y,z);
        x = x + h;
      }
      y = y + h;
    }
    z = z + h;
  }

  hh = h*h;
  nxyz=nxy*nz;
  for (ii=0; ii < nxyz; ii++)					/*!*/
    f[ii]=f[ii]*hh;
} 
/*==========================================================================*/
void bnd3d(REAL *u, REAL *f, INT nx, INT ny, INT nz)
{
/* Sets boundart conditions. External function used: nomxyz() */
INT i,j,k,ii;

  for (j=1; j<=ny; j++)
    for (k=1; k<=nz; k++)
    {
      ii = nomxyz(1,j,k,nx,ny,nz)-1;				/*!*/
      u[ii] = 0.;
      f[ii] = 0.;
      ii = nomxyz(nx,j,k,nx,ny,nz)-1;				/*!*/
      u[ii] = 0.;
      f[ii] = 0.;
    }

  for (i=1; i<=nx; i++)
    for (k=1; k<=nz; k++)
    {
      ii = nomxyz(i,1,k,nx,ny,nz)-1;				/*!*/
      u[ii] = 0.;
      f[ii] = 0.;
      ii = nomxyz(i,ny,k,nx,ny,nz)-1; 				/*!*/
      u[ii] = 0.;
      f[ii] = 0.;
    }

 for (i=1; i<=nx; i++)
    for (j=1; j<=ny; j++)
    {
      ii = nomxyz(i,j,1,nx,ny,nz)-1;				/*!*/
      u[ii] = 0.;
      f[ii] = 0.;
      ii = nomxyz(i,j,nz,nx,ny,nz)-1;				/*!*/
      u[ii] = 0.;
      f[ii] = 0.;
    }
}
/*==========================================================================*/
/* 2D CASE                                                                  */
/*==========================================================================*/
REAL fxy(REAL x, REAL y)
{
  return 1.;
}
/*==========================================================================*/
void rhs2d(REAL *f, INT nx, INT ny)
{
/* Forms right hand side. External function used: fxy() */
INT i,j,ii,j0,nxy;
REAL h,hh,x,y;

  h=1./((double)(nx-1));
  y=0.;

  for (j=1; j<=ny; j++)
  { 
    x=0.;
    j0 = (j-1)*nx;
    for (i=0; i < nx; i++)					/*!*/
    { 
      ii = i + j0;
      f[ii] = fxy(x,y);
      x = x + h;
    }
    y = y + h;
  }

  hh = h*h;
  nxy=nx*ny;
  for (ii=0; ii < nxy; ii++)					/*!*/
    f[ii]=f[ii]*hh;
} 
/*==========================================================================*/
void bnd2d(REAL *u, REAL *f, INT nx, INT ny)
{
/* Sets boundart conditions. External function used: nomxyz() */
INT i,j,ii;

  for (j=1; j<=ny; j++)
  {
    ii = nomxy(1,j,nx,ny)-1;					/*!*/
    u[ii] = 0.;
    f[ii] = 0.;
    ii = nomxy(nx,j,nx,ny)-1;					/*!*/
    u[ii] = 0.;
    f[ii] = 0.;
  }

  for (i=1; i<=nx; i++)
  {
    ii = nomxy(i,1,nx,ny)-1;					/*!*/
    u[ii] = 0.;
    f[ii] = 0.;
    ii = nomxy(i,ny,nx,ny)-1;					/*!*/
    u[ii] = 0.;
    f[ii] = 0.;
  }
}
/*==========================================================================*/
