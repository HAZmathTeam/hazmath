/*! \file examples/mg_geometric/gmg_routines.h
 *
 *  Copyright 2007_HAZmath_. (Hu, Adler, Zikatanov).
 *
 * \brief Header with core routines for the MG solver for Poisson
 *        Poisson equation on lattice grid
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
/* ddot0                                                                    */
/*==========================================================================*/
static REAL ddot0(REAL *u, REAL *v, INT n)
{
REAL ddotr=0.;
INT k;
  for (k=0; k < n; k++)						/*!*/
    ddotr += u[k]*v[k];

  return ddotr;
}
/*==========================================================================*/
/* res02d                                                                   */
/*==========================================================================*/
void res02d(REAL *f, REAL *u, REAL *r0, INT nx, INT ny)
{
  /*
    ... COMPUTES THE RESIDUAL
    ...    r0 = f-A*u
  */
  REAL zero=0.;
  REAL four=4.;

  INT i,i0,i1,i2,i3,i4,j,j0,jm1,jp1,k,nsize;

  nsize = nx*ny;

  for (k=0; k < nsize; k++)					/*!*/
    r0[k] = zero;

  for (j=2; j <= ny-1; j++)
    {
      jp1 = j*nx;
      j0 = jp1-nx;
      jm1 = j0-nx;
      for (i=1; i < nx-1; i++)					/*!*/
	{
	  i0 = i   +  j0;
	  i1 = i+1 +  j0;
	  i2 = i   +  jp1;
	  i3 = i-1 +  j0;
	  i4 = i   +  jm1;
	  r0[i0] = f[i0] - u[i0]*four + u[i1] + u[i2] + u[i3] + u[i4];
	}
    }
}
/*==========================================================================*/
/* res0                                                                     */
/*==========================================================================*/
void res0(REAL *f, REAL *u, REAL *r0, INT nx, INT ny, INT nz)
{
/*
... COMPUTES THE RESIDUAL
...    r0 = f-A*u
*/
REAL zero=0.;
REAL six =6e+00;

INT i,i0,i1,i2,i3,i4,i5,i6,j,j0,jm1,jp1,k,k0,km1,kp1,nxy,nsize;

  nxy = nx*ny;
  nsize = nxy*nz;

  for (k=0; k < nsize; k++)					/*!*/
    r0[k] = zero;

  for (k=2; k <= nz-1; k++)
  {
    kp1 = k*nxy;
    k0 = kp1-nxy;
    km1 = k0-nxy;
    for (j=2; j <= ny-1; j++)
    {
      jp1 = j*nx;
      j0 = jp1-nx;
      jm1 = j0-nx;
      for (i=1; i < nx-1; i++)					/*!*/
      {
        i0 = i   +  j0    + k0;
        i1 = i+1 +  j0    + k0;
        i2 = i   +  jp1   + k0;
        i3 = i   +  j0    + kp1;
        i4 = i-1 +  j0    + k0;
        i5 = i   +  jm1   + k0;
        i6 = i   +  j0    + km1;
        r0[i0] = f[i0] - u[i0]*six + u[i1] + u[i2] + u[i3] + u[i4] + u[i5] + u[i6];
      }
    }
  }
}
/*==========================================================================*/
/* Restriction-Prolongation*/
/*==========================================================================*/
/* 
...RESTRICTION & PROLONGATION
*/
/*==========================================================================*/
/* p3d                                                                      */
/*==========================================================================*/
void p3d(REAL *u, REAL *pu, INT nfx, INT nfy, INT nfz, INT ncx, INT ncy, INT ncz)
{
/*
...  This subroutine performs P1 interpolation on the unit cube in 3D
...  The numbering convention convention is that we have
...  lexicographical numbering, first goes x, then y then z Thus the
...  global number of a point (i,j,k) will be iglobal = (k-1)*nxy +
...  (j-1)*nx + i. That also agrees with the fortran convention, so a
...  parameters can be U(NCX,NCY,NCZ) and PU(NFX,NFY,NFZ)
...
...  u(*) is input vector on the coarse grid of size (ncx*ncy*ncz)
...  pu(*) is output vector, p.w. linear interpolation is used and
...  pu(*) = pu(*) + P*u.  The reason being that in MG we have to
...  perform UNEW = UOLD + P*UCOARSE, and this prolongation works
...  without introducing an extra vector UNEW, it directly gives 
...
...  PU <-- PU+interpolation(U).
...  
...  The vector pu(*) is of size nfx*nfy*nfz
*/
REAL two=2.;
REAL half=0.5;

INT nfxy,ncxy,nctot,nbig,nftot;
INT i,i0,ic0,ic1,k,k0,k1,k00,k10,kc0,kc1,j,j0,j1,j00,j10,jc0,jc1;
INT lf1,lf2,lf3,lf4,lf5,lf6,lf7,lf8,lc1,lc2,lc3,lc4,lc5,lc6,lc7,lc8;

  nfxy = nfx*nfy;
  ncxy = ncx*ncy;
  nctot = ncxy*ncz;
  nbig = -2*nctot;
  nftot = nfx*nfy*nfz;

  for (i = 0; i < nftot; i++)					/*!*/
    pu[i] = two*pu[i];

  for (k = 1; k <= ncz; k++)
  {
    k0 = 2*k - 1;
    k1 = 2*k;
    k00 = (k0-1)*nfxy;
    k10 = k00 + nfxy;
    kc0 = (k-1)*ncxy;

    if(k == ncz)
      kc1 = nbig;
    else
      kc1 = kc0 + ncxy;

    for (j = 1; j <= ncy; j++)
    {
      j0 = 2*j - 1;
      j1 = 2*j;
      j00 = (j0-1)*nfx;
      j10 = j00 + nfx;
      jc0 = (j-1)*ncx;
      if(j == ncy)
        jc1 = nbig; 
      else
        jc1 = jc0 + ncx;

      for (i = 1; i <= ncx; i++)
      {
        i0 = 2*i - 2;						/*!*/
        ic0 = i-1;						/*!*/
        if(i == ncx)
          ic1 = nbig-1;						/*!*/
        else
          ic1 = ic0 + 1;
/*
... Calculate indecies of the neighbors
*/
        lf1 = i0 + j00 + k00;
        lf2 = lf1 + 1;
        lf3 = i0 + j10 + k00;
        lf4 = lf3 + 1;
        lf5 = i0 + j00 + k10;
        lf6 = lf5 + 1;
        lf7 = i0 + j10 + k10;
        lf8 = lf7 + 1;
/*
... These should give something < 0 if needed
*/
        lc1 = ic0 + jc0 + kc0;
        lc2 = ic1 + jc0 + kc0;
        lc3 = ic0 + jc1 + kc0;
        lc4 = ic1 + jc1 + kc0;
        lc5 = ic0 + jc0 + kc1;
        lc6 = ic1 + jc0 + kc1;
        lc7 = ic0 + jc1 + kc1;
        lc8 = ic1 + jc1 + kc1;
/*
... Interpolate then:
*/
        pu[lf1]  = pu[lf1] + u[lc1] + u[lc1];
        if(lc2 > 0 && lc2 <= nctot) 
          pu[lf2]  =  pu[lf2] + (u[lc1] + u[lc2]);
        if(lc3 > 0 && lc3 <= nctot) 
          pu[lf3]  =  pu[lf3] + (u[lc1] + u[lc3]);
        if(lc4 > 0 && lc4 <= nctot) 
          pu[lf4]  = pu[lf4] + (u[lc1] + u[lc4]);
        if(lc5 > 0 && lc5 <= nctot)
          pu[lf5]  = pu[lf5] + (u[lc1] + u[lc5]);
        if(lc6 > 0 && lc6 <= nctot)
          pu[lf6]  = pu[lf6] + (u[lc1] + u[lc6]);
        if(lc7 > 0 && lc7 <= nctot)
          pu[lf7]  = pu[lf7] + (u[lc1] + u[lc7]);
        if(lc8 > 0 && lc8 <= nctot) 
          pu[lf8]  = pu[lf8] + (u[lc1] + u[lc8]);
      }
    }
  }
  for (i = 0; i < nftot; i++)					/*!*/
    pu[i] = half*pu[i];
}
/*==========================================================================*/
/* pt3d                                                                     */
/*==========================================================================*/
void pt3d(REAL *u, REAL *ru, INT nfx, INT nfy, INT nfz, INT ncx, INT ncy, INT ncz)
{
/*
.... This computes 4 times P^t*x and then scales it appropriately
     implicit real*8(a-h,o-z), integer(i-n)
     parameter (zero=0d0,qrtr=0.25d0)
     dimension  ru(*),u(*)
*/
REAL zero=0.;
REAL qrtr=0.25;

INT ijk,nfxy,ncxy,nctot,nbig;
INT i,i0,ic0,ic1,k,k0,k1,k00,k10,kc0,kc1,j,j0,j1,j00,j10,jc0,jc1;
INT lf1,lf2,lf3,lf4,lf5,lf6,lf7,lf8,lc1,lc2,lc3,lc4,lc5,lc6,lc7,lc8;

  nfxy = nfx*nfy;
  ncxy = ncx*ncy;
  nctot = ncxy*ncz;
  nbig = -2*nctot;

  for (k = 0; k < nctot; k++)					/*!*/
    ru[k] = zero;

  for (k = 1; k <= ncz-1; k++)
  {
    k0 = 2*k - 1; 
    k1 = 2*k;
    k00 = (k0-1)*nfxy;
    k10 = k00 + nfxy;
    kc0 = (k-1)*ncxy;
    kc1 = kc0 + ncxy;

    for (j = 1; j <= ncy-1; j++)
    {
      j0 = 2*j - 1;
      j1 = 2*j; 
      j00 = (j0-1)*nfx;
      j10 = j00 + nfx;
      jc0 = (j-1)*ncx;
      jc1 = jc0 + ncx;

      for (i = 1; i <= ncx-1; i++)
      {
        i0 = 2*i - 2;						/*!*/
        ic0 = i-1;						/*!*/
        ic1 = ic0 + 1;
/*
... Calculate indecies of the neighbors
*/
        lf1 = i0 + j00 + k00;
        lf2 = lf1 + 1;
        lf3 = i0 + j10 + k00;
        lf4 = lf3 + 1;
        lf5 = i0 + j00 + k10;
        lf6 = lf5 + 1;
        lf7 = i0 + j10 + k10;
        lf8 = lf7 + 1;

        lc1 = ic0 + jc0 + kc0;
        lc2 = ic1 + jc0 + kc0;
        lc3 = ic0 + jc1 + kc0;
        lc4 = ic1 + jc1 + kc0;
        lc5 = ic0 + jc0 + kc1;
        lc6 = ic1 + jc0 + kc1;
        lc7 = ic0 + jc1 + kc1;
        lc8 = ic1 + jc1 + kc1;

        if(i > 1 && j > 1 && k > 1)
          ru[lc1] = ru[lc1] + u[lf1]+u[lf1]+(u[lf2]+u[lf3]+u[lf4]+u[lf5]+u[lf6]+u[lf7]+u[lf8]);
        if(j > 1 && k > 1)
          ru[lc2]  =  ru[lc2] + u[lf2];
        if(i > 1 && k > 1)
          ru[lc3]  =  ru[lc3] + u[lf3];
        if(k > 1)
          ru[lc4]  =  ru[lc4] + u[lf4];
        if(i > 1 && j > 1)
           ru[lc5]  =  ru[lc5] + u[lf5];
        if(j > 1)
           ru[lc6]  =  ru[lc6] + u[lf6];
        if(i > 1)
           ru[lc7]  =  ru[lc7] + u[lf7];
        ru[lc8]  =  ru[lc8] + u[lf8];
      }
    }
  }

  for (ijk = 0; ijk < nctot; ijk++)				/*!*/
     ru[ijk] = qrtr*ru[ijk];

}
/*==========================================================================*/
/* 2D case                                                                  */
/*==========================================================================*/
/* p2d                                                                      */
/*==========================================================================*/
void p2d(REAL *u, REAL *pu, INT nfx, INT nfy, INT ncx, INT ncy)
{
/*
...  This subroutine performs P1 interpolation on the unit square in 2D
...  The numbering convention convention is that we have
...  lexicographical numbering, first goes x, then y then z Thus the
...  global number of a point (i,j) will be iglobal = (j-1)*nx +
...  i. That also agrees with the fortran convention, so a parameters
...  can be U(NCX,NCY) and PU(NFX,NFY)
...  
...  u(*) is input vector on the coarse grid of size (ncx*ncy*ncz)
...  pu(*) is output vector, p.w. linear interpolation is used and
...  pu(*) = pu(*) + P*u.  The reason being that in MG we have to
...  perform UNEW = UOLD + P*UCOARSE, and this prolongation works
...  without introducing an extra vector UNEW, it directly gives
...  
...  PU <-- PU+interpolation(U).
...  
...  The vector pu(*) is of size nfx*nfy
*/
REAL two=2.;
REAL half=0.5;

 INT nfxy,ncxy,nctot,nbig,nftot;
INT i,i0,ic0,ic1,j,j0,j1,j00,j10,jc0,jc1;
INT lf1,lf2,lf3,lf4,lc1,lc2,lc3,lc4;

  nfxy = nfx*nfy;
  ncxy = ncx*ncy;
  nctot = ncxy;
  nbig = -2*nctot;
  nftot = nfx*nfy;

  for (i = 0; i < nftot; i++)					/*!*/
    pu[i] = two*pu[i];

  for (j = 1; j <= ncy; j++)
  {
    j0 = 2*j - 1;
    j1 = 2*j;
    j00 = (j0-1)*nfx;
    j10 = j00 + nfx;
    jc0 = (j-1)*ncx;

    if(j == ncy)
      jc1 = nbig;
    else
      jc1 = jc0 + ncx;

    for (i = 1; i <= ncx; i++)
    {
      i0 = 2*i - 2;						/*!*/
      ic0 = i-1;						/*!*/

      if(i == ncx)
        ic1 = nbig-1;						/*!*/
      else
        ic1 = ic0 + 1;
/*
...  Calculate indecies of the neighbors
*/
      lf1 = i0 + j00;
      lf2 = lf1 + 1;
      lf3 = i0 + j10; 
      lf4 = lf3 + 1;
/*... These should give something < 0 if needed */
      lc1 = ic0 + jc0;
      lc2 = ic1 + jc0;
      lc3 = ic0 + jc1;
      lc4 = ic1 + jc1;
/*... Interpolate then: */
      pu[lf1]  = pu[lf1] + u[lc1] + u[lc1];
      if(lc2 > 0 && lc2 <= nctot)
        pu[lf2]  =  pu[lf2] + (u[lc1] + u[lc2]);
      if(lc3 > 0 && lc3 <= nctot)
         pu[lf3] =  pu[lf3] + (u[lc1] + u[lc3]);
      if(lc4 > 0 && lc4 <= nctot)
         pu[lf4]  = pu[lf4] + (u[lc1] + u[lc4]);
    }
  }

  for (i = 0; i < nftot; i++)					/*!*/
    pu[i] = half*pu[i];
}
/*==========================================================================*/
/* pt2d                                                                     */
/*==========================================================================*/
void pt2d(REAL *u, REAL *ru, INT nfx, INT nfy, INT ncx, INT ncy)
{
/*
.... This computes 4 times P^t*x and then scales it appropriately
*/
REAL zero=0.;
REAL half=0.5;

INT ij,nfxy,ncxy,nctot,nbig;
INT i,i0,ic0,ic1,j,j0,j1,j00,j10,jc0,jc1,k;
INT lf1,lf2,lf3,lf4,lc1,lc2,lc3,lc4;

  nfxy = nfx*nfy;
  ncxy = ncx*ncy;
  nctot = ncxy;
  nbig = -2*nctot;

  for (k = 0; k < nctot; k++)					/*!*/
    ru[k] = zero;

  for (j = 1; j <= ncy-1; j++)
  {
    j0 = 2*j - 1;
    j1 = 2*j; 
    j00 = (j0-1)*nfx;
    j10 = j00 + nfx;
    jc0 = (j-1)*ncx;
    jc1 = jc0 + ncx;

    for (i = 1; i <= ncx-1; i++)
    {
      i0 = 2*i - 2;						/*!*/
      ic0 = i-1;						/*!*/
      ic1 = ic0 + 1;
/*
...  Calculate indecies of the neighbors
*/
      lf1 = i0 + j00;
      lf2 = lf1 + 1;
      lf3 = i0 + j10;
      lf4 = lf3 + 1;

      lc1 = ic0 + jc0;
      lc2 = ic1 + jc0;
      lc3 = ic0 + jc1;
      lc4 = ic1 + jc1;

      if(i > 1 && j > 1)
        ru[lc1] = ru[lc1] + u[lf1]+u[lf1]+(u[lf2]+u[lf3]+ u[lf4]);
      if(j > 1)
        ru[lc2]  =  ru[lc2] + u[lf2];
      if(i > 1)
        ru[lc3]  =  ru[lc3] + u[lf3];
      ru[lc4]  =  ru[lc4] + u[lf4];
    }
  }

  for (ij = 0; ij < ncxy; ij++)
    ru[ij] = half*ru[ij];
}
/* 
...TWO COLOR SMOOTHERS
*/
/*==========================================================================*/
/* swep3db (backward)                                                       */
/*==========================================================================*/
void swep3db(REAL *u, REAL *f,
            INT nbegx, INT nendx, INT nstepx,
            INT nbegy, INT nendy, INT nstepy,
            INT nbegz, INT nendz, INT nstepz,
            INT nx, INT ny, INT nz, REAL dd0)
{
INT nxy,k,kp1,k0,km1,j,j0,jp1,jm1,i,i0,i1,i2,i3,i4,i5,i6,nbx1,nex2;

 if(nbegz < nendz || nbegy < nendy || nbegx < nendx) return;

  nxy=nx*ny;
  nbx1=nbegx-1;
  nex2=nendx-2;
  for (k=nbegz; k >= nendz; k+=nstepz)
  {
    kp1= k*nxy;
    k0 = kp1-nxy;
    km1= k0-nxy;
    for (j = nbegy; j >= nendy; j+=nstepy)
    {
      jp1= j*nx;
      j0 = jp1-nx;
      jm1= j0-nx;
      for (i = nbx1; i > nex2; i+=nstepx)			/*!*/
      {
        i0 = i   +  j0    + k0;
        i1 = i+1 +  j0    + k0;
        i2 = i   +  jp1   + k0;
        i3 = i   +  j0    + kp1;
        i4 = i-1 +  j0    + k0;
        i5 = i   +  jm1   + k0;
        i6 = i   +  j0    + km1;
        u[i0] = (f[i0] +  u[i1] + u[i2] + u[i3] + u[i4] + u[i5] + u[i6])*dd0;
      }
    }
  }
}
/*==========================================================================*/
/* swep3df (forward)                                                        */
/*==========================================================================*/
void swep3df(REAL *u, REAL *f,
            INT nbegx, INT nendx, INT nstepx,
            INT nbegy, INT nendy, INT nstepy,
            INT nbegz, INT nendz, INT nstepz,
            INT nx, INT ny, INT nz, REAL dd0)
{
INT nxy,k,kp1,k0,km1,j,j0,jp1,jm1,i,i0,i1,i2,i3,i4,i5,i6;

 if(nbegz > nendz || nbegy > nendy || nbegx > nendx) return;

  nxy=nx*ny;
  for (k=nbegz; k <= nendz; k+=nstepz)
  {
    kp1= k*nxy;
    k0 = kp1-nxy;
    km1= k0-nxy;
    for (j = nbegy; j <= nendy; j+=nstepy)
    {
      jp1= j*nx;
      j0 = jp1-nx;
      jm1= j0-nx;
      for (i = nbegx-1; i < nendx; i+=nstepx)			/*!*/
      {
        i0 = i   +  j0    + k0;
        i1 = i+1 +  j0    + k0;
        i2 = i   +  jp1   + k0;
        i3 = i   +  j0    + kp1;
        i4 = i-1 +  j0    + k0;
        i5 = i   +  jm1   + k0;
        i6 = i   +  j0    + km1;
        u[i0] = (f[i0] +  u[i1] + u[i2] + u[i3] + u[i4] + u[i5] + u[i6])*dd0;
      }
    }
  }
}
/*==========================================================================*/
/*==========================================================================*/
/* rb0b3d                                                                   */
/*==========================================================================*/
void rb0b3d(REAL *u, REAL *f, INT nx, INT ny, INT nz, INT nsweeps, REAL dd0)
{
/*
n1 = n+1
... Backward GAUSS - SEIDEL. ON INPUT THE INITIAL GUESS IS IN U
... On OUT the new iterate is also in U
...  OUTPUT AND THE DIFFERENCE BETWEEN THE INITIAL GUESS AND
...  THIS ITERATE IS IN U
...  the ordering is as follows (e for even o for odd):
...  e-e-e
...  e-e-o
...  e-o-e
...  e-o-o
...  o-e-e
...  o-e-o
...  o-o-e
...  o-o-o
... There are 8 loops that look the same.
*/
INT nxy,nsize,nxe,nxo,nye,nyo,nze,nzo,n0e,n0o,nstep,isweep;

  nxy = nx*ny;
  nsize = nxy*nz;

  nxe=nx-1;
  nxo=nx-2;
  nye=ny-1;
  nyo=ny-2;
  nze=nz-1;
  nzo=nz-2;
  n0e=2;
  n0o=3;
  nstep=-2;

  for (isweep = 1; isweep <= nsweeps; isweep++)
  {
/*...  e-e-e (and going backwards) */
    swep3db(u,f,nxe,n0e,nstep,nye,n0e,nstep,nze,n0e,nstep,nx,ny,nz,dd0);
/*...  e-e-o */
    swep3db(u,f,nxe,n0e,nstep,nye,n0e,nstep,nzo,n0o,nstep,nx,ny,nz,dd0);
/*...  e-o-e */
    swep3db(u,f,nxe,n0e,nstep,nyo,n0o,nstep,nze,n0e,nstep,nx,ny,nz,dd0);
/*...  e-o-o */
    swep3db(u,f,nxe,n0e,nstep,nyo,n0o,nstep,nzo,n0o,nstep,nx,ny,nz,dd0);
/*...  o-e-e */
    swep3db(u,f,nxo,n0o,nstep,nye,n0e,nstep,nze,n0e,nstep,nx,ny,nz,dd0);
/*...  o-e-o */
    swep3db(u,f,nxo,n0o,nstep,nye,n0e,nstep,nzo,n0o,nstep,nx,ny,nz,dd0);
/*...  o-o-e */
    swep3db(u,f,nxo,n0o,nstep,nyo,n0o,nstep,nze,n0e,nstep,nx,ny,nz,dd0);
/*...  o-o-o */
    swep3db(u,f,nxo,n0o,nstep,nyo,n0o,nstep,nzo,n0o,nstep,nx,ny,nz,dd0);
  }
}
/*==========================================================================*/
/* rb0f3d                                                                   */
/*==========================================================================*/
void rb0f3d(REAL *u, REAL *f, INT nx, INT ny, INT nz, INT nsweeps, REAL dd0)
{
/*
n1 = n+1
... FORWARD GAUSS - SEIDEL. ON INPUT THE INITIAL GUESS IS IN U
... THE NEW ITERATE IS IN V ON
...  OUTPUT AND THE DIFFERENCE BETWEEN THE INITIAL GUESS AND
...  THIS ITERATE IS IN U
...  the ordering is as follows (e for even o for odd):
...  o-o-o
...  o-o-e
...  o-e-o
...  o-e-e
...  e-o-o
...  e-o-e
...  e-e-o
...  e-e-e
*/
INT nxy,nsize,nxe,nxo,nye,nyo,nze,nzo,n0e,n0o,nstep,isweep;

  nxy = nx*ny;
  nsize = nxy*nz;

  nxe=nx-1;
  nxo=nx-2;
  nye=ny-1;
  nyo=ny-2;
  nze=nz-1;
  nzo=nz-2;
  n0e=2;
  n0o=3;
  nstep=2;

  for (isweep = 1; isweep <= nsweeps; isweep++)
  {
/*...  o-o-o */
    swep3df(u,f,n0o,nxo,nstep,n0o,nyo,nstep,n0o,nzo,nstep,nx,ny,nz,dd0);
/*...  o-o-e */
    swep3df(u,f,n0o,nxo,nstep,n0o,nyo,nstep,n0e,nze,nstep,nx,ny,nz,dd0);
/*...  o-e-o */
    swep3df(u,f,n0o,nxo,nstep,n0e,nye,nstep,n0o,nzo,nstep,nx,ny,nz,dd0);
/*...  o-e-e */
    swep3df(u,f,n0o,nxo,nstep,n0e,nye,nstep,n0e,nze,nstep,nx,ny,nz,dd0);
/*...  e-o-o */
    swep3df(u,f,n0e,nxe,nstep,n0o,nyo,nstep,n0o,nzo,nstep,nx,ny,nz,dd0);
/*...  e-o-e */
    swep3df(u,f,n0e,nxe,nstep,n0o,nyo,nstep,n0e,nze,nstep,nx,ny,nz,dd0);
/*...  e-e-o */
    swep3df(u,f,n0e,nxe,nstep,n0e,nye,nstep,n0o,nzo,nstep,nx,ny,nz,dd0);
/*...  e-e-e */
    swep3df(u,f,n0e,nxe,nstep,n0e,nye,nstep,n0e,nze,nstep,nx,ny,nz,dd0);
  }
}
/*==========================================================================*/
/* swep2db (backward)                                                       */
/*==========================================================================*/
void swep2db(REAL *u, REAL *f,
            INT nbegx, INT nendx, INT nstepx,
            INT nbegy, INT nendy, INT nstepy,
            INT nx, INT ny, REAL dd0)
{
INT j,j0,jp1,jm1,i,i0,i1,i2,i3,i4,nbx1,nex2;

//fprintf(stdout,"\nB-nx=%i ny=%i nbegx=%i nbendx=%i %i %i*********", nx,ny,nbegx,nendx,nbegy,nendy); 
 if(nbegy < nendy || nbegx < nendx) 
   {
     //     fprintf(stdout,"B-returns..."); 
     return;
   }
 nbx1=nbegx-1;
 nex2=nendx-2;

 for (j = nbegy; j >= nendy; j+=nstepy) 
   {
     jp1 = j*nx;
     j0  = jp1-nx;
     jm1 = j0-nx;
     for (i = nbx1; i > nex2; i+=nstepx)			/*!*/
       {
	 i0 = i   +  j0;
	 i1 = i+1 +  j0;
	 i2 = i   +  jp1;
	 i3 = i-1 +  j0;
	 i4 = i   +  jm1; 
	 //	 fprintf(stdout,"\n%i %i %i %i %i", i0,i1,i2,i3,i4); 
	 u[i0] = (f[i0] +  u[i1] + u[i2] + u[i3] + u[i4])*dd0;
       }
   }
}
/*==========================================================================*/
/* swep2df (forward)                                                        */
/*==========================================================================*/
void swep2df(REAL *u, REAL *f,
            INT nbegx, INT nendx, INT nstepx,
            INT nbegy, INT nendy, INT nstepy,
            INT nx, INT ny, REAL dd0)
{
INT j,j0,jp1,jm1,i,i0,i1,i2,i3,i4;

//fprintf(stdout,"\nF---nbeg,nend nx=%i ny=%i %i %i %i %i*********", nx,ny,nbegx,nendx,nbegy,nendy); 
 if(nbegy > nendy || nbegx > nendx) {
   //   fprintf(stderr,"F-returns..."); 
   return;
 }
 /**/
 for (j = nbegy; j <= nendy; j+=nstepy) 
   {
     jp1 = j*nx;
     j0  = jp1-nx;
     jm1 = j0-nx;
     for (i = nbegx-1; i < nendx; i+=nstepx)			/*!*/
       {
	 i0 = i   +  j0;
	 i1 = i+1 +  j0;
	 i2 = i   +  jp1;
	 i3 = i-1 +  j0;
	 i4 = i   +  jm1; 
	 //	 fprintf(stdout,"\nF---- %i %i %i %i %i", i0,i1,i2,i3,i4); 
	 u[i0] = (f[i0] +  u[i1] + u[i2] + u[i3] + u[i4])*dd0;
       }
   }
}
/*==========================================================================*/
/* 2D case                                                                  */
/*==========================================================================*/
/* rb0b2d                                                                   */
/*==========================================================================*/
void rb0b2d(REAL *u, REAL *f, INT nx, INT ny, INT nsweeps, REAL dd0)
{
INT nsize,nxy,nxe,nxo,nye,nyo,n0e,n0o,nstep,isweep;

/*
n1 = n+1
... FORWARD. THE NEW ITERATE IS IN U
*/

  nsize = nx*ny;
  nxy = nsize;
/*
... BACKWARD Gauss-Seidel we  do the coarse ones in reverse order at the end. 
*/
  nxe=nx-1;
  nxo=nx-2;
  nye=ny-1;
  nyo=ny-2;
  n0e=2;
  n0o=3;
  nstep=-2;

  for (isweep = 1; isweep <= nsweeps; isweep++)
  {
/*... even-even going back (all have e at the end) */
    swep2db(u,f,nxe,n0e,nstep,nye,n0e,nstep,nx,ny,dd0);
/*... odd-even */
    swep2db(u,f,nxo,n0o,nstep,nye,n0e,nstep,nx,ny,dd0);
/*... even-odd */
    swep2db(u,f,nxe,n0e,nstep,nyo,n0o,nstep,nx,ny,dd0);
/*... coarse coarse (odd-odd) */
    swep2db(u,f,nxo,n0o,nstep,nyo,n0o,nstep,nx,ny,dd0);
   }
}
/*==========================================================================*/
/* rb0f2d                                                                   */
/*==========================================================================*/
void rb0f2d(REAL *u, REAL *f, INT nx, INT ny, INT nsweeps, REAL dd0)
{
INT nsize,nxy,nxe,nxo,nye,nyo,n0e,n0o,nstep,isweep;
/*
n1 = n+1
...  FORWARD GAUSS - SEIDEL. ON INPUT THE INITIAL GUESS IS IN U
...  THE NEW ITERATE IS IN U ON
*/

  nsize = nx*ny;
  nxy = nsize;
/*
... FORWARD Gauss-Seidel we first do coarse ones (odd-odd). 
*/
  nxe=nx-1;
  nxo=nx-2;
  nye=ny-1;
  nyo=ny-2;
  n0e=2;
  n0o=3;
  nstep=2;

  for (isweep = 1; isweep <= nsweeps; isweep++)
  {
/*... odd-odd (coarse nodes) and going forward */
    swep2df(u,f,n0o,nxo,nstep,n0o,nyo,nstep,nx,ny,dd0);
/*... even-odd */
    swep2df(u,f,n0e,nxe,nstep,n0o,nyo,nstep,nx,ny,dd0);
/*... odd-even */
    swep2df(u,f,n0o,nxo,nstep,n0e,nye,nstep,nx,ny,dd0);
/*... even-even */
    swep2df(u,f,n0e,nxe,nstep,n0e,nye,nstep,nx,ny,dd0);
  }
}
/*******************************************************************/
/*==========================================================================*/
/* 3D V-cycle                                                                  */
/*==========================================================================*/
void v_gmg3D(REAL dd, REAL *u, REAL *f, INT levels, INT *ke,\
		REAL *tk, REAL *rk, REAL *uk,\
		INT nx, INT ny, INT nz, INT nsweeps, REAL *errnrm)
{
  INT j,k,nt,nfx,nfy,nfz,nsizef;
  INT  ncx=0,ncy=0,ncz=0;
  REAL *ukker=NULL;
  REAL *ukkerc=NULL;
  REAL *rkker=NULL;
  REAL *rkkerc=NULL;
  
  k = levels;
  nfx = nx;
  nfy = ny;
  nfz = nz;
  nt = nx*ny*nz;
  ukker= uk + ke[levels - k]-nt;
  rkker=rk + ke[levels - k]-nt;
  for (k=levels-1; k >= 1; k--) {
    nsizef = nfx*nfy*nfz;
    ncx = (nfx+1)/2;
    ncy = (nfy+1)/2;
    ncz = (nfz+1)/2;

    if(k == (levels-1)) {
      rb0f3d(u,f,nfx,nfy,nfz,nsweeps,dd);
      res0(f,u,tk,nfx,nfy,nfz);
    } else {
      for (j = 0; j < nsizef; j++) *(ukker+j) = 0.;
      rb0f3d(ukker,rkker,nfx,nfy,nfz,nsweeps,dd);
      res0(rkker,ukker,tk,nfx,nfy,nfz);
    }
    rkkerc = rk + ke[levels - k]-nt;
    ukkerc = uk + ke[levels - k]-nt;
    pt3d(tk,rkkerc,nfx,nfy,nfz,ncx,ncy,ncz);
    nfx = ncx;
    nfy = ncy;
    nfz = ncz;
    rkker = rkkerc;
    ukker=ukkerc;
  }
  nsizef=nfx*nfy*nfz;
  for (j = 0; j < nsizef; j++) *(ukker+j) = 0.;
  rb0f3d(ukker,rkker,nfx,nfy,nfz,nsweeps,dd);

  for (k=1; k < levels; k++)
    {
      ukker=uk+ke[levels - k - 1]-nt;
      rkker=rk+ke[levels - k - 1]-nt;
      ukkerc=uk+ke[levels - k]-nt;
      nfx = 2*(ncx-1) + 1;
      nfy = 2*(ncy-1) + 1;
      nfz = 2*(ncz-1) + 1;

      if(k < (levels-1))
	{
	  p3d(ukkerc,ukker,nfx,nfy,nfz,ncx,ncy,ncz);
	  rb0b3d(ukker,rkker,nfx,nfy,nfz,nsweeps,dd);
	}
      else
	{
	  p3d(ukkerc,u,nfx,nfy,nfz,ncx,ncy,ncz);
	  rb0b3d(u,f,nfx,nfy,nfz,nsweeps,dd);
	  /* 3D one more computation of the norm of the residual */
	  res0(f,u,tk,nfx,nfy,nfz);
	  *errnrm=ddot0(tk,tk,nt);
	  *errnrm = sqrt(*errnrm);
	}

      ncx = nfx;
      ncy = nfy;
      ncz = nfz;
    }
}

/*==========================================================================*/
/* 3D full multigrid V-cycle                                                */
/*==========================================================================*/
void fmg_gmg3D(REAL dd, REAL *u, REAL *f, INT levels, INT *ke, \
		  REAL *tk, REAL *rk, REAL *uk, \
		  INT nx, INT ny, INT nz, INT nsweeps, INT nvcycles, REAL *errnrm)
{
  INT j,k,nt,nfx,nfy,nfz, nsizef;
  INT  ncx=0,ncy=0,ncz=0;
  INT ket[40];
  INT nxxx, nyyy, nzzz, ntk,lvlk;
  REAL *ukker=NULL;
  REAL *ukkerc=NULL;
  REAL *rkker=NULL;
  REAL *rkkerc=NULL;
	
  k = levels;

  nfx = nx;
  nfy = ny;
  nfz = nz;
  nt = nx*ny*nz;

  /* restrict the right hand side f from finest grid to the coarest grid */
  /*	
	rkker=rk + ke[levels - k]-nt;
	ukker=uk + ke[levels - k]-nt;
  */
  rkker=rk;
  ukker=uk;
  for (k=levels-1; k >= 1; k--) {
    nsizef = nfx*nfy*nfz;
    ncx = (nfx+1)/2;
    ncy = (nfy+1)/2;
    ncz = (nfz+1)/2;
    rkkerc = rk + ke[levels - k] - nt;
    ukkerc = uk + ke[levels - k] - nt;
    
    if(k == (levels-1)) {
      pt3d(f, rkker, nfx, nfy, nfz, ncx, ncy, ncz);
    } else {
      pt3d(rkker, rkkerc, nfx, nfy, nfz, ncx, ncy, ncz);
    }
    
    nfx = ncx;
    nfy = ncy;
    nfz = ncz;
    rkker = rkkerc;
    ukker = ukkerc;
  }
  nsizef = nfx*nfy*nfz;
  for (j = 0; j < nsizef; j++) *(ukker+j) = 0.;

  rb0f3d(ukker,rkker,nfx,nfy,nfz, nsweeps,dd);

  for (k=1; k < levels; k++)
    {
      ukker  = uk + ke[levels - k-1]-nt;
      rkker  = rk + ke[levels - k-1]-nt;
      ukkerc = uk + ke[levels - k] - nt;
      rkkerc = rk + ke[levels - k] - nt;

      nfx = 2*(ncx-1) + 1;
      nfy = 2*(ncy-1) + 1;
      nfz = 2*(ncz-1) + 1;
      lvlk=k+1;		
      if(k < levels-1)
	{
	  /*    fprintf(stderr,"\ndiff= %i" , ukkerc-ukker); */
	  p3d(ukkerc,ukker,nfx,nfy,nfz, ncx,ncy, ncz);
	  /* preform V-cycle on level k+1 */
	  ket[0] = 1;
	  nxxx = nfx;
	  nyyy = nfy;
	  nzzz = nfz;
	  for(j =1; j < lvlk; j++){
	    ntk = nxxx * nyyy * nzzz;
	    ket[j] = ket[j-1] + ntk;
	    nxxx = (nxxx - 1)/2 + 1;
	    nyyy = (nyyy - 1)/2 + 1;
	    nzzz = (nzzz - 1)/2 + 1;
	  }
	  for (j=0; j<nvcycles; j++) {
	    v_gmg3D(dd, ukker, rkker, lvlk, ket, tk, rkkerc, ukkerc,\
		       nfx, nfy, nfz, nsweeps, errnrm);
	  }
	
	} else {
	p3d(ukkerc,u,nfx,nfy,nfz, ncx,ncy,ncz);
	for (j=0; j<nvcycles; j++) {
	  v_gmg3D(dd, u, f, levels, ke, tk, rk, uk,\
		     nfx, nfy, nfz, nsweeps, errnrm);
	}
      }
      ncx = nfx;
      ncy = nfy;
      ncz = nfz;
    }
}
/*==========================================================================*/
/* lap3d                                                                    */
/*==========================================================================*/
void gmg3D(INT isfmg, INT levels, INT nx, INT ny, INT nz, INT nt,\
	      REAL *f, REAL *u, REAL *tk, REAL *uk, REAL *rk, INT *ke,\
	      REAL tol, INT nsweeps, INT nvcycles, INT mxitr, INT *iterv,\
	      REAL *errrel, REAL *reduce, REAL h)
{
  REAL err0, errnrm, errold, dd;
  INT iter,n,k;

  n=nx;
  dd = 1./6.;
  /* since we assume that we start with 0 init guess, this will be the norm of f */
  res0(f,u,tk,nx,ny,nz);
  err0=ddot0(tk,tk,nt);
  err0=sqrt(err0);
  if(err0 <= 1.e-10) err0 = 1e+00;
  if(isfmg) {
    fmg_gmg3D(dd, u, f, levels, ke, tk, uk, \
		 rk, nx, ny, nz, nsweeps, nvcycles, &errnrm);
    *iterv=-10;
    *errrel=errnrm/err0;
  } else {
/* CALCULATE RESIDUAL */
    res0(f,u,tk,nx,ny,nz);
    err0=ddot0(tk,tk,nt);
    err0 = sqrt(err0);
    if(err0 <= 1.e-10) err0 = 1e+00;
    errnrm = err0;
    errold = err0;
    iter=0;
    do {
      iter = iter + 1;
      errold = errnrm;
      v_gmg3D(dd,u,f,levels,ke,tk,uk,rk,nx,ny,nz,nsweeps,&errnrm);
    } while (iter < mxitr && errnrm/err0 > tol);

    *iterv=iter;
    *reduce=errnrm/errold;
    *errrel=errnrm/err0;
  }
  /* fprintf(stdout,"\nPrint 2\n"); */
  /* for (k=0; k<nt; k++) { */
  /*   fprintf(stdout,"%d %16.8e\n",k,u[k]); */
  /* } */
}
/*==========================================================================*/
/* 2D V-cycle                                                               */
/*==========================================================================*/
void v_gmg2D(REAL dd, REAL *u, REAL *f, INT levels, INT *ke,\
		REAL *tk, REAL *rk, REAL *uk,
		INT nx, INT ny, INT nsweeps, REAL *errnrm)
{
  INT j,k,nt,nfx,nfy,nsizef;
  INT ncx=0,ncy=0;
  REAL *ukker=NULL;
  REAL *ukkerc=NULL;
  REAL *rkker=NULL;
  REAL *rkkerc=NULL;

  k = levels;
  nfx = nx;
  nfy = ny;
  nt = nx*ny;
  ukker= uk + ke[levels - k]-nt;
  rkker=rk + ke[levels - k]-nt;
  for (k=levels-1; k >= 1; k--)
    {
      nsizef = nfx*nfy;
      ncx = (nfx+1)/2;
      ncy = (nfy+1)/2;

      if(k == (levels-1))
	{
	  rb0f2d(u,f,nfx,nfy,nsweeps,dd);
	  res02d(f,u,tk,nfx,nfy);
	}
      else
	{
	  for (j = 0; j < nsizef; j++) /* ??? */
	    *(ukker+j) = 0.;
	  rb0f2d(ukker,rkker,nfx,nfy,nsweeps,dd);
	  res02d(rkker,ukker,tk,nfx,nfy);
	}

      rkkerc = rk + ke[levels - k]-nt;
      ukkerc = uk + ke[levels - k]-nt;
      pt2d(tk,rkkerc,nfx,nfy,ncx,ncy);
      nfx = ncx;
      nfy = ncy;
      rkker = rkkerc;
      ukker=ukkerc;
    }

  for (j = 0; j < nfx*nfy; j++) /* ??? */
    *(ukker+j) = 0.;

  rb0f2d(ukker,rkker,nfx,nfy,nsweeps,dd);

  for (k=1; k < levels; k++)
    {
      ukker=uk+ke[levels - k-1]-nt;
      rkker=rk+ke[levels - k-1]-nt;
      ukkerc=uk+ke[levels - k]-nt;

      nfx = 2*(ncx-1) + 1;
      nfy = 2*(ncy-1) + 1;

      if(k < levels-1) {
	//    fprintf(stderr,"\ndiff= %i" , ukkerc-ukker);
	p2d(ukkerc,ukker,nfx,nfy,ncx,ncy);
	rb0b2d(ukker,rkker,nfx,nfy,nsweeps,dd);
      } else {
	//      fprintf(stderr,"\ndiffend= %i" , ukkerc-uk);
	p2d(ukkerc,u,nfx,nfy,ncx,ncy);
	rb0b2d(u,f,nfx,nfy,nsweeps,dd);
	/* 2D one more computation of the norm of the residual */
	res02d(f,u,tk,nfx,nfy);
	*errnrm=ddot0(tk,tk,nt);
	*errnrm = sqrt(*errnrm);
      }
      ncx = nfx;
      ncy = nfy;
    }
}
/*==========================================================================*/
/* 2D full multigrid V-cycle                                                */
/*==========================================================================*/
void fmg_gmg2D(REAL dd, REAL *u, REAL *f, INT levels, INT *ke,\
		  REAL *tk, REAL *rk, REAL *uk,\
		  INT nx, INT ny, INT nsweeps, INT nvcycles, REAL *errnrm)
{
  INT j,k,nt,nfx,nfy,nsizef;
  INT ncx=0,ncy=0;
  INT ket[40];
  INT nxxx, nyyy, kk, ntk,lvlk;
  REAL *ukker=NULL;
  REAL *ukkerc=NULL;
  REAL *rkker=NULL;
  REAL *rkkerc=NULL;

  k = levels;

  nfx = nx;
  nfy = ny;
  nt = nx*ny;
  /*
    rkker=rk + ke[levels - k]-nt;
    ukker=uk + ke[levels - k]-nt;
  */
  rkker=rk;
  ukker=uk;
  for (k=levels-1; k >= 1; k--) {
    nsizef = nfx*nfy;
    ncx = (nfx+1)/2;
    ncy = (nfy+1)/2;
    rkkerc = rk + ke[levels - k] - nt;
    ukkerc = uk + ke[levels - k] - nt;
    if(k == (levels-1)) {
      pt2d(f, rkker, nfx, nfy, ncx, ncy);
    } else {
      pt2d(rkker, rkkerc, nfx, nfy, ncx, ncy);
    }
    nfx = ncx;
    nfy = ncy;
    rkker = rkkerc;
    ukker = ukkerc;
    /*    fprintf(stderr," The difference:  %li\n ",rkker-rk); */
  }

  for (j = 0; j < nfx*nfy; j++)  *(ukker+j) = 0e+00;
  rb0f2d(ukker,rkker,nfx,nfy,nsweeps,dd);
  for (k=1; k < levels; k++) {
    ukker  = uk + ke[levels - k-1]-nt;
    rkker  = rk + ke[levels - k-1]-nt;
    ukkerc = uk + ke[levels - k] - nt;
    rkkerc = rk + ke[levels - k] - nt;

    nfx = 2*(ncx-1) + 1;
    nfy = 2*(ncy-1) + 1;
    lvlk=k+1;		
    if(k < levels-1) {
      //    fprintf(stderr,"\ndiff= %i" , ukkerc-ukker);
      p2d(ukkerc,ukker,nfx,nfy,ncx,ncy);
      //preform V-cycle on level k+1
      ket[0] = 1;
      nxxx = nfx;
      nyyy = nfy;
      for(kk =1; kk < k+1; kk++){
	ntk = nxxx * nyyy;
	ket[kk] = ket[kk-1] + ntk;
	nxxx = (nxxx - 1)/2 + 1;
	nyyy = (nyyy - 1)/2 + 1;
      }
      /*	 fprintf(stderr," The difference SECOND:  %li %g\n ",rkkerc-rk,dd);*/
      for (j=0; j< nvcycles; j++) {
	v_gmg2D(dd, ukker, rkker, lvlk, ket, tk, rkkerc, ukkerc,\
		   nfx, nfy, nsweeps, errnrm);
      }
    } else {
      p2d(ukkerc,u,nfx,nfy,ncx,ncy);
      for (j=0; j < nvcycles; j++) {
	v_gmg2D(dd, u, f, levels, ke, tk, rk, uk, \
		   nfx, nfy, nsweeps, errnrm);
      }
    }
    ncx = nfx;
    ncy = nfy;
  }
}
/*==========================================================================*/
/* lap2d                                                                    */
/*==========================================================================*/
void gmg2D(INT isfmg, INT levels, INT nx, INT ny, INT nt,\
	      REAL *f, REAL *u, REAL *tk, REAL *uk, REAL *rk, INT *ke,\
	      REAL tol, INT nsweeps, INT nvcycles, \
	      INT mxitr, INT *iterv, REAL *errrel, REAL *reduce, REAL h)
{
  REAL err0, errnrm, errold, dd;
  INT iter,n;

  dd = 0.25e00;
  n = nx;
  res02d(f,u,tk,nx,ny);
  err0=ddot0(tk,tk,nt);
  err0 = sqrt(err0);
  if(err0 <= 1.e-10) err0 = 1.;
  if(isfmg) {
    fmg_gmg2D(dd, u, f, levels, ke, tk, uk, rk, nx, ny, nsweeps, nvcycles, &errnrm);
    *errrel=errnrm/err0;
  } else {
    errnrm = err0;
    errold = err0;
    iter = 0;
    do {
      iter = iter + 1;
      errold = errnrm;
      v_gmg2D(dd,u,f,levels,ke,tk,uk,rk,nx,ny,nsweeps,&errnrm);
    } while (iter < mxitr && errnrm/err0 > tol);

    *iterv=iter;
    *reduce=errnrm/errold;
    *errrel=errnrm/err0;
  }
}
/****************************************************************/
static void print_info(FILE *fp,		\
		       INT isfmg,		\
		       INT lvlloc,		\
		       INT nx,			\
		       INT nspdim,		\
		       INT nsweeps,		\
		       REAL tol,		\
		       long int memtotal,		\
		       long int ntall)
{
  /* Print some info */
  fprintf(fp,"\n");
  fprintf(fp,"*Number of levels is set to  :  %7li\n",(long int)lvlloc);
  fprintf(fp,"*Num. points in one direction:  %7li\n",(long int)nx);
  fprintf(fp,"*Spatial dimension is set to :  %7li\n",(long int)nspdim);
  fprintf(fp,"*Num. of Smoothing steps     :  %7li\n",(long int)nsweeps);
  if(isfmg) {
    fprintf(fp,"\n**FMG FMG FMG**\n");
  } else {
    fprintf(fp,"\n**V-cycle iterations\n");
    fprintf(fp,"\n**Stopping criteria: |rk|/|r0| < %12.3g\n",tol);
  }
  fprintf(fp, "\n***** NUMBER OF UNKNOWNS=%li\n",ntall);
  fprintf(fp,    "\n**MEMORY NEEDED (DOUBLE):  (%12ld) ==> Bytes/dof%10.2f\n",(long int)memtotal, ((double )memtotal)/((double ) ntall)*sizeof(REAL));
}

/**********************************MEMORY*********************************/
void memgmg(INT isfmg, INT levels, INT nspdim, INT nsweeps,	\
	    INT *malls,INT *ke,INT nx,INT ny,INT nz,		\
	    INT *ntotal, REAL *h, REAL tol)
{
  /*  REAL x0 = 0e+00, y0 = 0e+00, z0 = 0e+00;*/
  REAL hx,hy,hz;
  long int nxx,nyy,nzz,ntk,ntall;
  long int klastt,mallsloc,memuf,memwork,memtotal;
  INT lvlloc;
  INT i;
  hx = 1e+00/( (REAL ) (nx-1));
  hy = 1e+00/( (REAL ) (ny-1));
  if(nz != 1 ) 
    hz = 1e+00/( (REAL ) (nz-1));
  else 
    hz = 1e+00;
  *h=hx;
  lvlloc = levels;
  long int *keloc=malloc((lvlloc+1)*sizeof(long int));

  nxx = (long int )( nx);
  nyy = (long int )( ny);
  nzz = (long int )( nz);
  
  ntall = nxx*nyy*nzz;
    //
  keloc[0] = 0;
  for (i = 0 ; i < lvlloc-1; i++) {
    ntk = nxx*nyy*nzz;
    keloc[i+1] = keloc[i] + ntk;
    nxx = (nxx - 1)/2+1;
    nyy = (nyy - 1)/2+1;
    nzz = (nzz - 1)/2+1;
  }
  klastt = keloc[lvlloc-1] + nxx*nyy*nzz;
  mallsloc=klastt - ntall;
  memuf = ntall+ntall;
  memwork = ntall  + mallsloc + mallsloc;
  memtotal = memuf + memwork;
  /*show info*/
  print_info(stdout,
	     isfmg,				\
	     lvlloc,				\
	     nx,				\
	     nspdim,				\
	     nsweeps,				\
	     tol,				\
	     memtotal,				\
	     ntall);
  /* Convert to integer 4 */
  /*--------------------------------------------------------*/
  *malls=((INT )mallsloc);
  *ntotal=((INT )ntall); 
  for (i = 0 ; i< lvlloc; i++) {
    ke[i] = ((INT ) keloc[i]);
  }
  /*--------------------------------------------------------*/
  free(keloc);
  return;
}
