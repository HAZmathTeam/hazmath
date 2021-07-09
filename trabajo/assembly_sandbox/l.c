/*! \file src/amr/laplace_p1.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note Solving Laplace equation with linear finite elements on sequence of grids. 
 *
 *  \note: modified by ltz on 20210531
 *
 */
/*********************************************************************/
#include "hazmath.h"
#include "meshes_inline.h"
/****************************************************************************/
int main(int argc, char *argv[])
{
  INT dim=2;// 2d,3d,4d... example
  INT jlevel,k;
  scomplex *sc;
  switch(dim){
  case 4:
    sc=mesh4d();    
    break;
  case 3:
    sc=mesh3d();    
    break;
  default:
    sc=mesh2d();
    uniformrefine2d(sc);
  }
  return 0;
}


