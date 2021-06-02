/*! \file src/amr/meshes_inline.h
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note Just 2d and 3d mesh examples entered statically by hand/
 *
 *  \note: modified by ltz on 20210601
 *
 */
//
scomplex *mesh2d()
{
  const INT  dim=2;  
  const INT ns=8,nv=9;
  const INT dim1=dim+1;
  INT el2v[72]={0,  1,  3,			\
                3,  1,  4,			\
                1,  2,  4,			\
                4,  2,  5,			\
                3,  4,  6,			\
                6,  4,  7,			\
                4,  5,  7,			\
                7,  5,  8 }  ;
  REAL x[18]={ 0.0e0,  0.0e0,			\
	       0.5e0,  0.0e0,			\
	       1.0e0,  0.0e0,			\
	       0.0e0,  0.5e0,			\
	       0.5e0,  0.5e0,			\
	       1.0e0,  0.5e0,			\
	       0.0e0,  1.0e0,			\
	       0.5e0,  1.0e0,			\
	       1.0e0,  1.0e0 } ;
  scomplex *sc=haz_scomplex_init(dim,ns, nv);
  memcpy(sc->x,x,nv*dim*sizeof(REAL));
  memcpy(sc->nodes,el2v,ns*dim1*sizeof(INT));
  return sc;
}
scomplex *mesh3d()
{
  const INT  dim=3;  
  const INT ns=48,nv=27;
  const INT dim1=dim+1;
  INT el2v[192]={0,10,9,13,				\
		 0,1,10,13,				\
		 0,4,1,13,				\
		 0,9,12,13,				\
		 0,12,3,13,				\
		 0,3,4,13,				\
		 1,11,10,14,				\
		 1,2,11,14,				\
		 1,5,2,14,				\
		 1,10,13,14,				\
		 1,13,4,14,				\
		 1,4,5,14,				\
		 3,13,12,16,				\
		 3,4,13,16,				\
		 3,7,4,16,				\
		 3,12,15,16,			\
		 3,15,6,16,				\
		 3,6,7,16,				\
		 4,14,13,17,			\
		 4,5,14,17,				\
		 4,8,5,17,				\
		 4,13,16,17,			\
		 4,16,7,17,				\
		 4,7,8,17,				\
		 9,19,18,22,			\
		 9,10,19,22,			\
		 9,13,10,22,			\
		 9,18,21,22,			\
		 9,21,12,22,			\
		 9,12,13,22,			\
		 10,20,19,23,			\
		 10,11,20,23,			\
		 10,14,11,23,			\
		 10,19,22,23,			\
		 10,22,13,23,			\
		 10,13,14,23,			\
		 12,22,21,25,			\
		 12,13,22,25,			\
		 12,16,13,25,			\
		 12,21,24,25,			\
		 12,24,15,25,			\
		 12,15,16,25,			\
		 13,23,22,26,			\
		 13,14,23,26,			\
		 13,17,14,26,			\
		 13,22,25,26,			\
		 13,25,16,26,			\
		 13,16,17,26};
    // coordinates of vertices.
  REAL x[81]={0.0e0,0.0e0,0.0e0,			\
	      0.5e0,0.0e0,0.0e0,			\
	      1.0e0,0.0e0,0.0e0,			\
	      0.0e0,0.5e0,0.0e0,			\
	      0.5e0,0.5e0,0.0e0,			\
	      1.0e0,0.5e0,0.0e0,			\
	      0.0e0,1.0e0,0.0e0,			\
	      0.5e0,1.0e0,0.0e0,			\
	      1.0e0,1.0e0,0.0e0,			\
	      0.0e0,0.0e0,0.5e0,			\
	      0.5e0,0.0e0,0.5e0,			\
	      1.0e0,0.0e0,0.5e0,			\
	      0.0e0,0.5e0,0.5e0,			\
	      0.5e0,0.5e0,0.5e0,			\
	      1.0e0,0.5e0,0.5e0,			\
	      0.0e0,1.0e0,0.5e0,			\
	      0.5e0,1.0e0,0.5e0,			\
	      1.0e0,1.0e0,0.5e0,			\
	      0.0e0,0.0e0,1.0e0,			\
	      0.5e0,0.0e0,1.0e0,			\
	      1.0e0,0.0e0,1.0e0,			\
	      0.0e0,0.5e0,1.0e0,			\
	      0.5e0,0.5e0,1.0e0,			\
	      1.0e0,0.5e0,1.0e0,			\
	      0.0e0,1.0e0,1.0e0,			\
	      0.5e0,1.0e0,1.0e0,			\
	      1.0e0,1.0e0,1.0e0};
  //
  scomplex *sc=haz_scomplex_init(dim,ns, nv);
  memcpy(sc->x,x,nv*dim*sizeof(REAL));
  memcpy(sc->nodes,el2v,ns*dim1*sizeof(INT));
  return sc;
}
