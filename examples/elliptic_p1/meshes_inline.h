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
  const INT ns=2,nv=4;
  const INT dim1=dim+1;
  // 2*(dim+1)=6
  INT el2v[6]={0,2,3,				\
	       0,1,3};  
  // 4*(dim)=8
  REAL x[8]={0.,0.,				\
	     1.,0.,				\
	     0.,1.,				\
	     1.,1.};
  scomplex *sc=haz_scomplex_init(dim,ns, nv,dim);
  memcpy(sc->x,x,nv*dim*sizeof(REAL));
  memcpy(sc->nodes,el2v,ns*dim1*sizeof(INT));
  return sc;
}
scomplex *mesh3d()
{
  const INT  dim=3;  
  const INT ns=6,nv=8;
  const INT dim1=dim+1;
  // 6*(dim+1)=24
  INT el2v[24]={0,4,6,7,			\
		0,4,5,7,			\
		0,2,6,7,			\
		0,2,3,7,			\
		0,1,5,7,			\
		0,1,3,7};
  // 8*(dim)=24  
  REAL x[24]={0.,0.,0.,			\
	      0.,0.,1.,			\
	      0.,1.,0.,			\
	      0.,1.,1.,			\
	      1.,0.,0.,			\
	      1.,0.,1.,			\
	      1.,1.,0.,			\
	      1.,1.,1.};
  scomplex *sc=haz_scomplex_init(dim,ns, nv,dim);
  memcpy(sc->x,x,nv*dim*sizeof(REAL));
  memcpy(sc->nodes,el2v,ns*dim1*sizeof(INT));
  return sc;
}
/*****************************************************/
scomplex *mesh4d()
{
  const INT  dim=4;  
  const INT ns=24,nv=16;
  const INT dim1=dim+1;
  // 24*(dim+1)=120
  INT el2v[120]={0,8,12,14,15,			\
		 0,8,12,13,15,			\
		 0,8,10,14,15,			\
		 0,8,10,11,15,			\
		 0,8,9,13,15,			\
		 0,8,9,11,15,			\
		 0,4,12,14,15,			\
		 0,4,12,13,15,			\
		 0,4,6,14,15,			\
		 0,4,6,7,15,			\
		 0,4,5,13,15,			\
		 0,4,5,7,15,			\
		 0,2,10,14,15,			\
		 0,2,10,11,15,			\
		 0,2,6,14,15,			\
		 0,2,6,7,15,			\
		 0,2,3,11,15,			\
		 0,2,3,7,15,			\
		 0,1,9,13,15,			\
		 0,1,9,11,15,			\
		 0,1,5,13,15,			\
		 0,1,5,7,15,			\
		 0,1,3,11,15,			\
		 0,1,3,7,15};
  
  // 16*(dim)=64  
  REAL x[64]={0.,0.,0.,0.,			\
	      0.,0.,0.,1.,			\
	      0.,0.,1.,0.,			\
	      0.,0.,1.,1.,			\
	      0.,1.,0.,0.,			\
	      0.,1.,0.,1.,			\
	      0.,1.,1.,0.,			\
	      0.,1.,1.,1.,			\
	      1.,0.,0.,0.,			\
	      1.,0.,0.,1.,			\
	      1.,0.,1.,0.,			\
	      1.,0.,1.,1.,			\
	      1.,1.,0.,0.,			\
	      1.,1.,0.,1.,			\
	      1.,1.,1.,0.,			\
	      1.,1.,1.,1.};      
  scomplex *sc=haz_scomplex_init(dim,ns, nv,dim);
  memcpy(sc->x,x,nv*dim*sizeof(REAL));
  memcpy(sc->nodes,el2v,ns*dim1*sizeof(INT));
  return sc;
}