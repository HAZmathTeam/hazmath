#include "hazmath.h"
/******************************************************************/
static SHORT chk_above(INT dim, INT node, \
                       REAL *x, REAL *xmass, REAL *sn, \
                       REAL factorial, \
                       INT face, INT *fv, void *wrk)
{
  /*
    checks if a node is above a simplex of dimension (dim-1), i.e. a face.
    work space: wrk should be at least dim*(dim+2) REALS and dim integers.
    returns 0 if the node is below the face (already in the convex hull) and 1 otherwise.
  */
  INT dim1 = dim+1,l,i,j,ln;
  REAL di;
  REAL *xnode=(REAL *)(x+node*dim);
  // REAL *bt=(REAL *)wrk;
  // REAL *piv=bt+dim*dim;
  // REAL *snloc=piv+dim;
  // INT *p = (INT *)(wrk+(dim*dim + dim + dim)*sizeof(REAL));
  // this below only needed if the point lies on the face. so we ignore it for now.
  // REAL vo;
  // for (j = 0;j<dim;j++){
  //   ln=fv[j]*dim;
  //   for(i=0;i<dim;i++)
  //     bt[j*dim+i] = x[ln+i]-xnode[i];
  // }
  // if(lufull(1, dim, &vo, bt,p,piv))
  //     return 1; // node lies on the face and the distance is zero and it is considered as being above.
  // // sn should be (-1,...,-1) and compute the gradient of lambda0=-grad(sum(lambda_k))
  di=0e0;
  for (l=0;l<dim;l++){
    // take the inner product of (x-x0).n for this face;
    di+=sn[l]*(xnode[l]-xmass[l]);
  }
  fprintf(stdout,"\n(node,face)=(%i,%i); disti=%e;",node,face,di);
  print_full_mat_int(1,dim,fv,"face_vertex");
  if(di> (-1e-16))
    return 1;// the node is above the face.
  else
    return 0;// the node is below the face.
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
void visible(INT dim, REAL factorial,     \
             INT root, INT node, iCSRmat *stos, \
             INT *nfaces, iCSRmat *fv, SHORT *maskf, \
             REAL *x, REAL **xmass, REAL **sn, void *wrk)
{
	/* node is the node whose visible faces give the horizon ridges */
  INT chk,lsize,i,j,k,l,jk,iaa,iab,ila,ilb,face;
	INT maskf_value=(node+1);
  ivector *vset=malloc(1*sizeof(ivector));
  vset->row=1;
  vset->val=calloc(vset->row,sizeof(INT));
	vset->val[0]=root;
  ila = 0 ; ilb = 1;
  maskf[root]=maskf_value;
  while(1) {
    vset->row = ilb;
    lsize = vset->row;
    for (i=ila;i<ilb;i++){
      face = vset->val[i];
      iaa = stos->IA[face];
      iab = stos->IA[face+1];
      for(jk=iaa;jk<iab;jk++) {
				j = stos->JA[jk];
        if(j<0) continue;
        // maskf_value here should be the node number + 1.
        // print_full_mat(dim,1,sn[j],"snj");
        // print_full_mat(dim,1,xmass[j],"xmassj");
        fprintf(stdout,"\nface=%d; mask=%d; node=%d; test_face=%d",face,maskf[j],node,j);
        chk=chk_above(dim,node,x,xmass[j],sn[j], \
                      factorial, \
                      j,(fv->JA+fv->IA[j]),wrk);
				if((maskf[j] != maskf_value)&&chk) {
					vset->val=realloc(vset->val,(vset->row+1)*sizeof(INT));
	  			vset->val[vset->row] = j;
					maskf[j] = maskf_value;
          fprintf(stdout,"\n%d <= %d",vset->row,lsize);
	  			vset->row++;
				} //if
      }  //for
    }   //for
    ila = ilb;
    ilb = vset->row;
    if(vset->row <= lsize) break;
  } //while
  print_full_mat_int(1,vset->row,vset->val,"visible");
  nfold=*nfaces;
  nfnew=nfold;
  // add faces to the list formed by the point and horizon ridges.
  for(jk=0;jk<vset->row;jk++){
    face=vset->val[jk];
    iaa=stos->IA[face];
    iab=stos->IA[face+1];
    for(l=iaa;l<iab;l++){
      i=l-iaa;
      nbr=stos->JA[l];
      if((nbr<0) || (maskf[abs(nbr)]!=maskf_value)){
        // this subsimplex is a ridge. we will add a new face
        nfnew++;
        maskf=realloc(maskf,nfnew*sizeof(SHORT));
        maskf[nfnew-1]=0;
        //
        nfadd=nfnew-nfold;
        fvloc->IA=realloc(fvloc->IA,(nfadd+1)*sizeof(INT));
        fvloc->nnz=fvloc->IA[nfadd-1]+dim;
        fvloc->IA[nfadd]=fvloc->nnz;
        fvloc->JA=realloc(fvloc->JA,fvloc->nnz*sizeof(INT));
        fvloc->val=realloc(fvloc->val,fvloc->nnz*sizeof(INT));
        //
        lsize=1;
        for(k=0;k<dim;k++){
          // find where in
          if(i==k) continue;
          // point on the ridge opposite to the vertex opposite to the subsimplex i.
          kpt=fv->JA[k+fv->IA[face]];
          fvloc[lsize]=kpt;
          lsize++;
        }
        fvloc[0]=node;
      }
    }
  }
//  memcpy((xf+lsize*dim),(x+kpt*dim),dim*sizeof(REAL));
//  memcpy(xf,(x+node*dim),dim*sizeof(REAL));
  return;
}
