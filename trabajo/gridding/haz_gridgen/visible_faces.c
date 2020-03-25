#include "hazmath.h"
/******************************************************************/
static SHORT chk_above(INT dim, REAL fact, INT face, INT *fv, \
                REAL *x, REAL *xmass, REAL *sn0, void *wrk)
{
  /*
    checks if a point is above a simplex.
    work space: wrk should be at least dim*(dim+2) REALS and dim integers.
    returns 0 if the point is below the simplex (already in the convex hull) and 1 otherwise.
  */
  // INT dim1 = dim+1,i,j,ln,ln1;
  // REAL vo,di;
  // REAL *bt=(REAL *)wrk;
  // REAL *piv=bt+dim*dim;
  // REAL *sn0=piv+dim;
  // INT *p = (INT *)(wrk+(dim*dim + dim + dim)*sizeof(REAL));
  // ln1=node*dim;
  // for (j = 0;j<dim;j++){
  //   ln=fv[j]*dim;
  //   for(i=0;i<dim;i++){
  //     bt[j*dim+i] = x[ln+i]-x[ln1+i];
  //   }
  // }
  // if(lufull(1, dim, &vo, bt,p,piv)) {
  //   dist[0]=0.;
  //   return 1; // this is a degenerate simplex, the distance is zero and it is considered above.
  // } else
  //   vo=vo/fact;
  // memset(sn0,0,dim1*dim*sizeof(REAL));
  // for (j = 0; j<dim;j++) sn0[j]=-1.; // inward
  // // sn should be (-1,...,-1) and compute the gradient of lambda0=-grad(sum(lambda_k))
  // solve_pivot(0,dim,bt,sn, p,piv);
  // for(j=0;j<nfaces;j++){
  //   // construct the tet:
  //   ln=j*dim;
  //   disti=0e0;
  //   for (l=0;l<dim;l++){
  //     // take the inner product of (x-x0).n for this face;
  //     disti+=snsn[ln+l]*(x[i*dim+l]-xmass[j][l]);
  //   }
  //   fprintf(stdout,"\n%i::%i;(node,face)=(%i,%i); disti=%e; dist0=%e",iptr,nnzp2f,i,j,disti,dist0);
  //   if(disti<dist0) continue;
  //   j0=j;
  //   dist0=disti;
  // }
  //
  //
  // dist[0]=0.;
  // for(i=0;i<dim;i++) dist[0]+=sn0[i]*sn[i];
  // if(dist[0]>0.)
  //   return 0;// the point is above the face.
  // else
    return 1;// the point is below the face.
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx*/
void visible(INT root, INT node, iCSRmat *stos, \
                iCSRmat *fv, SHORT *maskf, \
                INT dim, REAL fact,     \
                REAL *x, REAL **xmass, REAL **sn, void *wrk)
{
	/* node is the node whose visible faces give the horizon ridges */
  INT chk,nvis,lsize,i,j,jk,iaa,iab,ila,ilb,nlvl,face;
	INT *vis=NULL,maskf_value=(node+1);
	ila = 0 ; ilb = 1;
  maskf[root]=maskf_value;
	vis=(INT *)malloc(1*sizeof(INT));
	vis[0]=root;
  while(1) {
    nlvl = ilb;
    lsize = nlvl;
    for (i=ila;i<ilb;i++){
      face = vis[i];
      iaa = stos->IA[face];
      iab = stos->IA[face+1];
      for(jk=iaa;jk<iab;jk++) {
				j = stos->JA[jk];
        // maskf_value here should be the node number.
        chk=chk_above(dim,fact,node,fv->JA+fv->IA[j], \
                      x, xmass[j], sn[j], wrk);
				if((maskf[j] != maskf_value)&&chk) {
					vis=realloc(vis,(nlvl+1)*sizeof(INT));
	  			vis[nlvl] = j;
					maskf[j] = maskf_value;
          fprintf(stdout,"\n%d: %d <= %d",j,nlvl,lsize);
	  			nlvl++;
				} //if
      }  //for
    }   //for
    ila = ilb;
    ilb = nlvl;
    if(nlvl <= lsize) break;
  } //while
  nvis = lsize;
	print_full_mat_int(1,nvis,vis,"vis");
  return;
}
