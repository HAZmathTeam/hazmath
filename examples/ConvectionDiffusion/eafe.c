
#include "hazmath.h"

static REAL bernoulli(const REAL z)
{
  // returns B(z) = z/(exp(z)-1)
  double tolb=1e-12,zlarge=256e+0;  
  if (fabs(z) < tolb)
    return (1.-z*0.5); // around 0 this is the asymptotic;
  else if(z<zlarge)
    return (z/(exp(z)-1.));
  else //z>tlarge this is zero pretty much
    return 0.;
}

void LumpMassBndry(const trimesh mesh, const INT bndry_marker,
		   void (*vector_val)(REAL *, REAL *, REAL),	\
		   void (*scalar_val)(REAL *, REAL *, REAL),	\
		   dvector *dmass,  dvector *rhs) 
{
  //calculates boundary integral using lumped mass.
  /* For an equation in divergence form, the natural boundary condition is
   \grad u . n + (b . n)u=g; 
   and also for 
   \grad u . n = g; 
   which is a Robin type condition when the equation is in divergence form.  
  */
  if(bndry_marker > 0 && bndry_marker < 3) return
  INT i=-1, j=-1,k=-1,ifa=-1,ifb=-1,attri=-16;
  REAL sixth=1./6.,vol=1e10,bdotn=1e10,rhsi=1e10;
  /* mesh entities: number of faces, number boundary faces and so on */
  INT  nf=mesh.nface,nfb=mesh.nbface,dim=mesh.dim; 
  REAL ad[dim],xm[dim];
  iCSRmat *f2v=mesh.f_v; /* face to vertex map as iCSR */
  INT *fonb=mesh.f_bdry; /* boundary markers for every face */
  /* get areas, normal vectors and barycenters  of faces */
  REAL *fa=mesh.f_area, *fn=mesh.f_norm, *fm=mesh.f_mid;
  fprintf(stdout,"Num Of Bdr Elements =%i\n",nfb);
  for (i = 0; i < nf; i++)   {
    attri=fonb[i];
    if(attri != bndry_marker)continue;    
    vol=sixth*sqrt(fa[i]);
    for(k=0;k<dim;k++){
      xm[k]=fm[i*dim+k];
    }
    vector_val(xm,ad,0.0);
    bdotn=0.;
    for(k=0;k<dim;k++) bdotn += ad[k]*fn[i*dim+k];
    bdotn*=vol;
    // bnormal is b*normal=b\cdot normal; 
    scalar_val(&rhsi,xm,0.0);
    rhsi *= vol;
    for(k=ifa;k<ifb;k++){
      j=f2v->JA[k];
      //      dmass.val[j]-= bdotn;
      rhs.val[j] += rhsi;
    }
  }
  return;
}

void eafe(const trimesh mesh,					\
	  void (*scalar_val)(REAL *, REAL *, REAL),		\
	  void (*vector_val)(REAL *, REAL *, REAL),		\
	  dCSRmat A, dvector rhs)
{
  INT i,j,jk,iaa,iab,nv=mesh.nv,dim=mesh.dim;
  INT *ia=A.IA, *ja=A.JA ;
  REAL *a=A.val;
  coordinates *xyz=mesh.cv;
  dvector diag0 = dvec_create(nv);
  dvector dmass = dvec_create(nv);
  for (i=0;i<nv;i++) {dmass.val[i]=diag0.val[i]=0.;}
  REAL ad[dim], xm[dim],te[dim];
  REAL xi,yi,zi,xj,yj,zj,bte,alpe;
  for (i = 0; i < nv; i++) {
    xi=xyz->x[i];
    yi=xyz->y[i];
    zi=xyz->z[i];
    iaa=ia[i]-1;
    iab=ia[i+1]-1;
    for (jk=iaa; jk<iab; jk++){
      j=ja[jk]-1; 
      if(i != j){
	xj=xyz->x[j];
	yj=xyz->y[j];
	zj=xyz->z[j];
	// compute the advection field at the middle of the edge and
	// then the bernoulli function
	te[0]  = xi - xj;
	te[1]  = yi - yj;
	xm[0] = (xi + xj)*0.5e+0;
	xm[1] = (yi + yj)*0.5e+0;
	if(dim>2){
	  te[2]  = zi - zj;
	  xm[2] = (zi + zj)*0.5e+0;
	}
	vector_val(ad,xm,0.0);
	bte = ad[0]*te[0]+ ad[1]*te[1]+ ad[2]*te[2];
	scalar_val(&alpe,xm,0.0);
	//	  xmid.Print(std::cout,3);
	//	  std::cout << "alpe = "<<alpe<<"; bte="<< bte<<std::endl<<std::flush;      
	// alpe=a(xmid)\approx harmonic_average=|e|/(int_e 1/a);
	// should be computed by quadrature in general for 1/a(x).
	// a_{ij}=B(beta.t_e/harmonic_average)*harmonic_average*omega_e;
	//	  for (i,j):    B(beta.t_e/harmonic_average);   
	//	  for (j,i):    B(-beta.t_e/harmonic_average), the minus comes from
	//     is because t_e=i-j;
	a[jk] *= (alpe*(bernoulli(bte/alpe))); 
	// the diagonal is equal to the negative column sum;
	diag0.val[j]-=a[jk];
      }
    } 
  }
  // Another loop to set up the diagonal equal to whatever it needs to equal. 
  for (i = 0; i < nv; i++) 
    for (jk=ia[i]; jk<ia[i+1]; jk++){
      j = ja[jk]; 
      if(i == j) a[jk]=diag0.val[i]+dmass.val[i];
    }
  return;
}
