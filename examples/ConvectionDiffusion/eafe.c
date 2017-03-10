
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

static void poisson_coeff(REAL *val,REAL* x, REAL t) {
  // a(x)
  *val = 1.0;
  return;
}

void SchurProduct(const trimesh mesh,					\
		  void (*diffusion)(REAL *, REAL *, REAL),		\
		  void (*advection_vector)(REAL *, REAL *, REAL),	\
		  dCSRmat A, dvector dmass)
{
  INT i,j,jk,nv=mesh.nv,dim=mesh.dim;
  INT *ia=A.IA, *ja=A.JA ;
  REAL *a=A.val;
  coordinates *xyz=mesh.cv;
  dvector diag0;
  diag0.row=xyz->n;
  diag0.val=(REAL *)calloc(nv,sizeof(REAL));
  for (i=0;i<nv;i++) diag0.val[i]=0.;
  REAL advcoeff[dim], xmid[dim],te[dim];
  REAL xi,yi,zi,xj,yj,zj,bte,bern,alpe;
  for (i = 0; i < nv; i++) {
    xi=xyz->x[i];
    yi=xyz->y[i];
    zi=xyz->z[i];
    for (jk=ia[i]; jk<ia[i+1]; jk++){
      j=ja[jk]; 
      if(i != j){
	xj=xyz->x[j];
	yj=xyz->y[j];
	zj=xyz->z[j];
	// compute the advection field at the middle of the edge and
	// then the bernoulli function
	te[0]  = xi - xj;
	te[1]  = yi - yj;
	te[2]  = zi - zj;
	xmid[0] = (xi + xj)*0.5e+0;
	xmid[1] = (yi + yj)*0.5e+0;
	xmid[2] = (zi + zj)*0.5e+0;
	advection_vector(advcoeff,xmid,0.0);
	//	  bte = (beta*t_e); //c++ overloaded "*" s..t we can multiply these...
	//	  I do it the old way, unrolled it is faster
	bte = advcoeff[0]*te[0]+ advcoeff[1]*te[1]+ advcoeff[2]*te[2];
	///
	// diffusion coefficient for the flux J = a(x)\nabla u + \beta u;
	diffusion(&alpe,xmid,0.0);
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
  // Another loop to set up the diagonal equal to the negative column sum;
  for (i = 0; i < nv; i++) 
    for (jk=ia[i]; jk<ia[i+1]; jk++){
      j = ja[jk]; 
      if(i == j) a[jk]=diag0.val[i]+dmass.val[i];
    }
  return;
}
/**/
void LumpMassBndry(const trimesh mesh, const INT bndry_marker,
		   dvector dmass, dvector rhs) 
{
  //calculates boundary integral using lumped mass.
  /* For an equation in divergence form, the natural boundary condition is
   \grad u . n + (b . n)u=g; 
  */
  INT i, j,k,nve,attri,node=-16;
  REAL sixth=1./6.,voltri=1e10,distz=0e0;
  /* mesh entities */
  INT ns = mesh.nelm,nv=mesh.nv,dim=mesh.dim; /* ns=number of simplices of highest dimension*/
  coordinates *xyz = mesh.cv;
  REAL bnormal[dim],xyzb[dim];
  INT *f2v=mesh.f_v; /* face to vertex map as iCSR */
  INT *fonb=mesh.f_bdry; /* boundary markers for every face */
  INT *fa=mesh.f_area; /* areas of faces */
  INT nf=mesh.nface,nfb=mesh.nbface;
  fprintf(stdout,"Num Of Bdr Elements =%i\n",nfb);
  if(dim==2) {
    nf=nedge;
    nfb=mesh.nbedge;
    f2v=mesh.ed_v;
    fonb=mesh.ed_bdry;
    fa=mesh.ed_len;
  }
  for (i = 0; i < nf; i++)   {
    attri=fonb[i];
    if(attri != bndry_marker)continue;
    voltri=fa[i];
    voltri=sixth*pow(voltri,0.5);
    xyzb[0]=xyz.x[i];
    xyzb[1]=xyz.y[i];
    if(dim>2)
      xyzb[2]=xyz.z[i];
    //    ifstrat = f2v->IA[i]-1;
    //    ifend =   f2v->IA[i+1]-1;
    // advection(xyzb,bnormal);
    // bnormal is b*normal=b\cdot (0,0,1); 
    //dmass[node[j]]-=voltri * bnormal[2];
    ifa=f2v->IA[i]-1;
    ifb=f2v->IA[i+1]-1;
    rhsi = rhs_neumann(xyzb);
    for(k=ifa;k<ifb;k++){
      j=f2v->JA[k];
      rhs[j] += voltri * rhsi;
    }
  }
  return;
}
