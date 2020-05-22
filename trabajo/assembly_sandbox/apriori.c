/**************************************************************************/
#include "hazmath.h"
#include "biot_headers.h"
/*************************************************************************/ 
REAL err_simplex(simplex_data *splex,				\
		 local_vec **ueall,				\
		 INT nq1d,					\
		 qcoordinates *cqelm,				\
		 mesh_struct *mesh,REAL time,INT elm,		\
		 void *wrk, void *wrkt)
{
  // compute the energy norm of the error on a simplex. 
  // Loop indices
  INT i,j,ij,quad,face,dim=splex->dim;
  REAL sigmaij,beta=0.,lam,mu,alpha,eei;
  REAL sol[15],dsol[100],zsigma[100]; // 7 dim: u(dim),w(dim),p(1),  du(dim*dim) , dw(1),dp(dim)
  REAL ea,ec,ee,esigma,ew,ep,integral = 0.0;
  //
  local_vec *ue=ueall[0], *uet=ueall[1];// time dependence
  //////////////////////////////////////////////////////////
  REAL *val=(REAL *)wrk; // no time derivatives yet
  REAL *valbub  = val;
  REAL *vallin  = valbub  + dim;
  REAL *valw    = vallin  + dim;//2*dim
  REAL *valp    = valw    + dim;
  REAL *dval    = valp    + 1; // last one for p. total=4*dim+1
  REAL *sigmab  = dval    + dim*dim + dim*dim + 1 + dim; // dim*dim
  // for b
  // and dim
  // for
  // u[k],k=1:dim;
  REAL *sigmau  = sigmab  + dim*dim; // dim ^2 for sigma(b)
  REAL *depsb   = sigmau  + dim*dim; // dim^2 for sigma(u);
  REAL *dtrepsb = depsb   + dim; // div of a matrix=vector
  REAL *conduct = dtrepsb + dim; // 
  REAL *wrkend  = conduct + dim*dim;
  
  //////////////Repetition for the minus element //
  REAL *valt=(REAL *)wrkt; //minus element: dim for b; dim for u, dim for w; 1 for p0.
  REAL *valbubt  = valt;
  REAL *vallint  = valbubt  + dim;
  REAL *valwt    = vallint  + dim;//2*dim
  REAL *valpt    = valwt    + dim;
  REAL *dvalt    = valpt    + 1; // last one for p. total=4*dim+1
  REAL *sigmabt  = dvalt    + dim*dim + dim*dim + 1 + dim; // dim*dim
							   // for b
							   // and dim
							   // for
							   // u[k],k=1:dim;
  REAL *sigmaut  = sigmabt  + dim*dim; // dim ^2 for sigma(b)
  REAL *depsbt   = sigmaut  + dim*dim; // dim^2 for sigma(u);
  REAL *dtrepsbt = depsbt   + dim; // div of a matrix=vector
  REAL *conductt = dtrepsbt + dim; // 
  REAL *wrkendt  = conductt + dim*dim;
  // some other pointers which are fit in the ones above
  REAL *dbub=dval;
  REAL *dlin=dbub + dim*dim; // address of the gradient of the linear part. 
  REAL *divw=dlin + dim*dim; // address of div(w)
  REAL *gradp=divw+1; // address of grad[p]
  // time dept
  REAL *dbubt=dvalt;
  REAL *dlint=dbubt  + dim*dim; // address of the gradient of the linear part. 
  REAL *divwt=dlint  + dim*dim; // address of div(w)
  REAL *gradpt=divwt + 1; // address of grad[p]
  //
  /* REAL *pt=valt+dim+dim+dim; // address of pt */
  /* REAL *dlineart=dvalt + dim*dim; // for the time dependent: address of the linear part. */
  //
  //////////////////////////////////////////////  //
  REAL *wptr=NULL, *wptrt=NULL;
  // Quadrature Weights and Nodes
  REAL w;// quad weights
  REAL *qx = (REAL *) calloc(dim,sizeof(REAL));
  // Quadrature on elm
  INT *v_on_elm=splex->v;
  INT *dof_on_elm=splex->dofs;
  REAL divu=0.,divb=0.,divut=0.,divbt=0.,divsol=0.;
  ea=0.;
  ec=0.;
  ee=0.;
  zquad_elm(cqelm,splex,nq1d);
  for (quad=0;quad<cqelm->nq_per_elm;quad++) {
    // elements
    qx[0] = cqelm->x[quad];
    qx[1] = cqelm->y[quad];
    if(dim==3) qx[2] = cqelm->z[quad];
    //////////////////////////////////////////////////////////////////////////
    conductivity2D(conduct,qx,time,NULL); // Actually this gives the
					  // inverse of K
    get_mu(&mu,qx,time,NULL);
    get_lam(&lam,qx,time,NULL);
    get_alpha(&alpha,qx,time,NULL);
    //////////////////////////////////////////////////////////////////////////
    zblockfe_interp(val,ue->b,qx,splex);
    zblockfe_interp(valt,uet->b,qx,splex);
    zblockfe_dinterp(dval,ue->b,qx,splex);
    zblockfe_dinterp(dvalt,uet->b,qx,splex);
    ////////////////////////////////////////////////////////////
    //    fprintf(stdout, "\n(x1,x2)=(%e,%e),t=%e",qx[0],qx[1],time);
    true_sol2D(sol,qx,time,NULL);
    Dtrue_sol2D(dsol,qx,time,NULL);
    /* fprintf(stdout,"\nelm=%d; (time,x)=(%.5e,%.5e,%.5e): ",elm,time,qx[0],qx[1]); */
    /* for(i=0;i<dim;i++){ */
    /*   fprintf(stdout,"%.7e ",val[i]); */
    /* } */
    /* fprintf(stdout,"\n");		        */
    /* print_full_mat(dim,dim,dsol,"dsol"); */
    /* print_full_mat(dim,dim,dlin,"dlin"); */
    /* fprintf(stdout,"\n========="); */
    ////////////////////////////////////////////////////////////
    // all vectors are ready, symmetrize the gradients of bubble and u1,u2
    divu=0.,divb=0.,divut=0.,divbt=0.,divsol=0.;
    // this below computes div u and div ub and div ut and div ubt
    /* for(j=0;j<dim*dim;j++) */
    /*   fprintf(stdout,"\n222222duet(%d)=%.7e",j,dbubt[j]); */
    for(i=0;i<dim;i++){
      // trace(sigma(u))
      divu+=dlin[i*dim+i];
      // trace(sigma(b))
      divb+=dbub[i*dim+i];
      // time derivatives
      //      divut+=dlint[i*dim+i];
      //      divbt+=dbubt[i*dim+i];
    }
    for(i=0;i<dim;i++){
      divsol+=dsol[i*dim+i];// u1x+u1y.
    }
    memset(sigmau,0,dim*dim*sizeof(REAL));
    memset(sigmab,0,dim*dim*sizeof(REAL));
    memset(zsigma,0,dim*dim*sizeof(REAL));
    for(i=0;i<dim;i++){
      for(j=0;j<dim;j++){
	sigmau[i*dim+j]=2.*mu*0.5*(dlin[i*dim+j]+dlin[j*dim+i]);
	sigmab[i*dim+j]=2.*mu*0.5*(dbub[i*dim+j]+dbub[j*dim+i]);
	zsigma[i*dim+j]=2.*mu*0.5*(dsol[i*dim+j]+dsol[j*dim+i]);
      }
    }
    for(i=0;i<dim;i++){
      sigmau[i*dim+i]+=lam*divu;
      sigmab[i*dim+i]+=lam*divb;
      zsigma[i*dim+i]+=lam*divsol;
    }      
    // weight
    w = cqelm->w[quad];
    // bubble part:
    esigma=0.;
    for(i=0;i<dim;i++){
      for(j=0;j<dim;j++){
	sigmaij=zsigma[i*dim+j] - sigmau[i*dim+j] - sigmab[i*dim+j]; 
	esigma+=sigmaij*sigmaij;
      }
    }
    //    dsol[0:dim-1,0:dim-1] is grad u
    //    dsol[dim*dim] is div w
    //    edivw=(dsol[4] - divw[0])*(dsol[4] - divw[0]);
    //    fprintf(stdout, "\nrhsi=%.6e,%.5e,%.5e,%.5e,%.6e",w,rhs[4],divut,divbt,divw[0]); 
    //    fprintf(stdout, "\nrhsi=%.5e,%.6e",w,rhsi); 
    ew=0.;
    //w=sol[dim],...sol[dim+dim-1]
    for(i=0;i<dim;i++){
      //      fprintf(stdout, "\nsolw,valw=%.6e,%.6e",sol[dim+i],valw[i]);
      for(j=0;j<dim;j++){
	ij=i*dim+j;
	ew+=conduct[ij]*(sol[dim+j]-valw[j])*(sol[dim+i]-valw[i]);
      }
    }
    //
    //p=sol[dim+dim]
    ep=(sol[dim+dim]-valp[0])*(sol[dim+dim]-valp[0]);
    ///////////////////////////////////////
    ea += w*esigma;
    ee += w*ew;
    ec += w*ep;
  }
  integral = ea + ee + ec;
  /* fprintf(stdout,"\nelm=%d; xintegral=%e",elm,e2i);//xxxxxintegral); */
  //  fprintf(stdout,"\ncqelm->nq_per_elm=%d",cqelm->nq_per_elm); 
  if(qx) free(qx);
  return integral;
}
/****************************************************************************/
REAL *apriori_err(dvector *uh0,
		  dvector *uh1,						\
		  block_fespace *FE,					\
		  mesh_struct *mesh,					\
		  REAL time,						\
		  REAL dt) 
{
  // uh0 is the solution from the previous time step;
  // uh1 is the current solution;
  INT dim = mesh->dim;
  INT nq1d=5;
  INT i,j,k,l,jk,row;
  // allocate a mask array to contain indicator for which elements the
  // bulk estimator and face estimators
  REAL *errelm=(REAL *)calloc(mesh->nelm,sizeof(REAL));  
  INT elmp;
  simplex_data *sp=simplex_data_init(mesh->dim,FE,1);
  // uep is for elmp and uem is for elmm, so we allocate both.
  local_vec **uep = dof_data_init(sp);
  qcoordinates *cqelm = allocateqcoords(nq1d,1,dim);
  INT memwrk=16*(dim*dim + dim + 1);// needs for each of these are 12*dim+5*dim*dim+3
  void *wrkp = (void *)calloc(2*memwrk,sizeof(REAL));
  void *wrkm = wrkp + memwrk*sizeof(REAL);
  //
  // end of data_face piece.  
  for(elmp=0;elmp<mesh->nelm;elmp++) {
    simplex_data_update(elmp,sp,mesh,FE,1);
    dof_data_update(uep,elmp,sp,uh0,uh1,dt);
    errelm[elmp]=err_simplex(sp,uep,				\
			     nq1d,cqelm,			\
			     mesh,time,elmp,			\
			     wrkp,wrkm);
    //    fprintf(stdout,"\nelmp:%d: %e",elmp,errelm[elmp]);
  }
  free_qcoords(cqelm);
  free(cqelm);
  //
  simplex_data_free(sp);
  dof_data_free(uep);
  free(wrkp);// wrkm is part of this;  
  return errelm;
}
