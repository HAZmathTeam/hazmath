/**************************************************************************/
#include "hazmath.h"
#include "biot_headers.h"
/*************************************************************************/ 
REAL estimator_simplex(simplex_data *splex,			\
		       local_vec **ueall,			\
		       INT nq1d,				\
		       qcoordinates *cqelm_pm,			\
		       mesh_struct *mesh,REAL time,INT elm,	\
		       void *wrk, void *wrkt)
{
  // Loop indices
  INT i,j,ij,quad,face,dim=splex->dim;
  REAL beta=0.,lam,mu,alpha,eei,rhsi,exact[10];
  REAL e1i,e2i,e3i,e1,e2,e3,e31,integral = 0.0;
  //
  local_vec *ue=ueall[0], *uet=ueall[1];// time dependence
  //////////////////////////////////////////////////////////
  // here is the same as for plus and minus elements except we have rhs and time dependence. 
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
  REAL *rhs     = conduct + dim*dim;
  REAL *wrkend  = rhs + dim+dim+dim+1; // one extra dim
  
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
  REAL *rhst     = conductt + dim*dim;
  REAL *wrkendt  = rhst + dim+dim+dim+1; // one extra dim for t too. 
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
  zquad_elm(cqelm_pm,splex,nq1d);
  /* for (quad=0;quad<cqelm->nq_per_elm;quad++) { */
  /*   fprintf(stdout,							\ */
  /* 	    "\ncq=(%.6e,%.6e),zcq=(%.6e,%.6e):diff=(%.12e,%.12e)" ,	\ */
  /* 	    cqelm->x[quad],cqelm->y[quad],				\ */
  /* 	    cqelm_pm->x[quad],cqelm_pm->y[quad],		       	\ */
  /* 	    cqelm->x[quad]-cqelm_pm->x[quad],				\ */
  /* 	    cqelm->y[quad]-cqelm_pm->y[quad]); */
  /*     } */
  //  fprintf(stdout,"\nb_dofs");
  INT dof_per_elm = 0,np,ndofsi=0;
  /* for(i=0;i<splex->nspaces;i++) { */
  /*   fprintf(stdout,"\nspace=%d; ",i); */
  /*   for(j=ndofsi;j<ndofsi+splex->ndofs[i];j++){ */
  /*     fprintf(stdout,"u_comp[%d]=%.6e; ",splex->dofs[j],ue->b[j]); */
  /*   } */
  /*   ndofsi+=splex->ndofs[i]; */
  /* } */
  REAL hk=pow(splex->vol,1./((REAL )dim));
  e1i=0.;
  e2i=0.;
  e3i=0.;
  for (quad=0;quad<cqelm_pm->nq_per_elm;quad++) {    
    // elements
    qx[0] = cqelm_pm->x[quad];
    qx[1] = cqelm_pm->y[quad];
    if(dim==3) qx[2] = cqelm_pm->z[quad];
    //////////////////////////////////////////////////////////////////////////
    conductivity2D(conduct,qx,time,NULL); // Actually this gives the
					  // inverse of K
    get_mu(&mu,qx,time,NULL);
    get_lam(&lam,qx,time,NULL);
    get_alpha(&alpha,qx,time,NULL);
    //////////////////////////////////////////////////////////////////////////
    zblockfe_interp(val,ue->b,qx,splex);
    zblockfe_interp(valt,uet->b,qx,splex);
    //////////////////////////////////////////////////////////////////////////
    /* wptr=val; */
    /* print_full_mat(1,2,wptr,"xb"); */
    /* wptr+=dim; */
    /* print_full_mat(1,1,wptr,"u1"); */
    /* wptr++; */
    /* print_full_mat(1,1,wptr,"u2"); */
    /* wptr++; */
    /* print_full_mat(1,2,wptr,"w"); */
    /* wptr+=dim; */
    /* print_full_mat(1,1,wptr,"p"); */
    /* fprintf(stdout,"\n**********************************************"); */
    /////////////////////////////////////////////////////////////////////////
    /* for(j=0;j<splex->num_dofs;j++) */
    /*   fprintf(stdout,"\nuet(%d)=%.7e",j,uet->b[j]); */
    /* for(j=2*dim*dim;j<dim*dim+dim*dim+1;j++) */
    /*   fprintf(stdout,"\ndue(%d)=%.7e",j,dval[j]); */
    zblockfe_dinterp(dval,ue->b,qx,splex);
    zblockfe_dinterp(dvalt,uet->b,qx,splex);
    /* for(j=0;j<dim*dim;j++) */
    /*   fprintf(stdout,"\nduet(%d)=%.7e",j,dvalt[j]); */
    /////////////////////////////////////////////////////////////////////////
    /* fprintf(stdout,"\nx1=%f,x2=%f",qx[0],qx[1]); */
    /* wptr=dvalt; */
    /* print_full_mat(dim,dim,wptr,"xgrad_bt"); */
    /* wptr+=dim*dim; */
    /* print_full_mat(1,dim,wptr,"xgradu1t"); */
    /* wptr+=dim; */
    /* print_full_mat(1,dim,wptr,"xgradu2t"); */
    /* wptr+=dim; */
    /* print_full_mat(1,1,dval,"div wt"); */
    /* wptr++; */
    /* print_full_mat(1,dim,wptr,"grad p"); */
    /* wptr+=dim; */
    //////////////////////////////////////////////////////////////////////////    
    /* print_full_mat(1,dim*dim,dvalt,"1dvalt"); */
    rhs_func(rhs,qx,time,NULL);
    /* print_full_mat(1,dim*dim,dvalt,"dvalt"); */
    //////////////////////////////////////////////////////////////////////////
    /* wptr=rhs; */
    /* print_full_mat(1,2,wptr,"rhsu"); */
    /* wptr+=dim; */
    /* print_full_mat(1,2,wptr,"rhsw"); */
    /* wptr++; */
    /* print_full_mat(1,1,wptr,"rhsp"); */
    /////////////////////////////////////////////////////////////////////////
    /* zzbubble_hessian2d(depsb,dtrepsb,			\ */
    /* 		       u,qx,v_on_elm,dof_on_elm,mesh); */
    /* print_full_mat(1,2,depsb,"depsb"); */
    /* print_full_mat(1,2,dtrepsb,"dtrepsb"); */
    zbubble_hessian2d(depsb,dtrepsb,		\
    		      ue->b,qx,splex);
    //    print_full_mat(1,2,depsb,"depsb1");
    //    print_full_mat(1,2,dtrepsb,"dtrepsb1");
    // compute exact solution at qx..
    //    true_sol2D(exact,qx,time,NULL);
    //    print_full_mat(1,5,exact,"exact");
    //compute the approximate solution at qx
    ///////////////////////////////////////
    // all vectors are ready, symmetrize the gradients of bubble and u1,u2
    REAL divu=0.,divb=0.,divut=0.,divbt=0.;
    // this below computes div u and div ub and div ut and div ubt
    /* fprintf(stdout,"\nQ=(%.15e,%.15e)",qx[0],qx[1]); */
    /* fprintf(stdout,"\nbubDOF=[%.15e,%.15e,%.15e]",uet->b[0],uet->b[1],uet->b[2]); */
    /* for(j=0;j<dim;j++){ */
    /*   fprintf(stdout,"\ndbubt(%d,1:%d)=[%.7e,%.7e]",j+1,dim,dbubt[j*dim+0],dbubt[j*dim+1]); */
    /* } */
    /* for(j=0;j<dim;j++){ */
    /*   fprintf(stdout,"\ndlint(%d,1:%d)=[%.7e,%.7e]",j+1,dim,dlint[0*dim+j],dlint[1*dim+j]); */
    /* } */
    for(i=0;i<dim;i++){
      // trace(sigma(u))
      divu+=dlin[i*dim+i];
      // trace(sigma(b))
      divb+=dbub[i*dim+i];
      // time derivatives
      divut+=dlint[i*dim+i];
      divbt+=dbubt[i*dim+i];
    }
    memset(sigmau,0,dim*dim*sizeof(REAL));
    memset(sigmab,0,dim*dim*sizeof(REAL));
    for(i=0;i<dim;i++){
      for(j=0;j<dim;j++){
	sigmau[i*dim+j]=2.*mu*0.5*(dlin[i*dim+j]+dlin[j*dim+i]);
	sigmab[i*dim+j]=2.*mu*0.5*(dbub[i*dim+j]+dbub[j*dim+i]);
      }
    }
    for(i=0;i<dim;i++){
      sigmau[i*dim+i]+=lam*divu;
      sigmab[i*dim+i]+=lam*divb;
    }      
    // weight
    w = cqelm_pm->w[quad];
    // bubble part:
    e1=0.;
    for(i=0;i<dim;i++){
      rhsi=rhs[i] + 2.*mu*depsb[i] + lam*dtrepsb[i] - alpha*gradp[i]; //f+div(sigma(b))-alpha \grad p
      e1+=rhsi*rhsi;
    }
    //    rhsi=rhs[4]-beta*valpt[0] - divut - divbt - divw[0];
    rhsi=rhs[4] + beta*valpt[0] +  alpha*(divut + divbt) + divw[0];
    //    fprintf(stdout, "\nelm:%d; rhsi=%7e=sum([%.5e,%.5e,%.6e, %.6e])",elm,rhsi,rhs[4],divut,divbt,divw[0]); 
    //    fprintf(stdout, "\nrhsi=%.5e,%.6e",w,rhsi); 
    e2=rhsi*rhsi;
    e3=0.;
    // K^{-1}*w+grad p
    for(i=0;i<dim;i++){
      e31=gradp[i];
      for(j=0;j<dim;j++){
	ij=i*dim+j;
	e31+=(conduct[ij]*valw[j]);	  
      }
      e3+=e31*e31;
    }
    ///////////////////////////////////////
    e1i += w*e1;
    e2i += w*e2;
    e3i += w*e3;
  }
  integral = e1i*hk*hk + e3i*hk*hk + e2i;
  /* fprintf(stdout,"\nelm=%d; xintegral=%e",elm,e2i);//xxxxxintegral); */
  //  fprintf(stdout,"\ncqelm_pm->nq_per_elm=%d",cqelm_pm->nq_per_elm); 
  if(qx) free(qx);
  return integral;
}
/*************************************************************************/ 
REAL estimator_face(simplex_data *sp,				\
		    local_vec **uepall,				\
		    simplex_data *sm,				\
		    local_vec **uemall,				\
		    REAL *xfi, REAL *finrm, REAL fiarea,	\
		    INT nq1d,					\
		    qcoordinates *cqface,			\
		    REAL time,INT elmp, INT elmm,		\
		    INT fflag, void *wrkp, void *wrkm)
{
  //fflag is the face flag
  //wrkp and wrkm must have enough space in them for
  // 10*dim+5 + 5*dim*dim + 2 more dim for the tangential vector and the f1 below!
  // Loop indices
  INT i,j,ij,quad,face,dim=sp->dim;
  REAL lam,mu,alpha,eei,divb,divu,divut,divbt,rhsi,pjump,exact[10];
  REAL integral = 0.0,e1=0.,e3=0.,e31;
  //
  local_vec *uep=uepall[0], *uept=uepall[1], *uem=uemall[0],*uemt=uemall[1];
  //////////////////////////////////////////////////////////
  REAL *valp=(REAL *)wrkp; //plus element: dim for b; dim for u, dim for w; 1 for p0.
  REAL *valbubp  = valp;
  REAL *vallinp  = valbubp  + dim;
  REAL *valwp    = vallinp  + dim;//2*dim
  REAL *valpp    = valwp    + dim;
  REAL *dvalp    = valpp    + 1; // last one for p. total=4*dim+1
  REAL *sigmabp  = dvalp    + dim*dim + dim*dim + 1 + dim; // dim*dim
							   // for b
							   // and dim
							   // for
							   // u[k],k=1:dim;
  REAL *sigmaup  = sigmabp  + dim*dim; // dim ^2 for sigma(b)
  REAL *depsbp   = sigmaup  + dim*dim; // dim^2 for sigma(u);
  REAL *dtrepsbp = depsbp   + dim; // div of a matrix=vector
  REAL *conductp = dtrepsbp + dim; // 
  REAL *wrkendp  = conductp + dim*dim;
  
  //////////////Repetition for the minus element //
  REAL *valm=(REAL *)wrkm; //minus element: dim for b; dim for u, dim for w; 1 for p0.
  REAL *valbubm  = valm;
  REAL *vallinm  = valbubm  + dim;
  REAL *valwm    = vallinm  + dim;//2*dim
  REAL *valpm    = valwm    + dim;
  REAL *dvalm    = valpm    + 1; // last one for p. total=4*dim+1
  REAL *sigmabm  = dvalm    + dim*dim + dim*dim + 1 + dim; // dim*dim
							   // for b
							   // and dim
							   // for
							   // u[k],k=1:dim;
  REAL *sigmaum  = sigmabm  + dim*dim; // dim ^2 for sigma(b)
  REAL *depsbm   = sigmaum  + dim*dim; // dim^2 for sigma(u);
  REAL *dtrepsbm = depsbm   + dim; // div of a matrix=vector
  REAL *conductm = dtrepsbm + dim; // 
  REAL *wrkendm  = conductm + dim*dim; // 
  // some other pointers:
  //plus
  REAL *dbubp=dvalp;
  REAL *dlinp=dbubp + dim*dim; // address of the gradient of the linear part. 
  REAL *divwp=dlinp + dim + dim; // address of div(w)
  REAL *gradpp=divwp+1; // address of grad[p]
  // minus
  REAL *dbubm=dvalm;
  REAL *dlinm=dbubm + dim*dim; // address of the gradient of the linear part. 
  REAL *divwm=dlinm + dim + dim; // address of div(w)
  REAL *gradpm=divwm+1; // address of grad[p]
  //
  /* REAL *pt=valt+dim+dim+dim; // address of pt */
  /* REAL *dlineart=dvalt + dim*dim; // for the time dependent: address of the linear part. */
  //
  // Quadrature Weights and Nodes
  REAL w;// quad weights
  REAL *qx = (REAL *) calloc(dim,sizeof(REAL));
  //  2donly: turn in positive direction. 
  REAL *fitau    = wrkendp; // 
  REAL *f1       = fitau    + dim; // 
  wrkendp  = f1       + dim;
  fitau[0] = -finrm[1];
  fitau[1] =  finrm[0];
  // Quadrature on face;
  zquad_face(cqface,nq1d,dim,xfi,fiarea);
  REAL hf=pow(fiarea,1./((REAL )(dim-1)));
  e1=0.;
  e3=0.;
  for (quad=0;quad<cqface->nq_per_elm;quad++) {
    // elements
    qx[0] = cqface->x[quad];
    qx[1] = cqface->y[quad];
    if(dim==3) qx[2] = cqface->z[quad]; 
    //////////////////////////////////////////////////////////////////////////
    conductivity2D(conductp,sp->mid,time,NULL); // Actually this gives the
    conductivity2D(conductm,sm->mid,time,NULL); // inverse of K at the centers of both elements
    get_mu(&mu,qx,time,NULL);
    get_lam(&lam,qx,time,NULL);
    get_alpha(&alpha,qx,time,NULL);
    //////////////////////////////////////////////////////////////////////////
    zblockfe_interp(valp,uep->b,qx,sp);
    zblockfe_interp(valm,uem->b,qx,sm);
    /* if(elmm>=0){ */
    /*   print_full_mat(1,dim,valbubp,"valp"); */
    /*   print_full_mat(1,dim,valbubm,"valm"); */
    /* } */
    //////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    zblockfe_dinterp(dvalp,uep->b,qx,sp);
    /* if(elmm>=0){ */
    /*   print_full_mat(1,dim*dim,dbubp,"dvalp"); */
    /* } */
    zblockfe_dinterp(dvalm,uem->b,qx,sm);
    /* if(elmm>=0){ */
    /*   print_full_mat(1,dim*dim,dbubm,"dvalm"); */
    /* } */
    // evaluations are done. 
    REAL divup=0.,divum=0.,divbp=0.,divbm=0.;
    // this below computes div u and div ub and div ut and div ubt
    for(i=0;i<dim;i++){
      // trace(sigma(u))
      divup+=dlinp[i*dim+i];
      divum+=dlinm[i*dim+i];
      // trace(sigma(b))
      divbp+=dbubp[i*dim+i];
      divbm+=dbubm[i*dim+i];
    }
    memset(sigmaup,0,dim*dim*sizeof(REAL));
    memset(sigmaum,0,dim*dim*sizeof(REAL));
    memset(sigmabp,0,dim*dim*sizeof(REAL));
    memset(sigmabm,0,dim*dim*sizeof(REAL));
    for(i=0;i<dim;i++){
      for(j=0;j<dim;j++){
	sigmaup[i*dim+j]=2.*mu*0.5*(dlinp[i*dim+j]+dlinp[j*dim+i]);
	sigmaum[i*dim+j]=2.*mu*0.5*(dlinm[i*dim+j]+dlinm[j*dim+i]);
	sigmabp[i*dim+j]=2.*mu*0.5*(dbubp[i*dim+j]+dbubp[j*dim+i]);
	sigmabm[i*dim+j]=2.*mu*0.5*(dbubm[i*dim+j]+dbubm[j*dim+i]);
      }
    }
    for(i=0;i<dim;i++){
      sigmaup[i*dim+i]+=lam*divup;
      sigmaum[i*dim+i]+=lam*divum;
      sigmabp[i*dim+i]+=lam*divbp;
      sigmabm[i*dim+i]+=lam*divbm;
    }      
    //    print_full_mat(dim,dim,sigmabp,"SP");
    //    print_full_mat(dim,dim,sigmabm,"SM");
    //    fprintf(stdout,"\nXXXXXYYYYYelm=(%d,%d)",elmp,elmm);
    // compute the integral here....
    e1=0.;
    w = cqface->w[quad];
    pjump=valpp[0]-valpm[0];
    for(i=0;i<dim;i++){
      f1[i]=-finrm[i]*alpha*pjump;// [alpha p] * n_F
      for(j=0;j<dim;j++){
	ij=i*dim+j;
	f1[i]+=(sigmaup[ij] - sigmaum[ij] + sigmabp[ij]-sigmabm[ij])*finrm[j];
      }
      e1+=f1[i]*f1[i];
    }
    // K^{-1}*w*tf
    e31=0.;
    for(i=0;i<dim;i++){
      for(j=0;j<dim;j++){
	ij=i*dim+j;
	e31+=(conductp[ij]*valwp[j]-conductm[ij]*valwm[j])*fitau[i];	  
      }
    }
    e3=(e31*e31 + pjump*pjump);
    integral+=w*e1*hf + w*e3*hf;
  }
  //  fprintf(stdout,"\ncqface->nq_per_elm=%d",cqface->nq_per_elm); 
  /* print_full_mat(1,dim,fitau,"fitau"); */
  /* print_full_mat(1,dim,finrm,"finrm"); */
  /* fprintf(stdout,"\nelm=(%d,%d); integral=%e",elmp,elmm,integral); */
  /* fprintf(stdout,"\n---------------------------\n"); */
  /* fprintf(stdout,"\nelms=(%d,%d); face flag=%d",elmp,elmm,fflag); */
  /* print_full_mat(dim+1,dim,sp->f_norm,"normalsp"); */
  /* print_full_mat(dim+1,dim,sp->f_norm,"normalsm"); */
  //fprintf(stdout,"\nhf=%.9e",hf);
  if(qx) free(qx);
  if(elmm<0 && fflag==1)
    return 0.0;
  return integral;
}
/****************************************************************************/
void estimator_block(estimator *est,
		     dvector *uh0,					\
		     dvector *uh1,					\
		     block_fespace *FE,					\
		     mesh_struct *mesh,					\
		     void (*f)(REAL *,REAL *,REAL,void *),		\
		     REAL time,						\
		     REAL dt) 
{
  // uh0 is the solution from the previous time step;
  // uh1 is the current solution;
  INT dim = mesh->dim;
  INT nq1d=5;
  INT i,j,k,l,jk,row;
  // allocate a mask array to contain indicator for which elements the
  INT *mask=(INT *)calloc(mesh->nelm,sizeof(INT));
  // bulk estimator and face estimators
  // Need face to element map
  iCSRmat f_el;
  icsr_trans(mesh->el_f,&f_el); // we need this one face,element
  //  icsr_print_matlab(stdout,&f_el);
  INT elmp,elmm; // simplices
  simplex_data *sp=simplex_data_init(mesh->dim,FE,1);
  simplex_data *sm=simplex_data_init(mesh->dim,FE,1);
  // uep is for elmp and uem is for elmm, so we allocate both.
  local_vec **uep = dof_data_init(sp);
  local_vec **uem = dof_data_init(sm); 
  qcoordinates *cqelm_pm = allocateqcoords(nq1d,1,dim);
  qcoordinates *cqface = allocateqcoords_bdry(nq1d,1,dim,2); // 1 is 1
							     // region;
							     // 2 is a
							     // face
							     // (in 2d
							     // this
							     // is an
							     // edge)
  REAL *data_face=calloc(dim+dim*dim, sizeof(REAL)); //coords of normal and vertices of a face. 
  REAL *xfi=data_face; //coords of vertices on face i. 
  REAL *finrm=xfi+dim*dim; //coords of normal on face i
  REAL *data_face_end=finrm+dim; //end
  INT memwrk=16*(dim*dim + dim + 1);// needs for each of these are 12*dim+5*dim*dim+3
  void *wrkp = (void *)calloc(2*memwrk,sizeof(REAL));
  void *wrkm = wrkp + memwrk*sizeof(REAL);
  //
  // end of data_face piece.
  memset(mask,0,mesh->nelm);
  for(i=0;i<mesh->nface;i++) {
    for(jk=mesh->f_v->IA[i];jk<mesh->f_v->IA[i+1];jk++){
      j=jk-mesh->f_v->IA[i];
      k=mesh->f_v->JA[jk];
      xfi[j*dim+0]=mesh->cv->x[k];
      if(dim>1)
	xfi[j*dim+1]=mesh->cv->y[k];
      if(dim>2)
	xfi[j*dim+2]=mesh->cv->z[k];	  
    }
    finrm[0]=mesh->f_norm[i*dim+0];
    if(dim>1)
      finrm[1]=mesh->f_norm[i*dim+1];
    if(dim>2)
      finrm[2]=mesh->f_norm[i*dim+2];
    /* fprintf(stdout,"\nface=%d",i); */
    /* print_full_mat(1,dim,finrm,"finrm"); */
    /* print_full_mat(dim,dim,xf,"xf"); */
    /* for (jk=0;jk<cqface->nq_per_elm;jk++) { */
    /*   fprintf(stdout,							\ */
    /* 	      "\ncq=(%.6e,%.6e),zcq=(%.6e,%.6e):diff=(%.12e,%.12e)" ,	\ */
    /* 	      cqface->x[jk],cqface->y[jk],				\ */
    /* 	      zcqface->x[jk],zcqface->y[jk],				\ */
    /* 	      cqface->x[jk]-zcqface->x[jk],				\ */
    /* 	      cqface->y[jk]-zcqface->y[jk]); */
    /* } */
    // find the elements that intersect at this face:
    elmp = f_el.JA[f_el.IA[i]];
    if((f_el.IA[i+1]-f_el.IA[i])>=2)
      elmm = f_el.JA[f_el.IA[i+1]-1];      
    else
      elmm=-1;
    //tplus
    //    fprintf(stdout, "\nelmp:%d====BEGIN======",elmp);
    simplex_data_update(elmp,sp,mesh,FE,1);
    dof_data_update(uep,elmp,sp,uh0,uh1,dt);
    //    simplex_data_print(elmp, sp,uep);      
    if(elmm>=0){
      simplex_data_update(elmm,sm,mesh,FE,1);
      dof_data_update(uem,elmm,sm,uh0,uh1,dt);
    } else {
      memset(uem[0]->b,0,sm->num_dofs);
      memset(uem[1]->b,0,sm->num_dofs);
    }
    est->face[i]=estimator_face(sp,uep,				\
				sm,uem,				\
				xfi,finrm,mesh->f_area[i],	\
				nq1d,cqface,			\
				time,				\
				elmp,elmm,			\
				mesh->f_flag[i],wrkp,wrkm);
    if(!mask[elmp]){
      est->bulk[elmp]=estimator_simplex(sp,uep,			\
					nq1d,cqelm_pm,		\
					mesh,time,elmp,		\
					wrkp,wrkm);
      mask[elmp]=1;
      //      fprintf(stdout, "\nelmp:%d====END======",elmp);
    }
    if(elmm>=0){
      if(!mask[elmm]){
	//	fprintf(stdout, "\nelmm:%d====BEGIN======",elmm);
	simplex_data_update(elmm,sm,mesh,FE,1);
	dof_data_update(uem,elmm,sm,uh0,uh1,dt);
	//	simplex_data_print(elmm, sm, uem);
	est->bulk[elmm]=estimator_simplex(sm,uem,			\
					  nq1d,cqelm_pm,		\
					  mesh,time,elmm,		\
					  wrkp,wrkm);
	//	fprintf(stdout, "\nelmp:%d====END======",elmm);
	mask[elmm]=1;
      }// mask>0
    }//elmm>=0
  }
  free_qcoords(cqelm_pm);
  free(cqelm_pm);
  free_qcoords(cqface);
  free(cqface);
  //
  simplex_data_free(sp);
  dof_data_free(uep);
  simplex_data_free(sm);
  dof_data_free(uem);
  free(data_face);// wrkm is part of this;  
  free(wrkp);// wrkm is part of this;  
  return;
}
