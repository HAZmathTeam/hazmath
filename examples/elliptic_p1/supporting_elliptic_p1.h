/*! \file examples/elliptic_p1/elliptic_p1_supporting.h
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20190420.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note includes assembling, fe solution, plotting (projection for
 *        3d and 4d from the whole domain to part of the boundary)
 *
 *  \note: modified by ltz on 20210602
 *
 */
/**********************************************************************************/
/*!
 * \fn static void mapit(scomplex *sc,const REAL *vc)
 *
 * \brief Constructs a simplicial mesh in a polyhedral domain O in
 *        dim-dimensions (dim=sc->n). The domain is assumed to be
 *        isomorphic to the cube in dim-dimensions. Examples: it is a
 *        quadrilateral when d=2 and hexagonal when d=3. To avoid
 *        ambiguity, we order the vertices vc[] of O lexicographicaly
 *        to get vcp[] so that the j-th vertex in the ordered array,
 *        with coordinates vcp[j*dim--(j+1)*dim-1] is mapped to the
 *        vertex of the unit cube whose coordinates are the digits in
 *        the binary representation of j. Here j=[0,...,2^(dim)-1]. 
 *        
 *
 * \param sc    I: simplicial complex defining the FE grid.
 *
 * \param vc[] I:    A REAL array with coordinates of the vertices of the
 *                   domain. The vertex k is with coordinates
 *                   vc[k*dim--(k+1)*dim-1]. These could be given in
 *                   any order.
 *
 * \note Ludmil (20210807)
 */
/**********************************************************************************/
static void mapit(scomplex *sc,REAL *vc)
{
  if(!vc) return;
  /* maps a mesh on the unit cube in d-dimensions to a domain with vertices vc[] */
  INT dim=sc->n;
  INT i,j,kf,dim1=dim+1;
  cube2simp *c2s=cube2simplex(dim);
  REAL *vcp_xhat = (REAL *)calloc(dim*(c2s->nvcube+1),sizeof(REAL));
  REAL *vcp = vcp_xhat;
  REAL *xhat = vcp + c2s->nvcube*dim;
  // work with a copy:
  memcpy(vcp,vc,(dim*c2s->nvcube)*sizeof(REAL));
  INT *p = (INT *)calloc(c2s->nvcube,sizeof(INT));
  /*order vcp lexicographically because it will be mapped to the unit
    cube which has lexicographically ordered vertices.*/
  dlexsort(c2s->nvcube,dim,vcp,p);
  /* the permutation of vertices is not needed, so we free it */
  free(p);
  /* transpose vcp to get one vertex per column as we first transform x,
     then y then z and so on */
  r2c(c2s->nvcube,dim,sizeof(REAL),vcp);
  for(kf=0;kf<sc->nv;kf++){
    for(i=0;i<dim;i++) xhat[i]=sc->x[kf*dim+i];
    for(i=0;i<dim;i++) sc->x[kf*dim+i]=interp4(c2s,vcp+i*c2s->nvcube,xhat);
  }
  cube2simp_free(c2s);
  free(vcp_xhat);  
  return;
}
/**********************************************************************************/
/*!
 * \fn static dvector fe_sol(scomplex *sc,const REAL alpha,const REAL gamma)
 *
 * \brief assembles and solves the finite element discretization of a
 *        reaction diffusion equation discretized with P1 continuous
 *        elements in any spatial dimension > 1. The right hand side
 *        is f(x)=1.
 *
 * \param sc    I: simplicial complex defining the FE grid.
 *
 * \param alpha I:  the diffusion coefficient
 *
 * \param gamma I: the reaction coefficient
 *
 * \return The FE solution as dvector.
 *
 * \note Ludmil (20210807)
 */
/**********************************************************************************/
static dvector fe_sol(scomplex *sc,const REAL alpha,const REAL gamma)
{
  REAL fi;
  INT i,ij,j,solver_flag=-10,print_level=0;
  dCSRmat A,M;
  clock_t clk_assembly_start = clock(); // begin assembly timing;
  assemble_p1(sc,&A,&M);
  // Now we solve the Neumann problem (A+M) u = f;
  for(i=0;i<A.nnz;++i)
    A.val[i]=alpha*A.val[i] + gamma*M.val[i];
  /*We have Neumann problem, so no boundary codes*/
  memset(sc->bndry,0,sc->nv*sizeof(INT));
  // RHS set up: We take as a right hand side the function f(x)=1;
  // the right hand side is just M*ones(n,1);
  dvector f=dvec_create(A.row);
  dvector sol=f;
  dvec_set(f.row,&f,1e0);// this we can have as function value at the vertices in general. 
  dvector rhs=dvec_create(M.row);// this is the "algebraic" right hand side 
  for(i=0;i<M.row;++i){
    fi=0e0;
    for(ij=M.IA[i];ij<M.IA[i+1];++ij){
      j=M.JA[ij];
      fi+=M.val[ij]*f.val[j];
    }
    rhs.val[i]=fi;
  }
  //  dvector_print(stdout,&f);
  //  dvector_print(stdout,&rhs);
  clock_t clk_assembly_end = clock(); // End of timing for mesh and FE setup
  fprintf(stdout,"\n%%%%%%CPUtime(assembly) = %.3f sec\n",
	  (REAL ) (clk_assembly_end - clk_assembly_start)/CLOCKS_PER_SEC);
  if(print_level > 50)
    csr_print_matlab(stdout,&M);  
  dcsr_free(&M); // not needed
  /*SOLVER SOLVER*/
  // use the same as f for the solution;
  linear_itsolver_param linear_itparam;
  linear_itparam.linear_precond_type = PREC_AMG;
  //  linear_itparam.linear_precond_type = PREC_DIAG;
  //  linear_itparam.linear_precond_type = SOLVER_AMG;
  linear_itparam.linear_itsolver_type     = SOLVER_CG;
  linear_itparam.linear_stop_type         = STOP_REL_RES;
  // Solver parameters
  linear_itparam.linear_tol      = 1e-6;
  linear_itparam.linear_maxit    = 100;
  linear_itparam.linear_restart       = 100;
  linear_itparam.linear_print_level    = 7;
  /* Set parameters for algebriac multigrid methods */
  AMG_param amgparam;
  param_amg_init(&amgparam);

  // adjust AMG parameters if needed
  // General AMG parameters
  amgparam.AMG_type             = UA_AMG;
  amgparam.print_level          = 3;
  amgparam.maxit                = 1;
  amgparam.max_levels           = 10;
  amgparam.coarse_dof           = 2000;
  amgparam.cycle_type           = W_CYCLE;
  amgparam.smoother             = SMOOTHER_GS;
  amgparam.presmooth_iter       = 1;
  amgparam.postsmooth_iter      = 1;
  amgparam.coarse_solver        = SOLVER_UMFPACK;
  //amgparam.relaxation           = 1.0;
  //amgparam.polynomial_degree    = 2;
  //amgparam.coarse_scaling       = ON;
  //amgparam.amli_degree          = 2;
  //amgparam.amli_coef            = NULL;
  //amgparam.nl_amli_krylov_type  = SOLVER_VFGMRES;

  // Aggregation AMG parameters
  amgparam.aggregation_type     = VMB;
  amgparam.strong_coupled       = 0.0;
  amgparam.max_aggregation      = 100;

  //amgparam.tentative_smooth     = 0.67;
  //amgparam.smooth_filter        = OFF;

  // print AMG paraemters
  //  param_amg_print(&amgparam);

  /* Actual solve */
  memset(sol.val,0,sol.row*sizeof(REAL));
  switch(linear_itparam.linear_precond_type){
  case SOLVER_AMG:
    // Use AMG as iterative solver
    solver_flag = linear_solver_amg(&A, &rhs, &sol, &amgparam);
    break;
  case PREC_AMG:
    // Use Krylov Iterative Solver with AMG
    solver_flag = linear_solver_dcsr_krylov_amg(&A, &rhs, &sol, &linear_itparam, &amgparam);
    break;
  case  PREC_DIAG:
    solver_flag = linear_solver_dcsr_krylov_diag(&A, &rhs, &sol, &linear_itparam);
    break;
  default:
    solver_flag = linear_solver_dcsr_krylov_amg(&A, &rhs, &sol, &linear_itparam, &amgparam);
    break;
  }
  dcsr_free(&A);
  dvec_free(&rhs);
  return sol;// return solution
}
/**********************************************************************************/
/*!
 * \fn static INT proj_lower_dim(scomplex *dsc)
 *
 * \brief projection of part of a (d-1) dimensional simplicial complex
 *        in R(d) onto a bounded domain in R(d-1). For example,
 *        projecting the boundary (x_4=1) of the 4 dimensional cube
 *        onto the unit cube in 3D. Similar to stereographic
 *        projection of the 3D sphere to the 2D plane.
 *
 * \param dsc I: simplicial complex (usually the boundary of the
 *               d-dimensional cube).
 *
 * \return integer 0(success) or 1(error ocurred)
 *
 * \note Ludmil (20210807)
 */
/**********************************************************************************/
static INT proj_lower_dim(scomplex *dsc)
{
// dsc should be a boundary complex.
  if(dsc->nbig != (dsc->n+1)){
    fprintf(stdout,"\n%%%%******* ERROR: Wrong dimensions of the boundary simplicial complex\n");
    return 1;
  }
  REAL tol=1e-10;
  INT i,j,node,n=dsc->n,n1=dsc->n+1,ns=dsc->ns,nv=dsc->nv;  
  INT nbig=dsc->nbig,nbig1=dsc->nbig+1;
  // mark for removal all vertices with last coordinate WITHIN TOL OF (1) and all simplices attached to them;
  INT *indv=calloc(nv,sizeof(INT));
  for(i=0;i<nv;i++) indv[i]=-1;
  REAL xn;
  INT nvnew=0;
  for(i=0;i<nv;i++){
    xn=dsc->x[nbig*i+nbig-1];
    if(fabs(xn-1e0)<tol) continue;
    indv[i]=nvnew;
    nvnew++;
  }
  INT nsnew=0,remove;
  for(i=0;i<ns;i++){
    remove=0;
    for(j=0;j<n1;j++){
      node=dsc->nodes[i*n1+j];
      if(indv[node]<0){
	remove=1;
	break;
      }
    }
    if(remove) continue;
    for(j=0;j<n1;j++){
      node=dsc->nodes[i*n1+j];
      dsc->nodes[nsnew*n1+j]=indv[node];
    }
    dsc->flags[nsnew]=0;
    nsnew++;
  }
  for(i=0;i<nv;i++){
    node=indv[i];
    if(node>=0) {
      // copy coordinates after projection:
      xn=1e0/(1e0-dsc->x[nbig*i+nbig-1]);
      for(j=0;j<dsc->n;j++){
	dsc->x[node*n+j]=dsc->x[i*nbig+j]*xn;
      }
      dsc->bndry[node]=0;
    }
  }
  free(indv);
  // adjust sizes;
  dsc->nbig=n;
  dsc->ns=nsnew;
  dsc->nv=nvnew;
  dsc->nodes=realloc(dsc->nodes,dsc->ns*(dsc->n+1)*sizeof(INT));
  dsc->x=realloc(dsc->x,(dsc->nv*dsc->n)*sizeof(REAL));
  return 0;
}
/**********************************************************************************/
/*!
 * \fn static void draw_grids(const SHORT todraw,scomplex *sc,
 *                            dvector *sol)
 *
 * \brief writing grids to vtu. Outpput defined for 2,3,4.
 *
 * \param sc  I: simplicial complex defining the grid
 *
 * \param sol I: the numerical solution evaluated at the vertices of
 *               the mesh
 *
 * \note Ludmil (20210807)
 */
/**********************************************************************************/
static void draw_grids(const SHORT todraw,scomplex *sc, dvector *sol)
{ 
  if(!todraw){
    // no drawing needed, so return;
    return;
  }
  /* 
     WRITE THE OUTPUT vtu file for paraview: can be viewed with
     paraview 
  */
  /* 
   * if we want to draw the boundary (this is not supported at the
   * moment) dsc is a simplicial complex defining the boundary of the
   * grid. if (dsc .NOT. null) then projection of the part of the dsc
   * on a simplicial complex in R(d-1) is drawn.
   */
  /* 
     scomplex *dsc=malloc(sizeof(scomplex));
     dsc[0]=sc_bndry(sc); 
  */
  /**/
  INT idsc;
  switch(sc->n){
  case 5:
    fprintf(stdout,"\n%%%% **** NO PLOT: Dimension=%d is too large for plotting\n\n",sc->n);
    break;
  case 4:
    fprintf(stdout,"\n%%%% **** NO PLOT: Dimension=%d is too large for plotting\n\n",sc->n);
    /* if(dsc) { */
    /*   idsc = proj_lower_dim(dsc); */
    /*   vtkw("output/4d_to_3d.vtu",dsc,0,1.); */
    /* } else { */
    /*   fprintf(stdout,"\nNO PLOT: no booundary data"); */
    /* } */
    break;
  case 3:
    vtkw("output/3d.vtu",sc,0,1.);
    /* if(dsc) */
    /*   vtkw("output/3d_to_2d.vtu",dsc,0,1.); */
    break;
  default:
    vtkw("output/2d.vtu",sc,0,1.);
    /* if(dsc) */
    /*   vtkw("output/2d_to_1d.vtu",dsc,0,1.); */
  }
  /*  haz_scomplex_free(dsc);*/
  return;
}
