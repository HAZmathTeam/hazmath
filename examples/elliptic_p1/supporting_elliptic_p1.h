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
/******************************************************************/
static void scomplex_partial_free(scomplex **sc_in, const INT all)
{
  scomplex *sc=NULL;
  if(!sc_in[0])
    return;
  if(all && sc_in[0]){
    sc=sc_in[0];
    if(sc_in[0]->nodes) {free(sc_in[0]->nodes);sc_in[0]->nodes=NULL;}
    if(sc_in[0]->x) {free(sc_in[0]->x);sc_in[0]->x=NULL;}
    free(sc);sc=NULL;
    return;
  }
  sc=sc_in[0];
  if(sc->marked) {
    free(sc->marked);sc->marked=NULL;
  }
  if(sc->gen) {
    free(sc->gen);sc->gen=NULL;
  }
  if(sc->parent) {
    free(sc->parent);sc->parent=NULL;
  }
  if(sc->child0) {
    free(sc->child0);sc->child0=NULL;
  }
  if(sc->childn) {
    free(sc->childn);sc->childn=NULL;
  }
  if(sc->bndry) {
    free(sc->bndry);sc->bndry=NULL;
  }
  if(sc->flags) {
    free(sc->flags);sc->flags=NULL;
  }
  if(sc->vols) {
    free(sc->vols);sc->vols=NULL;
  }
  if(sc->nbr) {
    free(sc->nbr);sc->nbr=NULL;
  }
  if(sc->csys) {
    free(sc->csys);sc->csys=NULL;
  }
  if(sc->etree) {
    free(sc->etree);sc->etree=NULL;
  }
  if(sc->bndry_v) {
    icsr_free(sc->bndry_v);free(sc->bndry_v);sc->bndry_v=NULL;
  }
  if(sc->parent_v) {
    icsr_free(sc->parent_v);free(sc->parent_v);sc->parent_v=NULL;
  }
  if(sc->bfs) {
    icsr_free(sc->bfs);free(sc->bfs);sc->bfs=NULL;
  }
  return;
}
/******************************************************************/
static void scomplex_print_matlab(const char *fname,scomplex *sc)
{
  INT ns_max=100000;
  FILE *fp=fopen(fname,"w");
  if(sc->ns>ns_max){
    fprintf(fp,"\n%%%% Too large:elements=%d>%d\n",sc->ns,ns_max);
    fclose(fp);
    return;
  }
  INT i, j, k, in,in1,n1=sc->n+1,dim=sc->n,nv=sc->nv,ns=sc->ns;
  fprintf(fp,"\nt=[");
  for(j=0;j<n1;j++){
    for(i=0;i<ns-1;++i){
      fprintf(fp,"%d,",sc->nodes[i*n1+j]);
    }
    fprintf(fp,"%d;\n",sc->nodes[(ns-1)*n1+j]);
  }
  fprintf(fp,"];t=t';\n");
  fprintf(fp,"\nnbr=[");
  for(j=0;j<n1;j++){
    for(i=0;i<ns-1;i++){
      fprintf(fp,"%d,",sc->nbr[i*n1+j]);
    }
    fprintf(fp,"%d;\n",sc->nbr[(ns-1)*n1+j]);
  }
  fprintf(fp,"];nbr=nbr';\n");
  //
  fprintf(fp,"\nxp=[");
  for(j=0;j<dim;j++){
    for(i=0;i<nv-1;i++){
      fprintf(fp,"%.16e,",sc->x[i*dim+j]);
    }
    fprintf(fp,"%.16e;\n",sc->x[(nv-1)*dim+j]);
  }
  fprintf(fp,"];xp=xp';\n");
  fprintf(fp,"\nib=[");  
  for(i=0;i<nv-1;i++){
    fprintf(fp,"%d,",sc->bndry[i]);
  }
  fprintf(fp,"%d];ib=ib';\n",sc->bndry[nv-1]);
  fclose(fp);
  return;
}
void find_bndry_vertices(const INT dim,const INT ns,INT *je, INT *idir)
{
  /*    
	%% we call an intersection interior if two simplices from je
if they intersect in dim points.  in a simplicial mesh of
	%% dimension dim uses face_vertex=fv and vertex_simplex=vt to
	%% find the boundary faces, interior faces boundary vertices,
	%% and interior vertices,
  */
  INT dim1=dim+1,nv=-1,i,j,k,m;
  iCSRmat f2v;
  nv=je[0];
  // find the number of DOFs
  for(j=1;j<dim1*ns;++j)
    if(je[j]>nv) nv=je[j];
  nv++; // in C the max of je is nv-1 because in C we index from 0.
  // now we have the number of vertices: it is nv.
  INT *nbr=calloc(ns*dim1,sizeof(INT));

  // find neighboring elements and construct nbring list which is in
  // accordance with the numbering in je, i.e. the neighbor of el at
  // place (j) is opposite to the j-th vertex in je(el,:)  
  find_nbr(ns,nv,dim,je,nbr);  
  ////////////////////////////////////////////////////////////
  INT isn1,is,nbf,nnzbf;
  // count the boundary vertices.
  nbf=0;
  nnzbf=0;
  for(i=0;i<ns;i++){
    isn1=i*dim1;
    for(j=0;j<dim1;j++){
      is=nbr[isn1+j];
      if(is<0){
	nbf++;
	nnzbf+=dim;
      }
    }
  }
  //
  // now working on the boundary:
  //
  f2v.row=nbf;
  f2v.col=nv;
  f2v.nnz=dim*nbf;
  f2v.IA=calloc(nbf+1,sizeof(INT));
  f2v.JA=calloc(nbf*dim,sizeof(INT));
  f2v.val=NULL;
  /* 
     forming the boundary_face_2_vertex matrix uses that the
     neighboring list of elements is in accordance with the
  */
  //  fprintf(stdout,"\nnbf=%d; GUESS=%d,f2v_nnz=%d\n",nbf,f2v.nnz,nbf*dim);fflush(stdout);
  INT nbfnew=0;
  INT nnzf2v=0;
  f2v.IA[0]=nnzf2v;
  for(i=0;i<ns;i++){
    for(j=0;j<dim1;j++){
      if(nbr[i*dim1+j]<0) {
	for(m=0;m<dim1;m++){
	  if(m==j) continue;
	  f2v.JA[nnzf2v]=je[i*dim1+m];
	  nnzf2v++;
	}
	nbfnew++;
	f2v.IA[nbfnew]=nnzf2v;
      }
    }
  }
  f2v.nnz=nnzf2v;
  if(nbf!=nbfnew){
    fprintf(stderr,"\n%%***ERROR(1): num. bndry faces mismatch (nbf=%d .ne. nbfnew=%d) in %s",nbf,nbfnew,__FUNCTION__);
    exit(65);
  }
  f2v.IA[nbf]=f2v.nnz;
  f2v.JA=realloc(f2v.JA,f2v.nnz*sizeof(INT));
  ///////////////////////////////////////////////////////////
  f2v.row=nbf;
  f2v.col=nv;
  //  fprintf(stdout,"\nnbf=%d; nv=%d; f2v.nnz=%d\n",nbf,nv,f2v.nnz);fflush(stdout);
  memset(idir,0,nv*sizeof(INT));
  // for bndry vertices k, idir[k]=(nv+1)
  for(k=0;k<f2v.nnz;++k){
    idir[f2v.JA[k]]=nv+1;
  }
  //  icsr_write_icoo("fv.dat",&f2v);
  icsr_free(&f2v);
  free(nbr);
  return;
}
/**************************************************************************/
/******************************************************************************/
static void symb_assembly(INT ns, INT ndof, INT ndofloc, INT *je, INT **ia_out, INT **ja_out,INT *idir)
{
  /*
   *  ====================================================================
   * ndofloc are the ndofs per element
   *  Input:
   *    je is sc->nodes (see scomplex) ie is allocated here.
   *    ndof        - number of dofs
   *    ndofloc        - number of local dofs per simplex
   *    ns        - number of elements
   *    idir      - array of dimension n which contains the value  
   *               n+1 in the positions corresponding to Dirichlet 
   *               nodes and 0 elsewhere.
   *
   *  Output:
   *    **ia_out, **ja_out   - structure of the global matrix, allocated here;
   *
   * --------------------------------------------------------------------
   */
  INT nnz=-1,nnzi=-1,np=-1,i,j,k;
  INT ieta,ietb,iea,ieb,ip,kp;
  dCSRmat el,elt;
  el.row=ns;  el.col=ndof; el.nnz=ndofloc*ns;
  el.IA=calloc(ns+1,sizeof(INT));
  el.IA[0]=0;
  for(j=0;j<ns;++j){
    el.IA[j+1]=el.IA[j]+ndofloc;
  }
  el.JA=je; el.val=NULL;
  //
  elt.row=ndof;  elt.col=ns; elt.nnz=el.nnz;
  elt.IA=calloc(ndof+1,sizeof(INT));
  elt.JA=calloc(elt.nnz,sizeof(INT));
  elt.val=NULL;
  //
  dcsr_transz(&el,NULL,&elt);
  //
  INT *ie=el.IA, *iet=elt.IA, *jet=elt.JA;
  //
  nnz = 0;
  // the array idir should have ndof+1 at all Dirichlet ndofs. It can also have ndof, but this is a choice; 
  np=ndof+1;
  //allocate ia[]
  ia_out[0]=calloc(ndof+1,sizeof(INT));
  INT *ia=ia_out[0];
  //
  for(i=0;i<ndof;++i){
    nnzi = nnz;
    if(idir[i]!=np){
      ieta = iet[i];
      ietb = iet[i+1];
      for(ip=ieta;ip<ietb;++ip){
	j = jet[ip];
	iea = ie[j];
	ieb = ie[j+1];
	for(kp = iea;kp<ieb;++kp){
	  k = je[kp];
	  if (idir[k]>=i) continue;
	  nnz++;
	  idir[k] = i;
	}
      }    
    } else {
      // diagonal only;
      nnz++;
    }
    ia[i] = nnzi;
  }
  ia[ndof] = nnz;
  // init idir
  for(i=0;i<ndof;++i)
    if(idir[i]<np && idir[i] > (-1))
      idir[i]=-1;
  //
  /*SECOND RUN, allocate ja and put the column indices in ja*/
  ja_out[0]=calloc(ia[ndof],sizeof(INT));
  INT *ja=ja_out[0];
  nnz=0;
  for(i=0;i<ndof;++i){
    nnzi = nnz;
    if(idir[i]!=np){
      ieta = iet[i];
      ietb = iet[i+1];
      for(ip=ieta;ip<ietb;++ip){
	j = jet[ip];
	iea = ie[j];
	ieb = ie[j+1];
	for(kp = iea;kp<ieb;++kp){
	  k = je[kp];
	  if (idir[k]>=i) continue;
	  ja[nnz] = k;
	  nnz++;
	  idir[k] = i;
	}
      }    
    } else {
      // diagonal only;
      ja[nnz] = i;
      nnz++;
    }
    // this is already done;    ia[i] = nnzi;
  }
  if(ia[ndof]!=nnz){
    fprintf(stderr,"\n%%%% ERROR: NNZ mismatch in symbolic assembly in %s",__FUNCTION__);
    if(ia) {free(ia);ia=NULL;}// this sets ia_out[0] to NULL;
    if(ja) {free(ja);ja=NULL;}// this sets ja_out[0] to null;
    dcsr_free(&elt);
    if(el.IA) free(el.IA);
    exit(16);
  } 
  dcsr_free(&elt);
  if(el.IA) free(el.IA);
  return;
}
/*C====================================================================*/
static void num_assembly(INT ndof, INT *ia, INT *ja,	\
			 INT ndofloc, INT *jeloc,	\
			 REAL *a, REAL *rhs,		\
			 REAL *aloc,REAL *rhsloc,		\
			 INT *idir, INT *ip)
{
  /* ====================================================================
   *  Numerical assembly of an element matrix AE and vector BE into 
   *  the nodal assembly matrix A and right-hand vector B: General 
   *  case.
   *
   *  Input:  
   *    IA, JA - structure of matrix A in RRCU. The order of A is N.
   *    IDIR   - Array which identifies Dirichlet nodes; it contains
   *             n for a Dirichlet node, 0 otherwise.
   *    aloc, rhsloc - element matrix and vector to be  assembled in AN, and B.
   *
   *    JELOC    - local to global map
   *    ndofloc - number of DOFs per element
   *
   *  Output: 
   *    a - values of nonzeros of the global matrix.
   *    rhs - right-hand side (dvec)
   *
   *  Working space:
   *    IP     - of demension N, initialized to -1 before the loop over elements; IP is used here and 
   *             then reset to -1. 
   *
   *  Note:
   *    The prescribed values of the Dirichlet unknowns should be 
   *    stored in the corresponding positions of B. In other words if
   *    i is a Dirichlet node with prescribed value C_i, then, before
   *    using this algorithm, set IDIR(i) = n+1, B(i) = C_i.
   --------------------------------------------------------------------
  */
  INT kl=-1,m=-1,ll=-1,l=-1,i=-1,j=-1,k=-1,iaa,iab;
  kl=0;
  for(l=0;l<ndofloc;++l){
    i = jeloc[l];
    // check if we are at a dirichlet ndof. If we are, there is nothing to do. 
    if (idir[i]>ndof) continue;
    rhs[i] = rhs[i] + rhsloc[l];// add to rhs
    //kl=l*ndofloc;
    for(ll = 0;ll<ndofloc;++ll){
      k=kl+ll;// k=l*ndofloc+ll = (l,ll) locally is  (i,j) in a
      j = jeloc[ll];
      if(idir[j]>ndof){
	rhs[i] -= aloc[k]*rhs[j]; // rhs[i]=rhs[i]-a[i,j]*rhs[j] if j is dirichlet and i is not. 
      } else {
	ip[j] = k; // save the position of a(i,j)=aloc(l,ll)
      }
    }
    iaa = ia[i];
    iab = ia[i+1];
    m = 0;
    for(j=iaa;j<iab;++j){
      k = ip[ja[j]];
      if(k<0) continue; 
      a[j] = a[j] + aloc[k]; // add
      ip[ja[j]] = -1;// init
      m++;
      if(m==ndofloc) break; // break if we got ndofloc elements in the row of iaa.
    }
    kl += ndofloc; //kl=l*ndofloc, so this is the next one;
  }
  return;
}
/*************************************************************************/
static dvector fe_sol_no_dg(scomplex *sc,const REAL alpha,const REAL gamma)
{
  INT solver_flag=-10,print_level=0;
  clock_t clk_assembly_start = clock(); // begin assembly timing;
  INT i,j,k,idim1,jdim1;  // loop and working
  INT ns,nv,nnz; // num simplices, vertices, num nonzeroes
  REAL volume,fact; //mass matrix entries and dim factorial.  
  // for simplices: number of vertices per simplex. 
  INT dim=0,dim1=1;
  //  scomplex *sc=hazr("3d_example");
  fprintf(stdout,"ns=%d,nv=%d,dim=%d\n",sc->ns,sc->nv,sc->n); fflush(stdout);
  dim=sc->n;
  dim1=sc->n+1;  
  nv=sc->nv;// shorthand for num vertices. 
  ns=sc->ns; // shorthand for num simplices
  /*=====================================================*/
  fact=sc->factorial;
  /*=====================================================*/
  //
  // local mass matrix: it is a constant matrix times the volume of an element.
  REAL *mlocal=local_mm(dim);
  fprintf(stdout,"\nnum_simplices=%d ; num_vertices=%d",ns,nv);fflush(stdout);
  // to compute the volume and to grab the local coordinates of the vertices in the simplex we need some work space
  REAL *slocal=calloc(dim1*dim1,sizeof(REAL));// local stiffness matrix.
  REAL *grad=calloc(dim1*dim,sizeof(REAL));// local matrix with
  					   // gradients from which
					   // many things can be
					   // computed.
  REAL *xs=calloc(dim*dim1,sizeof(REAL));// for the local coordinates of vertices of a simplex;
  REAL *rhsloc=calloc(dim1,sizeof(REAL));// for the local rhs
  void *wrk=calloc(dim1*dim1,sizeof(REAL));// this is used in every simplex but is allocated only once.
  //
  INT ndof=nv,ndofloc=dim1;
  ivector idir=ivec_create(sc->nv);
  ivector ip=ivec_create(sc->nv);
  ivec_set(idir.row,&idir,-1);
  ivec_set(ip.row,&ip,-1);
  find_bndry_vertices(sc->n,sc->ns,sc->nodes,idir.val);
  ivec_set(idir.row,&idir,-1);  
  dCSRmat A;
  symb_assembly(ns, ndof, ndofloc, sc->nodes,&A.IA, &A.JA,idir.val);
  //
  A.row=nv; A.col=nv;  A.nnz=A.IA[nv];
  A.val=calloc(A.nnz,sizeof(REAL));
  memset(A.val,0,A.nnz*sizeof(REAL));
  //
  dvector rhs=dvec_create(nv);
  for(k=0;k<nv;++k){
    if(idir.val[k]>nv){
      rhs.val[k]=0e0;// this is the global rhs;
    }
  }
  REAL smjk;
  sc_vols(sc);// compute the volumes if not yet computed;
  for(i=0;i<ns;++i){
    idim1=i*dim1;
    // grab the vertex coordinates and compute the volume of the simplex.
    local_coords(dim,xs,&sc->nodes[idim1],sc->x);// sc->nodes[idim1:idim1+dim]
						 // point to global
						 // vertex numbers in
						 // element i
    // compute gradients
    grad_compute(dim, fact, xs, grad,wrk);
    //    print_full_mat(dim1,dim,grad,"grad");fflush(stdout);
    volume=sc->vols[i];
    //    sc->vols[i]=volume;
    //    fprintf(stdout,"\nvolume[%d]=%.5e",i,sc->vols[i]);fflush(stdout);
    // copute local stiffness matrix as grad*transpose(grad);
    local_sm(slocal,grad,dim,volume);
    for(j=0;j<dim1;j++){
      jdim1=j*dim1;
      rhsloc[j]=0.;
      for(k=0;k<dim1;k++){
	smjk=mlocal[jdim1+k]*volume;
	rhsloc[j]+=smjk*1e0;// here rhs(x)=1;
	slocal[jdim1+k]=alpha*slocal[jdim1+k] + gamma*smjk;
      }
    }
    num_assembly(ndof,A.IA,A.JA,			\
		 ndofloc, &sc->nodes[idim1],		\
		 A.val, rhs.val,			\
		 slocal,rhsloc,				\
		 idir.val,ip.val);
  }
  ivec_free(&ip);
  free(xs);
  free(rhsloc);
  free(slocal);
  free(grad);
  free(wrk);
  free(mlocal);
  ivec_free(&idir);
  clock_t clk_assembly_end = clock(); // End of timing for mesh and FE setup
  fprintf(stdout,"\n%%%%%%CPUtime(assembly) = %.3f sec\n",
	  (REAL ) (clk_assembly_end - clk_assembly_start)/CLOCKS_PER_SEC);
  /*SOLVER SOLVER*/
  dvector sol=dvec_create(nv);
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
