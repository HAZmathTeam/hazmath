#include "hazmath.h"
/***********************************************************************************************/
/**
 * \fn static void Schwarz_levels0 (INT inroot, dCSRmat *A, INT *mask, INT *nlvl,
 *                                 INT *iblock, INT *jblock, INT maxlev)
 *
 * \brief Form the level hierarchy of input root node
 *
 * \param inroot  Root node
 * \param A       Pointer to CSR matrix
 * \param mask    Pointer to flag array
 * \param nlvl    The number of levels to expand from root node
 * \param iblock  Pointer to vertices number of each level
 * \param jblock  Pointer to vertices of each level
 * \param maxlev  The maximal number of levels to expand from root node
 *
 * \note  This needs to be rewritten -- Xiaozhe
 *
 */
static void Schwarz_levels0 (INT inroot,
			     dCSRmat *A,
			     INT *mask,
			     INT *nlvl,
			     INT *iblock,
			     INT *jblock,
			     INT maxlev)
{
  // performs bfs starting in inroot and stops at distance maxlev from
  // inroot.
  INT *ia = A->IA;
  INT *ja = A->JA;
  INT nnz = A->nnz;
  INT i, j, lvl, lbegin, lvlend, nsize, node;
  INT jstrt, jstop, nbr, lvsize;

  // This is diagonal
  if (ia[inroot+1]-ia[inroot] <= 1) {
    lvl = 0;
    iblock[lvl] = 0;
    jblock[iblock[lvl]] = inroot;
    lvl ++;
    iblock[lvl] = 1;
  }
  else {
    // input node as root node (level 0)
    lvl = 0;
    jblock[0] = inroot;
    lvlend = 0;
    nsize  = 1;
    // mark root node
    mask[inroot] = 1;

    lvsize = nnz;

    // start to form the level hierarchy for root node(level1, level2, ... maxlev)
    while (lvsize > 0 && lvl < maxlev) {
      lbegin = lvlend;
      lvlend = nsize;
      iblock[lvl] = lbegin;
      lvl ++;
      for(i=lbegin; i<lvlend; ++i) {
	node = jblock[i];
	jstrt = ia[node]-1;
	jstop = ia[node+1]-1;
	for (j = jstrt; j<jstop; ++j) {
	  nbr = ja[j]-1;
	  if (mask[nbr] == 0) {
	    jblock[nsize] = nbr;
	    mask[nbr] = lvl;
	    nsize ++;
	  }
	}
      }
      lvsize = nsize - lvlend;
    }

    iblock[lvl] = nsize;

    // reset mask array
    for (i = 0; i< nsize; ++i) {
      node = jblock[i];
      mask[node] = 0;
    }
  }

  *nlvl = lvl;
}
/*! \file examples/Solver/Solver.c
 *
 *  Created by Xiaozhe Hu on 01/01/19.
 *  Copyright 2019_HAZMATH__. All rights reserved.
 *
 * \brief This program read in a matrix and a right hand side and solve it by certain linear solvers
 *
 * \note
 *
 */

void Schwarz_get_block_matrix0 (Schwarz_data *Schwarz,
				INT nblk,
				INT *iblock,
				INT *jblock,
				INT *mask)
{
  INT i, j, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;
  INT maxbs = 0, count, nnz;

  dCSRmat A = Schwarz->A;
  dCSRmat *blk = Schwarz->blk_data;

  INT  *ia  = A.IA;
  INT  *ja  = A.JA;
  REAL *val = A.val;

  // get maximal block size
  for (is=0; is<nblk; ++is) {
    ibl0 = iblock[is];
    ibl1 = iblock[is+1];
    nloc = ibl1-ibl0;
    maxbs = MAX(maxbs, nloc);
  }

  Schwarz->maxbs = maxbs;
  fprintf(stdout,"\nmax_blk_size=%d",maxbs);
  // allocate memory for each sub_block's right hand
  Schwarz->xloc1   = dvec_create(maxbs);
  Schwarz->rhsloc1 = dvec_create(maxbs);

  for (is=0; is<nblk; ++is) {
    ibl0 = iblock[is];
    ibl1 = iblock[is+1];
    nloc = ibl1-ibl0;
    count = 0;
    for (i=0; i<nloc; ++i ) {
      iblk = ibl0 + i;
      ki   = jblock[iblk];
      iaa  = ia[ki];
      iab  = ia[ki+1];
      count += iab - iaa;
      mask[ki] = i+1;
    }

    blk[is] = dcsr_create(nloc, nloc, count);
    blk[is].IA[0] = 0;
    nnz = 0;

    for (i=0; i<nloc; ++i) {
      iblk = ibl0 + i;
      ki = jblock[iblk];
      iaa = ia[ki];
      iab = ia[ki+1];
      for (kij = iaa; kij<iab; ++kij) {
	kj = ja[kij];
	j  = mask[kj];
	if(j != 0) {
	  blk[is].JA[nnz] = j-1;
	  blk[is].val[nnz] = val[kij];
	  nnz ++;
	}
      }
      blk[is].IA[i+1] = nnz;
    }

    blk[is].nnz = nnz;

    // zero the mask so that everyting is as it was
    for (i=0; i<nloc; ++i) {
      iblk = ibl0 + i;
      ki   = jblock[iblk];
      mask[ki] = 0;
    }
  }
}

/**
 * \fn INT Schwarz_setup0 (Schwarz_data *Schwarz, Schwarz_param *param)
 *
 * \brief Setup phase for the Schwarz methods
 *
 * \param Schwarz    Pointer to the Schwarz data
 * \param param      Type of the Schwarz method
 *
 * \return           SUCCESS if succeed
 *
 */
INT Schwarz_setup0 (Schwarz_data *Schwarz,
		    Schwarz_param *param)
{
  // information about A
  dCSRmat A = Schwarz->A;
  INT n = A.row;
  INT  block_solver = param->Schwarz_blksolver;
  INT  maxlev = param->Schwarz_maxlvl;
  Schwarz->swzparam = param;
  // local variables
  INT n1=n+1,i,ij,j,							\
    ibegin=0,								\
    inroot = -10, nsizei = -10, nsizeall = -10, nlvl = 0;
  INT *jbegin=NULL;
  //
  // these must be set previously; set to 0 if no blocks were defined
  // so far
  INT nblk=Schwarz->nblk;
  INT *iblock=NULL,*jblock=NULL;
  ivector MaxIndSet;
  //
  // data for Schwarz method
  INT *mask = NULL, *maxa = NULL;

  // return
  INT flag = SUCCESS;
  // allocate memory
  maxa    = (INT *)calloc(n1,sizeof(INT));
  mask    = (INT *)calloc(n1,sizeof(INT));
  // if we have already allocated some blocks, here we define some more. 
  fprintf(stdout,"\nSchwarz method(solver=%d ; MIS levels=%d)\n",	\
	  block_solver,				\
	  maxlev);fflush(stdout);
  fprintf(stdout,"\nSchwarz method(nblk=%d; ",nblk);fflush(stdout);
  if(nblk>0){
    iblock=Schwarz->iblock;
    fprintf(stdout," 1end=%d; nblk+1+n1=%d)\n",iblock[nblk],nblk+1+n1);fflush(stdout);
    iblock=(INT *)realloc(iblock,(nblk+1+n1)*sizeof(INT));
    fprintf(stdout," 2end=%d)\n",iblock[nblk]);fflush(stdout);
    jblock=Schwarz->jblock;
    jblock=(INT *)realloc(jblock,(iblock[nblk]+n1)*sizeof(INT));
    memset((iblock+nblk+1), 0, sizeof(INT)*n1);
    ibegin=nblk;
  }else{
    iblock=(INT *)calloc(n1,sizeof(INT));
    jblock=(INT *)calloc(n1,sizeof(INT));
    memset(iblock, 0, sizeof(INT)*n1);
    ibegin=0;
    fprintf(stdout,"\nend=NA");fflush(stdout);
  }  
  maxa    = (INT *)calloc(n1,sizeof(INT));
  mask    = (INT *)calloc(n1,sizeof(INT));
  memset(mask,   0, sizeof(INT)*n1);
  memset(maxa,   0, sizeof(INT)*n1);
  nsizeall=0;
  maxa[0]=0;
  MaxIndSet = sparse_MIS(&A);
  memset(maxa,   0, sizeof(INT)*n);
  /*-------------------------------------------*/
  // find the blocks
  /*-------------------------------------------*/
  // first pass: do a maxlev level sets out for each node from MIS,
  // or perhaps from the last block.
  for ( i = 0; i < MaxIndSet.row; i++ ) {
    inroot = MaxIndSet.val[i];
    Schwarz_levels0(inroot,&A,mask,&nlvl,maxa,jblock,maxlev);
    nsizei=maxa[nlvl];
    nsizeall+=nsizei;
  }
  /* We only calculated the size of this up to here. So we can reallocate jblock */
  jblock = (INT *)realloc(jblock,(iblock[nblk]+nsizeall+n)*sizeof(INT));
  // second pass: redo the same again, but this time we store in jblock
  maxa[0]=0;
  nsizeall=0;
  jbegin=jblock+iblock[nblk];
  for (i=0;i<MaxIndSet.row;i++) {
    inroot = MaxIndSet.val[i];
    Schwarz_levels0(inroot,&A,mask,&nlvl,maxa,jbegin,maxlev);
    nsizei=maxa[nlvl];
    iblock[ibegin+i+1]=iblock[ibegin+i]+nsizei;
    nsizeall+=nsizei;
    jbegin+=nsizei;
  }
  nblk += MaxIndSet.row;
  fprintf(stdout,"\nBlocks print: Number of blocks=%d\n",nblk);
  for(i=0;i<nblk;i++){
    fprintf(stdout,"\nblock=%d; nodes=(",i);
    for(ij=iblock[i];ij<iblock[i+1];ij++){
      fprintf(stdout," %d ",jblock[ij]);
    }
    fprintf(stdout,")");
  }
  fprintf(stdout,"\n");
  /*-------------------------------------------*/
  //  LU decomposition of blocks
  /*-------------------------------------------*/
  memset(mask, 0, sizeof(INT)*n1);
  Schwarz->blk_data = (dCSRmat*)calloc(nblk, sizeof(dCSRmat));
  Schwarz_get_block_matrix0(Schwarz, nblk, iblock, jblock, mask);
  // Setup for each block solver
  switch (block_solver) {
#if WITH_SUITESPARSE
  case SOLVER_UMFPACK: {
    /* use UMFPACK direct solver on each block */
    dCSRmat *blk = Schwarz->blk_data;
    void **numeric	= (void**)calloc(nblk, sizeof(void*));
    dCSRmat Ac_tran;
    //printf("number of blocks = %d\n",nblk);
    for (i=0; i<nblk; ++i) {
      Ac_tran = dcsr_create(blk[i].row, blk[i].col, blk[i].nnz);
      dcsr_transz(&blk[i], NULL, &Ac_tran);
      dcsr_cp(&Ac_tran, &blk[i]);
      //printf("size of block %d: nrow=%d, nnz=%d\n",i, blk[i].row, blk[i].nnz);
      numeric[i] = umfpack_factorize(&blk[i], 0);
    }
    Schwarz->numeric = numeric;
    dcsr_free(&Ac_tran);

    break;
  }
#endif

  default: {
    /* do nothing for iterative methods */
  }
  }
  
  /*-------------------------------------------*/
  //  return and reassign the pointers and the added blocks. 
  /*-------------------------------------------*/
  Schwarz->nblk   = nblk;
  Schwarz->iblock = iblock;
  Schwarz->jblock = jblock;
  Schwarz->mask   = mask;
  Schwarz->maxa   = maxa;
  Schwarz->Schwarz_type = param->Schwarz_type;
  Schwarz->blk_solver = param->Schwarz_blksolver;

  printf("Schwarz method setup is done! Find %d blocks\n",nblk);

  return flag;
}
/************* HAZMATH FUNCTIONS and INCLUDES ***************************/
INT main (int argc, char* argv[])
{

  printf("\n===========================================================================\n");
  printf("Beginning Program to solve a linear system.\n");
  printf("===========================================================================\n");

  /* matrix and right hand side */
  dCOOmat Acoo;
  dCSRmat A;
  dvector b;
  dvector x;
  printf("\n===========================================================================\n");
  printf("Reading the matrix, right hand side, and parameters\n");
  printf("===========================================================================\n");

  /* read the matrix and right hand side */
  //  dcoo_read_dcsr("A.dat", &A);
  //  dvector_read("b.dat", &b);
  /****************************READING (fractures)******************************************/
  INT i,j,k,ij,jk,ik;
  SHORT  ifile=1;
  char *fname=(char *)malloc(256*sizeof(char));  
  char *dirname=strdup("../fractures/LS/");  
  i=sprintf(fname,"%smatrix%1d.ijv",dirname,ifile);
  FILE *fin = fopen(fname,"r");  
  fprintf(stdout,"\n%s:%i\n",fname,i);fflush(stdout);
  i=fscanf(fin,"%i",&(Acoo.row));
  i=fscanf(fin,"%i",&(Acoo.col));
  i=fscanf(fin,"%i",&(Acoo.nnz));
  Acoo.rowind=calloc(Acoo.nnz,sizeof(INT));
  Acoo.colind=calloc(Acoo.nnz,sizeof(INT));
  Acoo.val=calloc(Acoo.nnz,sizeof(REAL));
  fprintf(stdout,"\nReading the matrix...nrow=%d, ncol=%d, nnz=%d",Acoo.row,Acoo.col,Acoo.nnz);fflush(stdout);
  for(i=0;i<Acoo.nnz;i++){
    fscanf(fin,"%i %i %lg",(Acoo.rowind+i),(Acoo.colind+i),(Acoo.val+i));
    Acoo.rowind[i]--;Acoo.colind[i]--;
    //      fprintf(stdout,"\n%i: %i %i %23.16e",i,Acoo->rowind[i],Acoo->colind[i],Acoo->val[i]);
  }
  fprintf(stdout,"... %d nonzeroes: DONE.\n",Acoo.nnz);fflush(stdout);
  fclose(fin);
  i=sprintf(fname,"%srhs%1d.dat",dirname,ifile);
  fprintf(stdout,"\n%s\n",fname);fflush(stdout);
  fin = HAZ_fopen(fname,"r");
  b.row = Acoo.row; b.val = calloc(b.row,sizeof(REAL));
  fprintf(stdout,"\nReading the rhs...");
  rvecd_(fin,b.val,&(b.row));
  fprintf(stdout,"... %d rows: DONE.\n",b.row);fflush(stdout);
  fclose(fin);
  i=sprintf(fname,"%smatrix_structure%1d.dat",dirname,ifile);
  fprintf(stdout,"\n%s\n",fname);fflush(stdout);
  fin = HAZ_fopen(fname,"r");
  INT nblk,ibl0,ibl1;
  i=fscanf(fin,"%i",&nblk);
  INT *iblock=(INT *)calloc(nblk+1,sizeof(INT)); 
  iblock[nblk]=Acoo.row;
  fprintf(stdout,"\nReading the matrix structure...");fflush(stdout);
  for(i=0;i<nblk;i++){
    fscanf(fin,"%i",(iblock+i));
    iblock[i]--;
  }
  fclose(fin);
  if(fname)free(fname);
  INT *jblock=(INT *)calloc(iblock[nblk],sizeof(INT)); 
  for(i=0;i<nblk;i++){
    ibl0=iblock[i];
    ibl1=iblock[i+1];
    for(ij=ibl0;ij<ibl1;ij++){
      jblock[ij]=ij;
    }
  }
  INT nbmax=1,flag=0;
  INT mbi,mbsize[20],imax[20];
  for(j=0;j<nbmax;j++){
    imax[j]=-1;mbsize[j]=0;
    for(i=1;i<=(nblk);i++){
      flag=0;
      for(k=0;k<j;k++){
	if((i-1)==imax[k]) {flag=1;break;}
      }
      if(flag) continue;
      mbi=(iblock[i]-iblock[i-1]);
      if(mbsize[j]<mbi){mbsize[j]=mbi;imax[j]=i-1;}
    }
    fprintf(stdout,"\n...pass=%d; max size=%d at block %d;",j,mbsize[j],imax[j]);fflush(stdout);
  }
  fprintf(stdout,"\nDONE.\n");
  // let us see if we can solve the largest block:
  INT *mask=(INT *)calloc(Acoo.row+1,sizeof(INT)); 
  A.row=Acoo.row;
  A.col=Acoo.col;
  dcoo_2_dcsr(&Acoo,&A);
  /*****************************************************************************************/
  /* set Parameters from Reading in Input File */
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  Schwarz_data *swz=malloc(sizeof(Schwarz_data));
  Schwarz_param *swzparam=malloc(sizeof(Schwarz_param));
  swz->A=A;
  swz->nblk=nblk;
  swz->iblock=iblock;
  swz->jblock=jblock;
  swzparam->Schwarz_blksolver=inparam.Schwarz_blksolver;
  swzparam->Schwarz_maxlvl=inparam.Schwarz_maxlvl;
  /*********************************/
  Schwarz_setup0(swz,swzparam);
  /*****************************************************************************************/
  /* set initial guess */
  dvec_alloc(A.row, &x);
  dvec_set(x.row, &x, 0.0);

  /* Set Solver Parameters */
  INT solver_flag=-20;

  /* Set parameters for linear iterative methods */
  linear_itsolver_param linear_itparam;
  param_linear_solver_set(&linear_itparam, &inparam);

  /* Set parameters for algebriac multigrid methods */
  AMG_param amgparam;
  param_amg_init(&amgparam);
  param_amg_set(&amgparam, &inparam);
  param_amg_print(&amgparam);

  printf("\n===========================================================================\n");
  printf("Solving the linear system \n");
  printf("===========================================================================\n");
  exit(129);

  // Use AMG as iterative solver
  if (linear_itparam.linear_itsolver_type == SOLVER_AMG){
    solver_flag = linear_solver_amg(&A, &b, &x, &amgparam);
  } else { // Use Krylov Iterative Solver
    // Diagonal preconditioner
    if (linear_itparam.linear_precond_type == PREC_DIAG) {
      solver_flag = linear_solver_dcsr_krylov_diag(&A, &b, &x, &linear_itparam);
    }
    // AMG preconditioner
    else if (linear_itparam.linear_precond_type == PREC_AMG){
      solver_flag = linear_solver_dcsr_krylov_amg(&A, &b, &x, &linear_itparam, &amgparam);
    }
    // No preconditoner
    else{
      solver_flag = linear_solver_dcsr_krylov(&A, &b, &x, &linear_itparam);
    }
  }

  // Clean up memory
  dcsr_free(&A);
  dvec_free(&b);

}	/* End of Program */
/*******************************************************************/
