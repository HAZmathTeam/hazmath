/*! \file src/solver/Schwarz_setup.c
 *
 *  Setup phase for the Schwarz methods
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/25/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note  Done cleanup for releasing -- Xiaozhe Hu 03/12/2017
 *
 */

#include "hazmath.h"

/**
 * \fn static void Schwarz_levels (INT inroot, dCSRmat *A, INT *mask, INT *nlvl,
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
 * \note modified (ltz) 20230131
 *
 *
 */
static void Schwarz_levels (INT inroot,
                            dCSRmat *A,
                            INT *mask,
                            INT *nlvl,
                            INT *iblock,
                            INT *jblock,
                            INT maxlev)
{
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
    mask[inroot]=0;
  } else {
    // input node as root node (level 0)
    lvl = 0;
    jblock[0] = inroot;
    lvlend = 0;
    nsize  = 1;
    // mark root node
    mask[inroot] = 0; //
    lvsize = nnz;      
    // start to form the level hierarchy for root node(level1, level2, ... maxlev)
    while (lvsize > 0 && lvl < maxlev) {
      lbegin = lvlend;
      lvlend = nsize;
      iblock[lvl] = lbegin;
      lvl ++;
      for(i=lbegin; i<lvlend; ++i) {
	node = jblock[i];
	jstrt = ia[node];
	jstop = ia[node+1];
	for (j = jstrt; j<jstop; ++j) {
	  nbr = ja[j];
	  if (mask[nbr] < 0) {
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
    for (i=0;i<nsize;++i) {
      node = jblock[i];
      mask[node] = -1;
    }
  }
  *nlvl = lvl;
  return;
}
/***********************************************************************************************/
/**
 * \fn static void Schwarz_get_block_matrix (Schwarz_data *Schwarz, INT nblk,
 *                                    INT *iblock, INT *jblock, INT *mask)
 *
 * \brief Forms blocks for a Schwarz smoother. If a dof is not in any
 *        block, then a block of size=1 is created. In this way each
 *        dof is in at least one block.
 *
 * \param Schwarz Pointer to the Schwarz data
 *
 * \note  Done cleanup for releasing -- Xiaozhe Hu 08/28/2021
 *
 * \note modified (ltz) 20230131
 *
 */
static void Schwarz_get_block_matrix(Schwarz_data *Schwarz)
{
  INT i, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;
  INT maxbs,maxbnnz,nnz;    
  dCSRmat A = Schwarz->A;
  INT *iblock = Schwarz->iblock;
  INT *jblock = Schwarz->jblock;
  INT *mask = Schwarz->mask;
  INT nblk=Schwarz->nblk;
  dCSRmat *blk=NULL;
  INT  *ia  = A.IA;
  INT  *ja  = A.JA;
  REAL *val = A.val;
  
  for (is=0; is<A.row; ++is)
    mask[is]=-1;
  // get maximal block size and max nnz per block (max nnz is used later so not to allocate and free during the LU decomposition of the blocks.
  maxbs=0;
  maxbnnz=0;
  for (is=0; is<nblk; ++is) {
    ibl0 = iblock[is];
    ibl1 = iblock[is+1];
    nloc = ibl1-ibl0;
    if(maxbs<nloc) maxbs=nloc; //maxbs = MAX(maxbs, nloc);
    for (i=0; i<nloc; ++i) {
      iblk = ibl0 + i;
      ki   = jblock[iblk];
      mask[ki] = i;
    }
    nnz = 0;
    for (i=0; i<nloc; ++i) {
      iblk = ibl0 + i;
      ki = jblock[iblk];
      iaa = ia[ki];
      iab = ia[ki+1];
      for (kij = iaa; kij<iab; ++kij) {
	kj = ja[kij];
	if(mask[kj] < 0) continue;
	nnz++;
      }
    }
    if(maxbnnz<nnz) maxbnnz=nnz;
    for (i=0; i<nloc; ++i) {
      iblk = ibl0 + i;
      ki   = jblock[iblk];
      mask[ki] = -1;
    }
  }
  Schwarz->maxbs = maxbs;
  Schwarz->maxbnnz = maxbnnz;
  // allocate memory for each sub_block's right hand side and local solution
  Schwarz->xloc1   = dvec_create(maxbs);
  Schwarz->rhsloc1 = dvec_create(maxbs);
  /////////////////////////////////////////////////////////////////
  INT stype=(INT )Schwarz->Schwarz_type;
  if(stype != SCHWARZ_FORWARD &&		\
     stype != SCHWARZ_BACKWARD &&		\
     stype != SCHWARZ_SYMMETRIC &&		\
     stype != SCHWARZ_FORWARD_LOCAL &&	\
     stype != SCHWARZ_BACKWARD_LOCAL &&	\
     stype != SCHWARZ_SYMMETRIC_LOCAL){
    // if we end up here: default is symmetric Schwarz smoother with global LU.
    stype=SCHWARZ_SYMMETRIC;
    Schwarz->Schwarz_type=(SHORT )stype;
  }
  //
  // check if we do global LU.
  if(stype == SCHWARZ_FORWARD ||		\
     stype == SCHWARZ_BACKWARD ||		\
     stype == SCHWARZ_SYMMETRIC){
    // allocate the space for local matrices.
    Schwarz->blk_data = (dCSRmat*)calloc(nblk, sizeof(dCSRmat)); 
    blk=Schwarz->blk_data;
    for (is=0; is<nblk; ++is) {
      ibl0 = iblock[is];
      ibl1 = iblock[is+1];
      nloc = ibl1-ibl0;
      for (i=0; i<nloc; ++i) {
	iblk = ibl0 + i;
	ki = jblock[iblk];
	mask[ki]=i;
      }
      nnz = 0;
      for (i=0; i<nloc; ++i) {
	iblk = ibl0 + i;
	ki = jblock[iblk];
	iaa = ia[ki];
	iab = ia[ki+1];
	for (kij = iaa; kij<iab; ++kij) {
	  kj = ja[kij];
	  if(mask[kj] < 0) continue;
	  nnz++;
	}
      }
      blk[is] = dcsr_create(nloc, nloc, nnz);
      blk[is].IA[0] = 0;
      nnz = 0;
      for (i=0; i<nloc; ++i) {
	iblk = ibl0 + i;
	ki = jblock[iblk];
	iaa = ia[ki];
	iab = ia[ki+1];
	for (kij = iaa; kij<iab; ++kij) {
	  kj = ja[kij];	  
	  if(mask[kj] < 0) continue;
	  blk[is].JA[nnz] = mask[kj]; 
	  blk[is].val[nnz] = val[kij];
	  nnz++;
	}
	blk[is].IA[i+1] = nnz;
      }	
      blk[is].nnz = nnz;	
      // zero the mask so that everyting is as it was
      for (i=0; i<nloc; ++i) {
	iblk = ibl0 + i;
	ki   = jblock[iblk];
	mask[ki] = -1;
      }
    }
  } else {
    // this is only local LU
    Schwarz->blk_data = (dCSRmat*)calloc(1, sizeof(dCSRmat));      
    Schwarz->blk_data[0]=dcsr_create(maxbs,maxbs,maxbnnz);
  }
  return;
}
/***********************************************************************************************/
/**
 * \fn INT Schwarz_setup (Schwarz_data *Schwarz, Schwarz_param *param, ivector *seeds_in)
 *
 * \brief Setup phase for the Schwarz methods
 *
 * \param Schwarz    Pointer to the Schwarz data
 * \param param      Type of the Schwarz method
 *
 * \return           SUCCESS if succeed
 *
 */
INT Schwarz_setup(Schwarz_data *Schwarz,
		  Schwarz_param *param,
		  ivector *seeds_in)
{
  // information about A
  dCSRmat A = Schwarz->A;
  INT n   = A.row;
  INT  block_solver = param->Schwarz_blksolver;
  Schwarz->Schwarz_type=param->Schwarz_type;
  INT  maxlev = ABS(param->Schwarz_maxlvl);
  Schwarz->swzparam = param;
  INT i;
  INT inroot = -10, nsizei = -10, nsizeall = -10, nlvl = 0;
  INT *jb=NULL;
  // data for Schwarz method
  INT nblk=-1;
  INT *iblock = NULL, *jblock = NULL, *mask=NULL;
  void **numeric=NULL;
  INT max_blk_size = 0;
  // return
  INT flag = SUCCESS;
  nsizeall=0;
  // allocate memory
  ivector *seeds = (ivector *)calloc(1, sizeof(ivector));
  mask    = (INT *)calloc(n,sizeof(INT));
  for(i=0;i<n;++i)
    mask[i]=-1;
  iblock  = (INT *)calloc((n+1),sizeof(INT));
  memset(iblock, 0, sizeof(INT)*(n+1));
  jblock  = (INT *)calloc(n,sizeof(INT));
  if(seeds_in){
    seeds=seeds_in;
  } else {
    if (param->Schwarz_maxlvl < 0){
      ivec_alloc(A.row, seeds);
      for (i=0; i<A.row; i++) seeds->val[i] = i;
    } else {
      seeds = sparse_MIS(&A,NULL);
    }    
  }
  iblock[0]=0;
  /*-------------------------------------------*/
  // find the blocks
  /*-------------------------------------------*/
  // first pass: do a maxlev level sets out for each node
  for (i=0; i<seeds->row; i++ ) {
    inroot = seeds->val[i];
    Schwarz_levels(inroot,&A,mask,&nlvl,iblock,jblock,maxlev);
    nsizei=iblock[nlvl];
    if(max_blk_size<nsizei) max_blk_size=nsizei;
    nsizeall+=nsizei;
  }
  /* We only calculated the size of this up to here. So we can reallocate jblock */
  //    jblock = (INT *)realloc(jblock,(nsizeall+n)*sizeof(INT));
  jblock = (INT *)realloc(jblock,nsizeall*sizeof(INT));
  // second pass: redo the same again, but this time we store in jblock
  iblock[0]=0;
  nsizeall=0;
  jb=jblock;
  if(max_blk_size<1) max_blk_size=1;
  INT *iwork=calloc(max_blk_size,sizeof(INT));
  for (i=0;i<seeds->row;i++) {
    inroot = seeds->val[i];
    Schwarz_levels(inroot,&A,mask,&nlvl,iwork,jb,maxlev);
    nsizei=iwork[nlvl];
    iblock[i+1]=iblock[i]+nsizei;
    nsizeall+=nsizei;
    jb+=nsizei;
  }
  free(iwork);iwork=NULL;
  /*
    if there are points not in the jblock, or, equivalently they are not included in any block, then we make each of these to be a block.    
  */
  nblk = seeds->row;
  /*  here we check if there are points that are not covered by any Schwarz block */
  INT j,ibl0,ibl1,nblk_new=nblk;
  //
  //  fprintf(stdout,"iblk[%d]-1=%d,jblk[iblk[%d]-1]=%d",nblk,iblock[nblk]-1,nblk,jblock[iblock[nblk]-1]);fflush(stdout);
  for(i=0;i<nblk;++i){
    ibl0=iblock[i];
    ibl1=iblock[i+1];
    for(j=ibl0;j<ibl1;++j){
      mask[jblock[j]]=i;
    }
  }
  nblk_new=nblk;
  for(i=0;i<n;++i){
    if(mask[i]<0){
      nblk_new++;
    }      
  }
  if(nblk_new>nblk){
    // additional blocks needed to cover all dofs. 
    fprintf(stdout,"\n\n%d nodes not included in any block\n",nblk_new-nblk);fflush(stdout);
    iblock=(INT *)realloc(iblock,(nblk_new+1)*sizeof(INT));
    fprintf(stdout, "Local size: %d %d\n", iblock[nblk]+nblk_new, nblk); fflush(stdout);
    jblock=(INT *)realloc(jblock,(iblock[nblk]+nblk_new)*sizeof(INT));
    nsizeall=iblock[nblk];
    for(i=0;i<n;++i){
      if(mask[i]<0){
	jblock[nsizeall]=i;
	nsizeall++;
	nblk++;
	iblock[nblk]=nsizeall;
      }
    }
  }
  Schwarz->nblk   = nblk;
  Schwarz->iblock = iblock;
  Schwarz->jblock = jblock;
  Schwarz->mask   = mask;
  Schwarz_get_block_matrix(Schwarz);
  dCSRmat *blk = Schwarz->blk_data;
  /*-----------------------------------------------------------------*/
  /* now check for what kind of Schwarz method we have and setup all */
  /*-----------------------------------------------------------------*/
  if(Schwarz->Schwarz_type==SCHWARZ_FORWARD ||	\
     Schwarz->Schwarz_type==SCHWARZ_BACKWARD ||	\
     Schwarz->Schwarz_type==SCHWARZ_SYMMETRIC){
    // Setup for each block solver
    switch (block_solver) {      
      //#if WITH_SUITESPARSE
    case SOLVER_UMFPACK: {
      /* use UMFPACK direct solver on each block; 
	 find the max block size; then
	 Store the blocks if Schwarz_type<10 */
      numeric	= (void**)calloc(nblk, sizeof(void*));
      dCSRmat blk_tran=dcsr_create(Schwarz->maxbs,Schwarz->maxbs,Schwarz->maxbnnz);
      for (i=0; i<nblk; ++i) {
	//      dcsr_create(&blk_tran,blk[i].row, blk[i].col, blk[i].nnz);
	dcsr_transz(&blk[i], NULL, &blk_tran);
	dcsr_cp(&blk_tran, &blk[i]);
	//	dcsr_free(&blk_tran);
	//printf("size of block %d: nrow=%d, nnz=%d\n",i, blk[i].row, blk[i].nnz);
	numeric[i] = hazmath_factorize(&blk[i], 0);
      }
      dcsr_free(&blk_tran);
      break;
    }
      //#endif	
    default:
    break;
      /* do nothing for iterative methods */
    }
  } else {
    numeric=(void **)calloc(1,sizeof(void *));
  }
  /*-------------------------------------------*/
  //  return
  /*-------------------------------------------*/
  Schwarz->numeric = numeric;      
  //Schwarz->Schwarz_type = param->Schwarz_type;// this was set earlier.
  Schwarz->blk_solver = block_solver;
  // clean
  if(!seeds_in){
    ivec_free(seeds);
    if (seeds) free(seeds);
    seeds=NULL;
  }
  fprintf(stdout,"\nSchwarz method setup is done! Find %lld blocks. Maxmium block size = %lld\n",  (long long )nblk,  (long long )max_blk_size);    
  return flag;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
