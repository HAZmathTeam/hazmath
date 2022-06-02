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

static void Schwarz_levels (INT, dCSRmat *, INT *, INT *, INT *, INT *, INT);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/***********************************************************************************************/
/**
 * \fn void Schwarz_get_block_matrix (Schwarz_data *Schwarz, INT nblk,
 *                                    INT *iblock, INT *jblock, INT *mask)
 *
 * \brief Form Schwarz partition data
 *
 * \param Schwarz Pointer to the Schwarz data
 * \param nblk    Number of partitions
 * \param iblock  Pointer to number of vertices on each level
 * \param jblock  Pointer to vertices of each level
 * \param mask    Pointer to flag array
 *
 * \note  This needs to be rewritten -- Xiaozhe
 *
 * \note  Done cleanup for releasing -- Xiaozhe Hu 08/28/2021
 *
 */
void Schwarz_get_block_matrix (Schwarz_data *Schwarz,
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
            mask[ki] = i+1;  // The +1 -Peter
        }

        //printf("is=%d\n", is);
        //printf("nloc=%d, count=%d\n", nloc, count);
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
                    blk[is].JA[nnz] = j-1; // The -1 corresponds with +1 above. -Peter
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

/***********************************************************************************************/
/**
 * \fn INT Schwarz_setup (Schwarz_data *Schwarz, Schwarz_param *param)
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
                  Schwarz_param *param)
{
    // information about A
    dCSRmat A = Schwarz->A;
    INT n   = A.row;

    INT  block_solver = param->Schwarz_blksolver;
    INT  maxlev = ABS(param->Schwarz_maxlvl);
    Schwarz->swzparam = param;

    printf("param->Schwarz_maxlvl = %d\n", param->Schwarz_maxlvl);

    // local variables
    INT i;
    INT inroot = -10, nsizei = -10, nsizeall = -10, nlvl = 0;
    INT *jb=NULL;
    // data for Schwarz method
    INT nblk;
    INT *iblock = NULL, *jblock = NULL, *mask = NULL, *maxa = NULL;
    INT max_blk_size = 0;

    // return
    INT flag = SUCCESS;

    // allocate memory
    maxa    = (INT *)calloc(n,sizeof(INT));
    mask    = (INT *)calloc(n,sizeof(INT));
    iblock  = (INT *)calloc(n,sizeof(INT));
    jblock  = (INT *)calloc(n,sizeof(INT));

    nsizeall=0;
    memset(mask,   0, sizeof(INT)*n);
    memset(iblock, 0, sizeof(INT)*n);
    memset(maxa,   0, sizeof(INT)*n);

    maxa[0]=0;

    // select root nodes.
    ivector *MaxIndSet;
    //if (param->Schwarz_maxlvl < 0){
        MaxIndSet = (ivector *)calloc(1, sizeof(ivector));
        ivec_alloc(n, MaxIndSet);
        for (i=0; i<A.row; i++) MaxIndSet->val[i] = i;
        //maxlev = 1;
    //}
    //else {
    //    MaxIndSet = sparse_MIS(&A,NULL);
    //}

    /*-------------------------------------------*/
    // find the blocks
    /*-------------------------------------------*/
    // first pass: do a maxlev level sets out for each node
    for ( i = 0; i < MaxIndSet->row; i++ ) {
        inroot = MaxIndSet->val[i];
        Schwarz_levels(inroot,&A,mask,&nlvl,maxa,jblock,maxlev);
        nsizei=maxa[nlvl];
        max_blk_size = MAX(max_blk_size, nsizei);
        nsizeall+=nsizei;
    }

    /* We only calculated the size of this up to here. So we can reallocate jblock */
    jblock = (INT *)realloc(jblock,(nsizeall+n)*sizeof(INT));

    // second pass: redo the same again, but this time we store in jblock
    maxa[0]=0;
    iblock[0]=0;
    nsizeall=0;
    jb=jblock;
    for (i=0;i<MaxIndSet->row;i++) {
        inroot = MaxIndSet->val[i];
        Schwarz_levels(inroot,&A,mask,&nlvl,maxa,jb,maxlev);
        nsizei=maxa[nlvl];
        iblock[i+1]=iblock[i]+nsizei;
        nsizeall+=nsizei;
        jb+=nsizei;
    }
    nblk = MaxIndSet->row;

    /*-------------------------------------------*/
    //  LU decomposition of blocks
    /*-------------------------------------------*/
    memset(mask, 0, sizeof(INT)*n);
    Schwarz->blk_data = (dCSRmat*)calloc(nblk, sizeof(dCSRmat));
    Schwarz_get_block_matrix(Schwarz, nblk, iblock, jblock, mask);

    // Setup for each block solver
    switch (block_solver) {

#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            /* use UMFPACK direct solver on each block */
            dCSRmat *blk = Schwarz->blk_data;
            void **numeric	= (void**)calloc(nblk, sizeof(void*));
            dCSRmat blk_tran;
            for (i=0; i<nblk; ++i) {
                dcsr_alloc(blk[i].row, blk[i].col, blk[i].nnz, &blk_tran);
                dcsr_transz(&blk[i], NULL, &blk_tran);
                dcsr_cp(&blk_tran, &blk[i]);
                dcsr_free(&blk_tran);
                //printf("size of block %d: nrow=%d, nnz=%d\n",i, blk[i].row, blk[i].nnz);
                numeric[i] = umfpack_factorize(&blk[i], 0);
            }
            Schwarz->numeric = numeric;

            break;
        }
#endif

        default: {
            /* do nothing for iterative methods */
        }
    }

    /*-------------------------------------------*/
    //  return
    /*-------------------------------------------*/
    Schwarz->nblk   = nblk;
    Schwarz->iblock = iblock;
    Schwarz->jblock = jblock;
    Schwarz->mask   = mask;
    Schwarz->maxa   = maxa;
    Schwarz->Schwarz_type = param->Schwarz_type;
    Schwarz->blk_solver = param->Schwarz_blksolver;

    printf("Schwarz method setup is done! Find %d blocks. Maxmium block size = %d\n",nblk, max_blk_size);

    // clean
    ivec_free(MaxIndSet);
    if (MaxIndSet) free(MaxIndSet);
    MaxIndSet = NULL;

    return flag;
}

/***********************************************************************************************/
/**
 * \fn INT Schwarz_setup_with_seeds(Schwarz_data *Schwarz, Schwarz_param *param, ivector *seeds)
 *
 * \brief Setup phase for the Schwarz methods with user provided seeds for constructing blocks
 *
 * \param Schwarz    Pointer to the Schwarz data
 * \param param      Type of the Schwarz method
 * \param seeds      Pointer to the seeds for each block
 *
 * \return           SUCCESS if succeed
 *
 * \author           Xiaozhe Hu
 *
 */
INT Schwarz_setup_with_seeds(Schwarz_data *Schwarz,
                             Schwarz_param *param,
                             ivector *seeds)
{
    // information about A
    dCSRmat A = Schwarz->A;
    INT n   = A.row;

    INT  block_solver = param->Schwarz_blksolver;
    INT  maxlev = ABS(param->Schwarz_maxlvl);
    Schwarz->swzparam = param;

    // local variables
    INT i;
    INT inroot = -10, nsizei = -10, nsizeall = -10, nlvl = 0;
    INT *jb=NULL;
    // data for Schwarz method
    INT nblk;
    INT *iblock = NULL, *jblock = NULL, *mask = NULL, *maxa = NULL;
    INT max_blk_size = 0;

    // return
    INT flag = SUCCESS;

    // allocate memory
    maxa    = (INT *)calloc(n,sizeof(INT));
    mask    = (INT *)calloc(n,sizeof(INT));
    iblock  = (INT *)calloc(n,sizeof(INT));
    jblock  = (INT *)calloc(n,sizeof(INT));

    nsizeall=0;
    memset(mask,   0, sizeof(INT)*n);
    memset(iblock, 0, sizeof(INT)*n);
    memset(maxa,   0, sizeof(INT)*n);

    maxa[0]=0;

    /*-------------------------------------------*/
    // find the blocks
    /*-------------------------------------------*/
    // first pass: do a maxlev level sets out for each node
    for (i=0; i<seeds->row; i++ ) {
        inroot = seeds->val[i];
        Schwarz_levels(inroot,&A,mask,&nlvl,maxa,jblock,maxlev);
        nsizei=maxa[nlvl];
        max_blk_size = MAX(max_blk_size, nsizei);
        nsizeall+=nsizei;
    }

    /* We only calculated the size of this up to here. So we can reallocate jblock */
    jblock = (INT *)realloc(jblock,(nsizeall+n)*sizeof(INT));

    // second pass: redo the same again, but this time we store in jblock
    maxa[0]=0;
    iblock[0]=0;
    nsizeall=0;
    jb=jblock;
    for (i=0;i<seeds->row;i++) {
        inroot = seeds->val[i];
        Schwarz_levels(inroot,&A,mask,&nlvl,maxa,jb,maxlev);
        nsizei=maxa[nlvl];
        iblock[i+1]=iblock[i]+nsizei;
        nsizeall+=nsizei;
        jb+=nsizei;
    }
    nblk = seeds->row;

    /*-------------------------------------------*/
    //  LU decomposition of blocks
    /*-------------------------------------------*/
    memset(mask, 0, sizeof(INT)*n);
    Schwarz->blk_data = (dCSRmat*)calloc(nblk, sizeof(dCSRmat));
    Schwarz_get_block_matrix(Schwarz, nblk, iblock, jblock, mask);

    // Setup for each block solver
    switch (block_solver) {

#if WITH_SUITESPARSE
        case SOLVER_UMFPACK: {
            /* use UMFPACK direct solver on each block */
            dCSRmat *blk = Schwarz->blk_data;
            void **numeric	= (void**)calloc(nblk, sizeof(void*));
            dCSRmat blk_tran;
            for (i=0; i<nblk; ++i) {
                dcsr_alloc(blk[i].row, blk[i].col, blk[i].nnz, &blk_tran);
                dcsr_transz(&blk[i], NULL, &blk_tran);
                dcsr_cp(&blk_tran, &blk[i]);
                dcsr_free(&blk_tran);
                //printf("size of block %d: nrow=%d, nnz=%d\n",i, blk[i].row, blk[i].nnz);
                numeric[i] = umfpack_factorize(&blk[i], 0);
            }
            Schwarz->numeric = numeric;

            break;
        }
#endif

        default: {
            /* do nothing for iterative methods */
        }
    }

    /*-------------------------------------------*/
    //  return
    /*-------------------------------------------*/
    Schwarz->nblk   = nblk;
    Schwarz->iblock = iblock;
    Schwarz->jblock = jblock;
    Schwarz->mask   = mask;
    Schwarz->maxa   = maxa;
    Schwarz->Schwarz_type = param->Schwarz_type;
    Schwarz->blk_solver = param->Schwarz_blksolver;

    printf("Schwarz method setup is done! Find %d blocks. Maxmium block size = %d\n",nblk, max_blk_size);

    return flag;
}


/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

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
    }
    else {
        // input node as root node (level 0)
        lvl = 0;
        jblock[0] = inroot;
        lvlend = 0;
        nsize  = 1;
        // mark root node
        mask[inroot] = 1; //??

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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
