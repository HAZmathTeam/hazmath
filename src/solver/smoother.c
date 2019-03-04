/*! \file src/solver/smoother.c
 *
 *  Smoothers
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/25/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note  Done cleanup for releasing -- Xiaozhe Hu 03/12/2017
 *
 *  \todo allow different ordering in smoothers -- Xiaozhe Hu
 *  \todo add polynomial smoother, ilu smoothers, and block smoothers -- Xiaozhe Hu
 *
 */

#include "hazmath.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/**
 * \fn void smoother_dcsr_jacobi (dvector *u, const INT i_1, const INT i_n,
 *                                     const INT s, dCSRmat *A, dvector *b, INT L)
 *
 * \brief Jacobi smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1    Starting index
 * \param i_n    Ending index
 * \param s      Increasing step
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 */
void smoother_dcsr_jacobi(dvector *u,
                          const INT i_1,
                          const INT i_n,
                          const INT s,
                          dCSRmat *A,
                          dvector *b,
                          INT L)
{
    const INT    N = ABS(i_n - i_1)+1;
    const INT   *ia=A->IA, *ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    REAL        w = 0.8;  //0.8
    // local variables
    INT i,j,k,begin_row,end_row;

    REAL *t = (REAL *)calloc(N,sizeof(REAL));
    REAL *d = (REAL *)calloc(N,sizeof(REAL));

    while (L--) {
        if (s>0) {
            for (i=i_1;i<=i_n;i+=s) {
                t[i]=bval[i];
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    if (i!=j) t[i]-=aj[k]*uval[j];
                    else d[i]=aj[k];
                }
            }

            for (i=i_1;i<=i_n;i+=s) {
                if (ABS(d[i])>SMALLREAL) uval[i]=(1-w)*uval[i]+ w*t[i]/d[i];
            }
        }
        else {
            for (i=i_1;i>=i_n;i+=s) {
                t[i]=bval[i];
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    if (i!=j) t[i]-=aj[k]*uval[j];
                    else d[i]=aj[k];
                }
            }

            for (i=i_1;i>=i_n;i+=s) {
                if (ABS(d[i])>SMALLREAL) uval[i]=(1-w)*uval[i]+ w*t[i]/d[i];
            }
        }

    } // end while

    free(t);
    free(d);

    return;
}

/**
 * \fn void smoother_dcsr_gs(dvector *u, const INT i_1, const INT i_n,
 *                                 const INT s, dCSRmat *A, dvector *b, INT L)
 *
 * \brief Gauss-Seidel smoother
 *
 * \param u    Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1  Starting index
 * \param i_n  Ending index
 * \param s    Increasing step
 * \param A    Pointer to dBSRmat: the coefficient matrix
 * \param b    Pointer to dvector: the right hand side
 * \param L    Number of iterations
 *
 */
void smoother_dcsr_gs(dvector *u,
                      const INT i_1,
                      const INT i_n,
                      const INT s,
                      dCSRmat *A,
                      dvector *b,
                      INT L)
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;

    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  t,d=0.0;

    if (s > 0) {
        while (L--) {
            for (i=i_1;i<=i_n;i+=s) {
                t = bval[i];
                begin_row=ia[i],end_row=ia[i+1];

#if DIAGONAL_PREF // diagonal first
                d=aj[begin_row];
                if (ABS(d)>SMALLREAL) {
                    for (k=begin_row+1;k<end_row;++k) {
                        j=ja[k];
                        t-=aj[k]*uval[j];
                    }
                    uval[i]=t/d;
                }
#else // general order
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    if (i!=j)
                        t-=aj[k]*uval[j];
                    else if (ABS(aj[k])>SMALLREAL) d=1.e+0/aj[k];
                }
                uval[i]=t*d;
#endif
            } // end for i
        } // end while

    } // if s
    else {

        while (L--) {
            for (i=i_1;i>=i_n;i+=s) {
                t=bval[i];
                begin_row=ia[i],end_row=ia[i+1];
#if DIAGONAL_PREF // diagonal first
                d=aj[begin_row];
                if (ABS(d)>SMALLREAL) {
                    for (k=begin_row+1;k<end_row;++k) {
                        j=ja[k];
                        t-=aj[k]*uval[j];
                    }
                    uval[i]=t/d;
                }
#else // general order
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    if (i!=j)
                        t-=aj[k]*uval[j];
                    else if (ABS(aj[k])>SMALLREAL) d=1.0/aj[k];
                }
                uval[i]=t*d;
#endif
            } // end for i
        } // end while
    } // end if
    return;
}


/**
 * \fn void smoother_dcsr_sgs(dvector *u, dCSRmat *A, dvector *b, INT L)
 *
 * \brief Symmetric Gauss-Seidel smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 *
 */
void smoother_dcsr_sgs(dvector *u,
                       dCSRmat *A,
                       dvector *b,
                       INT L)
{
    const INT    nm1=b->row-1;
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;

    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  t,d=0;

    while (L--) {
        // forward sweep
        for (i=0;i<=nm1;++i) {
            t=bval[i];
            begin_row=ia[i], end_row=ia[i+1];
            for (k=begin_row;k<end_row;++k) {
                j=ja[k];
                if (i!=j) t-=aj[k]*uval[j];
                else d=aj[k];
            } // end for k
            if (ABS(d)>SMALLREAL) uval[i]=t/d;
        } // end for i

        // backward sweep
        for (i=nm1-1;i>=0;--i) {
            t=bval[i];
            begin_row=ia[i], end_row=ia[i+1];
            for (k=begin_row;k<end_row;++k) {
                j=ja[k];
                if (i!=j) t-=aj[k]*uval[j];
                else d=aj[k];
            } // end for k
            if (ABS(d)>SMALLREAL) uval[i]=t/d;
        } // end for i

    } // end while

    return;
}

/**
 * \fn void smoother_dcsr_sor(dvector *u, const INT i_1, const INT i_n, const INT s,
 *                                 dCSRmat *A, dvector *b, INT L, const REAL w)
 *
 * \brief SOR smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1    Starting index
 * \param i_n    Ending index
 * \param s      Increasing step
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 * \param w      Over-relaxation weight
 *
 *
 */
void smoother_dcsr_sor(dvector *u,
                       const INT i_1,
                       const INT i_n,
                       const INT s,
                       dCSRmat *A,
                       dvector *b,
                       INT L,
                       const REAL w)
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;

    // local variables
    INT    i,j,k,begin_row,end_row;
    REAL   t, d=0;


    while (L--) {
        if (s>0) {
            for (i=i_1; i<=i_n; i+=s) {
                t=bval[i];
                begin_row=ia[i], end_row=ia[i+1];
                for (k=begin_row; k<end_row; ++k) {
                    j=ja[k];
                    if (i!=j)
                        t-=aj[k]*uval[j];
                    else
                        d=aj[k];
                }
                if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
            }
        }
        else {
            for (i=i_1;i>=i_n;i+=s) {
                t=bval[i];
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    if (i!=j)
                        t-=aj[k]*uval[j];
                    else
                        d=aj[k];
                }
                if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
            }
        }
    }  // end while

    return;
}

/**
 * \fn void smoother_dcsr_L1diag(dvector *u, const INT i_1, const INT i_n, const INT s,
 *                                     dCSRmat *A, dvector *b, INT L)
 *
 * \brief Diagonal scaling (using L1 norm) smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1    Starting index
 * \param i_n    Ending index
 * \param s      Increasing step
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 */
void smoother_dcsr_L1diag(dvector *u,
                          const INT i_1,
                          const INT i_n,
                          const INT s,
                          dCSRmat *A,
                          dvector *b,
                          INT L)
{
    const INT    N = ABS(i_n - i_1)+1;
    const INT   *ia=A->IA, *ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;

    // local variables
    INT   i,j,k,begin_row,end_row;

    // Checks should be outside of for; t,d can be allocated before calling!!! --Chensong
    REAL *t = (REAL *)calloc(N,sizeof(REAL));
    REAL *d = (REAL *)calloc(N,sizeof(REAL));

    while (L--) {
        if (s>0) {
            for (i=i_1;i<=i_n;i+=s) {
                t[i]=bval[i]; d[i]=0.0;
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    t[i]-=aj[k]*uval[j];
                    d[i]+=ABS(aj[k]);
                }
            }

            for (i=i_1;i<=i_n;i+=s) {
                if (ABS(d[i])>SMALLREAL) u->val[i]+=t[i]/d[i];
            }
        }
        else {
            for (i=i_1;i>=i_n;i+=s) {
                t[i]=bval[i];d[i]=0.0;
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    t[i]-=aj[k]*uval[j];
                    d[i]+=ABS(aj[k]);
                }
            }

            for (i=i_1;i>=i_n;i+=s) {
                if (ABS(d[i])>SMALLREAL) u->val[i]+=t[i]/d[i];
            }
        }

    } // end while

    free(t);
    free(d);

    return;
}

/**
 * \fn void smoother_dcsr_Schwarz_forward (Schwarz_data  *Schwarz,
 *                                         Schwarz_param *param,
 *                                         dvector *x, dvector *b)
 *
 * \brief Schwarz smoother: forward sweep
 *
 * \param Schwarz Pointer to the Schwarz data
 * \param param   Pointer to the Schwarz parameter
 * \param x       Pointer to solution vector
 * \param b       Pointer to right hand
 *
 * \note Needs improvment -- Xiaozhe
 */
void smoother_dcsr_Schwarz_forward (Schwarz_data  *Schwarz,
                                    Schwarz_param *param,
                                    dvector       *x,
                                    dvector       *b)
{
    INT i, j, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;

    // Schwarz partition
    INT  nblk = Schwarz->nblk;
    dCSRmat *blk = Schwarz->blk_data;
    INT  *iblock = Schwarz->iblock;
    INT  *jblock = Schwarz->jblock;
    INT  *mask   = Schwarz->mask;
    INT  block_solver = param->Schwarz_blksolver;


    // Schwarz data
    dCSRmat A = Schwarz->A;
    INT *ia = A.IA;
    INT *ja = A.JA;
    REAL *val = A.val;

    // Local solution and right hand vectors
    dvector rhs = Schwarz->rhsloc1;
    dvector u   = Schwarz->xloc1;

#if WITH_SUITESPARSE
    void **numeric = Schwarz->numeric;
#endif

    for (is=0; is<nblk; ++is) {
        // Form the right hand of eack block
        ibl0 = iblock[is];
        ibl1 = iblock[is+1];
        nloc = ibl1-ibl0;
        for (i=0; i<nloc; ++i ) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            mask[ki] = i+1;
        }

        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki = jblock[iblk];
            rhs.val[i] = b->val[ki];
            iaa = ia[ki]-1;
            iab = ia[ki+1]-1;
            for (kij = iaa; kij<iab; ++kij) {
                kj = ja[kij]-1;
                j  = mask[kj];
                if(j == 0) {
                    rhs.val[i] -= val[kij]*x->val[kj];
                }
            }
        }

        // Solve each block
        switch (block_solver) {

#if WITH_SUITESPARSE
            case SOLVER_UMFPACK: {
                /* use UMFPACK direct solver on each block */
                umfpack_solve(&blk[is], &rhs, &u, numeric[is], 0);
                break;
            }
#endif
            default:
                /* use iterative solver on each block */
                u.row = blk[is].row;
                rhs.row = blk[is].row;
                dvec_set(u.row, &u, 0);
                dcsr_pvgmres(&blk[is], &rhs, &u, NULL, 1e-8, 20, 20, 1, 0);
        }

        //zero the mask so that everyting is as it was
        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            mask[ki] = 0;
            x->val[ki] = u.val[i];
        }
    }
}

/**
 * \fn void smoother_dcsr_Schwarz_backward (Schwarz_data  *Schwarz,
 *                                          Schwarz_param *param,
 *                                          dvector *x, dvector *b)
 *
 * \brief Schwarz smoother: backward sweep
 *
 * \param Schwarz Pointer to the Schwarz data
 * \param param   Pointer to the Schwarz parameter
 * \param x       Pointer to solution vector
 * \param b       Pointer to right hand
 *
 * \note Needs improvment -- Xiaozhe
 */
void smoother_dcsr_Schwarz_backward (Schwarz_data *Schwarz,
                                     Schwarz_param *param,
                                     dvector *x,
                                     dvector *b)
{
    INT i, j, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;

    // Schwarz partition
    INT  nblk = Schwarz->nblk;
    dCSRmat *blk = Schwarz->blk_data;
    INT  *iblock = Schwarz->iblock;
    INT  *jblock = Schwarz->jblock;
    INT  *mask   = Schwarz->mask;
    INT  block_solver = param->Schwarz_blksolver;


    // Schwarz data
    dCSRmat A = Schwarz->A;
    INT *ia = A.IA;
    INT *ja = A.JA;
    REAL *val = A.val;

    // Local solution and right hand vectors
    dvector rhs = Schwarz->rhsloc1;
    dvector u   = Schwarz->xloc1;

#if WITH_SUITESPARSE
    void **numeric = Schwarz->numeric;
#endif

    for (is=nblk-1; is>=0; --is) {
        // Form the right hand of eack block
        ibl0 = iblock[is];
        ibl1 = iblock[is+1];
        nloc = ibl1-ibl0;
        for (i=0; i<nloc; ++i ) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            mask[ki] = i+1;
        }

        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki = jblock[iblk];
            rhs.val[i] = b->val[ki];
            iaa = ia[ki]-1;
            iab = ia[ki+1]-1;
            for (kij = iaa; kij<iab; ++kij) {
                kj = ja[kij]-1;
                j  = mask[kj];
                if(j == 0) {
                    rhs.val[i] -= val[kij]*x->val[kj];
                }
            }
        }

        // Solve each block
        switch (block_solver) {

#if WITH_SUITESPARSE
            case SOLVER_UMFPACK: {
                /* use UMFPACK direct solver on each block */
                umfpack_solve(&blk[is], &rhs, &u, numeric[is], 0);
                break;
            }
#endif
            default:
                /* use iterative solver on each block */
                rhs.row = blk[is].row;
                u.row   = blk[is].row;
                dvec_set(u.row, &u, 0);
                dcsr_pvgmres (&blk[is], &rhs, &u, NULL, 1e-8, 20, 20, 1, 0);
        }

        //zero the mask so that everyting is as it was
        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki   = jblock[iblk];
            mask[ki] = 0;
            x->val[ki] = u.val[i];
        }
    }
}

/**
 * \fn void smoother_bdcsr_jacobi (dvector *u, const INT s, block_dCSRmat *A, dvector *b, INT L)
 *
 * \brief Jacobi smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param s      Increasing step
 * \param A      Pointer to block_dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 */
void smoother_bdcsr_jacobi(dvector *u,
                          const INT s,
                          block_dCSRmat *A,
                          dvector *b,
                          INT L)
{
    // Sub-vectors
    dvector utemp;
    dvector btemp;
    // Smooth
    INT i, istart;
    INT row;
    istart = 0;
    for(i=0; i<A->brow; i++){
        row = A->blocks[i+i*A->brow]->row;
        // Get sub-vectors
        utemp.row = row;
        utemp.val = u->val+istart;
        btemp.row = row;
        btemp.val = b->val+istart;
        // Call jacobi on specific block
        smoother_dcsr_jacobi(&utemp,0,row-1,s,A->blocks[i+i*A->brow],&btemp,L);
        // Move to next block
        istart += row;
    }

}

/**
 * \fn void smoother_block_setup( MG_blk_data *mgl, AMG_param *param)
 *
 * \brief Setup of block smoothers
 *
 * \param mgl       Pointer to MG_blk_data
 * \param param     Pointer to AMG_param
 *
 */

void smoother_block_setup( MG_blk_data *bmgl, AMG_param *param)
{
    // TODO: instead of param, may want a variable in bmgl that tracks params for each block
    // (since we don't want to create the Schwarz for each block, only those that will use it);
    INT blk, lvl, i;
    INT brow = bmgl[0].A.brow;
    SHORT max_levels = param->max_levels;
    fespace FE_fake;
    for(i=0;i<max_levels;i++){
      bmgl[i].mgl = (AMG_data **)calloc(brow,sizeof(AMG_data *));
    }
    for(i=1;i<max_levels;i++){
      bmgl[i].A_diag = (dCSRmat *)calloc(3,sizeof(dCSRmat));
    }

    Schwarz_param swzparam;

    for( blk=0; blk<brow; blk++ ){
      printf("Start %d\n",blk);
      // Initialize AMG for diagonal blocks
      bmgl[0].mgl[blk] = amg_data_create(max_levels);

      dcsr_alloc( bmgl[0].A_diag[blk].row, bmgl[0].A_diag[blk].row, bmgl[0].A_diag[blk].nnz, &bmgl[0].mgl[blk][0].A);
      dcsr_cp( &(bmgl[0].A_diag[blk]), &bmgl[0].mgl[blk][0].A );

      FE_fake.dirichlet = bmgl[0].dirichlet_blk[blk];
      FE_fake.ndof = bmgl[0].A_diag[blk].row;
      dcsr_shift( &bmgl[0].mgl[blk][0].A,  1);
      printf("Eliminating dirichlet BC for A_diag\n");
      eliminate_DirichletBC(NULL, &FE_fake , &bmgl[0].fine_level_mesh, NULL, &(bmgl[0].mgl[blk][0].A),0.0);
      printf("Eliminated dirichlet BC for A_diag\n");
      dcsr_shift( &bmgl[0].mgl[blk][0].A, -1);

      // Initialize Schwarz parameters
      bmgl->Schwarz_levels = param->Schwarz_levels;
      //if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_mmsize = param->Schwarz_mmsize;
        swzparam.Schwarz_maxlvl = param->Schwarz_maxlvl;
        swzparam.Schwarz_type   = param->Schwarz_type;
        //swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
        swzparam.Schwarz_blksolver = 32;
      //}

      // Fill in levels
      for( lvl=0; lvl<max_levels-1; lvl++ ){
        printf("\tLEVEL %d FILLING\n",lvl);
        // Copy block mgl R and P into mgl
        bmgl[0].mgl[blk][lvl].R.row = bmgl[lvl].R.blocks[blk+blk*brow]->row;
        bmgl[0].mgl[blk][lvl].R.col = bmgl[lvl].R.blocks[blk+blk*brow]->col;
        bmgl[0].mgl[blk][lvl].R.nnz = bmgl[lvl].R.blocks[blk+blk*brow]->nnz;
        bmgl[0].mgl[blk][lvl].R.val = bmgl[lvl].R.blocks[blk+blk*brow]->val;
        bmgl[0].mgl[blk][lvl].R.IA = bmgl[lvl].R.blocks[blk+blk*brow]->IA;
        bmgl[0].mgl[blk][lvl].R.JA = bmgl[lvl].R.blocks[blk+blk*brow]->JA;

        bmgl[0].mgl[blk][lvl].P.row = bmgl[lvl].P.blocks[blk+blk*brow]->row;
        bmgl[0].mgl[blk][lvl].P.col = bmgl[lvl].P.blocks[blk+blk*brow]->col;
        bmgl[0].mgl[blk][lvl].P.nnz = bmgl[lvl].P.blocks[blk+blk*brow]->nnz;
        bmgl[0].mgl[blk][lvl].P.val = bmgl[lvl].P.blocks[blk+blk*brow]->val;
        bmgl[0].mgl[blk][lvl].P.IA = bmgl[lvl].P.blocks[blk+blk*brow]->IA;
        bmgl[0].mgl[blk][lvl].P.JA = bmgl[lvl].P.blocks[blk+blk*brow]->JA;

        // Compute RAP for mgl
        //dcsr_rap( &bmgl[0].mgl[blk][lvl].R, &bmgl[0].mgl[blk][lvl].A, &bmgl[0].mgl[blk][lvl].P, &bmgl[0].mgl[blk][lvl+1].A);
        dcsr_rap( &bmgl[0].mgl[blk][lvl].R, &bmgl[lvl].A_diag[blk], &bmgl[0].mgl[blk][lvl].P, &bmgl[lvl+1].A_diag[blk]);

        dcsr_alloc( bmgl[lvl+1].A_diag[blk].row, bmgl[lvl+1].A_diag[blk].row, bmgl[lvl+1].A_diag[blk].nnz, &bmgl[0].mgl[blk][lvl+1].A);
        dcsr_cp( &(bmgl[lvl+1].A_diag[blk]), &bmgl[0].mgl[blk][lvl+1].A );
        FE_fake.dirichlet = bmgl[lvl+1].dirichlet_blk[blk];
        FE_fake.ndof      = bmgl[lvl+1].A_diag[blk].row;
        dcsr_shift( &bmgl[0].mgl[blk][lvl+1].A,  1);
        eliminate_DirichletBC(NULL, &FE_fake , &bmgl[lvl+1].fine_level_mesh, NULL, &(bmgl[0].mgl[blk][lvl+1].A),0.0);
        dcsr_shift( &bmgl[0].mgl[blk][lvl+1].A, -1);



        // setup total level number and current level
        bmgl[0].mgl[blk][lvl].num_levels = max_levels;
        bmgl[0].mgl[blk][lvl].cycle_type = param->cycle_type;
        bmgl[0].mgl[blk][lvl].b = dvec_create(bmgl[0].mgl[blk][lvl].A.row);
        bmgl[0].mgl[blk][lvl].x = dvec_create(bmgl[0].mgl[blk][lvl].A.row);
        bmgl[0].mgl[blk][lvl].w = dvec_create(2*bmgl[0].mgl[blk][lvl].A.row);

        /*-- Setup Schwarz smoother if necessary --*/
        //if( lvl < param->Schwarz_levels && blk == 1 ) { // && use_Schwarz[blk] == 1) or something similar
        if( blk == 1 ) { // && use_Schwarz[blk] == 1) or something similar
          printf("\tcalling Schwarz setup\n");
          bmgl[0].mgl[blk][lvl].Schwarz.A = dcsr_sympat( &bmgl[0].mgl[blk][lvl].A );
          dcsr_shift(&(bmgl[0].mgl[blk][lvl].Schwarz.A),1);
          Schwarz_setup_geometric( &bmgl[0].mgl[blk][lvl].Schwarz, &swzparam, bmgl[lvl].fine_level_mesh);
        }
      }

      // Set up Coarse level solve 
      if( lvl == max_levels-1 ){//TODO: is this correct level?
        printf("Coarse Level Solver Stuff...\n");
        switch ( param->coarse_solver ){
#if WITH_SUITESPARSE
          case SOLVER_UMFPACK: {
              dCSRmat Ac_tran;
              dcsr_trans(&bmgl[0].mgl[blk][lvl].A, &Ac_tran);
              dcsr_cp(&Ac_tran, &bmgl[0].mgl[blk][lvl].A);
              dcsr_free(&Ac_tran);
              bmgl[0].mgl[blk][lvl].Numeric = umfpack_factorize(&bmgl[0].mgl[blk][lvl].A, 0);
              break;}
#endif
          default:
              // Do nothing!
              break;
        }
      }
      // Propigate mgl across bmgl levels
      for( lvl=1; lvl<max_levels-1; lvl++){
        bmgl[lvl].mgl[blk] = &(bmgl[0].mgl[blk][lvl]);// maybe??
        printf("\tLVL %d propigated across bmgl levels\n",lvl);
      }
    }// blk
    return;
}

/**
 * \fn void smoother_block_biot_3field( MG_blk_data *mgl, AMG_param *param)
 *
 * \brief Setup of block smoothers
 *
 * \param mgl       Pointer to MG_blk_data
 * \param param     Pointer to AMG_param
 *
 */
void smoother_block_biot_3field( const INT lvl, MG_blk_data *bmgl, AMG_param *param, INT pre_post)
{
    printf("\tBeginning Block smoother, lvl=%d\n",lvl);
    param->presmooth_iter = 2;
    printf("\t\t| using %2d iterations\n",param->presmooth_iter);
    // Initialize
    INT n0, n1, n2, ntot;
    // Sub-vectors
    dvector x0, x1, x2;
    dvector b0, b1, b2;

    n0 = bmgl[lvl].mgl[0][0].A.row;
    n1 = bmgl[lvl].mgl[1][0].A.row;
    n2 = bmgl[lvl].mgl[2][0].A.row;
    ntot = n0+n1+n2;

    x0.row = n0;
    x0.val = bmgl[lvl].x.val;
    b0.row = n0;
    b0.val = bmgl[lvl].b.val;

    x1.row = n1;
    x1.val = &(bmgl[lvl].x.val[n0]);
    b1.row = n1;
    b1.val = &(bmgl[lvl].b.val[n0]);

    x2.row = n2;
    x2.val = &(bmgl[lvl].x.val[n0+n1]);
    b2.row = n2;
    b2.val = &(bmgl[lvl].b.val[n0+n1]);

    // Block 0: P1 + Bubble
    smoother_dcsr_sgs(&x0, &(bmgl[lvl].mgl[0][0].A), &b0, param->presmooth_iter);
    //directsolve_UMF(&(bmgl[lvl].mgl[0][0].A), &b0, &x0, 10);
    printf("\tBlock 0 done...\n");

    // Block 1: RT0
//    Schwarz_param swzparam;
//    swzparam.Schwarz_blksolver = bmgl[lvl].mgl[1][0].Schwarz.blk_solver;
//    if(pre_post == 1){
//      smoother_dcsr_Schwarz_forward( &(bmgl[lvl].mgl[1][0].Schwarz), &swzparam, &x1, &b1);
//      smoother_dcsr_Schwarz_backward(&(bmgl[lvl].mgl[1][0].Schwarz), &swzparam, &x1, &b1);
//    } else if (pre_post == 2){
//      smoother_dcsr_Schwarz_backward(&(bmgl[lvl].mgl[1][0].Schwarz), &swzparam, &x1, &b1);
//      smoother_dcsr_Schwarz_forward( &(bmgl[lvl].mgl[1][0].Schwarz), &swzparam, &x1, &b1);
//    }
    directsolve_UMF(&(bmgl[lvl].mgl[1][0].A), &b1, &x1, 10);
    printf("\tBlock 1 done...\n");

    // Block 2: P0
    smoother_dcsr_sgs(&x2, &(bmgl[lvl].mgl[2][0].A), &b2, param->presmooth_iter);
    //directsolve_UMF(&(bmgl[lvl].mgl[2][0].A), &b2, &x2, 10);
    printf("\tBlock 2 done...\n");
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
