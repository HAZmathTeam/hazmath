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
 * \fn void smoother_bdcsr_bsr (dvector *u, const INT s, block_dCSRmat *A, dvector *b, INT L)
 *
 * \brief BSR
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param s      Increasing step
 * \param A      Pointer to block_dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 * [ alpha C  B^T ] [v] = [d]
 * [   B      0   ] [q] = [e]
 *
 */
void smoother_bdcsr_bsr(dvector *u,
                        dvector *b,
                        REAL alpha,
                        block_dCSRmat *A,
                        dCSRmat *C,
                        dCSRmat *B,
                        dCSRmat *M,
                        INT L)
{
printf("Beginning BSR\n");
    // Local variables
    INT i;
    INT n0, n1;
    dvector r;
    dvector rhs;
    dvector temp;
    dvector v;
    dvector q;
    dvector d;
    dvector e;

    dcsr_axm(B, -1.0); /////////------------ TODO: biot not symmetric. need to deal with this.
    dCSRmat BT;
    dcsr_trans(B, &BT);

    dCSRmat Cinv;
    dcsr_alloc( C->row, C->col, C->nnz, &Cinv);
    dcsr_cp( C, &Cinv);
    for(i=0; i<Cinv.row; i++){ // Assuming C is diag
      if( ABS( Cinv.val[i] ) > SMALLREAL ) Cinv.val[i] = 1.0 / Cinv.val[i];
    }

    dCSRmat S; // S = BC^{-1}B^T
    dCSRmat Stemp;
    dcsr_rap( B, &Cinv, &BT, &S);
    if( M ){
      dcsr_alloc( S.row, S.col, S.nnz, &Stemp);
      dcsr_cp( &S, &Stemp);
      dcsr_free( &S );
      dcsr_add( &Stemp, 1.0, M, alpha, &S); // ( S - alpha M ) Check sign here
    }

    // fill local var
    n0 = BT.row;
    n1 = B->row;

    // solver loop
    dvec_alloc( b->row, &r);
    dvec_cp( b, &r);
    bdcsr_aAxpy( -1.0, A, u->val, r.val); // r = b - Au

    d.row = n0;
    d.val = r.val;
    e.row = n1;
    e.val = r.val + n0;

    // (B C^{-1} B^T) q = B C^{-1} d - alpha e
    dvec_alloc( e.row, &rhs);
    //dvec_axpy( -alpha, &e, &rhs);
    dvec_axpy( alpha, &e, &rhs); ///////---------
    dvec_alloc( Cinv.row, &temp);
    dcsr_aAxpy( 1.0, &Cinv, d.val, temp.val );
    dcsr_aAxpy( 1.0, B, temp.val, rhs.val );

    dvec_alloc( n1, &q);
    directsolve_UMF( &S, &rhs, &q, 2); // SOLVE
    //dcsr_pvgmres(&S, &rhs, &q, NULL, 1e-8, 20, 20, 1, 0);
    //smoother_dcsr_sgs(&q, &S, &rhs, 2);

    // v = 1/alpha C^{-1} (d - B^T q)
    dvec_alloc( n0, &v);
    dcsr_aAxpy( -1.0, &BT, q.val, d.val );
    dcsr_aAxpy( 1.0/alpha, &Cinv, d.val, v.val);

    // update u = u + [v, q]
    for( i=0; i<n0; i++ ){
      u->val[i] += v.val[i];
    }
    for( i=0; i<n1; i++ ){
      u->val[i+n0] += q.val[i];
      //u->val[i+n0] -= q.val[i]; ///------------?
    }


    // free
    dcsr_free(&BT);
    dcsr_free(&S);
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

      bmgl[0].mgl[blk][0].num_levels = max_levels;
      bmgl[0].mgl[blk][0].cycle_type = param->cycle_type;
      bmgl[0].mgl[blk][0].b = dvec_create(bmgl[0].mgl[blk][0].A.row);
      bmgl[0].mgl[blk][0].x = dvec_create(bmgl[0].mgl[blk][0].A.row);
      bmgl[0].mgl[blk][0].w = dvec_create(2*bmgl[0].mgl[blk][0].A.row);

      // Remove dirichlet boundaries
      FE_fake.dirichlet = bmgl[0].dirichlet_blk[blk];
      FE_fake.ndof = bmgl[0].A_diag[blk].row;
      dcsr_shift( &bmgl[0].mgl[blk][0].A,  1);
      printf("Eliminating dirichlet BC for A_diag\n");
      //eliminate_DirichletBC(NULL, &FE_fake , &bmgl[0].fine_level_mesh, NULL, &(bmgl[0].mgl[blk][0].A),0.0);
      eliminate_DirichletBC(NULL, &FE_fake , bmgl[0].fine_level_mesh, NULL, &(bmgl[0].mgl[blk][0].A),0.0);
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

        printf("RAP on A_diag\n");
        dcsr_alloc( bmgl[lvl+1].A_diag[blk].row, bmgl[lvl+1].A_diag[blk].row, bmgl[lvl+1].A_diag[blk].nnz, &bmgl[0].mgl[blk][lvl+1].A);
        dcsr_cp( &(bmgl[lvl+1].A_diag[blk]), &bmgl[0].mgl[blk][lvl+1].A );
        FE_fake.dirichlet = bmgl[lvl+1].dirichlet_blk[blk];
        FE_fake.ndof      = bmgl[lvl+1].A_diag[blk].row;
        dcsr_shift( &bmgl[0].mgl[blk][lvl+1].A,  1);
        //eliminate_DirichletBC(NULL, &FE_fake , &bmgl[lvl+1].fine_level_mesh, NULL, &(bmgl[0].mgl[blk][lvl+1].A),0.0);
        eliminate_DirichletBC(NULL, &FE_fake , bmgl[lvl+1].fine_level_mesh, NULL, &(bmgl[0].mgl[blk][lvl+1].A),0.0);
        dcsr_shift( &bmgl[0].mgl[blk][lvl+1].A, -1);
        printf("ELIM on A_diag fine\n");



        // setup total level number and current level
        bmgl[0].mgl[blk][lvl+1].num_levels = max_levels;
        bmgl[0].mgl[blk][lvl+1].cycle_type = param->cycle_type;
        bmgl[0].mgl[blk][lvl+1].b = dvec_create(bmgl[0].mgl[blk][lvl+1].A.row);
        bmgl[0].mgl[blk][lvl+1].x = dvec_create(bmgl[0].mgl[blk][lvl+1].A.row);
        bmgl[0].mgl[blk][lvl+1].w = dvec_create(2*bmgl[0].mgl[blk][lvl+1].A.row);

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
    //printf("\tBeginning Block smoother, lvl=%d\n",lvl);
    param->presmooth_iter = 5;
    //printf("\t\t| using %2d iterations\n",param->presmooth_iter);
    // Initialize
    INT n0, n1, n2;
    INT triType = 0;
    // Sub-vectors
    dvector x0, x1, x2;
    dvector b0, b1, b2;

    n0 = bmgl[lvl].mgl[0][0].A.row;
    n1 = bmgl[lvl].mgl[1][0].A.row;
    n2 = bmgl[lvl].mgl[2][0].A.row;

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

/*
    // Block 0: P1 + Bubble
//    smoother_dcsr_sgs(&x0, &(bmgl[lvl].mgl[0][0].A), &b0, param->presmooth_iter);
    directsolve_UMF(&(bmgl[lvl].mgl[0][0].A), &b0, &x0, 10);
//    bmgl[lvl].mgl[0]->b.row=n0; array_cp(n0, b0.val, bmgl[lvl].mgl[0]->b.val); // residual is an input
//    bmgl[lvl].mgl[0]->x.row=n0; dvec_set(n0, &bmgl[lvl].mgl[0]->x,0.0);
//    param->AMG_type = -1;
//    param->smoother = SMOOTHER_SGS;
//    INT i;
//    for(i=0;i<2;++i) mgcycle(bmgl[lvl].mgl[0], param);
//    array_cp(n0, bmgl[lvl].mgl[0]->x.val, x0.val);
    //printf("\tBlock 0 done...\n");

    // r1 = r1 - A3*z0
    if(triType==1){
      if (bmgl[lvl].A.blocks[3] != NULL)
        dcsr_aAxpy(-1.0, bmgl[lvl].A.blocks[3], x0.val, b1.val);
    }

    // Block 1: RT0
//    Schwarz_param swzparam;
//    swzparam.Schwarz_blksolver = bmgl[lvl].mgl[1][0].Schwarz.blk_solver;
//    if(pre_post == 1){
//      smoother_dcsr_Schwarz_forward( &(bmgl[lvl].mgl[1][0].Schwarz), &swzparam, &x1, &b1);
////      smoother_dcsr_Schwarz_backward(&(bmgl[lvl].mgl[1][0].Schwarz), &swzparam, &x1, &b1);
//    } else if (pre_post == 2){
//      smoother_dcsr_Schwarz_backward(&(bmgl[lvl].mgl[1][0].Schwarz), &swzparam, &x1, &b1);
////      smoother_dcsr_Schwarz_forward( &(bmgl[lvl].mgl[1][0].Schwarz), &swzparam, &x1, &b1);
//    }
    directsolve_UMF(&(bmgl[lvl].mgl[1][0].A), &b1, &x1, 10);
    //printf("\tBlock 1 done...\n");

    if(triType==1){
      // r2 = r2 - A6*z0 - A7*z1
      if (bmgl[lvl].A.blocks[6] != NULL)
          dcsr_aAxpy(-1.0, bmgl[lvl].A.blocks[6], x0.val, b2.val);
      if (bmgl[lvl].A.blocks[7] != NULL)
          dcsr_aAxpy(-1.0, bmgl[lvl].A.blocks[7], x1.val, b2.val);
    }

    // Block 2: P0
//    smoother_dcsr_sgs(&x2, &(bmgl[lvl].mgl[2][0].A), &b2, param->presmooth_iter);
    directsolve_UMF(&(bmgl[lvl].mgl[2][0].A), &b2, &x2, 10);
    //printf("\tBlock 2 done...\n");
*/

//block_dCSRmat Atemp;
//bdcsr_alloc( 2, 2, &Atemp);
//Atemp.blocks[0] = &bmgl[lvl].mgl[0][0].A;
//Atemp.blocks[1] = NULL;
//Atemp.blocks[2] = NULL;
//Atemp.blocks[3] = &bmgl[lvl].mgl[1][0].A;
//dCSRmat A = bdcsr_2_dcsr(&Atemp);

// BSR
dCSRmat A = bdcsr_subblk_2_dcsr ( &bmgl[lvl].A, 0, 1, 0, 1);
dCSRmat C = dcsr_create_identity_matrix( n0+n1, 0);
dvector diagvec;
dcsr_getdiag( A.row, &A, &diagvec);
C.val = diagvec.val;

dCSRmat B = bdcsr_subblk_2_dcsr ( &bmgl[lvl].A, 2, 2, 0, 1);

smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, 0.5, &bmgl[lvl].A, &C, &B, bmgl[lvl].A.blocks[8], 1);
//smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, 3.2484, &bmgl[lvl].A, &C, &B, bmgl[lvl].A.blocks[8], 1);
//smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, 1.8119, &bmgl[lvl].A, &C, &B, bmgl[lvl].A.blocks[8], 1);
//smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, 2, &bmgl[lvl].A, &C, &B, bmgl[lvl].A.blocks[8], 1);
//smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, 2.0, &bmgl[lvl].A, &C, &B, &bmgl[lvl].mgl[2][0].A, 1);
//smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, 1.0, &bmgl[lvl].A, &C, &B, NULL, 1);



    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
