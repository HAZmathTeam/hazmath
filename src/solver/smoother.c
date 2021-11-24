/*! \file src/solver/smoother.c
 *
 *  Smoothers
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 12/25/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note  Done cleanup for releasing -- Xiaozhe Hu 03/12/2017 & 08/28/2021
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
            mask[ki] = i+1;// TODO: zero-one fix?
        }

        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki = jblock[iblk];
            rhs.val[i] = b->val[ki];
            iaa = ia[ki];//-1; // TODO: zero-one fix?
            iab = ia[ki+1];//-1; // TODO: zero-one fix?
            for (kij = iaa; kij<iab; ++kij) {
                kj = ja[kij];//-1; // TODO: zero-one fix?
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
 * \fn void smoother_dcsr_Schwarz_forward_additive (Schwarz_data  *Schwarz,
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
void smoother_dcsr_Schwarz_forward_additive (Schwarz_data  *Schwarz,
                                    Schwarz_param *param,
                                    dvector       *x,
                                    dvector       *b,
                                    REAL       w)
{
    INT i, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;

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
    // Local solution storage
    dvector averaging_factor = dvec_create( x->row );
    dvector xout = dvec_create( x->row );//TODO: need to allocate

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
            mask[ki] = i+1;// TODO: zero-one fix?
        }// for i<nloc

        for (i=0; i<nloc; ++i) {
            iblk = ibl0 + i;
            ki = jblock[iblk];
            rhs.val[i] = b->val[ki];
            iaa = ia[ki];//-1; // TODO: zero-one fix?
            iab = ia[ki+1];//-1; // TODO: zero-one fix?
            for (kij = iaa; kij<iab; ++kij) {
                kj = ja[kij];//-1; // TODO: zero-one fix?
                //j  = mask[kj];
                //if(j == 0) {
                    rhs.val[i] -= val[kij]*x->val[kj];
                //}
            }
        }// for i<nloc

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
            ///////////////////////////////////////////// Option 1
            // averaging factor computed here, need loop later to sum over and apply.
            xout.val[ki] += u.val[i];//TODO: changed here
            averaging_factor.val[ki] += 1;
            ///////////////////////////////////////////// Option 2
            // averaging factor is precomputed. Averaging is done in the sum
        }
    }
    // Copy xout into x
    // Lazy way for now (memcpy or something is probably better)
    for (i=0; i<x->row; i++){
      if(averaging_factor.val[i] > 0){
        x->val[i] += w*xout.val[i]/averaging_factor.val[i];
      }
    }
    dvec_free(&xout);
    dvec_free(&averaging_factor);
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
            iaa = ia[ki];//-1;
            iab = ia[ki+1];//-1;
            for (kij = iaa; kij<iab; ++kij) {
                kj = ja[kij];//-1;
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
 * \fn void smoother_dcsr_Schwarz_backward_additive (Schwarz_data  *Schwarz,
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
void smoother_dcsr_Schwarz_backward_additive (Schwarz_data *Schwarz,
                                     Schwarz_param *param,
                                     dvector *x,
                                     dvector *b,
                                     REAL w)
{
    INT i, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;

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
    // Local solution storage
    dvector averaging_factor = dvec_create( x->row );
    dvector xout = dvec_create( x->row );//TODO: need to allocate

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
            iaa = ia[ki];//-1;
            iab = ia[ki+1];//-1;
            for (kij = iaa; kij<iab; ++kij) {
                kj = ja[kij];//-1;
                //j  = mask[kj];
                //if(j == 0) {
                    rhs.val[i] -= val[kij]*x->val[kj];
                //}
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
            // Option 1
            xout.val[ki] += u.val[i];
            averaging_factor.val[ki] += 1.0;
        }
    }
    // Copy xout into x
    printf("Addings\n");
    // Lazy way for now (memcpy or something is probably better)
    for (i=0; i<x->row; i++){
      if(averaging_factor.val[i] > 0){
        x->val[i] += w*xout.val[i]/averaging_factor.val[i];
      }
    }
    dvec_free(&xout);
    dvec_free(&averaging_factor);
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

/************************************************************************************************/
/**
 * \fn void smoother_dcsr_gs_graph_eigen(dvector *u, dCSRmat *A, const INT i_1, const INT i_n, const INT s, INT nsmooth, INT num_eigen)
 *
 * \brief Symmetric Gauss-Seidel smoother for graph Laplacian (used for eigenvalue solver)
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: has to be a zero vector
 * \param i_1  Starting index
 * \param i_n  Ending index
 * \param s    Increasing step
 * \param nsmooth      Number of iterations
 * \param num_eigen    Number of eigenvectors
 *
 *
 */
void smoother_dcsr_gs_graph_eigen(dvector *u,
                                  dCSRmat *A,
                                  dvector *b,
                                  const INT i_1,
                                  const INT i_n,
                                  const INT s,
                                  const INT nsmooth,
                                  const INT num_eigen)
{

  // local variables
  INT n = A->row;

  INT i, j;
  dvector x;
  x.row = n;

  // space for QR
  REAL *Q = (REAL *)calloc(n*num_eigen, sizeof(REAL));
  REAL *R = (REAL *)calloc(num_eigen*num_eigen, sizeof(REAL));

  // zero right hand side for eigenvalue problem
  //dvector b = dvec_create(n);
  //dvec_set(n, &b, 0.0);

  // main loop
  for (i=0; i<nsmooth; i++)
  {

    // loop over approximate eigenvectors
    for (j=0; j<num_eigen; j++)
    {

      // assign x
      x.val = &(u->val[j*n]);

      // GS smoothing
      smoother_dcsr_gs(&x, i_1, i_n, s, A, b, 1);

      // orthogonalize to the constant vector
      dvec_orthog_const(&x);

    }

    // Use QR decomposition to find x (orthonormal basis)
    // orthogonalize eigenvectors
    ddense_qr(n, num_eigen, u->val, Q, R);

    // copy orthonormal basis Q to x
    array_cp(n*num_eigen, Q, u->val);

  }

  // free
  free(Q);
  free(R);

}

/************************************************************************************************/
/**
 * \fn void smoother_dcsr_sgs_graph_eigen(dvector *u, dCSRmat *A, INT nsmooth, INT num_eigen)
 *
 * \brief Symmetric Gauss-Seidel smoother for graph Laplacian (used for eigenvalue solver)
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: has to be a zero vector
 * \param nsmooth      Number of iterations
 * \param num_eigen    Number of eigenvectors
 *
 *
 */
void smoother_dcsr_sgs_graph_eigen(dvector *u, dCSRmat *A, dvector *b, const INT nsmooth, const INT num_eigen)
{

  // local variables
  INT n = A->row;

  INT i, j;
  dvector x;
  x.row = n;

  // space for QR
  REAL *Q = (REAL *)calloc(n*num_eigen, sizeof(REAL));
  REAL *R = (REAL *)calloc(num_eigen*num_eigen, sizeof(REAL));

  // zero right hand side for eigenvalue problem
  //dvector b = dvec_create(n);
  //dvec_set(n, &b, 0.0);

  // main loop
  for (i=0; i<nsmooth; i++)
  {

    // loop over approximate eigenvectors
    for (j=0; j<num_eigen; j++)
    {

      // assign x
      x.val = &(u->val[j*n]);

      // SGS smoothing
      smoother_dcsr_sgs(&x, A, b, 1);

      // orthogonalize to the constant vector
      dvec_orthog_const(&x);

    }

    // Use QR decomposition to find x (orthonormal basis)
    // orthogonalize eigenvectors
    ddense_qr(n, num_eigen, u->val, Q, R);

    // copy orthonormal basis Q to x
    array_cp(n*num_eigen, Q, u->val);

  }

  // free
  free(Q);
  free(R);

  //dvec_free(&b);

}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
