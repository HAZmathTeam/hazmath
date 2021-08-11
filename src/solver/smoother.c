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
            //printf("%f\t\tu.row = %d, i=%d, ki=%d, ___rhs %f\n",u.val[i],u.row,i,ki,rhs.val[i]);
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
    //INT i, j, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;
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
            //xout.val[ki] += u.val[i] * averaging_factor.val[ki];
            //printf("%f\t\tu.row = %d, i=%d, ki=%d, ___rhs %f\n",u.val[i],u.row,i,ki,rhs.val[i]);
        }
    }
    // Copy xout into x
    printf("Addings\n");
    // Lazy way for now (memcpy or something is probably better)
    for (i=0; i<x->row; i++){
      // Using Option 1
      if(averaging_factor.val[i] > 0){
        x->val[i] += w*xout.val[i]/averaging_factor.val[i];
      }
      // Using Option 2
      //x->val[i] = xout.val[i];
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
    //INT i, j, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;
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
            //if(ki==0){ printf("ki=0 | i=%d | patch=%d\n",i,is);}
        }
    }
    // Copy xout into x
    printf("Addings\n");
    // Lazy way for now (memcpy or something is probably better)
    for (i=0; i<x->row; i++){
      // Using Option 1
      //printf("%f %f\n",x->val[i],xout.val[i]/averaging_factor.val[i]);
      if(averaging_factor.val[i] > 0){
        x->val[i] += w*xout.val[i]/averaging_factor.val[i];
      }
      //printf("\t AF[%d] = %f\n",i,averaging_factor.val[i]);
      // Using Option 2
      //x->val[i] = xout.val[i];
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
 * [   B      M   ] [q] = [e]
 *
 */
void smoother_bdcsr_bsr(dvector *u,
                        dvector *b,
                        REAL alpha,
                        REAL w,
                        block_dCSRmat *A,
                        dCSRmat *C,
                        dCSRmat *B,
                        dCSRmat *BT,
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

    dCSRmat Cinv;
    dcsr_alloc( C->row, C->col, C->nnz, &Cinv);
    dcsr_cp( C, &Cinv);
    for(i=0; i<Cinv.row; i++){ // Assuming C is diag
      if( ABS( Cinv.val[i] ) > SMALLREAL ) Cinv.val[i] = 1.0 / Cinv.val[i];
    }

    dCSRmat S; // S = -BC^{-1}B^T
    dCSRmat Stemp;
    dcsr_rap( B, &Cinv, BT, &S);
    dcsr_axm(&S, -1.0);
    if( M ){
      dcsr_alloc( S.row, S.col, S.nnz, &Stemp);
      dcsr_cp( &S, &Stemp);
      dcsr_free( &S );
      dcsr_add( &Stemp, 1.0, M, alpha, &S); // ( S - alpha M ) Check sign here
    }

    // fill local var
    n0 = BT->row;
    n1 = B->row;

    // solver loop
    dvec_alloc( b->row, &r);
    dvec_cp( b, &r);
    bdcsr_aAxpy( -1.0, A, u->val, r.val); // r = b - Au
    REAL absres0 = dvec_norm2(&r);

    d.row = n0;
    d.val = r.val;
    e.row = n1;
    e.val = r.val + n0;

    // (B C^{-1} B^T) q = B C^{-1} d - alpha e
    dvec_alloc( e.row, &rhs);
    dvec_axpy( alpha, &e, &rhs);// rhs = alpha*e
    dvec_alloc( Cinv.row, &temp);
    dcsr_aAxpy( 1.0, &Cinv, d.val, temp.val );// temp = Cinv*d
    dcsr_aAxpy( -1.0, B, temp.val, rhs.val );// rhs = alpha*e - B*Cinv*d

    dvec_alloc( n1, &q);
    directsolve_UMF( &S, &rhs, &q, 0); // SOLVE
    //dcsr_pvgmres(&S, &rhs, &q, NULL, 1e-8, 20, 20, 1, 0);
    //smoother_dcsr_sgs(&q, &S, &rhs, 2);

    // v = 1/alpha C^{-1} (d - B^T q)
    dvec_alloc( n0, &v);
    dcsr_aAxpy( -1.0, BT, q.val, d.val );// d - BT*q
    dcsr_aAxpy( 1.0/alpha, &Cinv, d.val, v.val);

    // update u = u + [v, q]
    for( i=0; i<n0; i++ ){
      u->val[i] += w * v.val[i];
    }
    for( i=0; i<n1; i++ ){
      u->val[i+n0] += w * q.val[i];
    }


    dvec_cp( b, &r);
    bdcsr_aAxpy( -1.0, A, u->val, r.val); // r = b - Au
    REAL absres  = dvec_norm2(&r);

    REAL factor  = absres/absres0;
    printf("SMOOTHER: %f\n",factor);

    // free
    //dcsr_free(&BT);
    dcsr_free(&S);
}

/**
 * \fn void smoother_bdcsr_bsr_biot3 (dvector *u, const INT s, block_dCSRmat *A, dvector *b, INT L)
 *
 * \brief BSR
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param s      Increasing step
 * \param A      Pointer to block_dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 * [   C   0    B   ] [v] = [d]
 * [   0   M    B   ] [w] = [e]
 * [   BT  BT   M   ] [q] = [f]
 *
 */
void smoother_bdcsr_bsr_biot3(dvector *u,
                              dvector *b,
                              REAL alpha,
                              REAL w,
                              block_dCSRmat *A,
                              INT L,
                              AMG_data* mgl_disp,
                              MG_blk_data *bmgl)
{
printf("Beginning BSR for Biot\n");
    // Local variables
    INT i;
    INT n0, n1, n2;
    n0 = A->blocks[0]->row;
    n1 = A->blocks[4]->row;
    n2 = A->blocks[8]->row;
    dvector r;
    dvector rhs;
    dvector temp;

    dvector v;
    dvector y;
    dvector q;

    dvector d;
    dvector e;
    dvector f;

    Schwarz_param swzparam;
    swzparam.Schwarz_blksolver = mgl_disp->Schwarz.blk_solver;

    INT Au_solve_TYPE = 2;
    // Problem Var...
    REAL M = 1e6;
    REAL nu = 0.49;
    REAL mu =  (3e4) / (1+2*nu);
    REAL lam = (3e4)*nu / ((1-2*nu)*(1+nu));
    REAL factor = (1.0) / (lam+(2*mu/2.0));
    REAL lamS = lam;

    // Approx Ainv with diag(A)^{-1}
    dCSRmat Ainv = dcsr_create_identity_matrix( A->blocks[0]->row, 0);
    dvector diagvec2;
    dcsr_getdiag( A->blocks[0]->row, A->blocks[0], &diagvec2);
    Ainv.val = diagvec2.val;
    for( i=0; i<Ainv.row; i++){ if( ABS(Ainv.val[i]) > SMALLREAL ) Ainv.val[i] = 1.0/Ainv.val[i];}
    // Approx Mwinv with diag(Mw)^{-1}
    dCSRmat Mwinv = dcsr_create_identity_matrix( A->blocks[4]->row, 0);
    dvector diagvec;
    dcsr_getdiag( A->blocks[4]->row, A->blocks[4], &diagvec);
    Mwinv.val = diagvec.val;
    for( i=0; i<Mwinv.row; i++){ if( ABS(Mwinv.val[i]) > SMALLREAL ) Mwinv.val[i] = 1.0/Mwinv.val[i];}
    // Form Schur complement. Using diag(Mw) and that BuT*Au^{-1}*Bu is spectrally equivalent to Mp
    dCSRmat S;
    dCSRmat Ap;
    dcsr_rap( A->blocks[7], &Mwinv, A->blocks[5], &Ap);
    dcsr_add( A->blocks[8], M*(factor+1.0/M), &Ap, -1.0, &S);// Mp is negative, probably...


    // Start
    dvec_alloc( b->row, &r);
    dvec_cp( b, &r);
//    printf("ThisIsDumb\n");
//    dCSRmat* blarg = A->blocks[0];
//    A->blocks[0] = bmgl[0].As->blocks[0];
    bdcsr_aAxpy( -1.0, A, u->val, r.val); // r = b - Au
//    A->blocks[0] = blarg;

    d.row = n0;
    d.val = r.val;
    e.row = n1;
    e.val = r.val + n0;
    f.row = n2;
    f.val = r.val + n0 + n1;

    // Get RHS
    dvec_alloc( f.row, &rhs);
    dvec_cp( &f, &rhs);

    dvec_alloc( e.row, &y);
    dcsr_aAxpy( 1.0, &Mwinv, e.val, y.val ); //y = Mwinv*e
    dcsr_aAxpy( -1.0, A->blocks[7], y.val, rhs.val );// rhs = rhs - B*Mwinv*e

    dvector rhs_stokes, sol_stokes, rhs_elast, sol_elast;
    //linear_itsolver_param linear_itparam;

    dvec_alloc( d.row, &temp);
    dvec_cp( &d, &temp);
    switch (Au_solve_TYPE) {
        case 0: // diag
          dvec_set( d.row, &d, 0.0);////////////////////////////
          dcsr_aAxpy( 1.0, &Ainv, temp.val, d.val);/////////////
          break;
        case 1: // GS
          dvec_set( d.row, &d, 0.0);////////////////////////////
          smoother_dcsr_sgs(&d, A->blocks[0], &temp, 1);
          //smoother_dcsr_gs(&d,0,A->blocks[0]->row,1, A->blocks[0], &temp, 2);
          break;
        case 2: // Schwarz
          dvec_set( d.row, &d, 0.0);////////////////////////////
          if(L==1){
            //smoother_dcsr_Schwarz_backward( &(mgl_disp->Schwarz), &swzparam, &d, &temp);
            //smoother_dcsr_Schwarz_forward( &(mgl_disp->Schwarz), &swzparam, &d, &temp);
            smoother_dcsr_Schwarz_forward_additive( &(mgl_disp->Schwarz), &swzparam, &d, &temp, 1.0);
          } else {
            //smoother_dcsr_Schwarz_forward( &(mgl_disp->Schwarz), &swzparam, &d, &temp);
            //smoother_dcsr_Schwarz_backward( &(mgl_disp->Schwarz), &swzparam, &d, &temp);
            smoother_dcsr_Schwarz_backward_additive( &(mgl_disp->Schwarz), &swzparam, &d, &temp, 1.0);
          }
          //dcsr_pvfgmres( A->blocks[0], &temp, &d, NULL, 1e-3, 1000, 1000, 1, 1);
          break;
        case 3: // Strange Thing
          dvec_alloc(n0+n2, &rhs_stokes);
          dvec_alloc(n0+n2, &sol_stokes);
          dvec_alloc(n0, &rhs_elast);
          dvec_alloc(n0, &sol_elast);
          for( i=0; i<n0; i++){
            rhs_stokes.val[i] = temp.val[i];
            sol_stokes.val[i] = 0.0;//d.val[i];
            rhs_elast.val[i]  = temp.val[i];
            sol_elast.val[i]  = 1.0;//d.val[i];
          }
          for( i=n0; i<(n0+n2); i++){
            rhs_stokes.val[i] = 0.0;
            sol_stokes.val[i] = 1.0;
          }

          block_directsolve_UMF( bmgl->As, &rhs_stokes, &sol_stokes, 0);
          directsolve_UMF( bmgl->As->blocks[0], &rhs_elast, &sol_elast, 0);
          //directsolve_UMF( A->blocks[0], &temp, &d, 0); // SOLVE
          for( i=0; i<n0; i++ ){
            //d.val[i] = ( lamS/(lamS+1) * (factor)*sol_stokes.val[i] + 1.0/(lamS+1) * sol_elast.val[i] /(2*mu) )/(2*mu);
            d.val[i] = ( lamS/((2*mu)*(lamS+2*mu)) * sol_stokes.val[i] + 1.0/(lamS+2*mu) * sol_elast.val[i] );
          }
          break;
        default: //Direct
          directsolve_UMF( A->blocks[0], &temp, &d, 0); // SOLVE
          break;
    }
    dcsr_aAxpy( -1.0, A->blocks[6], d.val, rhs.val );//rhs = rhs - BT*q

    dvec_alloc( n2, &q);
    directsolve_UMF( &S, &rhs, &q, 0); // SOLVE
    //dvec_set( q.row, &q, 0.0);////////////////////////////
    //smoother_dcsr_sgs(&q, &S, &rhs, 1);//////////////////

    // propigate back
    dvec_alloc( d.row, &v);
    dvec_cp( &d, &v); // d = Au\r for displacement, so now v is as well (v=Au\d).

    dvec_alloc( d.row, &temp);
    dcsr_aAxpy( 1.0, A->blocks[2], q.val, temp.val );
    switch (Au_solve_TYPE) {
        case 0: // diag
          dvec_set( d.row, &d, 0.0);//////////////////////////////
          dcsr_aAxpy( 1.0, &Ainv, temp.val, d.val);///////////////
          break;
        case 1: // GS
          dvec_set( d.row, &d, 0.0);//////////////////////////////
          smoother_dcsr_sgs(&d, A->blocks[0], &temp, 1);
          //smoother_dcsr_gs(&d,0,A->blocks[0]->row,1, A->blocks[0], &temp, 2);
          break;
        case 2: // Schwarz
          dvec_set( d.row, &d, 0.0);//////////////////////////////
          if(L==1){
            //smoother_dcsr_Schwarz_forward( &(mgl_disp->Schwarz), &swzparam, &d, &temp);
            smoother_dcsr_Schwarz_forward_additive( &(mgl_disp->Schwarz), &swzparam, &d, &temp, 1.0);
          } else {
            //smoother_dcsr_Schwarz_backward( &(mgl_disp->Schwarz), &swzparam, &d, &temp);
            smoother_dcsr_Schwarz_backward_additive( &(mgl_disp->Schwarz), &swzparam, &d, &temp, 1.0);
          }
          break;
        case 3: // Strange Thing
          dvec_alloc(n0+n2, &rhs_stokes);
          dvec_alloc(n0+n2, &sol_stokes);
          dvec_alloc(n0, &rhs_elast);
          dvec_alloc(n0, &sol_elast);
          for( i=0; i<(n0); i++){
            rhs_stokes.val[i] = temp.val[i];
            sol_stokes.val[i] = d.val[i];
            rhs_elast.val[i] = temp.val[i];
            sol_elast.val[i] = d.val[i];
          }
          for( i=n0; i<(n0+n2); i++){
            rhs_stokes.val[i] = 0.0;
            sol_stokes.val[i] = 0.0;
          }
          block_directsolve_UMF( bmgl->As, &rhs_stokes, &sol_stokes, 0);
          directsolve_UMF( bmgl->As->blocks[0], &rhs_elast, &sol_elast, 0);
          for( i=0; i<n0; i++ ){
            //d.val[i] = lamS/(lamS+1) * sol_stokes.val[i] + 1.0/(lamS+1) * sol_elast.val[i];
            d.val[i] = ( (lamS/(2*mu))/(lamS+2*mu) * sol_stokes.val[i] + 1.0/(lamS+2*mu) * sol_elast.val[i] );
          }
          //directsolve_UMF( A->blocks[0], &temp, &d, 0); // SOLVE
          break;
        default: // Direct
          directsolve_UMF( A->blocks[0], &temp, &d, 0); // SOLVE
          break;
    }
    dvec_axpy( -1.0, &d, &v);// v = v - d

    // propigate back
    dvec_alloc( d.row, &temp);
    dcsr_aAxpy( 1.0, A->blocks[5], q.val, temp.val );
    dcsr_aAxpy( -1.0, &Mwinv, temp.val, y.val);

    // update u = u + [v, q]
    for( i=0; i<n0; i++ ){
      u->val[i] += w * v.val[i];
    }
    for( i=0; i<n1; i++ ){
      u->val[i+n0] += w * y.val[i];
    }
    for( i=0; i<n2; i++ ){
      u->val[i+n0+n1] += w * q.val[i];
    }

    return;
}

/**
 * \fn void smoother_bdcsr_uzawa (dvector *u, const INT s, block_dCSRmat *A, dvector *b, INT L)
 *
 * \brief UZAWA
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param s      Increasing step
 * \param A      Pointer to block_dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 * [ alpha C  B^T ] [v] = [d]
 * [   B      M   ] [q] = [e]
 *
 */
void smoother_bdcsr_uzawa(dvector *u,
                          dvector *b,
                          REAL alpha,
                          REAL w,
                          block_dCSRmat *A,
                          dCSRmat *C,
                          dCSRmat *B,
                          dCSRmat *M,
                          INT L)
{
printf("Beginning UZAWA\n");
    // Local variables
    INT i;
    INT n0, n1;
    dvector r;
    dvector rhs;
    dvector v;
    dvector q;
    dvector d;
    dvector e;

    dCSRmat BT;
    dcsr_trans(B, &BT);
    dcsr_axm(&BT, -1.0); // Not symmetric

    dCSRmat Cinv;
    dcsr_alloc( C->row, C->col, C->nnz, &Cinv);
    dcsr_cp( C, &Cinv);
    for(i=0; i<Cinv.row; i++){ // Assuming C is diag
      if( ABS( Cinv.val[i] ) > SMALLREAL ) Cinv.val[i] = 1.0 / (alpha*Cinv.val[i]);
    }

    dCSRmat S; // S = M - BC^{-1}B^T
    dCSRmat Stemp;
    dcsr_rap( B, &Cinv, &BT, &S);
    if( M ){
      dcsr_alloc( S.row, S.col, S.nnz, &Stemp);
      dcsr_cp( &S, &Stemp);
      dcsr_free( &S );
      dcsr_add( &Stemp, -1.0, M, 1.0, &S); // ( M - S )
    }

    // fill local var
    n0 = B->col;
    n1 = B->row;

    // solver loop
    dvec_alloc( b->row, &r);
    dvec_cp( b, &r);
    bdcsr_aAxpy( -1.0, A, u->val, r.val); // r = b - Au

    d.row = n0;
    d.val = r.val;
    e.row = n1;
    e.val = r.val + n0;

    // v = C^{-1} (d)
    dvec_alloc( n0, &v);
    dcsr_aAxpy( 1.0, &Cinv, d.val, v.val);
    // update u = u + [v, q]
    for( i=0; i<n0; i++ ){
      u->val[i] += w * v.val[i];
    }

    // (B C^{-1} B^T) q = e - B C^{-1} d
    dvec_alloc( e.row, &rhs);
    dvec_axpy( 1.0, &e, &rhs);
    dcsr_aAxpy( -1.0, B, v.val, rhs.val );

    dvec_alloc( n1, &q);
    directsolve_UMF( &S, &rhs, &q, 0); // SOLVE

    // update u = u + [v, q]
    for( i=0; i<n1; i++ ){
      u->val[i+n0] += w * q.val[i];
    }

    // free
    dcsr_free(&BT);
    dcsr_free(&S);
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
    //param->presmooth_iter = 5;
    //printf("\t\t| using %2d iterations\n",param->presmooth_iter);
    // Initialize
    INT n0, n1, n2;
    INT triType = 1;
    // Sub-vectors
//    dvector x0, x1, x2;
//    dvector b0, b1, b2;
    //dvector res = dvec_create(bmgl[lvl].b.row);
    dvector res, y;
    dvector r0, r1, r2;
    dvector y0, y1, y2;

    //n0 = bmgl[lvl].mgl[0][0].A.row;
    n0 = bmgl[lvl].A.blocks[0]->row;
    //n1 = bmgl[lvl].mgl[1][0].A.row;
    n1 = bmgl[lvl].A.blocks[4]->row;
    //n2 = bmgl[lvl].mgl[2][0].A.row;
    n2 = bmgl[lvl].A.blocks[8]->row;

//    x0.row = n0;
//    x0.val = bmgl[lvl].x.val;
//    b0.row = n0;
//    b0.val = bmgl[lvl].b.val;
//
//    x1.row = n1;
//    x1.val = &(bmgl[lvl].x.val[n0]);
//    b1.row = n1;
//    b1.val = &(bmgl[lvl].b.val[n0]);
//
//    x2.row = n2;
//    x2.val = &(bmgl[lvl].x.val[n0+n1]);
//    b2.row = n2;
//    b2.val = &(bmgl[lvl].b.val[n0+n1]);
//
//    dvector b_disp;

    if(0){
      printf("Block relaxation. Type = %d\n",triType);
      // Alloc for correction
      dvec_alloc(bmgl[lvl].x.row, &y);
      y0.row = n0; y0.val = y.val;
      y1.row = n1; y1.val = &(y.val[n0]);
      y2.row = n2; y2.val = &(y.val[n0+n1]);

      // Calculate Residual, then smooth on residual
      dvec_alloc(bmgl[lvl].b.row, &res);
      dvec_cp(&(bmgl[lvl].b),&res);
      bdcsr_aAxpy(-1.0, &bmgl[lvl].A, bmgl[lvl].x.val, res.val);
      r0.row = n0; r0.val = res.val;
      r1.row = n1; r1.val = &(res.val[n0]);
      r2.row = n2; r2.val = &(res.val[n0+n1]);

      // Block 0: P1 + Bubble
    //dvector rhs_stokes, sol_stokes, rhs_elast, sol_elast;
    //REAL nu = 0.49;
    //REAL mu =  (3e4) / (1+2*nu);
    //REAL lam = (3e4)*nu / ((1-2*nu)*(1+nu));
    //REAL factor = (1.0) / (lam+(2*mu/2.0));
    //REAL lamS = lam;
    //INT i;
    //      dvec_alloc(n0+n2, &rhs_stokes);
    //      dvec_alloc(n0+n2, &sol_stokes);
    //      dvec_alloc(n0, &rhs_elast);
    //      dvec_alloc(n0, &sol_elast);
    //      for( i=0; i<(n0); i++){
    //        rhs_stokes.val[i] = r0.val[i];
    //        sol_stokes.val[i] = 0.0;
    //        rhs_elast.val[i] = r0.val[i];
    //        sol_elast.val[i] = 0.0;
    //      }
    //      for( i=n0; i<(n0+n2); i++){
    //        rhs_stokes.val[i] = 0.0;
    //        sol_stokes.val[i] = 0.0;
    //      }
    //      block_directsolve_UMF( bmgl[lvl].As, &rhs_stokes, &sol_stokes, 0);
    //      directsolve_UMF( bmgl[lvl].As->blocks[0], &rhs_elast, &sol_elast, 0);
    //      for( i=0; i<n0; i++ ){
    //        y0.val[i] = ( (lamS/(2*mu))/(lamS+2*mu) * sol_stokes.val[i] + 1.0/(lamS+2*mu) * sol_elast.val[i] );
    //      }
//    smoother_dcsr_sgs(&x0, &(bmgl[lvl].mgl[0][0].A), &b0, param->presmooth_iter);
      directsolve_UMF(&(bmgl[lvl].mgl[0][0].A), &r0, &y0, 0);
//    bmgl[lvl].mgl[0]->b.row=n0; array_cp(n0, b0.val, bmgl[lvl].mgl[0]->b.val); // residual is an input
//    bmgl[lvl].mgl[0]->x.row=n0; dvec_set(n0, &bmgl[lvl].mgl[0]->x,0.0);
//    param->AMG_type = -1;
//    param->smoother = SMOOTHER_SGS;
//    INT i;
//    for(i=0;i<2;++i) mgcycle(bmgl[lvl].mgl[0], param);
//    array_cp(n0, bmgl[lvl].mgl[0]->x.val, x0.val);
    //printf("\tBlock 0 done...\n");

      if(triType==1){
        // r1 = r1 - A3*z0
        if (bmgl[lvl].A.blocks[3] != NULL)
          dcsr_aAxpy(-1.0, bmgl[lvl].A.blocks[3], y0.val, r1.val);
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
      directsolve_UMF(&(bmgl[lvl].mgl[1][0].A), &r1, &y1, 0);
      //printf("\tBlock 1 done...\n");

      if(triType==1){
        // r2 = r2 - A6*z0 - A7*z1
        if (bmgl[lvl].A.blocks[6] != NULL)
            dcsr_aAxpy(-1.0, bmgl[lvl].A.blocks[6], y0.val, r2.val);
        if (bmgl[lvl].A.blocks[7] != NULL)
            dcsr_aAxpy(-1.0, bmgl[lvl].A.blocks[7], y1.val, r2.val);
      }

      // Block 2: P0
//    smoother_dcsr_sgs(&x2, &(bmgl[lvl].mgl[2][0].A), &b2, param->presmooth_iter);
      directsolve_UMF(&(bmgl[lvl].mgl[2][0].A), &r2, &y2, 0);
      //printf("\tBlock 2 done...\n");

      // Update Solution
      //dvec_axpy(1.0, &y, &bmgl[lvl].x);
      dvec_axpy(0.99, &y, &bmgl[lvl].x);

    } else if (1) {
    smoother_bdcsr_bsr_biot3( &bmgl[lvl].x, &bmgl[lvl].b, param->BSR_alpha, param->BSR_omega, &bmgl[lvl].A, pre_post, &bmgl[lvl].mgl[0][0],&bmgl[lvl]);
    }else {

    // BSR
    dCSRmat A = bdcsr_subblk_2_dcsr ( &bmgl[lvl].A, 0, 1, 0, 1);
    dCSRmat C = dcsr_create_identity_matrix( n0+n1, 0);
    dvector diagvec;
    dcsr_getdiag( A.row, &A, &diagvec);
    C.val = diagvec.val;

    dCSRmat B = bdcsr_subblk_2_dcsr (  &bmgl[lvl].A, 2, 2, 0, 1);
    dCSRmat BT = bdcsr_subblk_2_dcsr ( &bmgl[lvl].A, 0, 1, 2, 2);

    //smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, 2.58, 1.79, &bmgl[lvl].A, &C, &B, &BT, bmgl[lvl].A.blocks[8], 1);
    //smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, 5.06, 1.99, &bmgl[lvl].A, &C, &B, &BT, bmgl[lvl].A.blocks[8], 1);
    smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, param->BSR_alpha, param->BSR_omega, &bmgl[lvl].A, &C, &B, &BT, bmgl[lvl].A.blocks[8], 1);
    //smoother_bdcsr_uzawa( &bmgl[lvl].x, &bmgl[lvl].b, 1.8118, 1.0550, &bmgl[lvl].A, &C, &B, bmgl[lvl].A.blocks[8], 1);
    }



    return;
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

  //dvec_free(&b);

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


/************************************************************************************************/
/**
 * \fn void smoother_block_elasticity( MG_blk_data *mgl, AMG_param *param)
 *
 * \brief Setup of block smoothers
 *
 * \param mgl       Pointer to MG_blk_data
 * \param param     Pointer to AMG_param
 *
 */
void smoother_block_elasticity( const INT lvl, MG_blk_data *bmgl, AMG_param *param, INT pre_post)
{
    INT n0;//, n1;

    n0 = bmgl[lvl].A.blocks[1]->row;

    // BSR
    dvector diagvec;
    dcsr_getdiag( bmgl[lvl].A.blocks[0]->row, bmgl[lvl].A.blocks[0], &diagvec);
    dCSRmat C = dcsr_create_identity_matrix( n0, 0);
    C.val = diagvec.val;

    smoother_bdcsr_bsr( &bmgl[lvl].x, &bmgl[lvl].b, param->BSR_alpha, param->BSR_omega,
                        &bmgl[lvl].A,
                        &C,
                        bmgl[lvl].A.blocks[2],
                        bmgl[lvl].A.blocks[1],
                        NULL, pre_post);

//    ////////////////////////////////////////////////
//    // Precond elasticity with stokes solve
//    // /////////////////////////////////////////////
//    INT i;
//    REAL lam = 0;
//
//    dvector rhs_stokes;
//    dvector rhs_elast;
//    dvector sol_stokes;
//    dvector sol_elast;
//    dvector res;
//
//    dvec_alloc(bmgl[lvl].b.row, &res);
//    dvec_cp(&(bmgl[lvl].b),&res);
//
//    bdcsr_aAxpy(-1.0, &bmgl[lvl].A, bmgl[lvl].x.val, res.val);
//    printf("Full Norm Pre: %e\n", dvec_norm2(&res) );
//
//    dvec_alloc(n0+n1, &rhs_stokes);
//    dvec_alloc(n0+n1, &sol_stokes);
//    dvec_alloc(n0, &rhs_elast);
//    dvec_alloc(n0, &sol_elast);
//    for( i=0; i<(n0); i++){
//      rhs_stokes.val[i] = res.val[i];
//      sol_stokes.val[i] = 0.0;
//      rhs_elast.val[i] = res.val[i];
//      sol_elast.val[i] = 0.0;
//    }
//    for( i=n0; i<(n0+n1); i++){
//      rhs_stokes.val[i] = 0.0;//res.val[i];
//      sol_stokes.val[i] = 0.0;
//    }
//
//    block_directsolve_UMF( bmgl[lvl].As, &rhs_stokes, &sol_stokes, 0);
////    bdcsr_aAxpy(-1.0, bmgl[lvl].As, sol_stokes.val, rhs_stokes.val);
////    printf("Stokes Norm: %e\n", dvec_norm2(&rhs_stokes) );
//
//    //directsolve_UMF( bmgl[lvl].As->blocks[0], &rhs_elast, &sol_elast, 0);
//    //directsolve_UMF( bmgl[lvl].A.blocks[0], &rhs_elast, &sol_elast, 0);
//    //smoother_dcsr_sgs(&sol_elast, bmgl[lvl].As->blocks[0], &rhs_elast, 100);
//    smoother_dcsr_sgs(&sol_elast, bmgl[lvl].A.blocks[0], &rhs_elast, 2);
////    dcsr_aAxpy(-1.0, bmgl[lvl].As->blocks[0], sol_elast.val, rhs_elast.val);
////    printf("Elast Norm: %e\n", dvec_norm2(&rhs_elast) );
//
//    for( i=0; i<n0; i++ ){
//      bmgl[lvl].x.val[i] += param->BSR_omega*( lam/(lam+1.0) * sol_stokes.val[i] + 1.0/(lam+1.0) * sol_elast.val[i]);
//    }
//
////    for( i=0; i<n1; i++){//Jacobi?
////      bmgl[lvl].x.val[n0+i] = bmgl[lvl].x.val[n0+i] + 0.3*res.val[n0+i]/bmgl[lvl].A.blocks[3]->val[i];
////    }
//
//    dvec_cp(&(bmgl[lvl].b),&res);
//    bdcsr_aAxpy(-1.0, &bmgl[lvl].A, bmgl[lvl].x.val, res.val);
//    printf("Full Norm Post: %e\n", dvec_norm2(&res) );
//    printf("--------------------\n");

    return;
}
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
