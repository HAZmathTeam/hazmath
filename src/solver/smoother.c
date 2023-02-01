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
    REAL        w = 1.0;  //0.8
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

/* #if DIAGONAL_PREF // diagonal first */
/*                 d=aj[begin_row]; */
/*                 if (ABS(d)>SMALLREAL) { */
/*                     for (k=begin_row+1;k<end_row;++k) { */
/*                         j=ja[k]; */
/*                         t-=aj[k]*uval[j]; */
/*                     } */
/*                     uval[i]=t/d; */
/*                 } */
/* #else // general order */
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    if (i!=j)
                        t-=aj[k]*uval[j];
                    else if (ABS(aj[k])>SMALLREAL) d=1.e+0/aj[k];
                }
                uval[i]=t*d;
/* #endif */
            } // end for i
        } // end while

    } // if s
    else {

        while (L--) {
            for (i=i_1;i>=i_n;i+=s) {
                t=bval[i];
                begin_row=ia[i],end_row=ia[i+1];
/* #if DIAGONAL_PREF // diagonal first */
/*                 d=aj[begin_row]; */
/*                 if (ABS(d)>SMALLREAL) { */
/*                     for (k=begin_row+1;k<end_row;++k) { */
/*                         j=ja[k]; */
/*                         t-=aj[k]*uval[j]; */
/*                     } */
/*                     uval[i]=t/d; */
/*                 } */
/* #else // general order */
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    if (i!=j)
                        t-=aj[k]*uval[j];
                    else if (ABS(aj[k])>SMALLREAL) d=1.0/aj[k];
                }
                uval[i]=t*d;
/* #endif */
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

    //#if WITH_SUITESPARSE
    void **numeric = Schwarz->numeric;
    //#endif

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

	  //#if WITH_SUITESPARSE
            case SOLVER_UMFPACK: {
                /* use UMFPACK direct solver on each block */
                hazmath_solve(&blk[is], &rhs, &u, numeric[is], 0);
                break;
            }
	      //#endif
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

    //#if WITH_SUITESPARSE
    void **numeric = Schwarz->numeric;
    //#endif

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

	  //#if WITH_SUITESPARSE
            case SOLVER_UMFPACK: {
                /* use UMFPACK direct solver on each block */
                hazmath_solve(&blk[is], &rhs, &u, numeric[is], 0);
                break;
            }
	      //#endif
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

    //#if WITH_SUITESPARSE
    void **numeric = Schwarz->numeric;
    //#endif

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

	  //#if WITH_SUITESPARSE
            case SOLVER_UMFPACK: {
                /* use UMFPACK direct solver on each block */
                hazmath_solve(&blk[is], &rhs, &u, numeric[is], 0);
                break;
            }
	      //#endif
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

    //#if WITH_SUITESPARSE
    void **numeric = Schwarz->numeric;
    //#endif

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

	  //#if WITH_SUITESPARSE
            case SOLVER_UMFPACK: {
                /* use UMFPACK direct solver on each block */
                hazmath_solve(&blk[is], &rhs, &u, numeric[is], 0);
                break;
            }
	      //#endif
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

/************************************************************************************************/
/**
 * \fn void smoother_dbsr_jacobi(dBSRmat *A, dvector *b, dvector *u,
 *                               REAL *diaginv)
 *
 * \brief Jacobi relaxation
 *
 * \param A        Pointer to dBSRmat: the coefficient matrix
 * \param b        Pointer to dvector: the right hand side
 * \param u        Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param diaginv  Inverses for all the diagonal blocks of A
 *
 */
void smoother_dbsr_jacobi(dBSRmat *A,
                          dvector *b,
                          dvector *u,
                          REAL    *diaginv)
{
    // members of A
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    size = ROW*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;
    REAL         *val = A->val;

    // SHORT nthreads = 1;

    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;

    // auxiliary array
    REAL *b_tmp = NULL;

    // local variables
    INT i,j,k;
    INT pb;

    // b_tmp = b_val
    b_tmp = (REAL *)calloc(size, sizeof(REAL));
    memcpy(b_tmp, b_val, size*sizeof(REAL));

    // No need to assign the smoothing order since the result doesn't depend on it
    if (nb == 1) {
        for (i = 0; i < ROW; ++i) {
            for (k = IA[i]; k < IA[i+1]; ++k) {
                j = JA[k];
                if (j != i)
                    b_tmp[i] -= val[k]*u_val[j];
            }
        }
        for (i = 0; i < ROW; ++i) {
            u_val[i] = b_tmp[i]*diaginv[i];
        }

        free(b_tmp); b_tmp = NULL;
    }
    else if (nb > 1) {
        for (i = 0; i < ROW; ++i) {
            pb = i*nb;
            for (k = IA[i]; k < IA[i+1]; ++k) {
                j = JA[k];
                if (j != i)
                    ddense_ymAx(val+k*nb2, u_val+j*nb, b_tmp+pb, nb);
            }
        }

        for (i = 0; i < ROW; ++i) {
            pb = i*nb;
            ddense_mxv(diaginv+nb2*i, b_tmp+pb, u_val+pb, nb);
        }

        free(b_tmp); b_tmp = NULL;
    }
    else {
        printf("### HAZMATH ERROR: nb is illegal! [%s:%d]\n", __FILE__, __LINE__);
        check_error(ERROR_NUM_BLOCKS, __FUNCTION__);
    }

}

/**
 * \fn void smoother_dbsr_gs_ascend(dBSRmat *A, dvector *b, dvector *u,
 *                                  REAL *diaginv)
 *
 * \brief Gauss-Seidel relaxation in the ascending order
 *
 * \param A  Pointer to dBSRmat: the coefficient matrix
 * \param b  Pointer to dvector: the right hand side
 * \param u  Pointer to dvector: the unknowns (IN: initial guess, OUT: approximation)
 * \param diaginv  Inverses for all the diagonal blocks of A
 *
 */
void smoother_dbsr_gs_ascend(dBSRmat *A,
                             dvector *b,
                             dvector *u,
                             REAL    *diaginv)
{
    // members of A
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;
    REAL         *val = A->val;

    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;

    // local variables
    INT   i,j,k;
    INT   pb;
    REAL  rhs = 0.0;

    if (nb == 1) {
        for (i = 0; i < ROW; ++i) {
            rhs = b_val[i];
            for (k = IA[i]; k < IA[i+1]; ++k) {
                j = JA[k];
                if (j != i)
                    rhs -= val[k]*u_val[j];
            }
            u_val[i] = rhs*diaginv[i];
        }
    }
    else if (nb > 1) {
        REAL *b_tmp = (REAL *)calloc(nb, sizeof(REAL));

        for (i = 0; i < ROW; ++i) {
            pb = i*nb;
            memcpy(b_tmp, b_val+pb, nb*sizeof(REAL));
            for (k = IA[i]; k < IA[i+1]; ++k) {
                j = JA[k];
                if (j != i)
                    ddense_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
            }
            ddense_mxv(diaginv+nb2*i, b_tmp, u_val+pb, nb);
        }

        free(b_tmp); b_tmp = NULL;
    }
    else {
        printf("### HAZMATH ERROR: nb is illegal! [%s:%d]\n", __FILE__, __LINE__);
        check_error(ERROR_NUM_BLOCKS, __FUNCTION__);
    }

}

/**
 * \fn void smoother_dbsr_gs_descend (dBSRmat *A, dvector *b, dvector *u,
 *                                         REAL *diaginv)
 *
 * \brief Gauss-Seidel relaxation in the descending order
 *
 * \param A  Pointer to dBSRmat: the coefficient matrix
 * \param b  Pointer to dvector: the right hand side
 * \param u  Pointer to dvector: the unknowns (IN: initial guess, OUT: approximation)
 * \param diaginv  Inverses for all the diagonal blocks of A
 *
 */
void smoother_dbsr_gs_descend(dBSRmat *A,
                              dvector *b,
                              dvector *u,
                              REAL    *diaginv )
{
    // members of A
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;
    REAL         *val = A->val;

    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;

    // local variables
    INT i,j,k;
    INT pb;
    REAL rhs = 0.0;

    if (nb == 1) {
        for (i = ROW-1; i >= 0; i--) {
            rhs = b_val[i];
            for (k = IA[i]; k < IA[i+1]; ++k) {
                j = JA[k];
                if (j != i)
                    rhs -= val[k]*u_val[j];
            }
            u_val[i] = rhs*diaginv[i];
        }
    }
    else if (nb > 1) {
        REAL *b_tmp = (REAL *)calloc(nb, sizeof(REAL));

        for (i = ROW-1; i >= 0; i--) {
            pb = i*nb;
            memcpy(b_tmp, b_val+pb, nb*sizeof(REAL));
            for (k = IA[i]; k < IA[i+1]; ++k) {
                j = JA[k];
                if (j != i)
                    ddense_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
            }
            ddense_mxv(diaginv+nb2*i, b_tmp, u_val+pb, nb);
        }

        free(b_tmp); b_tmp = NULL;
    }
    else {
        printf("### HAZAMTH ERROR: nb is illegal! [%s:%d]\n", __FILE__, __LINE__);
        check_error(ERROR_NUM_BLOCKS, __FUNCTION__);
    }

}

/************************************************************************************************/
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
// void smoother_bdcsr_jacobi(dvector *u,
//                           const INT s,
//                           block_dCSRmat *A,
//                           dvector *b,
//                           INT L)
// {
//     // Sub-vectors
//     dvector utemp;
//     dvector btemp;
//     // Smooth
//     INT i, istart;
//     INT row;
//     istart = 0;
//     for(i=0; i<A->brow; i++){
//         row = A->blocks[i+i*A->brow]->row;
//         // Get sub-vectors
//         utemp.row = row;
//         utemp.val = u->val+istart;
//         btemp.row = row;
//         btemp.val = b->val+istart;
//         // Call jacobi on specific block
//         smoother_dcsr_jacobi(&utemp,0,row-1,s,A->blocks[i+i*A->brow],&btemp,L);
//         // Move to next block
//         istart += row;
//     }
//
// }

/************************************************************************************************/
/**
 * \fn void smoother_bdcsr_jacobi_jacobi(dvector *u, block_dCSRmat *A, dvector *b, dCSRmat *A_diag)
 *
 * \brief block Jacobi smoother and each block uses jacobi (using specified diagonal blocks)
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_jacobi_jacobi(dvector *u,
                                  block_dCSRmat *A,
                                  dvector *b,
                                  dCSRmat *A_diag)
{
    // Sub-vectors
    dvector utemp;
    dvector btemp;

    // Smooth
    INT i, istart;
    INT row;
    istart = 0;
    for(i=0; i<A->brow; i++){

        if (A_diag == NULL){
            row = A->blocks[i*A->brow+i]->row;
        }
        else{
            row = A_diag[i].row;
        }

        // Get sub-vectors
        utemp.row = row;
        utemp.val = u->val+istart;
        btemp.row = row;
        btemp.val = b->val+istart;

        // Call jacobi on specific block
        if (A_diag == NULL){
            smoother_dcsr_jacobi(&utemp, 0, row-1, 1, A->blocks[i*A->brow+i], &btemp, 1);
        }
        else {
            smoother_dcsr_jacobi(&utemp, 0, row-1, 1, &A_diag[i], &btemp, 1);
        }
        //directsolve_HAZ(&A_diag[i], &btemp, &utemp, 3);

        // Move to next block
        istart += row;
    }

}

/************************************************************************************************/
/**
 * \fn void smoother_bdcsr_jacobi_fgs(dvector *u, block_dCSRmat *A, dvector *b, dCSRmat *A_diag)
 *
 * \brief block Jacobi smoother and each block uses forward Gauss-Seidel (using specified diagonal blocks)
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_jacobi_fgs(dvector *u,
                               block_dCSRmat *A,
                               dvector *b,
                               dCSRmat *A_diag)
{
    // Sub-vectors
    dvector utemp;
    dvector btemp;

    // Smooth
    INT i, istart;
    INT row;
    istart = 0;
    for(i=0; i<A->brow; i++){

        if (A_diag == NULL){
            row = A->blocks[i*A->brow+i]->row;
        }
        else{
            row = A_diag[i].row;
        }

        // Get sub-vectors
        utemp.row = row;
        utemp.val = u->val+istart;
        btemp.row = row;
        btemp.val = b->val+istart;

        // Call gs on specific block
        if (A_diag == NULL){
            smoother_dcsr_gs(&utemp, 0, row-1, 1, A->blocks[i*A->brow+i], &btemp, 1);
        }
        else{
            smoother_dcsr_gs(&utemp, 0, row-1, 1, &A_diag[i], &btemp, 1);
        }

        // Move to next block
        istart += row;
    }

}

/************************************************************************************************/
/**
 * \fn void smoother_bdcsr_jacobi_bgs(dvector *u, block_dCSRmat *A, dvector *b, dCSRmat *A_diag)
 *
 * \brief block Jacobi smoother and each block uses backward Gauss-Seidel (using specified diagonal blocks)
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_jacobi_bgs(dvector *u,
                               block_dCSRmat *A,
                               dvector *b,
                               dCSRmat *A_diag)
{
    // Sub-vectors
    dvector utemp;
    dvector btemp;

    // Smooth
    INT i, istart;
    INT row;
    istart = 0;
    for(i=0; i<A->brow; i++){

        if (A_diag == NULL){
            row = A->blocks[i*A->brow+i]->row;
        }
        else{
            row = A_diag[i].row;
        }

        // Get sub-vectors
        utemp.row = row;
        utemp.val = u->val+istart;
        btemp.row = row;
        btemp.val = b->val+istart;

        // Call gs on specific block
        if (A_diag == NULL){
            smoother_dcsr_gs(&utemp, row-1, 0, -1, A->blocks[i*A->brow+i], &btemp, 1);
        }
        else{
            smoother_dcsr_gs(&utemp, row-1, 0, -1, &A_diag[i], &btemp, 1);
        }

        // Move to next block
        istart += row;
    }

}

/************************************************************************************************/
/**
 * \fn void smoother_bdcsr_jacobi_sgs(dvector *u, block_dCSRmat *A, dvector *b, dCSRmat *A_diag)
 *
 * \brief block Jacobi smoother and each block uses symmetric Gauss-Seidel (using specified diagonal blocks)
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_jacobi_sgs(dvector *u,
                               block_dCSRmat *A,
                               dvector *b,
                               dCSRmat *A_diag)
{
    // Sub-vectors
    dvector utemp;
    dvector btemp;

    // Smooth
    INT i, istart;
    INT row;
    istart = 0;
    for(i=0; i<A->brow; i++){

        if (A_diag == NULL){
            row = A->blocks[i*A->brow+i]->row;
        }
        else{
            row = A_diag[i].row;
        }

        // Get sub-vectors
        utemp.row = row;
        utemp.val = u->val+istart;
        btemp.row = row;
        btemp.val = b->val+istart;

        // Call gs on specific block
        if (A_diag == NULL){
            smoother_dcsr_sgs(&utemp, A->blocks[i*A->brow+i], &btemp, 1);
        }
        else {
            smoother_dcsr_sgs(&utemp, &A_diag[i], &btemp, 1);
        }

        // Move to next block
        istart += row;
    }

}

/************************************************************************************************/
/**
 * \fn void smoother_bdcsr_fgs_fgs(dvector *u, block_dCSRmat *A, dvector *b, dCSRmat *A_diag, REAL *work)
 *
 * \brief block Forward Gauss-Seidel smoother and each block uses forward Gauss-Seidel (using specified diagonal blocks)
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_fgs_fgs(dvector *u,
                            block_dCSRmat *A,
                            dvector *b,
                            dCSRmat *A_diag,
                            REAL *work)
{
    // Sub-vectors
    dvector utemp;
    dvector btemp;

    // save right hand side b in the workspace
    array_cp(b->row, b->val, work);

    // Smooth
    INT i, j, istart, jstart;
    INT row;
    istart = 0;
    for(i=0; i<A->brow; i++){

        if (A_diag == NULL){
            row = A->blocks[i*A->brow+i]->row;
        }
        else{
            row = A_diag[i].row;
        }

        // Get sub-vectors
        utemp.row = row;
        utemp.val = u->val+istart;
        btemp.row = row;
        btemp.val = work+istart;

        // Call gs on specific block
        if (A_diag == NULL){
            smoother_dcsr_gs(&utemp, 0, row-1, 1, A->blocks[i*A->brow+i], &btemp, 1);
        }
        else {
            smoother_dcsr_gs(&utemp, 0, row-1, 1, &A_diag[i], &btemp, 1);
        }
        //directsolve_HAZ(&A_diag[i], &btemp, &utemp, 3);

        // Move to next block
        istart += row;

        // update right hand side (stored in workspace)
        jstart = istart;
        for (j=i+1; j<A->brow; j++)
        {
            // set starting place
            dcsr_aAxpy(-1.0, A->blocks[j*A->bcol+i], utemp.val, work+jstart);

            // move jstart
            if (A_diag == NULL){
                jstart += A->blocks[j*A->brow+j]->row;
            }
            else{
                jstart += A_diag[j].row;
            }
        }
    }

}

/**
 * \fn void smoother_bdcsr_fgs_sgs(dvector *u, block_dCSRmat *A, dvector *b, dCSRmat *A_diag, REAL *work)
 *
 * \brief block Forward Gauss-Seidel smoother and each block uses symmetric Gauss-Seidel (using specified diagonal blocks)
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_fgs_sgs(dvector *u,
                            block_dCSRmat *A,
                            dvector *b,
                            dCSRmat *A_diag,
                            REAL *work)
{
    // Sub-vectors
    dvector utemp;
    dvector btemp;

    // save right hand side b in the workspace
    array_cp(b->row, b->val, work);

    // Smooth
    INT i, j, istart, jstart;
    INT row;
    istart = 0;
    for(i=0; i<A->brow; i++){

        if (A_diag == NULL){
            row = A->blocks[i*A->brow+i]->row;
        }
        else{
            row = A_diag[i].row;
        }

        // Get sub-vectors
        utemp.row = row;
        utemp.val = u->val+istart;
        btemp.row = row;
        btemp.val = work+istart;

        // Call gs on specific block
        if (A_diag == NULL){
            //smoother_dcsr_gs(&utemp, 0, row-1, 1, A->blocks[i*A->brow+i], &btemp, 1);
            smoother_dcsr_sgs(&utemp, A->blocks[i*A->brow+i], &btemp, 1);
        }
        else {
            //smoother_dcsr_gs(&utemp, 0, row-1, 1, &A_diag[i], &btemp, 1);
            smoother_dcsr_sgs(&utemp, &A_diag[i], &btemp, 1);
        }
        //directsolve_HAZ(&A_diag[i], &btemp, &utemp, 3);

        // Move to next block
        istart += row;

        // update right hand side (stored in workspace)
        jstart = istart;
        for (j=i+1; j<A->brow; j++)
        {
            // set starting place
            dcsr_aAxpy(-1.0, A->blocks[j*A->bcol+i], utemp.val, work+jstart);

            // move jstart
            if (A_diag == NULL){
                jstart += A->blocks[j*A->brow+j]->row;
            }
            else{
                jstart += A_diag[j].row;
            }
        }
    }

}

/************************************************************************************************/
/**
 * \fn void smoother_bdcsr_bgs_bgs(dvector *u, block_dCSRmat *A, dvector *b, dCSRmat *A_diag, REAL *work)
 *
 * \brief block backward Gauss-Seidel smoother and each block uses backward Gauss-Seidel (using specified diagonal blocks)
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_bgs_bgs(dvector *u,
                            block_dCSRmat *A,
                            dvector *b,
                            dCSRmat *A_diag,
                            REAL *work)
{
    // Sub-vectors
    dvector utemp;
    dvector btemp;

    // save right hand side b in the workspace
    array_cp(b->row, b->val, work);

    // Smooth
    INT i, j, istart, jstart;
    INT row;
    istart = u->row;
    for(i=A->brow-1; i>-1; i--){

        if (A_diag == NULL){
            row = A->blocks[i*A->brow+i]->row;
        }
        else{
            row = A_diag[i].row;
        }
        istart -= row;

        // Get sub-vectors
        utemp.row = row;
        utemp.val = u->val+istart;
        btemp.row = row;
        btemp.val = work+istart;

        // Call gs on specific block
        if (A_diag == NULL){
            smoother_dcsr_gs(&utemp, row-1, 0, -1, A->blocks[i*A->brow+i], &btemp, 1);
        }
        else {
            smoother_dcsr_gs(&utemp, row-1, 0, -1, &A_diag[i], &btemp, 1);
        }
        //directsolve_HAZ(&A_diag[i], &btemp, &utemp, 3);

        // update right hand side (stored in workspace)
        jstart = istart;
        for (j=i-1; j>-1; j--)
        {
            // move jstart
            if (A_diag == NULL){
                jstart -= A->blocks[j*A->brow+j]->row;
            }
            else{
                jstart -= A_diag[j].row;
            }

            // set starting place
            dcsr_aAxpy(-1.0, A->blocks[j*A->bcol+i], utemp.val, work+jstart);
        }
    }

}

/************************************************************************************************/
/**
 * \fn void smoother_bdcsr_bgs_sgs(dvector *u, block_dCSRmat *A, dvector *b, dCSRmat *A_diag, REAL *work)
 *
 * \brief block backward Gauss-Seidel smoother and each block uses symmetric Gauss-Seidel (using specified diagonal blocks)
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_bgs_sgs(dvector *u,
                            block_dCSRmat *A,
                            dvector *b,
                            dCSRmat *A_diag,
                            REAL *work)
{
    // Sub-vectors
    dvector utemp;
    dvector btemp;

    // save right hand side b in the workspace
    array_cp(b->row, b->val, work);

    // Smooth
    INT i, j, istart, jstart;
    INT row;
    istart = u->row;
    for(i=A->brow-1; i>-1; i--){

        if (A_diag == NULL){
            row = A->blocks[i*A->brow+i]->row;
        }
        else{
            row = A_diag[i].row;
        }
        istart -= row;

        // Get sub-vectors
        utemp.row = row;
        utemp.val = u->val+istart;
        btemp.row = row;
        btemp.val = work+istart;

        // Call gs on specific block
        if (A_diag == NULL){
            //smoother_dcsr_gs(&utemp, row-1, 0, -1, A->blocks[i*A->brow+i], &btemp, 1);
            smoother_dcsr_sgs(&utemp, A->blocks[i*A->brow+i], &btemp, 1);
        }
        else {
            //smoother_dcsr_gs(&utemp, row-1, 0, -1, &A_diag[i], &btemp, 1);
            smoother_dcsr_sgs(&utemp, &A_diag[i], &btemp, 1);
        }
        //directsolve_HAZ(&A_diag[i], &btemp, &utemp, 3);

        // update right hand side (stored in workspace)
        jstart = istart;
        for (j=i-1; j>-1; j--)
        {
            // move jstart
            if (A_diag == NULL){
                jstart -= A->blocks[j*A->brow+j]->row;
            }
            else{
                jstart -= A_diag[j].row;
            }

            // set starting place
            dcsr_aAxpy(-1.0, A->blocks[j*A->bcol+i], utemp.val, work+jstart);
        }
    }

}



/************************************************************************************************/
// Smoothers for interface problems
/************************************************************************************************/
/**
 * \fn void smoother_bdcsr_metric_additive(dvector *u, block_dCSRmat *A, dvector *b,
 *                                         dCSRmat *A_diag, REAL *work, dCSRmat *interface_dof,
 *                                         block_dCSRmat *A_gamma)
 *
 * \brief block additive smoother for metric AMG
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_metric_additive(dvector *u,
                                    block_dCSRmat *A,
                                    dvector *b,
                                    dCSRmat *A_diag,
                                    REAL *work,
                                    dCSRmat *interface_dof,
                                    block_dCSRmat * A_gamma)
{
    // local variables
    dvector r;
    dvector r_gamma;
    dvector e_gamma;
    dvector u_gamma;

    INT i;

    INT brow = A->brow;

    //--------------------------------
    // smoother on the interface part
    //--------------------------------
    // save right hand side b in the workspace
    array_cp(b->row, b->val, work + b->row);

    //printf("done copy b\n");

    // compute the overall residual
    r.row = b->row;
    r.val = work + b->row;
    bdcsr_aAxpy(-1.0, A, u->val, r.val);

    //printf("done computing r\n");

    // get the residual for the interface part
    r_gamma.row = brow*interface_dof->row;
    r_gamma.val = r.val + r.row;
    for (i=0; i<interface_dof->row; i++) r_gamma.val[i] = r.val[interface_dof->JA[i]];
    for (i=0; i<interface_dof->row; i++) r_gamma.val[interface_dof->row+i] = r.val[A->blocks[0]->row+i];

    //printf("done get r\n");

    // set initial to be zero
    e_gamma.row = brow*interface_dof->row;
    e_gamma.val = r_gamma.val + r_gamma.row;
    dvec_set(e_gamma.row, &e_gamma, 0.0);

    // save solution
    u_gamma.row = brow*interface_dof->row;
    u_gamma.val = e_gamma.val + e_gamma.row;
    for (i=0; i<interface_dof->row; i++) u_gamma.val[i] = u->val[interface_dof->JA[i]];
    for (i=0; i<interface_dof->row; i++) u_gamma.val[interface_dof->row+i] = u->val[A->blocks[0]->row+i];

    //printf("done set e\n");

    // solve for the interface part
    block_directsolve_HAZ(A_gamma, &r_gamma, &e_gamma, 0);
    //smoother_bdcsr_fgs_fgs(&e_gamma, A_gamma, &r_gamma, NULL, work);
    //smoother_bdcsr_bgs_bgs(&e_gamma, A_gamma, &r_gamma, NULL, work);
    //smoother_bdcsr_jacobi_jacobi(&e_gamma, A_gamma, &r_gamma, NULL);

    //printf("done solve interface\n");

    //--------------------------------------------
    // Gauss-Seidel smoother for the whole matrix
    //--------------------------------------------
    //smoother_bdcsr_fgs_fgs(u, A, b, A_diag, work);
    //smoother_bdcsr_bgs_bgs(u, A, b, A_diag, work);
    smoother_bdcsr_jacobi_jacobi(u, A, b, A_diag);

    //printf("done smooth whole matrix\n");

    //--------------------------------------------
    // update the solution
    //--------------------------------------------
    //for (i=0; i<interface_dof->row; i++) u->val[interface_dof->JA[i]] += e_gamma.val[i];
    //for (i=0; i<interface_dof->row; i++) u->val[A->blocks[0]->row+i] += e_gamma.val[interface_dof->row+i];
    for (i=0; i<interface_dof->row; i++) u->val[interface_dof->JA[i]] = u_gamma.val[i] + e_gamma.val[i];
    for (i=0; i<interface_dof->row; i++) u->val[A->blocks[0]->row+i] = u_gamma.val[interface_dof->row+i] + e_gamma.val[interface_dof->row+i];


    //printf("done update u\n");

}

/**
 * \fn void smoother_bdcsr_metric_additive_bsr(dvector *u, block_dCSRmat *A, dvector *b,
 *                                         dCSRmat *A_diag, REAL *work, dCSRmat *interface_dof,
 *                                         dBSRmat *A_gamma, REAL    *diaginv)
 *
 * \brief block additive smoother for metric AMG (interface part in dBSRmat)
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_metric_additive_bsr(dvector *u,
                                        block_dCSRmat *A,
                                        dvector *b,
                                        dCSRmat *A_diag,
                                        REAL *work,
                                        dCSRmat *interface_dof,
                                        dBSRmat * A_gamma,
                                        REAL    *diaginv)
{
    // local variables
    dvector r;
    dvector r_gamma;
    dvector e_gamma;
    dvector u_gamma;

    INT i;

    INT brow = A->brow;

    //--------------------------------
    // smoother on the interface part
    //--------------------------------
    // save right hand side b in the workspace
    array_cp(b->row, b->val, work + b->row);

    //printf("done copy b\n");

    // compute the overall residual
    r.row = b->row;
    r.val = work + b->row;
    bdcsr_aAxpy(-1.0, A, u->val, r.val);

    //printf("done computing r\n");

    // get the residual for the interface part
    r_gamma.row = brow*interface_dof->row;
    r_gamma.val = r.val + r.row;
    for (i=0; i<interface_dof->row; i++) r_gamma.val[i*brow] = r.val[interface_dof->JA[i]];
    for (i=0; i<interface_dof->row; i++) r_gamma.val[i*brow+1] = r.val[A->blocks[0]->row+i];

    //printf("done get r\n");

    // set initial to be zero
    e_gamma.row = brow*interface_dof->row;
    e_gamma.val = r_gamma.val + r_gamma.row;
    dvec_set(e_gamma.row, &e_gamma, 0.0);

    // save solution
    u_gamma.row = brow*interface_dof->row;
    u_gamma.val = e_gamma.val + e_gamma.row;
    for (i=0; i<interface_dof->row; i++) u_gamma.val[i*brow] = u->val[interface_dof->JA[i]];
    for (i=0; i<interface_dof->row; i++) u_gamma.val[i*brow+1] = u->val[A->blocks[0]->row+i];

    //printf("done set e\n");

    // solve for the interface part
    //printf("ROW = %d\n", A_gamma->ROW);
    //smoother_dbsr_gs_ascend(A_gamma, &r_gamma, &e_gamma, diaginv);
    //smoother_dbsr_gs_descend(A_gamma, &r_gamma, &e_gamma, diaginv);
    //smoother_dbsr_jacobi(A_gamma, &r_gamma, &e_gamma, diaginv);
    //printf("done solve interface\n");

    //--------------------------------------------
    // Gauss-Seidel smoother for the whole matrix
    //--------------------------------------------
    //smoother_bdcsr_fgs_fgs(u, A, b, A_diag, work);
    //smoother_bdcsr_bgs_bgs(u, A, b, A_diag, work);
    smoother_bdcsr_jacobi_jacobi(u, A, b, A_diag);

    //printf("done smooth whole matrix\n");

    //--------------------------------------------
    // update the solution
    //--------------------------------------------
    //for (i=0; i<interface_dof->row; i++) u->val[interface_dof->JA[i]] += e_gamma.val[i*brow];
    //for (i=0; i<interface_dof->row; i++) u->val[A->blocks[0]->row+i] += e_gamma.val[i*brow+1];
    for (i=0; i<interface_dof->row; i++) u->val[interface_dof->JA[i]] = u_gamma.val[i*brow] + e_gamma.val[i*brow];
    for (i=0; i<interface_dof->row; i++) u->val[A->blocks[0]->row+i] = u_gamma.val[i*brow+1] + e_gamma.val[i*brow+1];

    //printf("done update u\n");

}

/**
 * \fn void smoother_bdcsr_metric_multiplicative_omega_gamma(dvector *u, block_dCSRmat *A, dvector *b,
 *                                         dCSRmat *A_diag, REAL *work, dCSRmat *interface_dof,
 *                                         block_dCSRmat *A_gamma)
 *
 * \brief block additive smoother for metric AMG
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_metric_multiplicative_omega_gamma(dvector *u,
                                                      block_dCSRmat *A,
                                                      dvector *b,
                                                      dCSRmat *A_diag,
                                                      REAL *work,
                                                      dCSRmat *interface_dof,
                                                      block_dCSRmat * A_gamma)
{
    // local variables
    dvector r;
    dvector r_gamma;
    dvector e_gamma;

    INT i;

    INT brow = A->brow;

    //--------------------------------------------
    // Gauss-Seidel smoother for the whole matrix
    //--------------------------------------------
    smoother_bdcsr_fgs_fgs(u, A, b, A_diag, work);
    //smoother_bdcsr_bgs_bgs(u, A, b, A_diag, work);

    //printf("done smooth whole matrix\n");

    //--------------------------------
    // smoother on the interface part
    //--------------------------------
    // save right hand side b in the workspace
    array_cp(b->row, b->val, work + b->row);

    //printf("done copy b\n");

    // compute the overall residual
    r.row = b->row;
    r.val = work + b->row;
    bdcsr_aAxpy(-1.0, A, u->val, r.val);

    //printf("done computing r\n");

    // get the residual for the interface part
    r_gamma.row = brow*interface_dof->row;
    r_gamma.val = r.val + r.row;
    for (i=0; i<interface_dof->row; i++) r_gamma.val[i] = r.val[interface_dof->JA[i]];
    for (i=0; i<interface_dof->row; i++) r_gamma.val[interface_dof->row+i] = r.val[A->blocks[0]->row+i];

    //printf("done get r\n");

    // set initial to be zero
    e_gamma.row = brow*interface_dof->row;
    e_gamma.val = r_gamma.val + r_gamma.row;
    dvec_set(e_gamma.row, &e_gamma, 0.0);

    //printf("done set e\n");

    // solve for the interface part
    //block_directsolve_HAZ(A_gamma, &r_gamma, &e_gamma, 3);
    smoother_bdcsr_fgs_fgs(&e_gamma, A_gamma, &r_gamma, NULL, work);
    //smoother_bdcsr_bgs_bgs(&e_gamma, A_gamma, &r_gamma, NULL, work);

    //printf("done solve interface\n");

    //--------------------------------------------
    // update the solution
    //--------------------------------------------
    for (i=0; i<interface_dof->row; i++) u->val[interface_dof->JA[i]] += e_gamma.val[i];
    for (i=0; i<interface_dof->row; i++) u->val[A->blocks[0]->row+i] += e_gamma.val[interface_dof->row+i];

    //printf("done update u\n");

}

/**
 * \fn void smoother_bdcsr_metric_multiplicative_omega_gamma_bsr(dvector *u, block_dCSRmat *A, dvector *b,
 *                                         dCSRmat *A_diag, REAL *work, dCSRmat *interface_dof,
 *                                         dBSRmat *A_gamma, REAL * diaginv)
 *
 * \brief block additive smoother for metric AMG
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_metric_multiplicative_omega_gamma_bsr(dvector *u,
                                                          block_dCSRmat *A,
                                                          dvector *b,
                                                          dCSRmat *A_diag,
                                                          REAL *work,
                                                          dCSRmat *interface_dof,
                                                          dBSRmat * A_gamma,
                                                          REAL *diaginv)
{
    // local variables
    dvector r;
    dvector r_gamma;
    dvector e_gamma;

    INT i;

    INT brow = A->brow;

    //--------------------------------------------
    // Gauss-Seidel smoother for the whole matrix
    //--------------------------------------------
    smoother_bdcsr_fgs_fgs(u, A, b, A_diag, work);
    //smoother_bdcsr_bgs_bgs(u, A, b, A_diag, work);

    //printf("done smooth whole matrix\n");

    //--------------------------------
    // smoother on the interface part
    //--------------------------------
    // save right hand side b in the workspace
    array_cp(b->row, b->val, work + b->row);

    //printf("done copy b\n");

    // compute the overall residual
    r.row = b->row;
    r.val = work + b->row;
    bdcsr_aAxpy(-1.0, A, u->val, r.val);

    //printf("done computing r\n");

    // get the residual for the interface part
    r_gamma.row = brow*interface_dof->row;
    r_gamma.val = r.val + r.row;
    for (i=0; i<interface_dof->row; i++) r_gamma.val[i*brow] = r.val[interface_dof->JA[i]];
    for (i=0; i<interface_dof->row; i++) r_gamma.val[i*brow+1] = r.val[A->blocks[0]->row+i];

    //printf("done get r\n");

    // set initial to be zero
    e_gamma.row = brow*interface_dof->row;
    e_gamma.val = r_gamma.val + r_gamma.row;
    dvec_set(e_gamma.row, &e_gamma, 0.0);

    //printf("done set e\n");

    // solve for the interface part
    smoother_dbsr_gs_ascend(A_gamma, &r_gamma, &e_gamma, diaginv);
    //smoother_dbsr_gs_descend(A_gamma, &r_gamma, &e_gamma, diaginv);

    //printf("done solve interface\n");

    //--------------------------------------------
    // update the solution
    //--------------------------------------------
    for (i=0; i<interface_dof->row; i++) u->val[interface_dof->JA[i]] += e_gamma.val[i*brow];
    for (i=0; i<interface_dof->row; i++) u->val[A->blocks[0]->row+i] += e_gamma.val[i*brow+1];

    //printf("done update u\n");

}


/**
 * \fn void smoother_bdcsr_metric_multiplicative_gamma_omega(dvector *u, block_dCSRmat *A, dvector *b,
 *                                         dCSRmat *A_diag, REAL *work, dCSRmat *interface_dof,
 *                                         block_dCSRmat *A_gamma)
 *
 * \brief block additive smoother for metric AMG
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_metric_multiplicative_gamma_omega(dvector *u,
                                        block_dCSRmat *A,
                                          dvector *b,
                                          dCSRmat *A_diag,
                                          REAL *work,
                                          dCSRmat *interface_dof,
                                          block_dCSRmat * A_gamma)
{
    // local variables
    dvector r;
    dvector r_gamma;
    dvector e_gamma;

    INT i;

    INT brow = A->brow;

    //--------------------------------
    // smoother on the interface part
    //--------------------------------
    // save right hand side b in the workspace
    array_cp(b->row, b->val, work + b->row);

    //printf("done copy b\n");

    // compute the overall residual
    r.row = b->row;
    r.val = work + b->row;
    bdcsr_aAxpy(-1.0, A, u->val, r.val);

    //printf("done computing r\n");

    // get the residual for the interface part
    r_gamma.row = brow*interface_dof->row;
    r_gamma.val = r.val + r.row;
    for (i=0; i<interface_dof->row; i++) r_gamma.val[i] = r.val[interface_dof->JA[i]];
    for (i=0; i<interface_dof->row; i++) r_gamma.val[interface_dof->row+i] = r.val[A->blocks[0]->row+i];

    //printf("done get r\n");

    // set initial to be zero
    e_gamma.row = brow*interface_dof->row;
    e_gamma.val = r_gamma.val + r_gamma.row;
    dvec_set(e_gamma.row, &e_gamma, 0.0);

    //printf("done set e\n");

    // solve for the interface part
    //block_directsolve_HAZ(A_gamma, &r_gamma, &e_gamma, 3);
    //smoother_bdcsr_fgs_fgs(&e_gamma, A_gamma, &r_gamma, NULL, work);
    smoother_bdcsr_bgs_bgs(&e_gamma, A_gamma, &r_gamma, NULL, work);

    //printf("done solve interface\n");

    //--------------------------------------------
    // update the solution
    //--------------------------------------------
    for (i=0; i<interface_dof->row; i++) u->val[interface_dof->JA[i]] += e_gamma.val[i];
    for (i=0; i<interface_dof->row; i++) u->val[A->blocks[0]->row+i] += e_gamma.val[interface_dof->row+i];

    //printf("done update u\n");

    //--------------------------------------------
    // Gauss-Seidel smoother for the whole matrix
    //--------------------------------------------
    //smoother_bdcsr_fgs_fgs(u, A, b, A_diag, work);
    smoother_bdcsr_bgs_bgs(u, A, b, A_diag, work);

    //printf("done smooth whole matrix\n");

}

/**
 * \fn void smoother_bdcsr_metric_multiplicative_gamma_omega_bsr(dvector *u, block_dCSRmat *A, dvector *b,
 *                                         dCSRmat *A_diag, REAL *work, dCSRmat *interface_dof,
 *                                         dBSRmat *A_gamma, REAL    *diaginv)
 *
 * \brief block additive smoother for metric AMG
 *
 * \param u       Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A       Pointer to block_dBSRmat: the coefficient matrix
 * \param b       Pointer to dvector: the right hand side
 * \param A_diag  Pointer to specified diagonal blocks (e.g., approximation of Schur complements)
 *
 * \author  Xiaozhe Hu
 *
 */
void smoother_bdcsr_metric_multiplicative_gamma_omega_bsr(dvector *u,
                                                          block_dCSRmat *A,
                                                          dvector *b,
                                                          dCSRmat *A_diag,
                                                          REAL *work,
                                                          dCSRmat *interface_dof,
                                                          dBSRmat * A_gamma,
                                                          REAL    *diaginv)
{
    // local variables
    dvector r;
    dvector r_gamma;
    dvector e_gamma;

    INT i;

    INT brow = A->brow;

    //--------------------------------
    // smoother on the interface part
    //--------------------------------
    // save right hand side b in the workspace
    array_cp(b->row, b->val, work + b->row);

    //printf("done copy b\n");

    // compute the overall residual
    r.row = b->row;
    r.val = work + b->row;
    bdcsr_aAxpy(-1.0, A, u->val, r.val);

    //printf("done computing r\n");

    // get the residual for the interface part
    r_gamma.row = brow*interface_dof->row;
    r_gamma.val = r.val + r.row;
    for (i=0; i<interface_dof->row; i++) r_gamma.val[i*brow] = r.val[interface_dof->JA[i]];
    for (i=0; i<interface_dof->row; i++) r_gamma.val[i*brow+1] = r.val[A->blocks[0]->row+i];

    //printf("done get r\n");

    // set initial to be zero
    e_gamma.row = brow*interface_dof->row;
    e_gamma.val = r_gamma.val + r_gamma.row;
    dvec_set(e_gamma.row, &e_gamma, 0.0);

    //printf("done set e\n");

    // solve for the interface part
    //smoother_dbsr_gs_ascend(A_gamma, &r_gamma, &e_gamma, diaginv);
    smoother_dbsr_gs_descend(A_gamma, &r_gamma, &e_gamma, diaginv);

    //printf("done solve interface\n");

    //--------------------------------------------
    // update the solution
    //--------------------------------------------
    for (i=0; i<interface_dof->row; i++) u->val[interface_dof->JA[i]] += e_gamma.val[i*brow];
    for (i=0; i<interface_dof->row; i++) u->val[A->blocks[0]->row+i] += e_gamma.val[i*brow+1];

    //printf("done update u\n");

    //--------------------------------------------
    // Gauss-Seidel smoother for the whole matrix
    //--------------------------------------------
    //smoother_bdcsr_fgs_fgs(u, A, b, A_diag, work);
    smoother_bdcsr_bgs_bgs(u, A, b, A_diag, work);

    //printf("done smooth whole matrix\n");

}





/************************************************************************************************/
// Smoothers for Eigenvalue prolems
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

/**
 * \fn void smoother_dcsr_Schwarz (const Schwarz_data *Schwarz, 
 *                                 dvector *x_in,
 *                                 const dvector *b_in,
 *                                 const INT maxiter)
 *
 * \brief Schwarz smoothers: forward, backward, symmetric, local LU or
 *                           global LU.
 *
 * \param Schwarz Pointer to data about the Schwarz method, eg, blocks and:
 *                If Schwarz->Schwarz_type=1  (forward iter w global LU)
 *                If Schwarz->Schwarz_type=2  (backward iter w global LU)
 *                If Schwarz->Schwarz_type=3  (symmetric iter w global LU)
 *                If Schwarz->Schwarz_type=11 (forward iter w local LU)
 *                If Schwarz->Schwarz_type=12 (backward iter w local LU)
 *                If Schwarz->Schwarz_type=13 (symmetric iter w local LU)
 * \param x_in    Pointer to solution vector
 * \param b_in    Pointer to right hand
 * \param maxiter maximum number of iterations. 
 *
 * \note Needs improvment -- Xiaozhe
 *
 * \note Modified (ltz 20230131)
 */
void smoother_dcsr_Schwarz(const Schwarz_data *Schwarz,		\
			   dvector *x_in,			\
			   const dvector *b_in,			\
			   const INT maxiter)
{
  INT nzloc,nloc,ij,j,iter,jk,i,iloc,ibl0,ibl1;
  INT kblk,nblk_start,nblk_end, step;
  // local
  dCSRmat A=Schwarz->A;
  REAL *x=x_in->val,*b=b_in->val;
  INT stype=(INT )Schwarz->Schwarz_type;
  INT *mask=Schwarz->mask;
  INT *iblock=Schwarz->iblock;
  INT *jblock=Schwarz->jblock;
  INT block_solver=Schwarz->blk_solver;
  dCSRmat *blk=Schwarz->blk_data;
  //forward method always:)  
  nblk_start=0;
  nblk_end=Schwarz->nblk-1;
  step=1;
  dvector rhsloc=Schwarz->rhsloc1;
  dvector xloc=Schwarz->xloc1;  
  void **numeric = Schwarz->numeric; //
  dCSRmat Aloc=blk[0];// this is not needed if there is global LU.
  if(stype == SCHWARZ_BACKWARD_LOCAL || stype == SCHWARZ_BACKWARD){
    //swap;
    i=nblk_start;
    nblk_start=nblk_end;
    nblk_end=i;
    step=-step;
  }
  //
  // loop: a counter for symmetric Schwarz which loops twice over the blocks: forward and backwards.
  INT loop; 
  /*  
      fprintf(stdout,"\nstype=%d,%d,%d,%d::%d,%d,%d\n",			
      stype,SCHWARZ_FORWARD,SCHWARZ_BACKWARD,SCHWARZ_SYMMETRIC,		
      SCHWARZ_FORWARD_LOCAL,SCHWARZ_BACKWARD_LOCAL,SCHWARZ_SYMMETRIC_LOCAL);
  */
  if(stype==SCHWARZ_FORWARD ||		\
     stype==SCHWARZ_BACKWARD ||		\
     stype==SCHWARZ_SYMMETRIC){
    for(iter=0;iter<maxiter; iter++){
      loop=2;
      while(loop>0){
	kblk=nblk_start;
	while(kblk!=nblk_end+step){
	  //	  fprintf(stdout,"\nkblk(lu_global)=%d (nblk_end=%d)",kblk,nblk_end);fflush(stdout);
	  ibl0 = iblock[kblk];
	  ibl1 = iblock[kblk+1];
	  iloc=0;
	  nloc=ibl1-ibl0;
	  for(jk = ibl0;jk<ibl1;++jk){
	    i = jblock[jk];
	    mask[i] = iloc;
	    rhsloc.val[iloc] = b[i];
	    iloc++;
	  }
	  /* fprintf(stdout,"\nkblk(%d)=",kblk); */
	  /* for(jk=ibl0;jk<ibl1;++jk){ */
	  /*   i=jblock[jk]; //node number */
	  /*   fprintf(stdout," %d ",i);fflush(stdout); */
	  /* } */
	  // Solve each block
	  iloc=0;
	  for(jk=ibl0;jk<ibl1;++jk){
	    i=jblock[jk]; //node number
	    for(ij=A.IA[i];ij<A.IA[i+1];++ij){
	      j=A.JA[ij];
	      if(mask[j]<0) {
		rhsloc.val[iloc] -= A.val[ij]*x[j];
	      }
	    }
	    iloc++;
	  }
	  switch (block_solver) {
	  case SOLVER_UMFPACK: {
	    /* use UMFPACK direct solver on each block */
	    hazmath_solve(&blk[kblk], &rhsloc, &xloc, numeric[kblk], 0);
	    break;
	  }
	    //#endif
	  default:
	    /* use iterative solver on each block */
	    xloc.row = Aloc.row;
	    rhsloc.row = Aloc.row;
	    memset(xloc.val,0,xloc.row*sizeof(REAL));
	    dcsr_pvgmres(&Aloc, &rhsloc, &xloc, NULL, 1e-8, 20, 20, 1, 0);
	  }
	  iloc=0;
	  for(jk=ibl0;jk<ibl1;++jk){
	    i = jblock[jk];
	    mask[i] = -1;
	    x[i] = xloc.val[iloc];
	    iloc++;
	  }
	  /*...  DONE  */
	  kblk+=step;
	}
	if(stype==SCHWARZ_FORWARD || stype==SCHWARZ_BACKWARD){
	  break;
	} else { //?stype==SCHWARZ_SYMMETRIC
	  jk=nblk_end; // swap start and end;
	  nblk_end=nblk_start;
	  nblk_start=jk;
	  step=-step;// change step;
	  loop--; // decrease loop so once we pass twice through this we exit. 
	}
      }
    }
  } else {
    dCSRmat blk_tran=dcsr_create(Schwarz->maxbs,Schwarz->maxbs,Schwarz->maxbnnz);
    for(iter=0;iter<maxiter; iter++){
      loop=2;
      while(loop>0){
	// this loop may go forward and backward if needed
	kblk=nblk_start;
	while(kblk!=(nblk_end+step)){
	  //	  fprintf(stdout,"\nkblk(lu_local)=%d (nblk_end=%d)",kblk,nblk_end);fflush(stdout);
	  ibl0 = iblock[kblk];
	  ibl1 = iblock[kblk+1];
	  iloc=0;
	  nloc=ibl1-ibl0;
	  for(jk = ibl0;jk<ibl1;++jk){
	    i = jblock[jk];
	    mask[i] = iloc;
	    rhsloc.val[iloc] = b[i];
	    iloc++;
	  }
	  Aloc.row=nloc;
	  Aloc.col=nloc;
	  nzloc=0;
	  iloc=0;
	  Aloc.IA[iloc]=nzloc;
	  for(jk=ibl0;jk<ibl1;++jk){
	    i=jblock[jk]; //node number
	    for(ij=A.IA[i];ij<A.IA[i+1];++ij){
	      j=A.JA[ij];
	      if(mask[j]<0) {
		rhsloc.val[iloc] -= A.val[ij]*x[j];
	      }else {
		Aloc.JA[nzloc]=mask[j];
		Aloc.val[nzloc]=A.val[ij];
		nzloc++;
	      }
	    }
	    iloc++;
	    Aloc.IA[iloc]=nzloc;
	  }
	  Aloc.nnz=Aloc.IA[nloc];
	  /*... A=LU*/
	  switch(block_solver) {	
	  case SOLVER_UMFPACK: 
	    /* use direct solver on each block locally */
	    dcsr_transz(&Aloc, NULL, &blk_tran);
	    dcsr_cp(&blk_tran, &Aloc);
	    //	dcsr_free(&blk_tran);
	    //printf("size of block %d: nrow=%d, nnz=%d\n",i, blk[i].row, blk[i].nnz);
	    numeric[0] = hazmath_factorize(&Aloc, 0);
	    hazmath_solve(&Aloc, &rhsloc, &xloc, numeric[0], 0);
	    break;      
	  default:
	    /* use iterative solver on each block */
	    xloc.row = Aloc.row;
	    rhsloc.row = Aloc.row;
	    memset(xloc.val,0,xloc.row*sizeof(REAL));
	    dcsr_pvgmres(&Aloc, &rhsloc, &xloc, NULL, 1e-8, 20, 20, 1, 0);
	  }
	  iloc=0;
	  for(jk=ibl0;jk<ibl1;++jk){
	    i = jblock[jk];
	    mask[i] = -1;
	    x[i] = xloc.val[iloc];
	    iloc++;
	  }
	  /*...  DONE  */
	  kblk+=step;
	}
	if(stype==SCHWARZ_FORWARD_LOCAL || stype==SCHWARZ_BACKWARD_LOCAL){
	  break;
	} else { //?stype==SCHWARZ_SYMMETRIC_LOCAL
	  jk=nblk_end; // swap start and end;
	  nblk_end=nblk_start;
	  nblk_start=jk;
	  step=-step;// change step;
	  loop--; // decrease loop so once we pass twice through this we exit. 
	}
      }
    }
    dcsr_free(&blk_tran);
  }
  return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
