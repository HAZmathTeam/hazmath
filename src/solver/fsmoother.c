/*! \file src/solver/fsmoother.c
 *
 *  Fractional Smoothers
 *
 *  Created by Ana Budisa 2020-05-08 from original smoother.c (by HAZ people)
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *
 *
 */

#include "hazmath.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/**
 * \fn void smoother_dcsr_fjacobi (dvector *u, const INT i_1, const INT i_n,
 *                                 const INT s, dCSRmat *A, dvector *b, dCSRmat *M, const REAL p, const REAL w, INT L)
 *
 * \brief Fractional Jacobi smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1    Starting index
 * \param i_n    Ending index
 * \param s      Increasing step
 * \param A      Pointer to dCSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param M      Pointer to dCSRmat: the mass matrix
 * \param p      Fractional power/exponent
 * \param w      Relaxation parameter
 * \param L      Number of iterations
 *
 */
void smoother_dcsr_fjacobi(dvector *u,
                           const INT i_1,
                           const INT i_n,
                           const INT s,
                           dCSRmat *A,
                           dvector *b,
                           dCSRmat *M,
                           const REAL p,
                           const REAL w,
                           INT L)
{
    const INT    N = ABS(i_n - i_1)+1;
    const INT   *ia=A->IA, *ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    //REAL        w = 1.0;  //0.8
    // local variables
    INT i,j,k,begin_row,end_row;

    REAL *r = (REAL *)calloc(N,sizeof(REAL)); // b - \sum Aij * uj
    REAL *d = (REAL *)calloc(N,sizeof(REAL)); // jacobi part

    dvector Mdiag_1mp = dvec_create(M->row); // diag of mass matrix

    dcsr_getdiag_pow(0, 1-p, M, &Mdiag_1mp); // get M_ii^(1-p)

    while (L--) {
        if (s>0) {
            for (i=i_1;i<=i_n;i+=s) {
                r[i]=bval[i];
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    r[i]-=aj[k]*uval[j];
                    if (i==j) d[i]=pow(aj[k], p)*Mdiag_1mp.val[i]; //Aii^p * Mii^(1-p)
                }
            }

            for (i=i_1;i<=i_n;i+=s) {
                if (ABS(d[i])>SMALLREAL) uval[i]+= w*r[i]/d[i]; // u^k+1 = u^k + w * d^-1 * r^k
            }
        }
        else {
            for (i=i_1;i>=i_n;i+=s) {
                r[i]=bval[i];
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    r[i]-=aj[k]*uval[j];
                    if (i==j) d[i]=pow(aj[k], p)*Mdiag_1mp.val[i]; //Aii^p * Mii^(1-p)
                }
            }

            for (i=i_1;i>=i_n;i+=s) {
                if (ABS(d[i])>SMALLREAL) uval[i]+= w*r[i]/d[i]; // u^k+1 = u^k + w * d^-1 * r^k
            }
        }

    } // end while

    free(r);
    free(d);
    dvec_free(&Mdiag_1mp);

    return;
}

/**
 * \fn void smoother_dcsr_fgs(dvector *u, const INT i_1, const INT i_n,
 *                            const INT s, dCSRmat *A, dvector *b, dCSRmat *M,
 *                            const REAL p, INT L)
 *
 * \brief Fractional Gauss-Seidel smoother
 *
 * \param u    Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1  Starting index
 * \param i_n  Ending index
 * \param s    Increasing step
 * \param A    Pointer to dBSRmat: the coefficient matrix
 * \param b    Pointer to dvector: the right hand side
 * \param M    Pointer to dCSRmat: the mass matrix
 * \param p    Fractional power/exponent
 * \param L    Number of iterations
 *
 */
void smoother_dcsr_fgs(dvector *u,
                       const INT i_1,
                       const INT i_n,
                       const INT s,
                       dCSRmat *A,
                       dvector *b,
                       dCSRmat *M,
                       const REAL p,
                       INT L)
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;

    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  t,d=0.0;

    dvector Mdiag_1mp = dvec_create(M->row); // diag of mass matrix
    dcsr_getdiag_pow(0, 1-p, M, &Mdiag_1mp); // get M_ii^(1-p)

    if (s > 0) {
        while (L--) {
            for (i=i_1;i<=i_n;i+=s) {
                t = bval[i];
                begin_row=ia[i],end_row=ia[i+1];

#if DIAGONAL_PREF // diagonal first
                d=pow(aj[begin_row], p)*Mdiag_1mp.val[i]; // Aii^p * Mii^(1-p)
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
                    else if (ABS(aj[k])>SMALLREAL) d=pow(aj[k], p)*Mdiag_1mp.val[i]; // Aii^p * Mii^(1-p)
                }
                uval[i]=t/d;
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
                d=pow(aj[begin_row], p)*Mdiag_1mp.val[i]; // Aii^p * Mii^(1-p)
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
                    else if (ABS(aj[k])>SMALLREAL)  d=pow(aj[k], p)*Mdiag_1mp.val[i]; // Aii^p * Mii^(1-p)
                }
                uval[i]=t/d;
#endif
            } // end for i
        } // end while
    } // end if
    return;
}


/**
 * \fn void smoother_dcsr_sgs(dvector *u, dCSRmat *A, dvector *b, dCSRmat *M, REAL p, INT L)
 *
 * \brief Fractional Symmetric Gauss-Seidel smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param M      Pointer to dCSRmat: the mass matrix
 * \param p      Fractional power/exponent
 * \param L      Number of iterations
 *
 *
 */
void smoother_dcsr_fsgs(dvector *u,
                        dCSRmat *A,
                        dvector *b,
                        dCSRmat *M,
                        const REAL p,
                        INT L)
{
    const INT    nm1=b->row-1;
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;

    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  t,d=0;

    dvector Mdiag_1mp = dvec_create(M->row); // diag of mass matrix
    dcsr_getdiag_pow(0, 1-p, M, &Mdiag_1mp); // get M_ii^(1-p)

    while (L--) {
        // forward sweep
        for (i=0;i<=nm1;++i) {
            t=bval[i];
            begin_row=ia[i], end_row=ia[i+1];
            for (k=begin_row;k<end_row;++k) {
                j=ja[k];
                if (i!=j) t-=aj[k]*uval[j];
                else d=pow(aj[k], p)*Mdiag_1mp.val[i]; // Aii^p * Mii^(1-p)
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
                else d=pow(aj[k], p)*Mdiag_1mp.val[i]; // Aii^p * Mii^(1-p)
            } // end for k
            if (ABS(d)>SMALLREAL) uval[i]=t/d;
        } // end for i

    } // end while

    return;
}

