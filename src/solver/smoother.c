/*! \file src/solver/smoother.c
 *
 *  Smoothers
 *
 *  Created by James Adler and Xiaozhe Hu on 12/25/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 */

#include "hazmath.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/


void smoother_dcsr_jacobi (dvector *u,
                                const INT i_1,
                                const INT i_n,
                                const INT s,
                                dCSRmat *A,
                                dvector *b,
                                INT L)
{
    /**
     * \fn void smoother_dcsr_jacobi (dvector *u, const INT i_1, const INT i_n,
     *                                     const INT s, dCSRmat *A, dvector *b, INT L)
     *
     * \brief Jacobi method as a smoother
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
                //if (ABS(d[i])>SMALLREAL) uval[i]= w*t[i]/d[i];
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
                //if (ABS(d[i])>SMALLREAL) uval[i]=t[i]/d[i];
                if (ABS(d[i])>SMALLREAL) uval[i]=(1-w)*uval[i]+ w*t[i]/d[i];
            }
        }
        
    } // end while
    
    free(t);
    free(d);
    
    return;
}


void smoother_dcsr_gs (dvector *u,
                            const INT i_1,
                            const INT i_n,
                            const INT s,
                            dCSRmat *A,
                            dvector *b,
                            INT L)
{
    /**
     * \fn void smoother_dcsr_gs (dvector *u, const INT i_1, const INT i_n,
     *                                 const INT s, dCSRmat *A, dvector *b, INT L)
     *
     * \brief Gauss-Seidel method as a smoother
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


void smoother_dcsr_gs_cf (dvector *u,
                               dCSRmat *A,
                               dvector *b,
                               INT L,
                               INT *mark,
                               const INT order)
{
    /**
     * \fn void smoother_dcsr_gs_cf (dvector *u, dCSRmat *A, dvector *b, INT L,
     *                                    INT *mark, const INT order)
     *
     * \brief Gauss-Seidel smoother with C/F ordering for Au=b
     *
     * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
     * \param A      Pointer to dBSRmat: the coefficient matrix
     * \param b      Pointer to dvector: the right hand side
     * \param L      Number of iterations
     * \param mark   C/F marker array
     * \param order  C/F ordering: -1: F-first; 1: C-first
     *
     */
    
    const INT    nrow = b->row; // number of rows
    const INT   *ia = A->IA, *ja = A->JA;
    const REAL  *aj = A->val, *bval = b->val;
    REAL        *uval = u->val;
    
    INT i,j,k,begin_row,end_row;
    REAL t,d=0.0;
    
    // F-point first, C-point second
    if (order == FPFIRST) {
        
        while (L--) {
            for (i = 0; i < nrow; i ++) {
                if (mark[i] != 1) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 01/17/2013
                    d = aj[begin_row];
                    for (k = begin_row+1; k < end_row; k ++) {
                        j = ja[k];
                        t -= aj[k]*uval[j];
                    } // end for k
#else
                    for (k = begin_row; k < end_row; k ++) {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j];
                        else d = aj[k];
                    } // end for k
#endif // end if DIAG_PREF
                    if (ABS(d) > SMALLREAL) uval[i] = t/d;
                }
            } // end for i
            
            for (i = 0; i < nrow; i ++) {
                if (mark[i] == 1) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 01/17/2013
                    d = aj[begin_row];
                    for (k = begin_row+1; k < end_row; k ++) {
                        j = ja[k];
                        t -= aj[k]*uval[j];
                    } // end for k
#else
                    for (k = begin_row; k < end_row; k ++) {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j];
                        else d = aj[k];
                    } // end for k
#endif // end if DIAG_PREF
                    if (ABS(d) > SMALLREAL) uval[i] = t/d;
                }
            } // end for i
        } // end while
        
    }
    
    // C-point first, F-point second
    else {
        
        while (L--) {
            for (i = 0; i < nrow; i ++)  {
                if (mark[i] == 1) {
                    t = bval[i];
                    begin_row = ia[i],end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 09/22/2012
                    d = aj[begin_row];
                    for (k = begin_row+1; k < end_row; k ++) {
                        j = ja[k];
                        t -= aj[k]*uval[j];
                    } // end for k
#else
                    for (k = begin_row; k < end_row; k ++) {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j];
                        else d = aj[k];
                    } // end for k
#endif // end if DIAG_PREF
                    if (ABS(d) > SMALLREAL) uval[i] = t/d;
                }
            } // end for i
        
            for (i = 0; i < nrow; i ++) {
                if (mark[i] != 1) {
                    t = bval[i];
                    begin_row = ia[i],end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 09/22/2012
                    d = aj[begin_row];
                    for (k = begin_row+1; k < end_row; k ++) {
                        j = ja[k];
                        t -= aj[k]*uval[j];
                    } // end for k
#else
                    for (k = begin_row; k < end_row; k ++) {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j];
                        else d = aj[k];
                    } // end for k
#endif // end if DIAG_PREF
                    if (ABS(d) > SMALLREAL) uval[i] = t/d;
                }
            } // end for i
        } // end while
        
    } // end if order
    
    return;
}


void smoother_dcsr_sgs (dvector *u,
                             dCSRmat *A,
                             dvector *b,
                             INT L)
{
    /**
     * \fn void smoother_dcsr_sgs (dvector *u, dCSRmat *A, dvector *b, INT L)
     *
     * \brief Symmetric Gauss-Seidel method as a smoother
     *
     * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
     * \param A      Pointer to dBSRmat: the coefficient matrix
     * \param b      Pointer to dvector: the right hand side
     * \param L      Number of iterations
     *
     * \author Xiaozhe Hu
     * \date   10/26/2010
     *
     */
    
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

void smoother_dcsr_sor (dvector *u,
                             const INT i_1,
                             const INT i_n,
                             const INT s,
                             dCSRmat *A,
                             dvector *b,
                             INT L,
                             const REAL w)
{
    /**
     * \fn void smoother_dcsr_sor (dvector *u, const INT i_1, const INT i_n, const INT s,
     *                                  dCSRmat *A, dvector *b, INT L, const REAL w)
     *
     * \brief SOR method as a smoother
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
     * \author Xiaozhe Hu
     * \date   10/26/2010
     *
     */
    
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


void smoother_dcsr_sor_cf (dvector *u,
                                dCSRmat *A,
                                dvector *b,
                                INT L,
                                const REAL w,
                                INT *mark,
                                const INT order )
{
    /**
     * \fn void smoother_dcsr_sor_cf (dvector *u, dCSRmat *A, dvector *b, INT L,
     *                                     const REAL w, INT *mark, const INT order)
     *
     * \brief SOR smoother with C/F ordering for Au=b
     *
     * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
     * \param A      Pointer to dBSRmat: the coefficient matrix
     * \param b      Pointer to dvector: the right hand side
     * \param L      Number of iterations
     * \param w      Over-relaxation weight
     * \param mark   C/F marker array
     * \param order  C/F ordering: -1: F-first; 1: C-first
     *
     */
    const INT    nrow = b->row; // number of rows
    const INT   *ia = A->IA, *ja=A->JA;
    const REAL  *aj = A->val,*bval=b->val;
    REAL        *uval = u->val;
    
    // local variables
    INT    i,j,k,begin_row,end_row;
    REAL   t,d=0.0;
    
    // F-point first
    if (order == -1) {
        while (L--) {
            for (i = 0; i < nrow; i ++) {
                if (mark[i] == 0 || mark[i] == 2) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j];
                        else d = aj[k];
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
            } // end for i
            
            for (i = 0; i < nrow; i ++) {
                if (mark[i] == 1) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j];
                        else d = aj[k];
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
            } // end for i
        } // end while
    }
    else {
        while (L--) {
            for (i = 0; i < nrow; i ++) {
                if (mark[i] == 1) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j];
                        else d = aj[k];
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
            } // end for i
            
            for (i = 0; i < nrow; i ++) {
                if (mark[i] != 1) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j];
                        else d = aj[k];
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
            } // end for i

        } // end while
    }
    
    return;
}



void smoother_dcsr_L1diag (dvector *u,
                                const INT i_1,
                                const INT i_n,
                                const INT s,
                                dCSRmat *A,
                                dvector *b,
                                INT L)
{
    /**
     * \fn void smoother_dcsr_L1diag (dvector *u, const INT i_1, const INT i_n, const INT s,
     *                                     dCSRmat *A, dvector *b, INT L)
     *
     * \brief Diagonal scaling (using L1 norm) as a smoother
     *
     * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
     * \param i_1    Starting index
     * \param i_n    Ending index
     * \param s      Increasing step
     * \param A      Pointer to dBSRmat: the coefficient matrix
     * \param b      Pointer to dvector: the right hand side
     * \param L      Number of iterations
     *
     * \author Xiaozhe Hu, James Brannick
     * \date   01/26/2011
     *
     */
    
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


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
