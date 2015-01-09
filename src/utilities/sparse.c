//
//  sparse.c
//  
//
//  Created by Hu, Xiaozhe on 1/9/15.
//
//

#include <math.h>
#include <time.h>

#include "sparse.h"
#include "vec.h"


dCSRmat dcsr_create (const INT m,
                     const INT n,
                     const INT nnz)
{
    dCSRmat A;
    
    if ( m > 0 ) {
        A.IA = (INT *)calloc(m+1, sizeof(INT));
    }
    else {
        A.IA = NULL;
    }
    
    if ( n > 0 ) {
        A.JA = (INT *)calloc(nnz, sizeof(INT));
    }
    else {
        A.JA = NULL;
    }
    
    if ( nnz > 0 ) {
        A.val = (REAL *)calloc(nnz, sizeof(REAL));
    }
    else {
        A.val = NULL;
    }
    
    A.row = m; A.col = n; A.nnz = nnz;
    
    return A;
}


iCSRmat icsr_create (const INT m,
                     const INT n,
                     const INT nnz)
{
    iCSRmat A;
    
    if ( m > 0 ) {
        A.IA = (INT *)calloc(m+1, sizeof(INT));
    }
    else {
        A.IA = NULL;
    }
    
    if ( n > 0 ) {
        A.JA = (INT *)calloc(nnz, sizeof(INT));
    }
    else {
        A.JA = NULL;
    }
    
    if ( nnz > 0 ) {
        A.val = (INT *)calloc(nnz, sizeof(INT));
    }
    else {
        A.val = NULL;
    }
    
    A.row = m; A.col = n; A.nnz = nnz;
    
    return A;
}


void dcsr_free (dCSRmat *A)
{
    if ( A == NULL ) return;
    
    free(A->IA);  A->IA  = NULL;
    free(A->JA);  A->JA  = NULL;
    free(A->val); A->val = NULL;
}


void icsr_free (iCSRmat *A)
{
    if ( A == NULL ) return;
    
    free(A->IA);  A->IA  = NULL;
    free(A->JA);  A->JA  = NULL;
    free(A->val); A->val = NULL;
}

INT dcsr_trans (dCSRmat *A,
                dCSRmat *AT)
{
    const INT n=A->row, m=A->col, nnz=A->nnz;
    
    // Local variables
    INT i,j,k,p;
    
    AT->row=m;
    AT->col=n;
    AT->nnz=nnz;
    
    AT->IA=(INT*)calloc(m+1,sizeof(INT));
    
    AT->JA=(INT*)calloc(nnz,sizeof(INT));
    
    if (A->val) {
        AT->val=(REAL*)calloc(nnz,sizeof(REAL));
        
    }
    else { AT->val=NULL; }
    
    // first pass: find the Number of nonzeros in the first m-1 columns of A
    // Note: these Numbers are stored in the array AT.IA from 1 to m-1
    
    // fasp_iarray_set(m+1, AT->IA, 0);
    memset(AT->IA, 0, sizeof(INT)*(m+1));
    
    for (j=0;j<nnz;++j) {
        i=A->JA[j]; // column Number of A = row Number of A'
        if (i<m-1) AT->IA[i+2]++;
    }
    
    for (i=2;i<=m;++i) AT->IA[i]+=AT->IA[i-1];
    
    // second pass: form A'
    if (A->val) {
        for (i=0;i<n;++i) {
            INT ibegin=A->IA[i], iend=A->IA[i+1];
            for (p=ibegin;p<iend;p++) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->val[k]=A->val[p];
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
    }
    else {
        for (i=0;i<n;++i) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
            for (p=ibegin;p<iend1;p++) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->IA[j]=k+1;
            } // end for p
        } // end of i
    } // end if
    
    return 0;
}

void icsr_trans (iCSRmat *A,
                 iCSRmat *AT)
{
    const INT n=A->row, m=A->col, nnz=A->nnz, m1=m-1;
    
    // Local variables
    INT i,j,k,p;
    INT ibegin, iend;
    
#if DEBUG_MODE
    printf("### DEBUG: m=%d, n=%d, nnz=%d\n",m,n,nnz);
#endif
    
    AT->row=m; AT->col=n; AT->nnz=nnz;
    
    AT->IA=(INT*)calloc(m+1,sizeof(INT));
    
    AT->JA=(INT*)calloc(nnz,sizeof(INT));
    
    if (A->val) {
        AT->val=(INT*)calloc(nnz,sizeof(INT));
    }
    else {
        AT->val=NULL;
    }
    
    // first pass: find the Number of nonzeros in the first m-1 columns of A
    // Note: these Numbers are stored in the array AT.IA from 1 to m-1
    memset(AT->IA, 0, sizeof(INT)*(m+1));
    
    for (j=0;j<nnz;++j) {
        i=A->JA[j]; // column Number of A = row Number of A'
        if (i<m1) AT->IA[i+2]++;
    }
    
    for (i=2;i<=m;++i) AT->IA[i]+=AT->IA[i-1];
    
    // second pass: form A'
    if (A->val != NULL) {
        for (i=0;i<n;++i) {
            ibegin=A->IA[i], iend=A->IA[i+1];
            for (p=ibegin;p<iend;p++) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->val[k]=A->val[p];
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
    }
    else {
        for (i=0;i<n;++i) {
            ibegin=A->IA[i], iend=A->IA[i+1];
            for (p=ibegin;p<iend;p++) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
    } // end if
}

void dcsr_mxm (dCSRmat *A,
               dCSRmat *B,
               dCSRmat *C)
{
    INT i,j,k,l,count;
    
    INT *JD = (INT *)calloc(B->col,sizeof(INT));
    
    C->row=A->row;
    C->col=B->col;
    C->val = NULL;
    C->JA  = NULL;
    C->IA  = (INT*)calloc(C->row+1,sizeof(INT));
    
    for (i=0;i<B->col;++i) JD[i]=-1;
    
    // step 1: Find first the structure IA of C
    for(i=0;i<C->row;++i) {
        count=0;
        
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<count;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
                
                if (l==count) {
                    JD[count]=B->JA[j];
                    count++;
                }
            }
        }
        C->IA[i+1]=count;
        for (j=0;j<count;++j) {
            JD[j]=-1;
        }
    }
    
    for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];
    
    // step 2: Find the structure JA of C
    INT countJD;
    
    C->JA=(INT*)calloc(C->IA[C->row],sizeof(INT));
    
    for (i=0;i<C->row;++i) {
        countJD=0;
        count=C->IA[i];
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<countJD;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
                
                if (l==countJD) {
                    C->JA[count]=B->JA[j];
                    JD[countJD]=B->JA[j];
                    count++;
                    countJD++;
                }
            }
        }
        
        //for (j=0;j<countJD;++j) JD[j]=-1;
        memset(JD, -1, sizeof(INT)*countJD);
    }
    
    free(JD);
    
    // step 3: Find the structure A of C
    C->val=(REAL*)calloc(C->IA[C->row],sizeof(REAL));
    
    for (i=0;i<C->row;++i) {
        for (j=C->IA[i];j<C->IA[i+1];++j) {
            C->val[j]=0;
            for (k=A->IA[i];k<A->IA[i+1];++k) {
                for (l=B->IA[A->JA[k]];l<B->IA[A->JA[k]+1];l++) {
                    if (B->JA[l]==C->JA[j]) {
                        C->val[j]+=A->val[k]*B->val[l];
                    } // end if
                } // end for l
            } // end for k
        } // end for j
    }    // end for i
    
    C->nnz = C->IA[C->row]-C->IA[0];
    
}
