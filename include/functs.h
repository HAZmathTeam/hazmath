/*******************************************************************/  
/* This header file was automatically generated with "make fheaders".   */
/* WARNING: DO NOT EDIT!!!                               */  
/*******************************************************************/  
#include "macro.h"
#include "grid.h"
#include "quad.h"
#include "sparse.h"
#include "vec.h"
void creategrid(FILE *gfid,INT dim,trimesh* mesh) ;
iCSRmat convert_elmnode(INT *element_node,INT nelm,INT nvert,INT nve) ;
void get_nedge(INT *nedge,iCSRmat el_v,INT nv,INT nelm,INT v_per_elm) ;
iCSRmat get_edge_v(INT nedge,iCSRmat el_v,INT nv,INT nelm,INT v_per_elm) ;
void isboundary_v(INT nv,INT *bdry_v,INT *v_bdry,INT nbedge,INT *nbv) ;
void isboundary_ed(iCSRmat ed_v,INT nedge,INT nbedge,INT *bdry_v,INT *ed_bdry) ;
void isboundary_ed3D(iCSRmat ed_v,INT nedge,coordinates cv,INT *nbedge,INT *v_bdry,INT *ed_bdry) ;
void rveci_(FILE *fp, INT *vec, INT *nn)       ;
void rvecd_(FILE *fp,  REAL *vec, INT *nn);
dCSRmat dcsr_create (const INT m,
                     const INT n,
                     const INT nnz);
iCSRmat icsr_create (const INT m,
                     const INT n,
                     const INT nnz);
void dcsr_free (dCSRmat *A);
void icsr_free (iCSRmat *A);
INT dcsr_trans (dCSRmat *A,
                dCSRmat *AT);
void icsr_trans (iCSRmat *A,
                 iCSRmat *AT);
void dcsr_mxm (dCSRmat *A,
               dCSRmat *B,
               dCSRmat *C);
void icsr_mxm (iCSRmat *A,
               iCSRmat *B,
               iCSRmat *C);
void icsr_mxm_symb (iCSRmat *A,
                    iCSRmat *B,
                    iCSRmat *C);
void icsr_mxm_symb_max (iCSRmat *A,
                        iCSRmat *B,
                        iCSRmat *C,
                        INT multmax);

/* End of header file */
