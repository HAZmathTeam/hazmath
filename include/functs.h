/*******************************************************************/  
/* This header file was automatically generated with "make fheaders".   */
/* WARNING: DO NOT EDIT!!!                               */  
/*******************************************************************/  
#include "macro.h"
#include "grid.h"
#include "sparse.h"
#include "vec.h"
#include "fem.h"
#include "solver.h"
void assemble_DuDv_global(dCSRmat* A,dvector *b,fespace *FE,trimesh *mesh,qcoordinates *cq,void (*rhs)(REAL *,REAL *,REAL),void (*bc)(REAL *,REAL *,REAL),void (*coeff)(REAL *,REAL *,REAL),REAL time) ;
void assemble_DuDv_local(REAL* ALoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT elm,void (*coeff)(REAL *,REAL *,REAL),REAL time) ;
void FEM_RHS_Local(REAL* bLoc,fespace *FE,trimesh *mesh,qcoordinates *cq,INT *dof_on_elm,INT elm,void (*rhs)(REAL *,REAL *,REAL),REAL time) ;
void stiffG_nnz(dCSRmat *A, fespace *FE) ;
void stiffG_cols(dCSRmat *A, fespace *FE) ;
void LocaltoGlobal(INT *dof_on_elm,fespace *FE,dvector *b,dCSRmat *A,REAL *ALoc,REAL *bLoc) ;
void PX_H1_basis(REAL *p,REAL *dpx,REAL *dpy,REAL *dpz,REAL x,REAL y,REAL z,INT *dof,INT porder,trimesh *mesh) ;
void quad_tri_2D_2der(REAL *p,REAL *dpx,REAL *dpy,REAL *dpxx,REAL *dpyy,REAL *dpxy,REAL x,REAL y,REAL z,INT *dof,INT porder,trimesh *mesh) ;
void initialize_fespace(fespace *FE) ;
void create_fespace(fespace *FE,trimesh* mesh,INT FEtype);
void free_fespace(fespace* FE);
void get_P2(fespace* FE,trimesh* mesh) ;
void dump_el_dof(FILE* fid,iCSRmat *el_dof) ;
void dump_fespace(fespace *FE) ;
struct qcoordinates *allocateqcoords(INT nq1d,INT nelm,INT mydim);
void free_qcoords(qcoordinates* A);
qcoordinates* get_quadrature(trimesh *mesh,INT nq1d) ;
void quad_elm(qcoordinates *cqelm,trimesh *mesh,INT nq1d,INT elm) ;
qcoordinates* get_quadrature_edge(trimesh *mesh,INT nq1d) ;
void quad_edge(qcoordinates *cqedge,trimesh *mesh,INT nq1d,INT edge) ;
void dump_qcoords(qcoordinates *q) ;
void quad1d(REAL *gaussp, REAL *gaussc, INT ng1d);
void triquad_(REAL *gp, REAL *gc, INT ng1d);
void tetquad_(REAL *gp, REAL *gc, INT ng1d);
void creategrid(FILE *gfid,INT dim,INT nholes,trimesh* mesh) ;
void initialize_mesh(trimesh* mesh) ;
iCSRmat convert_elmnode(INT *element_node,INT nelm,INT nv,INT nve) ;
void get_nedge(INT* nedge, iCSRmat el_v) ;
iCSRmat get_edge_v(INT nedge,iCSRmat el_v) ;
void isboundary_v(INT nv,INT *bdry_v,INT *v_bdry,INT nbedge,INT *nbv) ;
void isboundary_ed(iCSRmat ed_v,INT nedge,INT nbedge,INT *bdry_v,INT *ed_bdry) ;
void isboundary_ed3D(iCSRmat ed_v,INT nedge,coordinates *cv,INT *nbedge,INT *v_bdry,INT *ed_bdry) ;
iCSRmat get_el_ed(iCSRmat el_v,iCSRmat ed_v) ;
struct coordinates *allocatecoords(INT ndof,INT mydim);
void free_coords(coordinates* A);
void edge_stats_all(REAL *ed_len,REAL *ed_tau,REAL *ed_mid,coordinates *cv,iCSRmat ed_v,INT dim) ;
void get_face_ordering(INT el_order,INT dim,INT f_order,INT *fel_order);
void get_face_maps(iCSRmat el_v,INT el_order,INT nface,INT dim,INT f_order,iCSRmat *el_f,INT *f_bdry,INT *nbface,iCSRmat *f_v,INT *fel_order);
void find_facenumber(iCSRmat el_v,INT elm,INT* nd,INT dim,INT *f_num)          ;
void face_stats(REAL *f_area,REAL *f_mid,REAL *f_norm,trimesh *mesh) ;
void get_el_mid(REAL *el_mid,iCSRmat el_v,coordinates *cv,INT dim) ;
void get_el_vol(REAL *el_vol,iCSRmat el_v,coordinates *cv,INT dim,INT v_per_elm) ;
void free_mesh(trimesh* mesh);
void get_incidence_row(INT row,iCSRmat *fem_map,INT* thisrow);
void dump_coords(FILE* fid,coordinates *c) ;
INT dcsr_pcg (dCSRmat *A,
              dvector *b,
              dvector *u,
              precond *pc,
              const REAL tol,
              const INT MaxIt,
              const SHORT stop_type,
              const SHORT prtlvl);
void array_null (REAL *x);
void array_set (const INT n,
                REAL *x,
                const REAL val);
void iarray_set (const INT n,
                 INT *x,
                 const INT val);
void array_cp (const INT n,
               REAL *x,
               REAL *y);
void iarray_cp (const INT n, 
                INT *x,
                INT *y);
void array_ax (const INT n,
               const REAL a,
               REAL *x);
void array_axpy (const INT n,
                 const REAL a,
                 REAL *x,
                 REAL *y);
void array_axpyz (const INT n,
                  const REAL a,
                  REAL *x,
                  REAL *y,
                  REAL *z);
void array_axpby (const INT n,
                  const REAL a,
                  REAL *x,
                  const REAL b,
                  REAL *y);
REAL array_dotprod (const INT n,
                    const REAL * x,
                    const REAL * y);
REAL array_norm1 (const INT n,
                  const REAL * x);
REAL array_norm2 (const INT n,
                  const REAL * x);
REAL array_norminf (const INT n,
                    const REAL * x);
void rveci_(FILE *fp, INT *vec, INT *nn)       ;
void rvecd_(FILE *fp,  REAL *vec, INT *nn);
void baddimension()          ;
void getinput(char* gridfile,REAL* fpar,INT* ipar);
void iarray_print(INT *vec, INT n   );
void array_print(REAL *vec, INT n   );
void csr_print_matlab(FILE* fid,dCSRmat *A);
void dvector_print(FILE* fid,dvector *b);
void print_itsolver_info (const INT ptrlvl,
                   const INT stop_type,
                   const INT iter,
                   const REAL relres,
                   const REAL absres,
                   const REAL factor);
void print_cputime (const char *message,
                    const REAL cputime);
void print_message (const INT ptrlvl,
                    const char *message);
dCSRmat dcsr_create (const INT m,
                     const INT n,
                     const INT nnz);
void dcsr_alloc (const INT m,
                      const INT n,
                      const INT nnz,
                      dCSRmat *A);
iCSRmat icsr_create (const INT m,
                     const INT n,
                     const INT nnz);
void dcsr_free (dCSRmat *A);
void icsr_free (iCSRmat *A);
void dcsr_null (dCSRmat *A);
void icsr_null (iCSRmat *A);
dCSRmat dcsr_perm (dCSRmat *A,
                   INT *P);
void icsr_cp (iCSRmat *A,
              iCSRmat *B);
void dcsr_cp (dCSRmat *A,
              dCSRmat *B);
INT dcsr_trans (dCSRmat *A,
                dCSRmat *AT);
void dcsr_trans_1 (dCSRmat *A,
                   dCSRmat *AT);
void icsr_trans (iCSRmat *A,
                 iCSRmat *AT);
void dcsr_shift (dCSRmat *A,
                 INT offset);
void icsr_trans_1 (iCSRmat *A,
                 iCSRmat *AT);
INT dcsr_add (dCSRmat *A,
              const REAL alpha,
              dCSRmat *B,
              const REAL beta,
              dCSRmat *C);
void dcsr_axm (dCSRmat *A,
               const REAL alpha);
void dcsr_mxv (dCSRmat *A,
                REAL *x,
                REAL *y);
void dcsr_mxv_agg (dCSRmat *A,
                   REAL *x,
                   REAL *y);
void dcsr_aAxpy (const REAL alpha,
                 dCSRmat *A,
                 REAL *x,
                 REAL *y);
void dcsr_aAxpy_agg (const REAL alpha,
                     dCSRmat *A,
                     REAL *x,
                     REAL *y);
REAL dcsr_vmv (dCSRmat *A,
               REAL *x,
               REAL *y);
void dcsr_mxm (dCSRmat *A,
               dCSRmat *B,
               dCSRmat *C);
void dcsr_mxm_1 (dCSRmat *A,
                 dCSRmat *B,
                 dCSRmat *C);
void icsr_mxm (iCSRmat *A,
               iCSRmat *B,
               iCSRmat *C);
void icsr_mxm_1 (iCSRmat *A,
               iCSRmat *B,
               iCSRmat *C);
void icsr_mxm_symb (iCSRmat *A,
                    iCSRmat *B,
                    iCSRmat *C);
void icsr_mxm_symb_1 (iCSRmat *A,
                 iCSRmat *B,
                 iCSRmat *C);
void icsr_mxm_symb_max (iCSRmat *A,
                        iCSRmat *B,
                        iCSRmat *C,
                        INT multmax);
void icsr_mxm_symb_max_1 (iCSRmat *A,
                        iCSRmat *B,
                        iCSRmat *C,
                        INT multmax);
INT dvec_isnan (dvector *u);
dvector dvec_create (const INT m);
ivector ivec_create (const INT m);
void dvec_alloc (const INT m,
                 dvector *u);
void ivec_alloc (const INT m,
                 ivector *u);
void dvec_free (dvector *u);
void ivec_free (ivector *u);
void dvec_null (dvector *x);
void dvec_rand (const INT n,
                dvector *x);
void dvec_set (INT n,
               dvector *x,
               REAL val);
void ivec_set (const INT m,
               ivector *u);
void dvec_cp (dvector *x,
                   dvector *y);
REAL dvec_maxdiff (dvector *x,
                   dvector *y);
void dvec_symdiagscale (dvector *b,
                             dvector *diag);
void dvec_axpy (const REAL a,
                dvector *x,
                dvector *y);
void dvec_axpyz(const REAL a,
                dvector *x,
                dvector *y,
                dvector *z);
REAL dvec_dotprod (dvector *x,
                   dvector *y);
REAL dvec_relerr (dvector *x,
                  dvector *y);
REAL dvec_norm1 (dvector *x);
REAL dvec_norm2 (dvector *x);
REAL dvec_norminf (dvector *x);

/* End of header file */
