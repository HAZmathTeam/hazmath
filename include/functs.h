/*******************************************************************/  
/* This header file was automatically generated with "make fheaders".   */
/* WARNING: DO NOT EDIT!!!                               */  
/*******************************************************************/  
#include "macro.h"
#include "grid.h"
#include "sparse.h"
#include "vec.h"
#include "fem.h"
void PX_H1_basis(REAL *p,REAL *dpx,REAL *dpy,REAL *dpz,REAL x,REAL y,REAL z,INT *dof,INT porder,trimesh mesh) ;
void quad_tri_2D_2der(REAL *p,REAL *dpx,REAL *dpy,REAL *dpxx,REAL *dpyy,REAL *dpxy,REAL x,REAL y,REAL z,INT *dof,INT porder,trimesh mesh) ;
void create_fespace(fespace *FE,trimesh* mesh,INT FEtype);
void free_fespace(fespace* FE);
void get_P2(fespace* FE,trimesh* mesh) ;
void dump_el_dof(FILE* fid,iCSRmat *el_dof) ;
void dump_fespace(fespace *FE) ;
void allocateqcoords(qcoordinates *A,INT nq1d,INT nelm,INT mydim);
void free_qcoords(qcoordinates* A);
qcoordinates get_quadrature(trimesh *mesh,INT nq1d) ;
void quad_elm(qcoordinates *cqelm,trimesh *mesh,INT nq1d,INT elm) ;
qcoordinates get_quadrature_edge(trimesh *mesh,INT nq1d) ;
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
void isboundary_ed3D(iCSRmat ed_v,INT nedge,coordinates cv,INT *nbedge,INT *v_bdry,INT *ed_bdry) ;
iCSRmat get_el_ed(iCSRmat el_v,iCSRmat ed_v) ;
void allocatecoords(coordinates *A,INT ndof,INT mydim);
void free_coords(coordinates* A);
void edge_stats_all(REAL *ed_len,REAL *ed_tau,REAL *ed_mid,coordinates cv,iCSRmat ed_v,INT dim) ;
void get_face_ordering(INT el_order,INT dim,INT f_order,INT *fel_order);
void get_face_maps(iCSRmat el_v,INT el_order,INT nface,INT dim,INT f_order,iCSRmat *el_f,INT *f_bdry,INT *nbface,iCSRmat *f_v,INT *fel_order);
void find_facenumber(iCSRmat el_v,INT* fel_order,INT elm,INT* nd,INT dim,INT *f_num)          ;
void face_stats(REAL *f_area,REAL *f_mid,REAL *f_norm,INT* fel_order,trimesh mesh) ;
void get_el_mid(REAL *el_mid,iCSRmat el_v,coordinates cv,INT dim) ;
void get_el_vol(REAL *el_vol,iCSRmat el_v,coordinates cv,INT dim,INT v_per_elm) ;
void free_mesh(trimesh* mesh);
void get_incidence_row(INT row,iCSRmat *fem_map,INT* thisrow);
void dump_coords(FILE* fid,coordinates *c) ;
void rveci_(FILE *fp, INT *vec, INT *nn)       ;
void rvecd_(FILE *fp,  REAL *vec, INT *nn);
void baddimension()          ;
void getinput(char* gridfile,REAL* fpar,INT* ipar);
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
void icsr_trans_1 (iCSRmat *A,
                 iCSRmat *AT);
void dcsr_mxm (dCSRmat *A,
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
void icsr_mxm_symb_max (iCSRmat *A,
                        iCSRmat *B,
                        iCSRmat *C,
                        INT multmax);
void icsr_mxm_symb_max_1 (iCSRmat *A,
                        iCSRmat *B,
                        iCSRmat *C,
                        INT multmax);

/* End of header file */
