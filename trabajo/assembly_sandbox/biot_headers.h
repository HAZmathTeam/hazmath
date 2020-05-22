/*STRUCTS*/
/******************************************************************************************************/
// local structure to hold vectors and other things locally on an
// element and ungeneralize whatever is general here.
typedef struct simplex_data{
  INT dim; //vertices in an element
  INT ns; //elements in an element; this must be 1;
  INT nv; //vertices in an element
  INT ne;//edges in an element;
  INT nf;
  INT nspaces;// how many FE spaces;
  INT num_dofs;// how dofs on a simplexnum_dofs

  INT *fetypes; // types of elements (could be only one type);
  INT *ndofs;// size is the same as fespaces
  REAL vol;// volume
  REAL *x; // coordinates of element vertices used to map the
  INT *v; // vertices in element
  INT *f;// faces in an element;
  INT *fv;// face to vertex map
  INT *dofs;
  INT *gdofs;
  //! normal vector on face
  REAL* f_norm;
  REAL *mid;// barycenter
  REAL* f_area; // areas of faces.
  REAL* f_mid;
  //edge data
  INT *ed;// edges in an element
  INT *edv;//edge_vertex
  INT *fed;//face_edge
  REAL *ed_len;
  REAL *ed_tau;
  REAL *ed_mid;
  //! indicates a flag for face such as whether a face is on boundary
  INT *f_flag;
  //! indicates a flag for element such as in what domain the element is.
  INT flag;
  REAL *p;
  REAL *dp;
  REAL *ddp;
  void *wrk;
} simplex_data;
/****************************************************************************/
typedef struct local_vec {
  REAL *b; // bubbles
  REAL *u; // displacements
  REAL *w; // Darcy velocity
  REAL *p; // pressure
} local_vec; /**/
/*************************************************************************/
typedef struct estimator{
  REAL *bulk;
  REAL *face;
} estimator;
/*************************************************************************/
/*FUNCTIONS*/
/*************************************************************************/
/* Coefficients and such */
void true_sol2D(REAL *val,REAL* x, REAL time,void *param);
void Dtrue_sol2D(REAL *val,REAL* x, REAL time,void *param);
void rhs_func(REAL *val,REAL* x, REAL time,void *param);
void conductivity2D(REAL *val,REAL* x, REAL time,void *param);
void get_mu(REAL *val,REAL* x, REAL time,void *param);
void get_lam(REAL *val,REAL* x, REAL time,void *param);
void get_alpha(REAL *val,REAL* x, REAL time,void *param);
/*************************************************************************/
INT calc_ndofs(INT fetype,INT dim,					\
	       const INT nv,const INT ne,const INT nf,const INT ns);
simplex_data *simplex_data_init(INT dim,void *fe_in, const INT isit_block);
void simplex_data_update(INT snum,simplex_data *s, mesh_struct *mesh,	\
			 void *fe_in,const INT isit_block);
local_vec **dof_data_init(simplex_data *s);
void dof_data_update(local_vec **ue,					\
		     INT snum, simplex_data *s,				\
		     dvector *uh0, dvector *uh1, const REAL dt);
void simplex_data_print(INT snum, simplex_data *s, local_vec **ue);
void simplex_data_free(simplex_data *s);
void dof_data_free(local_vec **ue);
void zquad_elm(qcoordinates *cqelm,simplex_data *splex,INT nq1d);
void zquad_face(qcoordinates *cqbdry,INT nq1d, INT dim, REAL *xf, REAL farea);
/*interpolations*/
void zpx_h1_basis(REAL *p,REAL *dp,REAL *x,INT  dim, REAL *xel,INT  porder);
void zrt_basis(REAL *phi,REAL *dphi,REAL *x,simplex_data *splex,	\
	       REAL *wrk);
void zbubble_face_basis(REAL *phi, REAL *dphi, REAL *x,
			simplex_data *splex, REAL *wrk);
INT zget_fem_basis(REAL *phi,REAL *dphi,REAL *x,		\
		   INT  fetype,simplex_data *splex,REAL *wrk);
void zfe_interp(REAL* val,REAL *u,REAL* x,			\
		INT fetype,simplex_data *splex, REAL *wrk);
void zfe_dinterp(REAL* val,REAL *u,REAL *x,
		 INT fetype, simplex_data *splex, REAL *wrk);
void zblockfe_interp(REAL* val,REAL *u,REAL* x,simplex_data *splex);
void zblockfe_dinterp(REAL* val,REAL *u,REAL* x,
		      simplex_data *splex);
void zbubble_hessian2d(REAL *depsp,REAL *dtrepsp,		\
		       REAL *u, REAL *x, simplex_data *splex);
/*EO_HEADER*/
