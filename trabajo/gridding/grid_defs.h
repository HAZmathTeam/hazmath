#ifndef REAL
#define REAL double
#endif

#ifndef INT
#define INT int
#endif


#ifndef FILENAMESIZE
#define FILENAMESIZE  1024
#endif

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE  0
#endif

#ifndef nula
#define nula  -1
#endif

#ifndef uint
#define uint unsigned int
#endif

void haz_neigh(INT ns,INT nv,INT n,INT *sv,INT *stos);

FILE *HAZ_fopen1( const char *fname, const char *mode );

/********************************************************************/
subscomplex *haz_subscomplex_init(scomplex *sc);
INT haz_refine_simplex(scomplex *sc, const INT is, const INT it);
void faces_attr(subscomplex *subsc);
void faces_cnt(subscomplex *subsc);
void vol_simplex(INT dim, REAL fact, REAL *xf, REAL *volt, void *wrk);
void area_face(INT dim, REAL fact, REAL *xf, REAL *sn,	\
	       REAL *areas,REAL *volt,			\
	       void *wrk);
void longest(scomplex *sc);
//void bfstree(INT it, scomplex *sc,INT *wrk);
void refining(INT ref_levels, scomplex *sc, INT nstar, REAL *xstar);
unsigned int reflect2(INT n, INT is, INT it,				\
		      INT* sv1, INT *sv2, INT* stos1, INT* stos2,	\
		      INT visited, INT *wrk);

void zinterp(unigrid *ug,REAL *x, const INT nvert, INT *mask, const INT maskvalue);

void vert_layer(INT nlayers,scomplex *sc,INT *mask, INT *nzptl);

void swne(INT n, INT nv, REAL *xo, REAL *xn, REAL *x);
void mapcoords(REAL *xp, INT np, INT nbig,			\
	       INT n, REAL *xo0,REAL *xn0,REAL *xo1,REAL *xn1);
void scalec(const INT flag, INT n, INT nv, REAL *xo, REAL *x, REAL *scale);
void scaleg(INT n, INT nv, REAL *xo, REAL *xn, REAL *x, REAL *scale);

INT features_r(INT dimorig, INT use_features, features *feat, REAL vfill);
INT features_w(features *feat, REAL *extra);

/* marking */
void markall(scomplex *sc,const int amark);
void markstar(scomplex *sc, dvector *w);
unsigned int markeql(scomplex *sc, dvector *w);

/* refinement */
void n_refine(INT ref_type, INT ref_levels, scomplex *sc,	\
	      dvector *errors,					\
	      void (*solving)(INT, scomplex *, void *),		\
	      void (*estimating)(INT , scomplex *, void *),	\
	      void (*marking)(INT , scomplex *, void *));
unsigned int reflect2(INT n, INT is, INT it,				\
		      INT* sv1, INT *sv2, INT* stos1, INT* stos2,	\
		      INT visited, INT *wrk);
void refining(INT ref_levels, scomplex *sc, INT nstar, REAL *xstar);
void scfinalize(scomplex *sc);
/************************************************************************/
macrocomplex *set_mmesh(input_grid *g0,					\
			cube2simp *c2s,					\
			INT *wrk);
/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
void scomplex_merge1(const INT nvall,		\
		     const INT nsall,		\
		     macrocomplex *mc,		\
		     scomplex **sc0,		\
		     cube2simp *c2s);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
INT locate1(INT *b,				\
	    INT *a,INT n,			\
	    INT *a2,INT n2,INT m2);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void macrocomplex_free(macrocomplex *mc);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void fix_grid(macrocomplex *mc,		\
	      scomplex **scin,			\
	      cube2simp *c2s,			\
	      input_grid *g0);
/*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/
scomplex *generate_grid(input_grid *g0);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
void scomplex_merge(scomplex **sc0,			\
		    const INT nsall, const INT nvall,	\
		    const INT cc, const INT bndry_cc,	\
		    input_grid *g0,cube2simp *c2s);
