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
void bfstree(INT it, scomplex *sc,INT *wrk);
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

/*******************************************************************/
/*INPUT GRID*/
/***********************************************************/
void input_grid_free(input_grid *g);
void input_grid_print(input_grid *g);
char **splits(char *s, const char *d, INT *num);
void read_data(char *data_coordsystems,		\
	       char *data_vertices,		\
	       char *data_edges,		\
	       input_grid *g);
void get_out(const char *pattern, size_t le);
char *make_string_from_file(FILE *the_file, size_t *length_string);
/********************************************************************/
char *get_substring(const char *pattern,	\
		    size_t *length_substring,	\
		    char *the_string);
input_grid *parse_input_grid(const char *input_file_grid);
/* NOT USED void set_filenames (INT nlayers, const char *geofile) */
/* { */
/*   INT i; */
/*   char *layer[nlayers]; */
/*   FILE *fp = fopen(geofile,"r"); */
/*   layer[0]=strdup("Basement_Rock"); */
/*   layer[1]=strdup("MiddleLower_Siwalik"); */
/*   layer[2]=strdup("Upper_Siwalik"); */
/*   layer[3]=strdup("Surface_Elevation"); */
/*   for(i = 0;i<nlayers;i++){ */
/*     fprintf(stdout," Layer: %s\n",layer[i]); */
/*     if(layer[i]) free(layer[i]); */
/*   } */
/*   return; */
/* } */



