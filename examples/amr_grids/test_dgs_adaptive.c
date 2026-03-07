/*
 * test_dgs_adaptive.c
 *
 * Tests DGS (Diening-Gehring-Storn) adaptive refinement in 2D-5D.
 * Refines simplices containing at least one of 3 random points
 * inside the unit cube.
 *
 * Usage: ./test_dgs_adaptive.ex <dim> <nref> <seed>
 */
#include "hazmath.h"

/* Mark all leaf simplices that contain at least one of the given points */
static ivector mark_containing(scomplex *sctop,
                               INT npts, REAL *pts)
{
  INT dim = sctop->n, n1 = dim + 1;
  INT ns = sctop->ns, i, k;
  ivector marked;
  marked.row = ns;
  marked.val = (INT *)calloc(ns, sizeof(INT));
  for (i = 0; i < ns; i++) {
    INT *el = sctop->nodes + i * n1;
    for (k = 0; k < npts; k++) {
      if (!xins(dim, el, sctop->x, pts + k * dim)) {
        marked.val[i] = 1;
        break;
      }
    }
  }
  return marked;
}

int main(int argc, char *argv[])
{
  INT dim = 2, nref = 8, seed = 42;
  if (argc > 1) dim = atoi(argv[1]);
  if (argc > 2) nref = atoi(argv[2]);
  if (argc > 3) seed = atoi(argv[3]);

  /* Generate 3 random points in (0,1)^dim (mesh_cube_init uses (0,1)^d) */
  srand(seed);
  INT npts = 3;
  REAL *pts = (REAL *)calloc(npts * dim, sizeof(REAL));
  INT i, j;
  fprintf(stdout, "\n%% dim=%d, nref=%d, seed=%d\n", (int)dim, (int)nref, (int)seed);
  fprintf(stdout, "%% Points in (0,1)^%d:\n", (int)dim);
  for (i = 0; i < npts; i++) {
    fprintf(stdout, "%%   P%d = (", (int)i);
    for (j = 0; j < dim; j++) {
      pts[i * dim + j] = 0.05 + 0.9 * ((REAL)rand() / (REAL)RAND_MAX);
      fprintf(stdout, "%8.4f%s", pts[i * dim + j], (j < dim - 1) ? ", " : "");
    }
    fprintf(stdout, ")\n");
  }

  /* Build initial mesh on (0,1)^dim using mesh_cube_init with ref_type=20 */
  INT ndiv = 1;
  if (dim == 2) ndiv = 2; /* finer initial mesh for 2D */
  scomplex **sc_all = mesh_cube_init(dim, ndiv, 20);
  scomplex *sc = sc_all[0];
  sc->ref_type = 20; /* DGS initialization */

  fprintf(stdout, "%% Initial mesh: ns=%d, nv=%d\n", (int)sc->ns, (int)sc->nv);

  /* Adaptive refinement loop */
  scomplex *sctop = NULL;
  ivector marked;
  for (i = 0; i < nref; i++) {
    sctop = scfinest(sc);
    marked = mark_containing(sctop, npts, pts);
    INT kmarked = 0;
    for (j = 0; j < marked.row; j++)
      if (marked.val[j]) kmarked++;
    fprintf(stdout, "%% level %2d: total_ns=%8d, leaves=%8d, marked=%6d\n",
            (int)i, (int)sc->ns, (int)sctop->ns, (int)kmarked);
    if (!kmarked) {
      fprintf(stdout, "%% No simplices marked, stopping.\n");
      ivec_free(&marked);
      haz_scomplex_free(sctop);
      break;
    }
    refine(1, sc, &marked);
    ivec_free(&marked);
    haz_scomplex_free(sctop);
  }

  /* Final stats and conformity check */
  sctop = scfinest(sc);
  fprintf(stdout, "%% Final mesh: total_ns=%d, leaves=%d, nv=%d\n",
          (int)sc->ns, (int)sctop->ns, (int)sc->nv);
  INT nerr = sc_conformity_check(sctop);
  if (nerr) {
    fprintf(stderr, "%% FAIL: mesh is non-conforming (%d bad facets)\n", (int)nerr);
  }
  haz_scomplex_free(sctop);

  /* Write VTU if dim <= 3 */
  if (dim <= 3) {
    /* (debug removed) */
    scfinalize(sc, (INT)1);
    {
      INT nerr2 = sc_conformity_check(sc);
      if (nerr2)
        fprintf(stderr, "%% FAIL: finalized mesh is non-conforming (%d bad facets)\n", (int)nerr2);
      else
        fprintf(stdout, "%% Finalized mesh conformity check PASSED (ns=%d, nv=%d)\n", (int)sc->ns, (int)sc->nv);
    }
    char fname[256];
    snprintf(fname, sizeof(fname), "output/dgs_adaptive_%dd.vtu", (int)dim);
    vtu_data vdata;
    vtu_data_init(sc, &vdata);
    vtkw(fname, &vdata);
    vtu_data_free(&vdata);
    fprintf(stdout, "%% VTU written to %s\n", fname);
  }

  /* Write points to a CSV file for the visualization script */
  {
    char fname[256];
    snprintf(fname, sizeof(fname), "output/dgs_points_%dd.csv", (int)dim);
    FILE *fp = fopen(fname, "w");
    for (i = 0; i < npts; i++) {
      for (j = 0; j < dim; j++)
        fprintf(fp, "%.10e%s", pts[i * dim + j], (j < dim - 1) ? "," : "\n");
    }
    fclose(fp);
  }

  haz_scomplex_free(sc);
  free(sc_all);
  free(pts);
  return 0;
}
