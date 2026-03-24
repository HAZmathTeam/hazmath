/*! \file test_bndry_codes.c
 *
 *  Authors: HAZmath (https://hazmath.net)
 *           Created with the help of Claude (Anthropic)
 *
 *  Tests boundary code propagation after refinement.
 *  Writes boundary faces as VTU with face codes as cell data.
 *
 *  Usage: ./test_bndry_codes.ex [input_file] [ref_type]
 *         ref_type: 0=DGS, 20=Bey, 22=marked-Bey-corner
 */
#include "hazmath.h"

/* Write a boundary scomplex (from sc_bndry) as VTU with flags as cell data */
static void write_bndry_vtu(const char *fname, scomplex *bsc)
{
  INT ns = bsc->ns, nv = bsc->nv, dim = bsc->dim, nbig = bsc->nbig;
  INT n1 = dim + 1;

  FILE *fp = fopen(fname, "w");
  if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); return; }

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
    (long long)nv, (long long)ns);

  /* Points */
  fprintf(fp, "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (INT i = 0; i < nv; i++) {
    for (INT d = 0; d < 3; d++) {
      if (d < nbig) fprintf(fp, "%.15e ", bsc->x[i * nbig + d]);
      else fprintf(fp, "0.0 ");
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n</Points>\n");

  /* Cells */
  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (INT f = 0; f < ns; f++) {
    for (INT j = 0; j < n1; j++)
      fprintf(fp, "%lld ", (long long)bsc->nodes[f * n1 + j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  for (INT f = 0; f < ns; f++) fprintf(fp, "%lld ", (long long)((f + 1) * n1));
  fprintf(fp, "\n</DataArray>\n");
  /* VTK types: 3=line(1D), 5=triangle(2D), 10=tetra(3D) */
  INT vtktype = (dim == 1) ? 3 : (dim == 2) ? 5 : 10;
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (INT f = 0; f < ns; f++) fprintf(fp, "%lld ", (long long)vtktype);
  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  /* Cell data: face code from flags */
  fprintf(fp, "<CellData Scalars=\"FaceCode\">\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"FaceCode\" format=\"ascii\">\n");
  for (INT f = 0; f < ns; f++)
    fprintf(fp, "%lld ", (long long)bsc->flags[f]);
  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "</CellData>\n");

  fprintf(fp, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
  fclose(fp);

  /* Print summary */
  fprintf(stdout, "  sc_bndry: %lld faces, %lld vertices, codes:", (long long)ns, (long long)nv);
  {
    INT codes[256]; INT counts[256]; INT nc = 0;
    for (INT f = 0; f < ns; f++) {
      INT c = bsc->flags[f];
      INT found = 0;
      for (INT k = 0; k < nc; k++) if (codes[k] == c) { found = 1; counts[k]++; break; }
      if (!found && nc < 256) { codes[nc] = c; counts[nc] = 1; nc++; }
    }
    for (INT c = 0; c < nc; c++)
      for (INT d = c+1; d < nc; d++)
        if (codes[c] > codes[d]) {
          INT t; t=codes[c]; codes[c]=codes[d]; codes[d]=t;
          t=counts[c]; counts[c]=counts[d]; counts[d]=t;
        }
    for (INT c = 0; c < nc; c++)
      fprintf(stdout, " %lld(%lld)", (long long)codes[c], (long long)counts[c]);
  }
  fprintf(stdout, "\n  Written to %s\n", fname);
}

/* Write all boundary faces (f_flag != 0 or boundary) as VTU from sc->fem */
static void write_bndry_vtu_from_fem(const char *fname, scomplex *sc)
{
  sc_fem *fem = sc->fem;
  iCSRmat *fv = fem->f_v;
  INT nface = fem->nface, dim = sc->dim, nbig = sc->nbig;
  INT n1f = dim; /* vertices per face */

  /* Collect boundary faces (f_flag != 0 or face has only 1 adjacent element) */
  iCSRmat f_el;
  icsr_trans(fem->el_f, &f_el);
  INT nbf = 0;
  for (INT f = 0; f < nface; f++) {
    INT deg = f_el.IA[f+1] - f_el.IA[f];
    if (fem->f_flag[f] || deg < 2) nbf++;
  }
  INT *bf_list = (INT *)malloc(nbf * sizeof(INT));
  INT ci = 0;
  for (INT f = 0; f < nface; f++) {
    INT deg = f_el.IA[f+1] - f_el.IA[f];
    if (fem->f_flag[f] || deg < 2) bf_list[ci++] = f;
  }
  icsr_free(&f_el);

  /* Renumber vertices */
  INT *gmap = (INT *)calloc(sc->nv, sizeof(INT));
  for (INT i = 0; i < sc->nv; i++) gmap[i] = -1;
  for (INT ci = 0; ci < nbf; ci++) {
    INT f = bf_list[ci];
    for (INT j = fv->IA[f]; j < fv->IA[f+1]; j++) gmap[fv->JA[j]] = 0;
  }
  INT nv_b = 0;
  for (INT i = 0; i < sc->nv; i++) if (gmap[i] >= 0) gmap[i] = nv_b++;

  FILE *fp = fopen(fname, "w");
  fprintf(fp, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", (long long)nv_b, (long long)nbf);
  fprintf(fp, "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (INT i = 0; i < sc->nv; i++) {
    if (gmap[i] < 0) continue;
    for (INT d = 0; d < 3; d++)
      fprintf(fp, "%.15e ", (d < nbig) ? sc->x[i*nbig+d] : 0.0);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n</Points>\n<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (INT ci = 0; ci < nbf; ci++) {
    INT f = bf_list[ci];
    for (INT j = fv->IA[f]; j < fv->IA[f+1]; j++) fprintf(fp, "%lld ", (long long)gmap[fv->JA[j]]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  for (INT ci = 0; ci < nbf; ci++) fprintf(fp, "%lld ", (long long)((ci+1)*n1f));
  fprintf(fp, "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  INT vtktype = (dim == 2) ? 3 : (dim == 3) ? 5 : 10;
  for (INT ci = 0; ci < nbf; ci++) fprintf(fp, "%lld ", (long long)vtktype);
  fprintf(fp, "\n</DataArray>\n</Cells>\n");
  fprintf(fp, "<CellData Scalars=\"FaceCode\">\n<DataArray type=\"Int64\" Name=\"FaceCode\" format=\"ascii\">\n");
  for (INT ci = 0; ci < nbf; ci++) fprintf(fp, "%lld ", (long long)fem->f_flag[bf_list[ci]]);
  fprintf(fp, "\n</DataArray>\n</CellData>\n");
  fprintf(fp, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
  fclose(fp);
  free(gmap);
  free(bf_list);
  fprintf(stdout, "  Written to %s\n", fname);
}

INT main(INT argc, char *argv[])
{
  const char *input_file = "input/3d_cube.input";
  INT ref_type = 20;
  if (argc > 1) input_file = argv[1];
  if (argc > 2) ref_type = atoi(argv[2]);

  config_z config = get_input(input_file);
  input_grid *g = input_grid_alloc();
  config2vars_amr(&config, g);
  free_config(&config);
  input_grid_set_output(g, input_file);

  scomplex **sc_all = generate_initial_grid(g);
  scomplex *sc = sc_all[0];
  INT dim = sc->dim, n1 = dim + 1;
  INT nref = 2;

  fprintf(stdout, "Input: %s  dim=%lld  ref_type=%lld  nref=%lld\n",
    input_file, (long long)dim, (long long)ref_type, (long long)nref);
  fprintf(stdout, "Initial: ns=%lld, nv=%lld\n", (long long)sc->ns, (long long)sc->nv);

  /* Refine (ref_type < 0 means no refinement) */
  if (ref_type < 0) {
    /* no refinement — just finalize */
    nref = 0;
  } else if (ref_type == 22) {
    /* Marked Bey near origin */
    REAL threshold = 1.0;
    for (INT lvl = 0; lvl < nref; lvl++) {
      ivector mark = ivec_create(sc->ns);
      INT nmarked = 0;
      for (INT k = 0; k < sc->ns; k++) {
        REAL dist2 = 0.0;
        for (INT dd = 0; dd < dim; dd++) {
          REAL c = 0.0;
          for (INT vv = 0; vv < n1; vv++)
            c += sc->x[sc->nbig * sc->nodes[n1 * k + vv] + dd];
          c /= (REAL)n1;
          dist2 += c * c;
        }
        mark.val[k] = (dist2 < threshold * threshold) ? 1 : 0;
        if (mark.val[k]) nmarked++;
      }
      fprintf(stdout, "  Level %lld: marked %lld / %lld\n",
        (long long)lvl, (long long)nmarked, (long long)sc->ns);
      uniformrefine_marked(sc, &mark);
      ivec_free(&mark);
      sc_vols(sc);
      threshold *= 0.6;
    }
    find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
  } else if (ref_type > 10) {
    /* Uniform Bey */
    for (INT lvl = 0; lvl < nref; lvl++) {
      uniformrefine(sc);
      sc_vols(sc);
    }
    find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
  } else {
    /* DGS bisection */
    refine(nref, sc, NULL);
  }

  /* Finalize and rebuild boundary data */
  scfinalize(sc, NULL, (INT)1);
  sc_vols(sc);
  fprintf(stdout, "Refined: ns=%lld, nv=%lld\n", (long long)sc->ns, (long long)sc->nv);

  /* Build FEM data (f_v, f_flag, coded faces) */
  if (dim <= 3) sc_build_fem_data(sc);

  /* Conformity check */
  INT nerr = sc_conformity_check(sc);
  if (nerr)
    fprintf(stderr, "FAIL: non-conforming (%lld errors)\n", (long long)nerr);
  else
    fprintf(stdout, "Conformity check PASSED\n");

  /* Use f_flag from sc->fem for face codes — no sc_bndry needed */
  if (dim <= 3 && sc->fem) {
    sc_fem *fem = sc->fem;
    /* Print coded-faces summary */
    fprintf(stdout, "  Coded faces: %lld / %lld total\n",
      (long long)fem->n_coded_faces, (long long)fem->nface);
    {
      INT codes[256]; INT counts[256]; INT nc = 0;
      for (INT ci = 0; ci < fem->n_coded_faces; ci++) {
        INT c = fem->f_flag[fem->coded_faces[ci]];
        INT found = 0;
        for (INT k = 0; k < nc; k++) if (codes[k] == c) { found = 1; counts[k]++; break; }
        if (!found && nc < 256) { codes[nc] = c; counts[nc] = 1; nc++; }
      }
      for (INT c = 0; c < nc; c++)
        for (INT d = c+1; d < nc; d++)
          if (codes[c] > codes[d]) {
            INT t; t=codes[c]; codes[c]=codes[d]; codes[d]=t;
            t=counts[c]; counts[c]=counts[d]; counts[d]=t;
          }
      fprintf(stdout, "  Face codes:");
      for (INT c = 0; c < nc; c++)
        fprintf(stdout, " %lld(%lld)", (long long)codes[c], (long long)counts[c]);
      fprintf(stdout, "\n");
      /* Count boundary vs interior */
      INT nbf = 0, nif = 0;
      for (INT ci = 0; ci < fem->n_coded_faces; ci++) {
        if (fem->coded_f_btype[ci]) nif++; else nbf++;
      }
      fprintf(stdout, "  Boundary: %lld, Interior: %lld\n", (long long)nbf, (long long)nif);
    }
    /* Write boundary faces VTU using f_v + f_flag */
    write_bndry_vtu_from_fem("output/bndry_refined.vtu", sc);
    /* Also write the volume mesh */
    {
      vtu_data vdata;
      vtu_data_init(sc, &vdata);
      sc_write_vtk("output/vol_refined.vtu", &vdata);
      vtu_data_free(&vdata);
    }
  }

  input_grid_free(g);
  haz_scomplex_free(sc);
  free(sc_all);
  fprintf(stdout, "Done.\n");
  return 0;
}
