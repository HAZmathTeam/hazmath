/*! \file examples/amr_grids/test_msh_roundtrip.c
 *
 *  \brief Demonstrates Gmsh .msh round-trip: generate a grid from an
 *         input file, write it to .msh, read it back, refine, and
 *         write the refined mesh to .msh and .vtu.
 *
 *  Usage: ./test_msh_roundtrip.ex [input_file]
 *         Default: input/3d_2cubes_edge_bcodes.input
 *
 *  Created 2026-03-19 (ltz with the help of Claude)
 */

#include "hazmath.h"
#include "supporting_amr_grids.h"

INT main(INT argc, char* argv[])
{
  const char* input_file = "input/3d_2cubes_edge.input";
  if (argc > 1) input_file = argv[1];

  printf("===========================================================\n");
  printf("  Gmsh .msh Round-Trip Test\n");
  printf("===========================================================\n\n");

  /* ---- Phase 1: Generate initial grid from input file ---- */
  printf("Phase 1: Generating initial grid from %s\n", input_file);
  config_z config = get_input(input_file);
  input_grid* g = input_grid_alloc();
  config2vars_amr(&config, g);
  free_config(&config);
  input_grid_set_output(g, input_file);

  scomplex** sc_all = generate_initial_grid(g);
  scomplex* sc = sc_all[0];
  scfinalize(sc, NULL, (INT)1);
  sc_vols(sc);

  printf("  Initial mesh: ns=%lld, nv=%lld, dim=%lld\n",
         (long long)sc->ns, (long long)sc->nv, (long long)sc->dim);

  /* Count boundary vertices */
  INT i, nbv = 0;
  for (i = 0; i < sc->nv; i++)
    if (sc->bndry[i] != 0) nbv++;
  printf("  Boundary vertices: %lld\n", (long long)nbv);
  printf("  Connected components: bulk=%lld, boundary=%lld\n",
         (long long)sc->cc, (long long)sc->bndry_cc);

  /* ---- Phase 2: Write to .msh ---- */
  const char* msh_file = "output/2cubes_original.msh";
  printf("\nPhase 2: Writing to %s\n", msh_file);
  sc_write_gmsh(msh_file, sc, 1);

  /* Also write VTU for visualization */
  {
    vtu_data vdata;
    vtu_data_init(sc, &vdata);
    sc_write_vtk("output/2cubes_original.vtu", &vdata);
    vtu_data_free(&vdata);
  }

  /* Free original mesh */
  haz_scomplex_free(sc);
  free(sc_all);
  input_grid_free(g);

  /* ---- Phase 3: Read back from .msh ---- */
  printf("\nPhase 3: Reading back from %s\n", msh_file);
  sc = sc_read_gmsh(msh_file);

  printf("  Read mesh: ns=%lld, nv=%lld, dim=%lld\n",
         (long long)sc->ns, (long long)sc->nv, (long long)sc->dim);

  /* Build connectivity for refinement */
  find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
  sc_vols(sc);

  /* Print boundary code summary */
  INT max_code = 0;
  for (i = 0; i < sc->nv; i++)
    if (sc->bndry[i] > max_code) max_code = sc->bndry[i];
  printf("  Boundary codes range: 0 to %lld\n", (long long)max_code);

  if (sc->bndry_f2v) {
    printf("  Boundary faces (from .msh): %lld\n",
           (long long)sc->bndry_f2v->row);
  }

  /* Conformity check */
  INT nerr = sc_conformity_check(sc);
  if (nerr)
    fprintf(stderr, "  FAIL: non-conforming (%lld errors)\n", (long long)nerr);

  /* ---- Phase 4: Refine twice ---- */
  printf("\nPhase 4: Refine (2 levels, newest vertex bisection)\n");
  INT dim = sc->dim;

  /* First refinement: uniform */
  refine(1, sc, NULL);
  /* Extract leaf mesh — this destroys the hierarchy, giving a flat mesh */
  scfinalize(sc, NULL, (INT)1);
  /* Reset level and generations so refine treats this as a fresh mesh */
  sc->level = 0;
  for (i = 0; i < sc->ns; i++) sc->gen[i] = 0;

  /* Second refinement: uniform on the leaf mesh */
  refine(1, sc, NULL);

  /* Finalize: compact to leaves only */
  scfinalize(sc, NULL, (INT)1);
  sc_vols(sc);

  printf("  Refined mesh: ns=%lld, nv=%lld\n",
         (long long)sc->ns, (long long)sc->nv);

  nbv = 0;
  for (i = 0; i < sc->nv; i++)
    if (sc->bndry[i] != 0) nbv++;
  printf("  Boundary vertices after refinement: %lld\n", (long long)nbv);

  /* Check no boundary vertex has zero code */
  INT zero_bv = 0;
  if (sc->bndry_v) {
    for (i = 0; i < sc->bndry_v->row && i < sc->nv; i++) {
      if (sc->bndry_v->IA[i+1] > sc->bndry_v->IA[i] && sc->bndry[i] == 0)
        zero_bv++;
    }
  }
  if (zero_bv)
    fprintf(stderr, "  WARNING: %lld boundary vertices have zero code\n", (long long)zero_bv);
  else
    printf("  All boundary vertices have nonzero codes.\n");

  /* Conformity check on refined mesh */
  nerr = sc_conformity_check(sc);
  if (nerr)
    fprintf(stderr, "  FAIL: refined mesh non-conforming (%lld errors)\n", (long long)nerr);

  /* ---- Phase 5: Write refined mesh ---- */
  printf("\nPhase 5: Writing refined mesh\n");
  sc_write_gmsh("output/2cubes_refined.msh", sc, 1);

  {
    vtu_data vdata;
    vtu_data_init(sc, &vdata);
    sc_write_vtk("output/2cubes_refined.vtu", &vdata);
    vtu_data_free(&vdata);
  }

  /* Also write MATLAB file */
  sc_matlab_write(sc, "output/2cubes_refined.m");

  printf("\nOutput files:\n");
  printf("  output/2cubes_original.msh  (initial mesh)\n");
  printf("  output/2cubes_original.vtu  (initial mesh, VTK)\n");
  printf("  output/2cubes_refined.msh   (refined mesh)\n");
  printf("  output/2cubes_refined.vtu   (refined mesh, VTK)\n");
  printf("  output/2cubes_refined.m     (refined mesh, MATLAB)\n");

  /* ---- Cleanup ---- */
  haz_scomplex_free(sc);

  printf("\nDone.\n");
  return 0;
}
