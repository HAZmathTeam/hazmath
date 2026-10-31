/*! \file examples/stokes_biot/stokes_br.c
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief Stokes equations discretized with the Bernardi-Raugel element
 *        (P1 velocity enriched with normal face bubbles) and P0 pressure:
 *
 *        -div(2 mu eps(u)) + grad(p) = f,   div(u) = 0   in (0,1)^d,  d = 2,3,
 *                                     u = 0              on the boundary.
 *
 *        For refinement levels 1, ..., refinement_levels (input.dat) the
 *        program builds a uniform mesh, assembles and solves the system and
 *        prints the errors against the manufactured solution in
 *        stokes_br_data.h together with the convergence rates.
 *
 * \note Usage: ./stokes_br.ex [mu=<viscosity>] [full] [solver=gmres|umfpack|compare]
 *                            [amg=rs|ua] [cycle=amli|w|v] [rhs=pr|std] [test=0|1|2]
 *       By default the bubbles are statically condensed (bubble-bubble block
 *       replaced by its diagonal); "full" solves the full BR/P0 system.
 *       gmres:   (default) GMRES with the block lower triangular AMG
 *                preconditioner (stokes_br_precond.h).  If GMRES with RS-AMG
 *                does not converge, GMRES with UA-AMG (AMLI-cycle) is tried,
 *                and if that fails as well, UMFPACK;
 *       umfpack: direct solve with the constraint int p = 0 appended;
 *       compare: both, the errors are those of the GMRES solution and the
 *                last column is the relative difference of the two solutions.
 *       amg:   AMG in the preconditioner, classical RS (V-cycle, default) or UA;
 *       cycle: cycle for UA: AMLI (default), W or V.
 *       rhs:   pr (default): source tested with Pi_h v (BDM1 interpolant),
 *              pressure robust (Li and Zikatanov, CMAME 393 (2022) 114815);
 *              std: source tested with v.
 *       test:  test problem in stokes_br_data.h: 0 (default) div u = 0,
 *              Dirichlet; 1 div u != 0, Dirichlet; 2 div u != 0, Dirichlet
 *              on x_1 = 0 and traction on the rest of the boundary.
 *
 */

/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
#include "stokes_br_data.h"
#include "stokes_br_system.h"
#include "stokes_br_precond.h"
/*********************************************************************************/

/****** OPTIONS ******************************************************************/
// Options that are not in input.dat; mu, test and rhs are the globals
// stokes_mu, stokes_test and stokes_br_pr
typedef struct {
  INT condense;  // 1 condensed | 0 full
  INT solver;    // 0 umfpack | 1 gmres | 2 compare
  INT amg;       // BR_AMG_RS | BR_AMG_UA
  INT cycle;     // UA only: V_CYCLE | W_CYCLE | AMLI_CYCLE
} br_options;

// One option in command-line form (mu=1e-6, full, ...); returns 0 if unknown
static INT br_option(const char* a, br_options* opt)
{
  if (!strncmp(a, "mu=", 3)) stokes_mu = atof(a + 3);
  else if (!strcmp(a, "full")) opt->condense = 0;
  else if (!strcmp(a, "condensed")) opt->condense = 1;
  else if (!strcmp(a, "solver=gmres")) opt->solver = 1;
  else if (!strcmp(a, "solver=umfpack")) opt->solver = 0;
  else if (!strcmp(a, "solver=compare")) opt->solver = 2;
  else if (!strcmp(a, "amg=rs")) opt->amg = BR_AMG_RS;
  else if (!strcmp(a, "amg=ua")) opt->amg = BR_AMG_UA;
  else if (!strcmp(a, "cycle=amli")) opt->cycle = AMLI_CYCLE;
  else if (!strcmp(a, "cycle=w")) opt->cycle = W_CYCLE;
  else if (!strcmp(a, "cycle=v")) opt->cycle = V_CYCLE;
  else if (!strcmp(a, "rhs=pr")) stokes_br_pr = 1;
  else if (!strcmp(a, "rhs=std")) stokes_br_pr = 0;
  else if (!strncmp(a, "test=", 5)) stokes_test = atoi(a + 5);
  else return 0;
  return 1;
}

// stokes_br.dat: lines "key = value", % comments; each line is mapped to
// the command-line form of the option
static void read_stokes_br_dat(const char* file, br_options* opt)
{
  const char* keys[5][2] = {{"viscosity", "mu="}, {"test_problem", "test="},
                            {"rhs", "rhs="}, {"system", ""}, {"solver", "solver="}};
  char line[512], key[128], eq[8], val[128], a[256];
  FILE* fp = fopen(file, "r");
  if (fp == NULL) {
    printf("### WARNING: %s not found; using the defaults\n", file);
    return;
  }
  while (fgets(line, sizeof(line), fp) != NULL) {
    if (sscanf(line, "%127s", key) != 1 || key[0] == '%') continue;
    if (sscanf(line, "%127s %7s %127s", key, eq, val) != 3 || strcmp(eq, "=")) {
      printf("### WARNING: %s: cannot read line: %s", file, line);
      continue;
    }
    INT i, found = 0;
    for (i = 0; i < 5 && !found; i++)
      if (!strcmp(key, keys[i][0])) {
        snprintf(a, sizeof(a), "%s%s", keys[i][1], val);
        found = br_option(a, opt);
      }
    if (!found) printf("### WARNING: %s: unknown key or value: %s = %s\n", file, key, val);
  }
  fclose(fp);
}

/****** MAIN DRIVER **************************************************************/
int main(int argc, char* argv[]) {

  INT k, lvl;

  // Set Parameters from Reading in Input File
  input_param inparam;
  param_input_init(&inparam);
  param_input("./input.dat", &inparam);

  INT dim = inparam.spatial_dim;
  INT nq1d = inparam.nquad;
  INT nlevels = inparam.refinement_levels;
  stokes_br_dim = dim;

  // Options: defaults, then input.dat (AMG_type, AMG_cycle_type), then
  // stokes_br.dat (BR-only options), then the command line
  br_options opt = {1, 1, BR_AMG_RS, AMLI_CYCLE};
  if (inparam.AMG_type == UA_AMG) opt.amg = BR_AMG_UA;
  else if (inparam.AMG_type != C_AMG)
    printf("### WARNING: AMG_type in input.dat must be C or UA; using C\n");
  if (inparam.AMG_cycle_type == V_CYCLE || inparam.AMG_cycle_type == W_CYCLE ||
      inparam.AMG_cycle_type == AMLI_CYCLE)
    opt.cycle = inparam.AMG_cycle_type;
  read_stokes_br_dat("./stokes_br.dat", &opt);
  for (k = 1; k < argc; k++)
    if (!br_option(argv[k], &opt)) printf("### WARNING: unknown option %s\n", argv[k]);
  INT condense = opt.condense, solver = opt.solver, amg = opt.amg, cycle = opt.cycle;
  stokes_br_condense = condense;
  const char* solver_name[3] = {"UMFPACK", "GMRES", "GMRES vs UMFPACK"};

  REAL err[4], err_prev[4] = {0.0, 0.0, 0.0, 0.0};
  const char* cycle_name[4] = {"", "V", "W", "AMLI"};
  printf("Bernardi-Raugel/P0 Stokes, %lldD, test %lld, mu = %g, %s system, %s rhs, %s, uniform refinement:\n",
         (long long)dim, (long long)stokes_test, stokes_mu, condense ? "condensed" : "full",
         stokes_br_pr ? "pressure-robust" : "standard", solver_name[solver]);
  if (solver >= 1)
    printf("preconditioner: block lower triangular, %s-AMG %s-cycle\n",
           amg == BR_AMG_RS ? "RS" : "UA", amg == BR_AMG_RS ? "V" : cycle_name[cycle]);
  printf("%3s %10s %11s %5s %11s %5s %11s %5s %11s %4s %8s %8s %8s %10s\n",
         "lvl", "ndof", "|u-uh|_L2", "rate", "|u-uh|_H1", "rate",
         "|p-ph|_L2", "rate", "max|divT|", "its", "t_asm", "t_gmres", "t_umf", "|xg-xu|");

  for (lvl = 1; lvl <= nlevels; lvl++) {

    // Uniform mesh on the unit cube; boundary flag 1 (Dirichlet)
    scomplex* sc = make_uniform_mesh(dim, lvl, inparam.refinement_type, inparam.boundary_codes);
    qcoordinates* cq = get_quadrature(sc, nq1d);

    // BR velocity = face bubbles + P1 in each component; P0 pressure
    fespace FE_b, FE_u[3], FE_p;
    create_fespace(&FE_b, sc, 61);
    for (k = 0; k < dim; k++) create_fespace(&FE_u[k], sc, 1);
    create_fespace(&FE_p, sc, 0);
    set_dirichlet_bdry(&FE_b, sc, 1, 16);
    for (k = 0; k < dim; k++) set_dirichlet_bdry(&FE_u[k], sc, 1, 16);
    set_dirichlet_bdry(&FE_p, sc, -1, -1);
    if (stokes_test == 2) {
      // Dirichlet only on x_1 = 0 (vertices and face midpoints); the mesh
      // codes of edge/corner vertices belong to one side only
      for (k = 0; k < dim; k++)
        for (INT i = 0; i < FE_u[k].ndof; i++)
          FE_u[k].dirichlet[i] = (FE_u[k].cdof->x[i] < 1e-12) ? 1 : 0;
      for (INT i = 0; i < FE_b.ndof; i++)
        FE_b.dirichlet[i] = (sc->fem->f_mid[i * dim] < 1e-12) ? 1 : 0;
    }

    INT nspaces = dim + 2;
    INT ndof = FE_b.ndof + FE_p.ndof;
    for (k = 0; k < dim; k++) ndof += FE_u[k].ndof;
    block_fespace FE;
    initialize_fesystem(&FE, nspaces, 2 * dim + 1, ndof, sc->fem->ns_leaf);
    FE.var_spaces[0] = &FE_b;
    for (k = 0; k < dim; k++) FE.var_spaces[1 + k] = &FE_u[k];
    FE.var_spaces[dim + 1] = &FE_p;
    set_dirichlet_bdry_block(&FE, sc);

    // Assemble and eliminate boundary conditions
    block_dCSRmat A;
    bdcsr_alloc(nspaces, nspaces, &A);
    dvector b; b.row = 0; b.val = NULL;
    clock_t t0 = clock();
    assemble_global_system(&A, &b, &FE, sc, cq, local_assembly_stokes_br, NULL, source_rhs, NULL, 0.0);
    REAL ta = (REAL)(clock() - t0) / CLOCKS_PER_SEC;
    eliminate_DirichletBC_blockFE_blockA(bc_zero, &FE, sc, &b, &A, 0.0);

    // Solve: the condensed (or full) system As x = bs
    dvector sol = dvec_create(ndof);
    block_dCSRmat Ac;
    dvector bc;
    REAL* Dinv = NULL;
    block_dCSRmat* As = &A;
    dvector* bs = &b;
    if (condense) {
      Dinv = (REAL*)calloc(FE_b.ndof, sizeof(REAL));
      condense_bubbles(&A, &b, &Ac, &bc, Dinv);
      As = &Ac;
      bs = &bc;
    }
    INT nsys = bs->row;
    dvector xg = dvec_create(nsys), xu = dvec_create(nsys);
    INT its = 0, flag = 0;
    REAL tg = 0.0, tu = 0.0, diff = 0.0;
    INT use_umf = (solver != 1);
    if (solver >= 1) {
      t0 = clock();
      its = solve_br_gmres(As, bs, sc, stokes_mu, amg, cycle, &xg);
      if (its < 0 && amg == BR_AMG_RS) {
        printf("    GMRES with RS-AMG failed (code %lld); trying GMRES with UA-AMG (AMLI-cycle)\n",
               (long long)its);
        its = solve_br_gmres(As, bs, sc, stokes_mu, BR_AMG_UA, AMLI_CYCLE, &xg);
      }
      tg = (REAL)(clock() - t0) / CLOCKS_PER_SEC;
      if (its < 0) {
        printf("    GMRES failed (code %lld); falling back to UMFPACK\n", (long long)its);
        use_umf = 1;
      }
    }
    if (use_umf) {
      t0 = clock();
      flag = solve_br_umfpack(As, bs, sc, &xu);  // extends As, bs in place
      tu = (REAL)(clock() - t0) / CLOCKS_PER_SEC;
      if (flag < 0) printf("### ERROR: UMFPACK failed with error code = %lld!\n", (long long)flag);
    }
    if (solver == 2) {
      REAL nu = 0.0;
      for (k = 0; k < nsys; k++) {
        diff += (xg.val[k] - xu.val[k]) * (xg.val[k] - xu.val[k]);
        nu += xu.val[k] * xu.val[k];
      }
      diff = sqrt(diff / nu);
    }
    dvector* xs = (solver >= 1 && its >= 0) ? &xg : &xu;
    if (condense) {
      recover_bubbles(&A, &b, Dinv, xs->val, sol.val);
      bdcsr_free(&Ac);
      dvec_free(&bc);
      free(Dinv);
    } else {
      for (k = 0; k < ndof; k++) sol.val[k] = xs->val[k];
    }
    dvec_free(&xg);
    dvec_free(&xu);

    // Errors and rates (h halves with each level)
    stokes_br_errors(err, sol.val, &FE, sc, cq);
    REAL r[3] = {0.0, 0.0, 0.0};
    if (lvl > 1)
      for (k = 0; k < 3; k++) r[k] = log(err_prev[k] / err[k]) / log(2.0);
    printf("%3lld %10lld %11.4e %5.2f %11.4e %5.2f %11.4e %5.2f %11.4e %4lld %8.2f %8.2f %8.2f %10.3e\n",
           (long long)lvl, (long long)ndof, err[0], r[0], err[1], r[1], err[2], r[2], err[3],
           (long long)its, ta, tg, tu, diff);
    fflush(stdout);
    for (k = 0; k < 4; k++) err_prev[k] = err[k];

    // Free
    bdcsr_free(&A);
    dvec_free(&b);
    dvec_free(&sol);
    free_fespace(&FE_b);
    for (k = 0; k < dim; k++) free_fespace(&FE_u[k]);
    free_fespace(&FE_p);
    free_blockfespace(&FE);
    free_qcoords(cq);
    free(cq);
    haz_scomplex_free(sc);
  }

  return 0;
}  /* End of Program */
/*********************************************************************************************/
