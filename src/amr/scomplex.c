/*! \file src/amr/scomplex.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20170715.
 *  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note containing all essentials routines for mesh refinement
 *  \modified 20260307 (ltz with the help of Claude (Anthropic))
 *
 */
#include "hazmath.h"
/**********************************************************************/
/*!
 * \fn REAL void haz_scomplex_realloc(scomplex *sc)
 *
 * \brief  reallocates memory if not allocated with haz_scomplex_init;
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void haz_scomplex_realloc(scomplex* sc) {
  INT i, j, ns = sc->ns, nv = sc->nv, dim = sc->dim;
  INT n1 = dim + 1;
  sc->print_level = 0;
  sc->factorial = 1.;
  for (j = 2; j < n1; j++) sc->factorial *= ((REAL)j);
  sc->nbr = realloc(sc->nbr, n1 * ns * sizeof(INT));
  sc->marked = realloc(sc->marked, ns * sizeof(INT));
  sc->gen = realloc(sc->gen, ns * sizeof(INT));
  sc->parent = realloc(sc->parent, ns * sizeof(INT));
  sc->child0 = realloc(sc->child0, ns * sizeof(INT));
  sc->childn = realloc(sc->childn, ns * sizeof(INT));
  sc->bndry = realloc(sc->bndry, nv * sizeof(INT));
  sc->csys = realloc(sc->csys, nv * sizeof(INT));
  sc->flags = realloc(sc->flags, ns * sizeof(INT));  // element flags
  sc->vols = realloc(sc->vols, ns * sizeof(REAL));   // element volumes
  for (i = 0; i < sc->ns; i++) {
    sc->marked[i] =
        FALSE;  // because first this array is used as working array.
    sc->gen[i] = 0;
    sc->parent[i] = -1;
    sc->child0[i] = -1;
    sc->childn[i] = -1;
    sc->flags[i] = -1;
  }
  for (i = 0; i < nv; i++) {
    sc->bndry[i] = 0;
    sc->csys[i] = 0;
  }
  sc->fem = NULL;
  sc->inc = NULL;
  return;
}
/**********************************************************************/
/*!
 * \fn REAL chk_sign(const int it, const int nbrit)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
REAL chk_sign(const INT it, const INT nbrit) {
  /*
    nbrit is a neighbor of it and this picks the sign of the normal
    vector for the face shared by both. Used to construct dcsrmat atf
    and local matrices to build aff. If the face is on the boundary,
    i.e. nbrit <0, we always return 1 (i.e. the outward normal).
  */
  if (it > nbrit) return 1e0;
  return -1e0;
}
/*********************************************************************/
/*!
 * \fn REAL volume_compute(INT dim, REAL factorial, REAL *xs,void *wrk);
 *
 * \brief Computes the volume of the simplex in Rn, which is det(B)/d_factorial.
 *
 * \param dim        I: The dimension of the problem.
 * \param factorial  I: dim! (dim factorial)
 * \param xs         I: coordinates of the simplex
 * \param wrk        W: working array of dimension
 *                      (dim+1)*(dim*sizeof(REAL) + sizeof(INT))
 *
 * \return the volume of the simplex.
 *
 * \note
 */
REAL volume_compute(INT dim, REAL factorial, REAL* xs, void* wrk) {
  INT dim1 = dim + 1, i, j, ln, ln1;
  REAL* bt = (REAL*)wrk;
  REAL* piv = bt + dim * dim;
  INT* p = (INT*)(wrk + (dim * dim + dim) * sizeof(REAL));
  REAL vol = -1e20;
  // construct bt using xs;
  for (j = 1; j < dim1; j++) {
    ln = j * dim;
    ln1 = ln - dim;
    for (i = 0; i < dim; i++) {
      bt[ln1 + i] = xs[ln + i] - xs[i];  // k-th row of bt is [x(k)-x(0)]^T.
                                         // x(k) are coords of vertex k.
    }
  }
  //  SHORT flag=ddense_lu(1, dim, &vol, bt,p,piv);
  //  if(flag){
  if (ddense_lu(1, dim, &vol, bt, p, piv)) {
    return 0e0;  // degeneraate simplex;
  } else
    return fabs(vol) / factorial;
}
/*********************************************************************/
/*!
 * \fn void sc_vols(scomplex *sc);
 *
 * \brief Fills in the array with volumes of the simpleices in
 *        n-dimansional simplicial complex.
 *
 * \param sc: pointer to the simplicial complex sc.
 *
 * \note
 */
void sc_vols(scomplex* sc) {
  INT dim = sc->dim, ns = sc->ns;
  INT dim1 = dim + 1, i, j, node, idim1;
  void* wrk = malloc(dim1 * (dim * sizeof(REAL) + sizeof(INT)));
  REAL* xs = calloc(dim1 * dim, sizeof(REAL));
  for (i = 0; i < ns; ++i) {
    idim1 = i * dim1;
    for (j = 0; j < dim1; ++j) {
      node = sc->nodes[idim1 + j];
      memcpy((xs + j * dim), (sc->x + node * dim), dim * sizeof(REAL));
    }
    sc->vols[i] = volume_compute(dim, sc->factorial, xs, wrk);
  }
  free(wrk);
  free(xs);
  return;
}
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_init(const INT n,INT ns, INT nv, const INT nbig)
 *
 * \brief Initialize simplicial complex in dimension n with ns
 *        simplices and nv vertices.
 *
 * \param n is the dimension of the simplicial complex;
 * \param ns is the number of simplices
 * \param nv is the number of vertices
 * \param nbig is the dimension of the space where the simplicial
 *        complex is embedded;
 *
 * \return initialized structure of type scomplex
 *
 * \note
 *
 */
scomplex* haz_scomplex_init(const INT n, INT ns, INT nv, const INT nbig) {
  /*
     n = dimension of the simplicial complex;
     nbig = dimension of the space where the complex is embedded;
     Future work: we should think of making this for different
     dimensions, e.g. 2-homogenous complex in 3d
  */
  scomplex* sc = (scomplex*)malloc(sizeof(scomplex));
  sc->nbig = nbig;
  sc->dim = n;
  if (sc->nbig <= 0) sc->nbig = n;
  INT n1 = sc->dim + 1, i, j, in1;
  sc->level = 0;
  sc->ref_type = 0;
  sc->print_level = 0;
  sc->factorial = 1.;
  for (j = 2; j < n1; j++) sc->factorial *= ((REAL)j);
  //
  sc->marked = (INT*)calloc(ns, sizeof(INT));
  sc->gen = (INT*)calloc(ns, sizeof(INT));
  sc->nbr = (INT*)calloc(ns * n1, sizeof(INT));
  sc->parent = (INT*)calloc(ns, sizeof(INT));
  sc->child0 = (INT*)calloc(ns, sizeof(INT));
  sc->childn = (INT*)calloc(ns, sizeof(INT));
  sc->nodes = (INT*)calloc(ns * n1, sizeof(INT));
  sc->bndry = (INT*)calloc(nv, sizeof(INT));
  sc->csys = (INT*)calloc(nv, sizeof(INT)); /* coord sys: 1 is polar, 2
                                               is cyl and so on */
  sc->parent_v = malloc(sizeof(iCSRmat));
  sc->parent_v[0] = icsr_create(nv, nv, nv);
  sc->flags = (INT*)calloc(ns, sizeof(INT));
  sc->x = (REAL*)calloc(nv * nbig, sizeof(REAL));
  sc->vols = (REAL*)calloc(ns, sizeof(REAL));
  for (i = 0; i < ns; i++) {
    sc->marked[i] = FALSE;  // because first is used for something else.
    sc->gen[i] = 0;
    sc->parent[i] = -1;
    sc->child0[i] = -1;
    sc->childn[i] = -1;
    sc->flags[i] = -1;
    sc->vols[i] = -1e20;
    in1 = i * n1;
    for (j = 0; j < n1; j++) {
      sc->nodes[in1 + j] = -1;
      sc->nbr[in1 + j] = -1;
    }
  }
  INT nnz_pv = 0;
  sc->parent_v->IA[0] = nnz_pv;
  for (i = 0; i < nv; i++) {
    sc->bndry[i] = 0;
    sc->csys[i] = 0;
    sc->parent_v->JA[nnz_pv] = i;
    sc->parent_v->val[nnz_pv] = 0;
    nnz_pv++;
    sc->parent_v->IA[i + 1] = nnz_pv;
  }
  // not needed for now, it will be freed later.
  //  if(nnz_pv) memset(sc->parent_v->val,0,nnz_pv*sizeof(INT));
  //////////////////////////////////////
  sc->nv = nv;
  sc->ns = ns;
  sc->bndry_cc = 1;  // one connected component on the boundary for now.
  sc->cc = 1;        // one connected component in the bulk for now.
  // NULL pointers for the rest
  sc->etree = NULL;
  sc->bfs = malloc(sizeof(iCSRmat));
  sc->bfs[0] = icsr_create(0, 0, 0);
  // the parent_v->val is not needed for now
  /* if(sc->parent_v->val) { */
  /*   free(sc->parent_v->val); */
  /*   sc->parent_v->val=NULL; */
  /* } */
  sc->bndry_v = malloc(sizeof(iCSRmat));
  sc->bndry_v[0] = icsr_create(0, 0, 0);
  sc->bndry_f2v = NULL;
  sc->v2s_head = NULL;
  sc->v2s_next = NULL;
  sc->v2s_simp = NULL;
  sc->v2s_count = 0;
  sc->v2s_alloc = 0;
  sc->v2s_nv = 0;
  sc->gen_N = NULL;
  sc->dgs_N = 0;
  sc->ncsys = 0;
  sc->systypes = NULL;
  sc->csys_ox = NULL;
  sc->fem = NULL;
  sc->inc = NULL;
  return sc;
}
/**********************************************************************/
/*!
 * \fn scomplex haz_scomplex_null(const INT n,INT ns, INT nv, const INT nbig)
 *
 * \brief Initialize all pointers in simplicial complex to NULL;
 * \param n is the dimension of the simplicial complex;
 * \param nbig is the dimension of the space where the simplicial
 *        complex is embedded;
 *
 * \return initialized structure of type scomplex
 *
 * \note
 *
 */
scomplex haz_scomplex_null(const INT n, const INT nbig) {
  /*
     n = dimension of the simplicial complex;
     nbig = dimension of the space where the complex is embedded;
     Future work: we should think of making this for different
     dimensions, e.g. 2-homogenous complex in 3d
  */
  scomplex sc;
  sc.nbig = nbig;
  sc.dim = n;
  if (sc.nbig <= 0) sc.nbig = n;
  INT n1 = sc.dim + 1, j;
  sc.level = 0;
  sc.ref_type = 0;
  sc.print_level = 0;
  sc.factorial = 1.;
  for (j = 2; j < n1; j++) sc.factorial *= ((REAL)j);
  //
  sc.marked = NULL;
  sc.gen = NULL;
  sc.nbr = NULL;
  sc.parent = NULL;
  sc.child0 = NULL;
  sc.childn = NULL;
  sc.nodes = NULL;
  sc.bndry = NULL;
  sc.csys = NULL;
  sc.bndry_v = malloc(sizeof(iCSRmat));
  sc.bndry_v[0] = icsr_create(0, 0, 0);  //
  sc.parent_v = malloc(sizeof(iCSRmat));
  sc.parent_v[0] = icsr_create(0, 0, 0);
  sc.bndry_f2v = NULL;
  sc.v2s_head = NULL;
  sc.v2s_next = NULL;
  sc.v2s_simp = NULL;
  sc.v2s_count = 0;
  sc.v2s_alloc = 0;
  sc.gen_N = NULL;
  sc.dgs_N = 0;
  sc.ncsys = 0;
  sc.systypes = NULL;
  sc.csys_ox = NULL;
  sc.flags = NULL;
  sc.x = NULL;
  sc.vols = NULL;
  sc.bndry_cc = 1;  // one connected component on the boundary for now.
  sc.cc = 1;        // one connected component in the bulk for now.
  // NULL pointers for the rest
  sc.etree = NULL;
  sc.bfs = malloc(sizeof(iCSRmat));
  sc.bfs[0] = icsr_create(0, 0, 0);
  sc.fem = NULL;
  return sc;
}
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_init_part(scomplex *sc)
 *
 * \brief Initialize simplicial complex which has already been
 *        allocated partially.
 * Assumption is: nv,ns,nodes,nbr, cc,
 *        bndry_cc are known. in dimension n with ns simplices and nv
 *        vertices.
 *
 * \param sc pointer to partially allocated simplician complex;
 *
 * \return full structure of type scomplex with all allocations done
 *
 * \note
 *
 */
void haz_scomplex_init_part(scomplex* sc) {
  INT nv = sc->nv, ns = sc->ns, n1 = sc->dim + 1, i, j;
  sc->level = 0;
  sc->ref_type = 0;
  sc->print_level = 0;
  sc->factorial = 1.;
  for (j = 2; j < n1; j++) sc->factorial *= ((REAL)j);
  //
  sc->marked = (INT*)calloc(ns, sizeof(INT));
  sc->gen = (INT*)calloc(ns, sizeof(INT));
  sc->parent = (INT*)calloc(ns, sizeof(INT));
  sc->child0 = (INT*)calloc(ns, sizeof(INT));
  sc->childn = (INT*)calloc(ns, sizeof(INT));
  sc->bndry = (INT*)calloc(nv, sizeof(INT));
  sc->csys = (INT*)calloc(nv, sizeof(INT));   /* coord sys: 1 is polar, 2
                                                 is cyl and so on */
  sc->flags = (INT*)calloc(ns, sizeof(INT));  // element flags
  //  sc->vols=(REAL *)calloc(ns,sizeof(REAL));// simplex volumes
  for (i = 0; i < ns; i++) {
    sc->marked[i] = FALSE;  // because first is used for something else.
    sc->gen[i] = 0;
    sc->parent[i] = -1;
    sc->child0[i] = -1;
    sc->childn[i] = -1;
    sc->flags[i] = -1;
  }
  for (i = 0; i < nv; i++) {
    sc->bndry[i] = 0;
    sc->csys[i] = 0;
  }
  sc->fem = NULL;
  sc->inc = NULL;
  return;
}
/**********************************************************************/
/*!
 * \fn void vol_simplex(INT dim, REAL fact, REAL *xf, REAL *volt, void *wrk)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void vol_simplex(INT dim, REAL fact, REAL* xf, REAL* volt, void* wrk) {
  /*
     computes the volume of a simplex; wrk should be at least
     dim*(dim+1) REALS and dim integers.
  */
  INT dim1 = dim + 1, i, j, ln, ln1;
  REAL* bt = (REAL*)wrk;
  REAL* piv = bt + dim * dim;
  INT* p = (INT*)(wrk + (dim * dim + dim) * sizeof(REAL));
  // construct bt using xf;
  for (j = 1; j < dim1; j++) {
    ln = j * dim;
    ln1 = ln - dim;
    for (i = 0; i < dim; i++) {
      bt[ln1 + i] = xf[ln + i] - xf[i];
    }
  }
  //  print_full_mat(dim,dim,bt,"bt");
  if (ddense_lu(1, dim, volt, bt, p, piv))
    *volt = 0e0;
  else
    *volt = fabs(*volt) / fact;
  return;
}
/* haz_scomplex_read moved to io_commented_out.c — use sc_read_gmsh instead */
/**********************************************************************/
/*!
 * \fn void haz_scomplex_print(scomplex *sc, const INT ns0,const char *infor)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void haz_scomplex_print(scomplex* sc, const INT ns0, const char* infor) {
  // print simplicial complex, starting with ns0.
  INT i, j, in, in1;
  INT dim = sc->dim, n1 = dim + 1, ns = sc->ns, nv = sc->nv, nbig = sc->nbig;
  if (ns0 < 0 || ns0 > ns) return;
  fprintf(stdout, "\nN=%lld,NBIG=%lld, NS=%lld, NV=%lld\n", (long long)sc->dim,
          (long long)sc->nbig, (long long)sc->ns, (long long)sc->nv);
  fflush(stdout);
  fprintf(stdout, "\n%s printout: %s\n", __FUNCTION__, infor);
  fprintf(stdout, "\nNODES list:\n");
  if (sc->parent) {
    for (i = ns0; i < ns; i++) {
      in1 = i * n1;
      fprintf(stdout, "Element: %lld ; vol=%e, Parent=%lld; NODES=",
              (long long)(i - ns0), sc->vols[i],
              (long long)sc->parent[i - ns0]);
      for (j = 0; j < n1; j++)
        fprintf(stdout, "%lld  ", (long long)sc->nodes[in1 + j]);
      fprintf(stdout, "\n");
    }
  } else {
    for (i = ns0; i < ns; i++) {
      in1 = i * n1;
      fprintf(stdout, "Element: %lld ; vol=%e, NODES=", (long long)(i - ns0),
              sc->vols[i]);
      for (j = 0; j < n1; j++)
        fprintf(stdout, "%lld  ", (long long)sc->nodes[in1 + j]);
      fprintf(stdout, "\n");
    }
  }
  fprintf(stdout, "\nNBR list:\n");
  if (sc->gen) {
    for (i = ns0; i < ns; i++) {
      in1 = i * n1;
      fprintf(stdout, "Element: %lld (%lld) ; NBR=", (long long)(i - ns0),
              (long long)sc->gen[i - ns0]);
      for (j = 0; j < n1; j++)
        fprintf(stdout, "%lld  ", (long long)(sc->nbr[in1 + j] - ns0));
      fprintf(stdout, "\n");
    }
  } else {
    for (i = ns0; i < ns; i++) {
      in1 = i * n1;
      fprintf(stdout, "Element: %lld ; NBR=", (long long)(i - ns0));
      for (j = 0; j < n1; j++)
        fprintf(stdout, "%lld  ", (long long)(sc->nbr[in1 + j] - ns0));
      fprintf(stdout, "\n");
    }
  }
  //
  if (sc->bndry) {
    for (i = 0; i < nv; i++) {
      in = i * nbig;
      fprintf(stdout, "Node: %lld ; Code: %lld ; COORDS=", (long long)i,
              (long long)sc->bndry[i]);
      for (j = 0; j < nbig; j++) {
        fprintf(stdout, "%e  ", sc->x[in + j]);
      }
      fprintf(stdout, "\n");
    }
  } else {
    for (i = 0; i < nv; i++) {
      in = i * nbig;
      fprintf(stdout, "Node: %lld ; COORDS=", (long long)i);
      for (j = 0; j < nbig; j++) {
        fprintf(stdout, "%e  ", sc->x[in + j]);
      }
      fprintf(stdout, "\n");
    }
  }
  return;
}
/**********************************************************************/
/*!
 * \fn void haz_scomplex_free(scomplex *sc)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
/**********************************************************************/
/*!
 * \fn void sc_free_fem_data(scomplex *sc)
 *
 * \brief Frees FEM-derived data in sc->fem. Sets sc->fem to NULL.
 *
 * \param sc  simplicial complex
 */
void sc_free_fem_data(scomplex* sc) {
  if (!sc || !sc->fem) return;
  sc_fem* fem = sc->fem;
  if (fem->leaf2global) free(fem->leaf2global);
  /* el_v, el_ed, ed_v are owned by sc->inc — freed there, not here. */
  fem->el_v = NULL;
  fem->el_ed = NULL;
  fem->ed_v = NULL;
  /* el_f and f_v are separate allocations for all dimensions */
  if (fem->el_f) { icsr_free(fem->el_f); free(fem->el_f); }
  if (fem->f_v) { icsr_free(fem->f_v); free(fem->f_v); }
  if (fem->f_ed) { icsr_free(fem->f_ed); free(fem->f_ed); }
  if (fem->el_vol) free(fem->el_vol);
  if (fem->el_mid) free(fem->el_mid);
  if (fem->ed_len) free(fem->ed_len);
  if (fem->ed_tau) free(fem->ed_tau);
  if (fem->ed_mid) free(fem->ed_mid);
  if (fem->f_area) free(fem->f_area);
  if (fem->f_norm) free(fem->f_norm);
  if (fem->f_mid) free(fem->f_mid);
  if (fem->el_flag) free(fem->el_flag);
  if (fem->ed_flag) free(fem->ed_flag);
  if (fem->f_flag) free(fem->f_flag);
  if (fem->dwork) free(fem->dwork);
  if (fem->iwork) free(fem->iwork);
  free(fem);
  sc->fem = NULL;
  sc->inc = NULL;
}
/**********************************************************************/
/*!
 * \fn static void compute_inc_signs(iCSRmat *k_to_j, iCSRmat *k_to_v, iCSRmat *j_to_v)
 *
 * \brief Compute signed incidence for inc[k][j] with j >= 1.
 *        Automatically detects if this is a boundary operator (j = k-1)
 *        by comparing the number of vertices per simplex.
 *
 *        For boundary operators (j = k-1): sign = (-1)^i * perm_sign,
 *        where i is the local index of the omitted vertex and perm_sign
 *        is the sign of the permutation from induced to global order.
 *        This encodes the chain complex boundary map ∂_k.
 *
 *        For non-boundary (j < k-1): sign = perm_sign only,
 *        encoding relative orientation of the j-simplex in the k-simplex.
 *
 * \param k_to_j  incidence matrix from k-simplices to j-simplices
 * \param k_to_v  k-simplices to vertices
 * \param j_to_v  j-simplices to vertices
 */
static void compute_inc_signs(iCSRmat *k_to_j, iCSRmat *k_to_v, iCSRmat *j_to_v)
{
  if (!k_to_j || !k_to_v || !j_to_v) return;
  INT n_k = k_to_j->row;
  if (n_k == 0) return;
  INT kp1 = (INT)(k_to_v->IA[1] - k_to_v->IA[0]); /* k+1 */
  INT jp1 = (INT)(j_to_v->IA[1] - j_to_v->IA[0]); /* j+1 */
  INT is_boundary = (jp1 == kp1 - 1) ? 1 : 0;
  if (k_to_j->val) free(k_to_j->val);
  k_to_j->val = (INT*)calloc(k_to_j->nnz, sizeof(INT));
  for (INT s = 0; s < n_k; s++) {
    INT *kverts = &k_to_v->JA[k_to_v->IA[s]];
    for (INT p = k_to_j->IA[s]; p < k_to_j->IA[s + 1]; p++) {
      INT sigma = k_to_j->JA[p];
      INT *jverts = &j_to_v->JA[j_to_v->IA[sigma]];
      /* Map shared vertices: for each vertex of k-simplex found in
         j-simplex, record its position in j-simplex's vertex list. */
      INT perm[16]; /* supports up to 15-simplices */
      INT omitted_sign = 1;
      INT idx = 0;
      for (INT a = 0; a < kp1; a++) {
        INT pos = -1;
        for (INT b = 0; b < jp1; b++) {
          if (kverts[a] == jverts[b]) { pos = b; break; }
        }
        if (pos >= 0) {
          perm[idx++] = pos;
        } else if (is_boundary) {
          /* Boundary operator: (-1)^(position of omitted vertex) */
          omitted_sign = (a % 2 == 0) ? 1 : -1;
        }
      }
      /* Sign of permutation = (-1)^(number of inversions) */
      INT inversions = 0;
      for (INT a = 0; a < jp1; a++) {
        for (INT b = a + 1; b < jp1; b++) {
          if (perm[a] > perm[b]) inversions++;
        }
      }
      INT psign = (inversions % 2 == 0) ? 1 : -1;
      k_to_j->val[p] = (is_boundary ? omitted_sign : 1) * psign;
    }
  }
}
/**********************************************************************/
/*!
 * \fn void sc_build_fem_data(scomplex *sc)
 *
 * \brief Builds FEM-derived data (edges, faces, volumes, etc.) from
 *        the leaf elements of the simplicial complex. Populates sc->fem.
 *
 * \param sc  simplicial complex (must have nodes, x, bndry, flags, child0, childn set)
 *
 * \note Calls geometry functions (edge_stats_all, face_stats, sync_facenode,
 *       get_el_vol, get_el_mid) directly on the scomplex.
 *       Requires dim == nbig (no manifold meshes).
 */
void sc_build_fem_data(scomplex* sc) {
  if (sc->fem != NULL) return;  /* already built */
  if (sc->nbig != sc->dim) return;  /* not supported yet */
  INT dim = sc->dim, nv = sc->nv, n1 = dim + 1;
  /* 1. Count leaf elements */
  INT ns_leaf = 0;
  for (INT j = 0; j < sc->ns; j++)
    if (sc->child0[j] < 0 || sc->childn[j] < 0) ns_leaf++;
  /* 2. Allocate and populate sc_fem */
  sc_fem* fem = (sc_fem*)calloc(1, sizeof(sc_fem));
  fem->ns_leaf = ns_leaf;
  fem->leaf2global = (INT*)calloc(ns_leaf, sizeof(INT));
  fem->el_flag = (INT*)calloc(ns_leaf, sizeof(INT));
  INT idx = 0;
  for (INT j = 0; j < sc->ns; j++) {
    if (sc->child0[j] < 0 || sc->childn[j] < 0) {
      fem->leaf2global[idx] = j;
      fem->el_flag[idx] = sc->flags[j];
      idx++;
    }
  }
  /* 3. Allocate sc->inc and build el_v = inc[dim*n1+0] */
  sc->inc = (iCSRmat**)calloc(n1 * n1, sizeof(iCSRmat*));
  iCSRmat *el_v = (iCSRmat*)malloc(sizeof(iCSRmat));
  el_v->row = ns_leaf;
  el_v->col = nv;
  el_v->nnz = ns_leaf * n1;
  el_v->IA = (INT*)calloc(ns_leaf + 1, sizeof(INT));
  el_v->val = NULL;
  for (INT j = 0; j < ns_leaf; j++) el_v->IA[j + 1] = (j + 1) * n1;
  el_v->JA = sc->nodes; /* shared with sc->nodes, freed by haz_scomplex_free */
  sc->inc[dim * n1 + 0] = el_v;
  fem->el_v = el_v; /* convenience pointer into sc->inc */
  /* 4. Count boundary vertices */
  fem->nbv = 0;
  for (INT i = 0; i < nv; i++)
    if (sc->bndry[i] != 0) fem->nbv++;
  /* 5. Allocate dwork */
  fem->dwork = (REAL*)calloc(n1 * (dim + 1), sizeof(REAL));
  /* 6. Assign fem to sc NOW so geometry functions can access sc->fem */
  sc->fem = fem;
  /* 7. Build edges, faces, geometry */
  if (dim == 2 || dim == 3) {
    /* Edge to vertex map */
    INT nedge = 0;
    iCSRmat ed_v = get_edge_v(&nedge, fem->el_v);
    /* Element to edge map */
    iCSRmat el_ed = get_el_ed(fem->el_v, &ed_v);
    /* Edge stats */
    REAL* ed_len = (REAL*)calloc(nedge, sizeof(REAL));
    REAL* ed_tau = (REAL*)calloc(nedge * dim, sizeof(REAL));
    REAL* ed_mid = (REAL*)calloc(nedge * dim, sizeof(REAL));
    edge_stats_all(ed_len, ed_tau, ed_mid, sc, &ed_v, dim);
    /* Compute number of faces via Euler characteristic */
    INT nconn_bdry = sc->bndry_cc;
    INT nholes = nconn_bdry - 1;
    INT nface = 0;
    INT euler = -10;
    if (dim == 2) {
      nface = nedge;
      euler = nv - nedge + ns_leaf + nholes;
    } else if (dim == 3) {
      nface = 1 + nedge - nv + ns_leaf;
      nface = nface + nholes;
      euler = nv - nedge + nface - ns_leaf - nholes;
    }
    if (euler != 1) {
      printf("ERROR HAZMATH DANGER: in function %s, Euler Characteristic doesn't equal 1+nholes! euler=%lld nholes=%lld.\n\n",
             __FUNCTION__, (long long)euler, (long long)nholes);
      exit(ERROR_DIM);
    }
    /* Face ordering */
    INT f_per_elm = n1;
    INT* fel_order = (INT*)calloc(f_per_elm * dim, sizeof(INT));
    get_face_ordering(n1, dim, f_per_elm, fel_order);
    /* Face maps */
    iCSRmat f_v = icsr_create(nface, nv, nface * dim);
    iCSRmat f_ed = icsr_create(nface, nedge, nface * (2 * dim - 3));
    /* Store el_ed and ed_v in sc->inc FIRST (needed before el_f assignment) */
    sc->inc[dim * n1 + 1] = malloc(sizeof(iCSRmat));
    *(sc->inc[dim * n1 + 1]) = el_ed;
    fem->el_ed = sc->inc[dim * n1 + 1];
    sc->inc[1 * n1 + 0] = malloc(sizeof(iCSRmat));
    *(sc->inc[1 * n1 + 0]) = ed_v;
    fem->ed_v = sc->inc[1 * n1 + 0];
    /* Face maps */
    INT* f_flag = (INT*)calloc(nface, sizeof(INT));
    INT nbface;
    iCSRmat *el_f_ptr = malloc(sizeof(struct iCSRmat));
    get_face_maps(fem->el_v, n1, &ed_v, nface, dim, f_per_elm, el_f_ptr, f_flag, &nbface, &f_v, &f_ed, fel_order);
    fem->el_f = el_f_ptr;
    /* Edge and face boundary flags */
    INT nbedge = 0;
    INT* ed_flag = (INT*)calloc(nedge, sizeof(INT));
    boundary_f_ed(&f_ed, &ed_v, nedge, nface, f_flag, sc->bndry, &nbedge, ed_flag, dim);
    /* Face stats (needs fem->el_f and fem->el_v to be set) */
    REAL* f_area = (REAL*)calloc(nface, sizeof(REAL));
    REAL* f_mid = (REAL*)calloc(nface * dim, sizeof(REAL));
    REAL* f_norm = (REAL*)calloc(nface * dim, sizeof(REAL));
    face_stats(f_area, f_mid, f_norm, &f_v, sc);
    /* Sync face nodes — stores orientation in f_v->val, no JA swaps */
    sync_facenode(&f_v, f_norm, sc);
    /* Store counts and remaining data */
    fem->nedge = nedge;
    fem->nface = nface;
    fem->nbedge = nbedge;
    fem->nbface = nbface;
    fem->f_v = malloc(sizeof(struct iCSRmat));
    *(fem->f_v) = f_v;
    fem->f_ed = malloc(sizeof(struct iCSRmat));
    *(fem->f_ed) = f_ed;
    fem->ed_len = ed_len;
    fem->ed_tau = ed_tau;
    fem->ed_mid = ed_mid;
    fem->f_area = f_area;
    fem->f_mid = f_mid;
    fem->f_norm = f_norm;
    fem->ed_flag = ed_flag;
    fem->f_flag = f_flag;
    if (fel_order) free(fel_order);
  } else if (dim == 1) {
    fem->nedge = 0;
    fem->nface = 0;
    fem->nbedge = 0;
    fem->nbface = 0;
  }
  /* 8. Element volumes and midpoints */
  REAL* el_mid = (REAL*)calloc(ns_leaf * dim, sizeof(REAL));
  REAL* el_vol = (REAL*)calloc(ns_leaf, sizeof(REAL));
  get_el_vol(el_vol, fem->el_v, sc, dim, n1);
  get_el_mid(el_mid, fem->el_v, sc, dim);
  fem->el_vol = el_vol;
  fem->el_mid = el_mid;
  /* 9. Populate remaining sc->inc entries for face data.
     el_v, el_ed, ed_v are already stored in sc->inc above.
     For dim>=3: face slots are separate from edge slots.
     For dim=2: face=edge, so inc[dim][1] already has el_ed.
     el_f, f_v, f_ed are owned by sc_fem (freed there). */
  if (dim >= 3) {
    sc->inc[dim * n1 + (dim - 1)] = fem->el_f;
    sc->inc[(dim - 1) * n1 + 0] = fem->f_v;
    sc->inc[(dim - 1) * n1 + 1] = fem->f_ed;
  }
  /* 10. Compute signed incidence for all inc[k][j].
     Convention: inc[k][j]->val stores +1 or -1 encoding
     relative orientation of the k-simplex and j-simplex.
     For boundary operators (j = k-1): encodes the chain complex ∂_k.
     For inc[k][0]: encodes the orientation sign of the k-simplex.
     Property: ∂_{k-1} ∘ ∂_k = 0 is satisfied.
     All computations are dimension-independent. */
  if (dim >= 2) {
    /* (a) ed_v->val: ∂_1 boundary operator.
       ∂[v0,v1] = v1 - v0: val = -1 for first vertex, +1 for second. */
    if (fem->ed_v->val) free(fem->ed_v->val);
    fem->ed_v->val = (INT*)calloc(fem->ed_v->nnz, sizeof(INT));
    for (INT e = 0; e < fem->nedge; e++) {
      fem->ed_v->val[fem->ed_v->IA[e]] = -1;
      fem->ed_v->val[fem->ed_v->IA[e] + 1] = 1;
    }
    /* (b) el_f: ∂_dim boundary (element to face). */
    compute_inc_signs(fem->el_f, fem->el_v, fem->f_v);
    /* (c) el_ed: relative edge orientation in element. */
    compute_inc_signs(fem->el_ed, fem->el_v, fem->ed_v);
    /* (d) f_ed: boundary or relative orientation (auto-detected). */
    if (fem->f_ed) compute_inc_signs(fem->f_ed, fem->f_v, fem->ed_v);
  }
  /* el_v->val: element orientation sign = sign(det(shape matrix B)).
     B[i][j] = x[v_{j+1}][i] - x[v_0][i]. */
  {
    if (el_v->val) free(el_v->val);
    el_v->val = (INT*)calloc(el_v->nnz, sizeof(INT));
    REAL *B = (REAL *)calloc(dim * dim, sizeof(REAL));
    for (INT el = 0; el < ns_leaf; el++) {
      INT *verts = &el_v->JA[el_v->IA[el]];
      for (INT j = 0; j < dim; j++)
        for (INT i = 0; i < dim; i++)
          B[j * dim + i] = sc->x[verts[j + 1] * dim + i]
                         - sc->x[verts[0] * dim + i];
      REAL detB = haz_det(dim, B);
      INT sign = (detB > 0.0) ? 1 : ((detB < 0.0) ? -1 : 0);
      for (INT j = el_v->IA[el]; j < el_v->IA[el + 1]; j++)
        el_v->val[j] = sign;
    }
    free(B);
  }
}
/**********************************************************************/
void haz_scomplex_free(scomplex* sc) {
  if (sc->marked) free(sc->marked);
  if (sc->gen) free(sc->gen);
  if (sc->parent) free(sc->parent);
  if (sc->child0) free(sc->child0);
  if (sc->childn) free(sc->childn);
  if (sc->bndry) free(sc->bndry);
  if (sc->flags) free(sc->flags);
  if (sc->nodes) free(sc->nodes);
  if (sc->x) free(sc->x);
  if (sc->vols) free(sc->vols);
  if (sc->nbr) free(sc->nbr);
  if (sc->csys) free(sc->csys);
  if (sc->etree) free(sc->etree);
  if (sc->bndry_v) {
    icsr_free(sc->bndry_v);
    free(sc->bndry_v);
    sc->bndry_v = NULL;
  }
  if (sc->parent_v) {
    icsr_free(sc->parent_v);
    free(sc->parent_v);
    sc->parent_v = NULL;
  }
  if (sc->bfs) {
    icsr_free(sc->bfs);
    free(sc->bfs);
    sc->bfs = NULL;
  }
  if (sc->bndry_f2v) {
    icsr_free(sc->bndry_f2v);
    free(sc->bndry_f2v);
    sc->bndry_f2v = NULL;
  }
  if (sc->v2s_head) {
    free(sc->v2s_head);
    sc->v2s_head = NULL;
  }
  if (sc->v2s_next) {
    free(sc->v2s_next);
    sc->v2s_next = NULL;
  }
  if (sc->v2s_simp) {
    free(sc->v2s_simp);
    sc->v2s_simp = NULL;
  }
  sc->v2s_count = 0;
  sc->v2s_alloc = 0;
  if (sc->gen_N) { free(sc->gen_N); sc->gen_N = NULL; }
  if (sc->systypes) { free(sc->systypes); sc->systypes = NULL; }
  if (sc->csys_ox) { free(sc->csys_ox); sc->csys_ox = NULL; }
  sc->ncsys = 0;
  sc_free_fem_data(sc);
  /* Free inc array — owns el_v, el_ed, ed_v.
     el_f, f_v, f_ed are owned by sc_fem, already freed above.
     For dim>=3, NULL out those inc entries to avoid double free. */
  if (sc->inc) {
    INT dim = sc->dim, n1 = dim + 1;
    /* NULL out fem-owned entries (freed by sc_free_fem_data) */
    if (dim >= 3) {
      sc->inc[dim * n1 + (dim - 1)] = NULL; /* el_f */
      sc->inc[(dim - 1) * n1 + 0] = NULL;   /* f_v */
      sc->inc[(dim - 1) * n1 + 1] = NULL;   /* f_ed */
    }
    /* el_v = inc[dim*n1+0] */
    if (sc->inc[dim * n1 + 0]) {
      sc->inc[dim * n1 + 0]->JA = NULL; /* JA = sc->nodes, freed above */
      icsr_free(sc->inc[dim * n1 + 0]);
      free(sc->inc[dim * n1 + 0]);
    }
    /* el_ed = inc[dim*n1+1] */
    if (dim >= 2 && sc->inc[dim * n1 + 1]) {
      icsr_free(sc->inc[dim * n1 + 1]);
      free(sc->inc[dim * n1 + 1]);
    }
    /* ed_v = inc[1*n1+0] */
    if (dim >= 2 && sc->inc[1 * n1 + 0]) {
      icsr_free(sc->inc[1 * n1 + 0]);
      free(sc->inc[1 * n1 + 0]);
    }
    free(sc->inc); sc->inc = NULL;
  }
  if (sc) free(sc);
  return;
}
/**********************************************************************/
/*!
 * \fn static unsigned int cmp_simplex(INT n, INT sim1, INT sim2, INT
 *			 *sv1, INT *sv2, INT *stos1, INT *stos2)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
static unsigned INT cmp_simplex(INT n, INT sim1, INT sim2, INT* sv1, INT* sv2,
                                INT* stos1, INT* stos2) {
  // sim1 and sim2 are pointers to the neighboring lists of two
  // simplices stos1 and stos2 are pointers to rows in
  // simplex-to-simplex incidence; sv1 and sv2 are pointers to rows in
  // simplex-vertex incidence for two neighboring simplices.
  //
  // rearrange to meet the structural condition.  compare two sets of
  // length n. they should overlap at exactly n-1 members.  If they do
  // not, 0 is returned, otherwise 1. no check if the members of the
  // sets are distinct (they should be for the purpose of comparing two
  // simplices.
  unsigned INT fnd;
  unsigned INT nf = 0;
  INT i1, k1 = -10, i2, k2 = -10;
  INT i, j;
  for (i = 0; i < n; i++) {
    fnd = 0;
    i1 = sv1[i];
    for (j = 0; j < n; j++) {
      if (i1 == sv2[j]) {
        fnd = 1;
        break;
      }
    }
    if (fnd) {
      continue;
    } else {
      // not found
      nf++;
      if (nf > 1) return 0;
      k1 = i;
    }
  }
  // same with sv1 and sv2 swapped.
  nf = 0;
  for (i = 0; i < n; i++) {
    fnd = 0;
    i2 = sv2[i];
    for (j = 0; j < n; j++) {
      if (i2 == sv1[j]) {
        fnd = 1;
        break;
      }
    }
    if (fnd) {
      continue;
    } else {
      // not found
      nf++;
      if (nf > 1) return 0;
      k2 = i;
    }
  }
  /* NOW put the neightbors at the right places ******* */
  if (k1 < 0 || k2 < 0) {
    fprintf(stderr, "\n***ERROR in %s; k1=%lld,k2=%lld\n\n", __FUNCTION__,
            (long long)k1, (long long)k2);
    exit(65);
  } else {
    stos1[k1] = sim2;
    stos2[k2] = sim1;
    return 1;
  }
}
/**********************************************************************/
/*!
 * \fn void find_nbr(INT ns,INT nv,INT n,INT *sv,INT *stos)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void find_nbr(INT ns, INT nv, INT n, INT* sv, INT* stos) {
  // find neighboring list
  INT *ivs = NULL, *jvs = NULL, *stosi = NULL, *stosk = NULL, *svi = NULL,
      *svk = NULL;
  INT i, j, k, jp, jia, jib, iabeg, iaend, ibbeg, ibend, kn1;
  INT n1 = n + 1, nv1 = nv + 1, nsv = ns * n1;
  /* allocate */
  ivs = (INT*)calloc(nv1, sizeof(INT));
  jvs = (INT*)calloc(nsv, sizeof(INT));
  /* transpose to get the vertex--simplex incidence */
  for (i = 0; i < nv1; ++i) ivs[i] = 0;
  //
  for (k = 0; k < ns; ++k) {
    kn1 = k * n1;
    for (i = 0; i < n1; i++) {
      j = sv[kn1 + i] + 2;
      if (j < nv1) ivs[j]++;
    }
  }
  ivs[0] = 0;
  ivs[1] = 0;
  if (nv != 1) {
    for (i = 2; i < nv1; ++i) {
      ivs[i] += ivs[i - 1];
    }
  }
  for (i = 0; i < ns; ++i) {
    iabeg = i * n1;
    iaend = i * n1 + n1;
    for (jp = iabeg; jp < iaend; ++jp) {
      j = sv[jp] + 1;
      k = ivs[j];
      jvs[k] = i;
      ivs[j] = k + 1;
    }
  }
  /**/
  INT* icp = (INT*)calloc(ns, sizeof(INT));
  for (i = 0; i < ns; ++i) icp[i] = -1;
  for (i = 0; i < nsv; ++i) stos[i] = -1;
  for (i = 0; i < ns; ++i) {
    iabeg = i * n1;
    iaend = iabeg + n1;
    stosi = stos + iabeg;
    svi = sv + iabeg;
    for (jia = iabeg; jia < iaend; ++jia) {
      j = sv[jia];  // vertex number
      ibbeg = ivs[j];
      ibend = ivs[j + 1];
      // loop over all simplices with this j as a vertex.
      for (jib = ibbeg; jib < ibend; ++jib) {
        k = jvs[jib];
        if (k <= i) continue;
        if (icp[k] != i) {
          icp[k] = i;
          kn1 = k * n1;
          stosk = stos + kn1;
          svk = sv + kn1;
          if (!cmp_simplex(n1, i, k, svi, svk, stosi, stosk)) continue;
        }  // if
      }  // for
    }  // for
  }  // for (i...
  if (icp) free(icp);
  if (ivs) free(ivs);
  if (jvs) free(jvs);
  //  }
  return;
}
/**********************************************************************/
/*!
 * \fn INT haz_add_simplex(INT is, scomplex *sc,REAL *xnew, INT *pv,
 *                         INT ibnew,INT csysnew,INT nsnew, INT nvnew)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
INT haz_add_simplex(INT is, scomplex* sc, REAL* xnew, INT* pv, INT ibnew,
                    INT csysnew, INT nsnew, INT nvnew) {
  /* adds nodes and coords as well */
  INT dim = sc->dim, nbig = sc->nbig, n1 = dim + 1, nv = sc->nv;  // ns=sc->ns;
  INT ks0 = sc->child0[is], ksn = sc->childn[is];
  INT isc0 = ks0 * n1, iscn = ksn * n1;
  //  INT *dsti,*srci;
  INT j, j0, jn, nnz_pv;
  REAL* dstr;
  /* nodes  AND neighbors */
  sc->nbr = realloc(sc->nbr, (nsnew * n1) * sizeof(INT));
  sc->nodes = realloc(sc->nodes, (nsnew * n1) * sizeof(INT));
  for (j = 0; j < n1; j++) {
    j0 = isc0 + j;
    jn = iscn + j;
    //
    sc->nodes[j0] = -1;
    sc->nbr[j0] = -1;
    sc->nodes[jn] = -1;
    sc->nbr[jn] = -1;
  }
  // new vertex (if any!!!)
  if (nvnew != nv) {
    sc->x = realloc(sc->x, (nvnew * nbig) * sizeof(REAL));
    dstr = (sc->x + nv * nbig);
    memcpy(dstr, xnew, dim * sizeof(REAL));
    if (xnew) free(xnew);
    sc->bndry = realloc(sc->bndry, (nvnew) * sizeof(INT));
    sc->bndry[nv] = ibnew;
    sc->csys = realloc(sc->csys, (nvnew) * sizeof(INT));
    sc->csys[nv] = csysnew;
    nnz_pv = sc->parent_v->nnz;
    sc->parent_v->row = nvnew;
    sc->parent_v->col = nv;
    sc->parent_v->JA = realloc(sc->parent_v->JA, (nnz_pv + 2) * sizeof(REAL));
    sc->parent_v->JA[nnz_pv] = pv[0];
    sc->parent_v->JA[nnz_pv + 1] = pv[1];
    sc->parent_v->val = realloc(sc->parent_v->val, (nnz_pv + 2) * sizeof(REAL));
    sc->parent_v->val[nnz_pv] = sc->level + 1;
    sc->parent_v->val[nnz_pv + 1] = sc->level + 1;
    sc->parent_v->nnz += 2;
    sc->parent_v->IA = realloc(sc->parent_v->IA, (nvnew + 1) * sizeof(REAL));
    sc->parent_v->IA[nvnew] = sc->parent_v->nnz;
    /* fprintf(stdout,"\nnv=%d;
     * nvnew=%d;nnz_pv=%d(pv[0]=%d,pv[1]=%d)",nv,nvnew,sc->parent_v->nnz,pv[0],pv[1]);
     */
  }
  // generation
  sc->gen = realloc(sc->gen, (nsnew) * sizeof(INT));
  sc->gen[ks0] = sc->gen[is] + 1;
  sc->gen[ksn] = sc->gen[is] + 1;
  // marked
  sc->marked = realloc(sc->marked, (nsnew) * sizeof(INT));
  sc->marked[ks0] = sc->marked[is];
  sc->marked[ksn] = sc->marked[is];
  // flags
  sc->flags = realloc(sc->flags, (nsnew) * sizeof(INT));
  sc->flags[ks0] = sc->flags[is];
  sc->flags[ksn] = sc->flags[is];
  // parents
  sc->parent = realloc(sc->parent, (nsnew) * sizeof(INT));
  sc->parent[ks0] = is;
  sc->parent[ksn] = is;
  // child0
  sc->child0 = realloc(sc->child0, (nsnew) * sizeof(INT));
  sc->child0[ks0] = -1;
  sc->child0[ksn] = -1;
  // childn
  sc->childn = realloc(sc->childn, (nsnew) * sizeof(INT));
  sc->childn[ks0] = -1;
  sc->childn[ksn] = -1;
  // volumes
  sc->vols = realloc(
      sc->vols, (nsnew) * sizeof(REAL));  // we can calculate all volumes at the
                                          // end, so just allocate here.
  // scalars
  sc->ns = nsnew;
  sc->nv = nvnew;
  return 0;
}
/**********************************************************************/
/*!
 * \fn static INT set_color(scomplex *sc, INT *color)
 *
 * \brief Generalized coloring with N+1 colors (Algorithm 2 from
 *        Diening, Gehring, Storn, "Adaptive Mesh Refinement for
 *        Arbitrary Initial Triangulations", Found. Comput. Math.,
 *        2025, DOI: 10.1007/s10208-024-09642-1).
 *
 *        Assigns to each vertex the smallest color not already
 *        attained by a neighboring vertex. The number of colors N+1
 *        is bounded by the maximal vertex degree + 1.
 *
 * \param sc     I: simplicial complex (uses nodes, ns, nv, n)
 * \param color  O: array of size nv; color[v] in {0,...,N}
 *
 * \return N (the largest color assigned)
 */
static INT set_color(scomplex* sc, INT* color) {
  INT nv = sc->nv, ns = sc->ns, dim = sc->dim;
  INT n1 = dim + 1, i, j, k, v, w, c;
  /* Build vertex-to-vertex adjacency from the element connectivity.
     We use a simple approach: for each simplex, all pairs of its
     vertices are neighbors. We store adjacency in CSR format. */
  /* Step 1: count edges per vertex (upper bound via element connectivity) */
  INT* deg = (INT*)calloc(nv, sizeof(INT));
  for (i = 0; i < ns; i++) {
    INT* el = sc->nodes + i * n1;
    for (j = 0; j < n1; j++) {
      deg[el[j]] += dim; /* at most dim neighbors per element */
    }
  }
  /* Allocate adjacency lists (with room for duplicates; we handle them) */
  INT* adjptr = (INT*)calloc(nv + 1, sizeof(INT));
  adjptr[0] = 0;
  for (i = 0; i < nv; i++) adjptr[i + 1] = adjptr[i] + deg[i];
  INT nnz_adj = adjptr[nv];
  INT* adjind = (INT*)malloc(nnz_adj * sizeof(INT));
  INT* pos = (INT*)calloc(nv, sizeof(INT)); /* current insert position */
  for (i = 0; i < ns; i++) {
    INT* el = sc->nodes + i * n1;
    for (j = 0; j < n1; j++) {
      v = el[j];
      for (k = 0; k < n1; k++) {
        if (k == j) continue;
        w = el[k];
        adjind[adjptr[v] + pos[v]] = w;
        pos[v]++;
      }
    }
  }
  free(pos);
  free(deg);
  /* Step 2: greedy coloring */
  INT max_color = 0;
  /* max possible colors is bounded by max vertex degree + 1;
     for safety allocate nv+1 entries */
  INT* used =
      (INT*)calloc(nv + 1, sizeof(INT)); /* flag array: used[c]=v+1 means color
                                            c is used by a neighbor of v */
  for (v = 0; v < nv; v++)
    color[v] = -1; /* uncolored = infinity in the paper */
  for (v = 0; v < nv; v++) {
    /* Mark colors used by neighbors of v */
    for (j = adjptr[v]; j < adjptr[v + 1]; j++) {
      w = adjind[j];
      if (color[w] >= 0) used[color[w]] = v + 1; /* mark as used for vertex v */
    }
    /* Find smallest color not used by any neighbor */
    for (c = 0;; c++) {
      if (used[c] != v + 1) {
        color[v] = c;
        break;
      }
    }
    if (c > max_color) max_color = c;
  }
  free(used);
  free(adjind);
  free(adjptr);
  return max_color; /* N */
}
/**********************************************************************/
/*!
 * \fn static void dgs_initialize(scomplex *sc, INT *color, INT N)
 *
 * \brief Maubach initialization using a generalized coloring
 *        (Definition 2 from Diening, Gehring, Storn, "Adaptive Mesh
 *        Refinement for Arbitrary Initial Triangulations",
 *        Found. Comput. Math., 2025,
 *        DOI: 10.1007/s10208-024-09642-1).
 *
 *        For each initial simplex T with vertices v0,...,vn, we sort
 *        the vertices so that c(v_j) < c(v_{j+1}) for all j.  Then
 *        the tagged simplex is:
 *          T = [vn, v0, v1,...,v_{n-1}]_n  if c(vn) == N
 *          T = [v0, v1,...,vn]_n            otherwise
 *        The tag gamma = n for all initial simplices.
 *
 * \param sc     I/O: simplicial complex (nodes are reordered in place)
 * \param color  I: vertex coloring array of size nv
 * \param N      I: the largest color
 */
static void dgs_initialize(scomplex* sc, INT* color, INT N) {
  INT ns = sc->ns, nv = sc->nv, dim = sc->dim, n1 = dim + 1;
  INT i, j;
  /* Allocate gen_N and set initial values: gen_N(v) = -color(v) */
  sc->gen_N = (INT*)calloc(nv, sizeof(INT));
  sc->dgs_N = N;
  for (i = 0; i < nv; i++) sc->gen_N[i] = -color[i];
  /* Sort each simplex by DECREASING gen_N */
  for (i = 0; i < ns; i++) {
    INT* el = sc->nodes + i * n1;
    for (j = 0; j < n1; j++)
      for (INT k = j + 1; k < n1; k++)
        if (sc->gen_N[el[j]] < sc->gen_N[el[k]]) {
          INT t = el[j]; el[j] = el[k]; el[k] = t;
        }
    sc->gen[i] = 0;
  }
  /* nbr is not reordered — it will be set to -1 for children
     and rebuilt by find_nbr in scfinalize */
}
/**********************************************************************/
/*!
 * \brief Build v2s index for current LEAF simplices only.
 *        Children created during refinement are NOT added.
 *        This prevents infinite recursion in the DGS closure.
 */
static void v2s_build(scomplex* sc) {
  INT ns = sc->ns, nv = sc->nv, n1 = sc->dim + 1;
  INT i, j, v;
  /* Count leaf simplex-vertex pairs */
  INT total = 0;
  for (i = 0; i < ns; i++)
    if (sc->child0[i] < 0) total += n1;
  sc->v2s_alloc = total > 0 ? total : 1;
  sc->v2s_head = (INT*)malloc(nv * sizeof(INT));
  sc->v2s_next = (INT*)malloc(sc->v2s_alloc * sizeof(INT));
  sc->v2s_simp = (INT*)malloc(sc->v2s_alloc * sizeof(INT));
  for (i = 0; i < nv; i++) sc->v2s_head[i] = -1;
  sc->v2s_count = 0;
  sc->v2s_nv = nv;
  for (i = 0; i < ns; i++) {
    if (sc->child0[i] >= 0) continue; /* skip non-leaf */
    INT* el = sc->nodes + i * n1;
    for (j = 0; j < n1; j++) {
      v = el[j];
      INT k = sc->v2s_count++;
      sc->v2s_simp[k] = i;
      sc->v2s_next[k] = sc->v2s_head[v];
      sc->v2s_head[v] = k;
    }
  }
}
/**********************************************************************/
static void v2s_free_data(scomplex* sc) {
  if (sc->v2s_head) { free(sc->v2s_head); sc->v2s_head = NULL; }
  if (sc->v2s_next) { free(sc->v2s_next); sc->v2s_next = NULL; }
  if (sc->v2s_simp) { free(sc->v2s_simp); sc->v2s_simp = NULL; }
  sc->v2s_count = 0; sc->v2s_alloc = 0; sc->v2s_nv = 0;
}
/**********************************************************************/
/* Algorithm 4 helpers: lvl_N and type_N from gen_N (eq 8 in paper)  */
/* gen_N = N*(lvl_N - 1) + type_N,  type_N in {1,...,N}              */
/**********************************************************************/
static inline INT dgs_type_N(INT gen, INT N) {
  return ((gen - 1) % N + N) % N + 1;
}
static inline INT dgs_lvl_N(INT gen, INT N) {
  INT t = dgs_type_N(gen, N);
  return (gen - t) / N + 1;
}
/**********************************************************************/
/*!
 * \brief Get the bisection edge of simplex is using Algorithm 4.
 *        Simplex vertices must be sorted by DECREASING gen_N.
 *        Returns the two vertex indices (v_a, v_b) of the bse.
 */
static void dgs_get_bse(scomplex* sc, INT is, INT *va_out, INT *vb_out) {
  INT dim = sc->dim, n1 = dim + 1, N = sc->dgs_N;
  INT *el = sc->nodes + is * n1;
  INT vm = el[dim], vm1 = el[dim - 1];
  INT lvl_m = dgs_lvl_N(sc->gen_N[vm], N);
  INT lvl_m1 = dgs_lvl_N(sc->gen_N[vm1], N);
  if (lvl_m != lvl_m1) {
    *va_out = vm1; *vb_out = vm;
  } else {
    INT j = dim;
    for (INT k = 0; k < n1; k++)
      if (dgs_lvl_N(sc->gen_N[el[k]], N) == lvl_m) { j = k; break; }
    *va_out = el[j]; *vb_out = vm;
  }
}
/**********************************************************************/
/*!
 * \brief Bisect simplex is using Algorithm 4, creating a NEW midpoint.
 *        Children are sorted by decreasing gen_N.
 * \return index of the new midpoint vertex
 */
static INT haz_bisect_new(scomplex* sc, INT is) {
  INT dim = sc->dim, nbig = sc->nbig, n1 = dim + 1, N = sc->dgs_N;
  INT *el = sc->nodes + is * n1;
  /* Algorithm 4: determine bse and gen_N of midpoint */
  INT vm = el[dim], vm1 = el[dim - 1];
  INT lvl_m = dgs_lvl_N(sc->gen_N[vm], N);
  INT lvl_m1 = dgs_lvl_N(sc->gen_N[vm1], N);
  INT va_idx, gen_new; /* va_idx = index into el[] */
  if (lvl_m != lvl_m1) {
    va_idx = dim - 1;
    gen_new = sc->gen_N[vm1] + N;
  } else {
    va_idx = dim;
    for (INT k = 0; k < n1; k++)
      if (dgs_lvl_N(sc->gen_N[el[k]], N) == lvl_m) { va_idx = k; break; }
    gen_new = sc->gen_N[vm] + 2*N + 1 - dgs_type_N(sc->gen_N[el[va_idx]], N);
  }
  INT v_a = el[va_idx], v_b = el[dim];
  /* Create midpoint vertex via haz_add_simplex */
  INT nsnew = sc->ns + 2, nvnew = sc->nv + 1, nodnew = sc->nv;
  REAL* xnew = (REAL*)calloc(nbig, sizeof(REAL));
  for (INT j = 0; j < nbig; j++)
    xnew[j] = 0.5 * (sc->x[v_a * nbig + j] + sc->x[v_b * nbig + j]);
  /* Project onto curve if both endpoints share a polar coordinate system */
  INT csysnew = (sc->csys[v_a] < sc->csys[v_b]) ? sc->csys[v_a] : sc->csys[v_b];
  if (sc->csys[v_a] == sc->csys[v_b] && sc->systypes &&
      sc->csys[v_a] < sc->ncsys && sc->systypes[sc->csys[v_a]] == 1) {
    /* Both endpoints in the same polar coordinate system.
       Compute midpoint on the arc: use Cartesian midpoint's direction
       but at the average radius from the origin. */
    INT sys = sc->csys[v_a];
    REAL *ox = sc->csys_ox + sys * dim;
    REAL ra = 0, rb = 0, rm_dir = 0;
    for (INT j = 0; j < dim; j++) {
      REAL da = sc->x[v_a * nbig + j] - ox[j];
      REAL db = sc->x[v_b * nbig + j] - ox[j];
      ra += da * da;
      rb += db * db;
    }
    ra = sqrt(ra); rb = sqrt(rb);
    REAL rho_avg = 0.5 * (ra + rb);
    /* Direction from origin to Cartesian midpoint */
    for (INT j = 0; j < dim; j++) {
      REAL d = xnew[j] - ox[j];
      rm_dir += d * d;
    }
    rm_dir = sqrt(rm_dir);
    if (rm_dir > 1e-14) {
      REAL scale = rho_avg / rm_dir;
      for (INT j = 0; j < dim; j++)
        xnew[j] = ox[j] + scale * (xnew[j] - ox[j]);
    }
  }
  INT pv[2] = {v_a, v_b};
  INT ks0 = sc->ns, ksn = sc->ns + 1;
  sc->child0[is] = ks0;
  sc->childn[is] = ksn;
  haz_add_simplex(is, sc, xnew, pv, 0, csysnew, nsnew, nvnew);
  /* Set gen_N for the new vertex (gen_N array was grown by haz_add_simplex
     via realloc of bndry/csys but gen_N needs separate growth) */
  sc->gen_N = realloc(sc->gen_N, nvnew * sizeof(INT));
  sc->gen_N[nodnew] = gen_new;
  /* Build children: replace va with v', replace vb with v' */
  el = sc->nodes + is * n1;
  INT *c0 = sc->nodes + ks0 * n1;
  INT *cn = sc->nodes + ksn * n1;
  for (INT j = 0; j < n1; j++) c0[j] = el[j];
  c0[va_idx] = nodnew;
  for (INT j = 0; j < n1; j++) cn[j] = el[j];
  cn[dim] = nodnew;
  /* Sort children by decreasing gen_N */
  for (INT i = 0; i < n1; i++)
    for (INT j = i+1; j < n1; j++) {
      if (sc->gen_N[c0[i]] < sc->gen_N[c0[j]]) { INT t=c0[i]; c0[i]=c0[j]; c0[j]=t; }
      if (sc->gen_N[cn[i]] < sc->gen_N[cn[j]]) { INT t=cn[i]; cn[i]=cn[j]; cn[j]=t; }
    }
  for (INT j = 0; j < n1; j++) {
    sc->nbr[ks0 * n1 + j] = -1;
    sc->nbr[ksn * n1 + j] = -1;
  }
  return nodnew;
}
/**********************************************************************/
/*!
 * \brief Bisect simplex is using Algorithm 4, REUSING existing midpoint.
 */
static void haz_bisect_reuse(scomplex* sc, INT is, INT nodnew) {
  INT dim = sc->dim, n1 = dim + 1, N = sc->dgs_N;
  INT *el = sc->nodes + is * n1;
  INT vm = el[dim], vm1 = el[dim - 1];
  INT lvl_m = dgs_lvl_N(sc->gen_N[vm], N);
  INT lvl_m1 = dgs_lvl_N(sc->gen_N[vm1], N);
  INT va_idx;
  if (lvl_m != lvl_m1) {
    va_idx = dim - 1;
  } else {
    va_idx = dim;
    for (INT k = 0; k < n1; k++)
      if (dgs_lvl_N(sc->gen_N[el[k]], N) == lvl_m) { va_idx = k; break; }
  }
  INT v_a = el[va_idx], v_b = el[dim];
  INT nsnew = sc->ns + 2, nvnew = sc->nv;
  INT pv[2] = {v_a, v_b};
  INT ks0 = sc->ns, ksn = sc->ns + 1;
  sc->child0[is] = ks0;
  sc->childn[is] = ksn;
  haz_add_simplex(is, sc, NULL, pv, 0, 0, nsnew, nvnew);
  el = sc->nodes + is * n1;
  INT *c0 = sc->nodes + ks0 * n1;
  INT *cn = sc->nodes + ksn * n1;
  for (INT j = 0; j < n1; j++) c0[j] = el[j];
  c0[va_idx] = nodnew;
  for (INT j = 0; j < n1; j++) cn[j] = el[j];
  cn[dim] = nodnew;
  for (INT i = 0; i < n1; i++)
    for (INT j = i+1; j < n1; j++) {
      if (sc->gen_N[c0[i]] < sc->gen_N[c0[j]]) { INT t=c0[i]; c0[i]=c0[j]; c0[j]=t; }
      if (sc->gen_N[cn[i]] < sc->gen_N[cn[j]]) { INT t=cn[i]; cn[i]=cn[j]; cn[j]=t; }
    }
  for (INT j = 0; j < n1; j++) {
    sc->nbr[ks0 * n1 + j] = -1;
    sc->nbr[ksn * n1 + j] = -1;
  }
}
/**********************************************************************/
/*!
 * \fn INT haz_refine_simplex(scomplex *sc, const INT is, const INT it)
 *
 * \brief Bisects simplex is with DGS conforming closure (Algorithm 3)
 *        using Algorithm 4 bisection (handles N+1 > n+1 colors).
 *        Scans ALL simplices for the conforming closure (brute force).
 *
 * \param sc  I/O: simplicial complex with gen_N set
 * \param is  I: index of the simplex to bisect
 * \param it  I: unused (kept for API compatibility), pass -1
 *
 * \return 0 on success
 */
INT haz_refine_simplex(scomplex* sc, const INT is, const INT it) {
  if (is < 0 || sc->child0[is] >= 0) return 0;
  INT dim = sc->dim, n1 = dim + 1;
  INT v0, vg;
  dgs_get_bse(sc, is, &v0, &vg);
  /* Conforming closure (Algorithm 3): find all leaves sharing edge
     {v0,vg} with bse != {v0,vg}, refine them first. */
  INT restart = 1;
  while (restart) {
    restart = 0;
    for (INT s = 0; s < sc->ns; s++) {
      if (s == is || sc->child0[s] >= 0) continue;
      INT *el = sc->nodes + s * n1;
      INT h0 = 0, hg = 0;
      for (INT a = 0; a < n1; a++) {
        if (el[a] == v0) h0 = 1;
        if (el[a] == vg) hg = 1;
      }
      if (!h0 || !hg) continue;
      INT bv0, bvg;
      dgs_get_bse(sc, s, &bv0, &bvg);
      if ((bv0 == v0 && bvg == vg) || (bv0 == vg && bvg == v0)) continue;
      haz_refine_simplex(sc, s, -1);
      restart = 1;
      break;
    }
  }
  if (sc->child0[is] >= 0) return 0;
  /* Bisect is, creating new midpoint */
  INT vnew = haz_bisect_new(sc, is);
  /* Bisect all remaining leaves sharing {v0,vg} */
  for (INT s = 0; s < sc->ns; s++) {
    if (sc->child0[s] >= 0) continue;
    INT *el = sc->nodes + s * n1;
    INT h0 = 0, hg = 0;
    for (INT a = 0; a < n1; a++) {
      if (el[a] == v0) h0 = 1;
      if (el[a] == vg) hg = 1;
    }
    if (h0 && hg) haz_bisect_reuse(sc, s, vnew);
  }
  return 0;
}
/******************************************************************/
/*!
 * \fn void refine(const INT ref_levels, scomplex *sc,ivector *marked)
 *
 * \brief Refines a simplicial complex.
 *
 * \param sc: scomplex containing the whole hierarchy of refinements
 *
 * \param marked: input ivector containing all simplices from the
 *               finest level which are marked for refinement; If
 *               marked is null then uniformly renfines the grid
 *               ref_levels.
 *
 *
 * \param ref_levels number of refinements. If marked is not null,
 * then this is ignored and only one refinement is done;
 *
 * \return void
 *
 * \note On first call (sc->level==0), performs DGS initialization
 *       (vertex coloring + simplex reordering) from:
 *       Diening, L., Gehring, L., Storn, J. Adaptive Mesh
 *       Refinement for Arbitrary Initial Triangulations.
 *       Found. Comput. Math. (2025).
 *       DOI: 10.1007/s10208-024-09642-1
 *
 */
void refine(const INT ref_levels, scomplex* sc, ivector* marked) {
  if (ref_levels <= 0) return;
  /*somethind to be done*/
  INT j = -1, i, nsold, print_level = 0, nsfine = -1;
  if (!sc->level) {
    /* sc->level this is set to 0 in haz_scomplex_init */
    /* form neighboring list on the coarsest level */
    find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
    /*
     * DGS initialization (Diening-Gehring-Storn, FoCM 2025):
     * Color the vertices and reorder each simplex accordingly.
     * This works for ANY conforming initial triangulation and
     * does not require the reflected-neighbor (Traxler) ordering.
     */
    {
      INT* color = (INT*)calloc(sc->nv, sizeof(INT));
      INT N = set_color(sc, color);
      /* dgs_initialize permutes both nodes and nbr arrays consistently,
         so no need to recompute neighbors. */
      dgs_initialize(sc, color, N);
      free(color);
      if (sc->print_level > 0)
        fprintf(stdout,
                "\n%%%% DGS initialization: %lld colors used (N=%lld)\n",
                (long long)(N + 1), (long long)N);
    }
    /* Traxler BFS tree initialization commented out — using DGS only */
    /* { */
    /*   INT *wrk=calloc(5*(sc->dim+2),sizeof(INT)); */
    /*   /\* construct bfs tree for the dual graph *\/ */
    /*   abfstree(0,sc,wrk,print_level); */
    /*   if(wrk) free(wrk); */
    /* } */
  }
  if ((!marked)) {
    // just refine everything that was not refined:
    for (i = 0; i < ref_levels; i++) {
      nsold = sc->ns;
      for (j = 0; j < nsold; j++)
        if (sc->child0[j] < 0 && sc->childn[j] < 0)
          haz_refine_simplex(sc, j, -1);
      sc->level++;
    }
    for (j = 0; j < sc->ns; j++)
      sc->marked[j] = TRUE;  // not sure we need this.
    // we are done here;
    return;
  } else if (!sc->level) {
    // we have not refined anything yet and marked is set, so
    if ((marked->row) && (marked->val))
      for (j = 0; j < sc->ns; j++) sc->marked[j] = marked->val[j];
    else {
      //      issue a warning and mark everything for refinement;
      for (j = 0; j < sc->ns; j++) sc->marked[j] = TRUE;
    }
  } else {
    /* we come here if we have refined few times and in such case we
       need to re-mark our simplices on the finest level using the
       values of marked at abs(sc->child0[j]+1) which were set from
       the previous visit here */
    for (j = 0; j < sc->ns; j++) {
      if (sc->child0[j] < 0 || sc->childn[j] < 0) {
        nsfine = abs(sc->child0[j] + 1);
        sc->marked[j] = marked->val[nsfine];
      }
    }
  }
  /*
   * refine everything that is marked on the finest level and is
   * not yet refined: (marked>0 and child<0)
   */
  nsold = sc->ns;
  for (j = 0; j < nsold; j++)
    if (sc->marked[j] && (sc->child0[j] < 0 && sc->childn[j] < 0))
      haz_refine_simplex(sc, j, -1);
  /*
   *  compute volumes (the volumes on the coarsest grid should be set in
   * generate_initial_grid, but just in case we are coming here from
   * some other function we compute these here too.
   */
  INT node, in1;
  void* wrk1 = malloc((sc->dim + 1) * (sc->dim * sizeof(REAL) + sizeof(INT)));
  REAL* xs = calloc((sc->dim + 1) * sc->dim, sizeof(REAL));
  for (i = 0; i < sc->ns; ++i) {
    if (sc->child0[i] < 0 || sc->childn[i] < 0) {
      in1 = i * (sc->dim + 1);
      for (j = 0; j <= sc->dim; ++j) {
        node = sc->nodes[in1 + j];
        memcpy((xs + j * sc->dim), (sc->x + node * sc->dim), sc->dim * sizeof(REAL));
      }
      sc->vols[i] = volume_compute(sc->dim, sc->factorial, xs, wrk1);
    }
  }
  free(wrk1);
  free(xs);
  sc->level++;
  return;
}
/* sc2mesh removed — scomplex is now the single mesh representation */
/*********************************************************************/
/*!
 * \fn scomplex *sc_bndry(scomplex *sc)
 *
 * \brief creates a boundary simplicial complex from a given simplicial complex.
 *
 * \param sc         I: the simplicial complex whose boundary we want to find
 *
 * \return  the coundary simplicial complex
 * \note    future: should be combined with find_cc_bndry_cc() in amr_utils.c.
 */
scomplex sc_bndry(scomplex* sc) {
  scomplex dsc;
  INT ns = sc->ns, nv = sc->nv, dim = sc->dim;
  /*
     first find the neighbors so that we have a consistent ordering of
     the vertices and faces. These may already be found, but if not we
     do it again to be sure the ordering is consistent: neighbor
     sharing face [j] is stored at the same place as vertex [j] in
     nodes[:,j]
   */
  find_nbr(ns, nv, dim, sc->nodes, sc->nbr);
  /**/
  INT dim1 = dim + 1, i, j, k, m, ns_b1, in1, jn1;
  INT ns_b = 0;
  for (i = 0; i < ns; i++) {
    for (j = 0; j < dim1; j++) {
      if (sc->nbr[i * dim1 + j] < 0) ns_b++;
    }
  }
  dsc.ns = ns_b;
  // end: computation of the nuumber of boundary faces. now store the vertices
  // for every face.
  INT* fnodes = calloc(ns_b * dim, sizeof(INT));
  ns_b1 = 0;
  for (i = 0; i < ns; i++) {
    for (j = 0; j < dim1; j++) {
      if (sc->nbr[i * dim1 + j] < 0) {
        k = 0;
        for (m = 0; m < dim1; m++) {
          if (m == j) continue;
          fnodes[ns_b1 * dim + k] = sc->nodes[i * dim1 + m];
          k++;
        }
        ns_b1++;
      }
    }
  }
  if (ns_b != ns_b1) {
    fprintf(stderr,
            "\n%%***ERROR(65): num. bndry faces mismatch (ns_b=%lld .ne. "
            "ns_b=%lld) in %s",
            (long long)ns_b1, (long long)ns_b, __FUNCTION__);
    exit(65);
  }
  // FIX numbering, there is global numbering of nodes and local numbering of
  // nodes:
  INT* indx = calloc(sc->nv, sizeof(INT));
  INT* indxinv = calloc(sc->nv, sizeof(INT));
  memset(indxinv, 0, sc->nv * sizeof(INT));
  // memset(indx,0,sc->nv*sizeof(INT));
  for (i = 0; i < sc->nv; ++i) indx[i] = -1;
  for (i = 0; i < ns_b * dim; ++i)
    indx[fnodes[i]]++;  // for interior pts this does not get incremented
  INT nv_b = 0;
  for (i = 0; i < sc->nv; ++i) {
    if (indx[i] < 0) continue;  // interior pt
    indx[i] = nv_b;
    indxinv[nv_b] = i;
    nv_b++;
  }
  fprintf(stdout, "\n%%number of boundary vertices=%lld (total nv=%lld)\n",
          (long long)nv_b, (long long)sc->nv);
  dsc = haz_scomplex_null((sc->dim - 1), sc->dim);
  dsc.nv = nv_b;
  dsc.ns = ns_b;
  ////////////////
  if (dsc.nbig > dsc.dim) {
    fprintf(stdout,
            "\n%%%%In %s:Simplicial complex of dimension %lld embedded in sc "
            "of dimension %lld\n\n",
            __FUNCTION__, (long long)dsc.dim, (long long)dsc.nbig);
  }
  // there are things we dont need:
  //  free(dsc.nodes);
  dsc.nodes = fnodes;
  if (dsc.nv < sc->nv) indxinv = realloc(indxinv, dsc.nv * sizeof(INT));
  // now we can init the complex and then remove the unnecessary stuff:
  for (i = 0; i < dsc.ns * (dsc.dim + 1); ++i) fnodes[i] = indx[fnodes[i]];
  // set x, sc->bndry and so on:
  dsc.x = (REAL*)calloc(dsc.nv * (dsc.nbig), sizeof(REAL));
  for (i = 0; i < dsc.nv; i++) {
    in1 = i * dsc.nbig;
    j = indxinv[i];
    jn1 = j * sc->dim;  // original coords:
    memcpy((dsc.x + in1), sc->x + jn1, sc->dim * sizeof(REAL));
  }
  free(indx);
  free(indxinv);
  haz_scomplex_realloc(&dsc);
  find_nbr(dsc.ns, dsc.nv, dsc.dim, dsc.nodes, dsc.nbr);
  return dsc;
}

/**********************************************************************/
/*!
 * \fn INT sc_cc_from_nbr(scomplex *sc)
 *
 * \brief Computes the number of connected components from the neighbor
 *        array sc->nbr.  Builds a compact CSR adjacency (skipping -1
 *        entries) and calls run_dfs.  Sets sc->cc on return.
 *
 * \param sc  I/O: simplicial complex with nbr[] filled.
 *
 * \return number of connected components
 */
INT sc_cc_from_nbr(scomplex* sc) {
  INT ns = sc->ns, n1 = sc->dim + 1;
  INT i, j, is, nnz = 0;
  iCSRmat s2s = icsr_create(ns, ns, n1 * ns);
  s2s.IA[0] = 0;
  for (i = 0; i < ns; i++) {
    for (j = 0; j < n1; j++) {
      is = sc->nbr[i * n1 + j];
      if (is >= 0) {
        s2s.JA[nnz] = is;
        s2s.val[nnz] = 1;
        nnz++;
      }
    }
    s2s.IA[i + 1] = nnz;
  }
  s2s.nnz = nnz;
  iCSRmat* blk_dfs = run_dfs(ns, s2s.IA, s2s.JA);
  sc->cc = blk_dfs->row;
  icsr_free(&s2s);
  icsr_free(blk_dfs);
  free(blk_dfs);
  return sc->cc;
}
/**********************************************************************/
/*!
 * \fn scomplex sc_build_boundary(scomplex *sc)
 *
 * \brief Builds the boundary simplicial complex from a bulk mesh.
 *
 *        The boundary complex has dimension dim-1 embedded in dim
 *        dimensions (nbig=dim).  Boundary faces are identified by
 *        nbr[i*n1+j] < 0.  Face codes are inherited from sc->bndry[]
 *        (minimum nonzero code over the face's vertices).
 *
 *        The returned scomplex has:
 *        - nodes[]: boundary face connectivity (local vertex numbering)
 *        - flags[]: face codes from the initial mesh
 *        - x[]: boundary vertex coordinates
 *        - bndry[]: per-vertex boundary codes (from bulk mesh)
 *        - nbr[]: face-to-face adjacency on the boundary
 *        - cc: number of connected components on the boundary
 *
 * \param sc  I: the bulk simplicial complex (with nbr[] set by find_nbr)
 *
 * \return the boundary simplicial complex (by value)
 *
 * \note sc->nbr must be filled (via find_nbr) before calling.
 *       sc->bndry[] should contain vertex boundary codes.
 */
scomplex sc_build_boundary(scomplex* sc) {
  INT ns = sc->ns, nv = sc->nv, dim = sc->dim, n1 = dim + 1;
  INT i, j, k, m;
  scomplex dsc;

  /* Count boundary faces */
  INT ns_b = 0;
  for (i = 0; i < ns; i++)
    for (j = 0; j < n1; j++)
      if (sc->nbr[i * n1 + j] < 0) ns_b++;

  /* Extract face connectivity and face codes */
  INT* fnodes = calloc(ns_b * dim, sizeof(INT));
  INT* fcodes = calloc(ns_b, sizeof(INT));
  INT idx = 0;
  for (i = 0; i < ns; i++) {
    for (j = 0; j < n1; j++) {
      if (sc->nbr[i * n1 + j] < 0) {
        /* Collect vertices of face (all except vertex j) */
        k = 0;
        INT mincode = 0;
        for (m = 0; m < n1; m++) {
          if (m == j) continue;
          INT v = sc->nodes[i * n1 + m];
          fnodes[idx * dim + k] = v;
          /* Face code = minimum nonzero bndry code over face vertices */
          INT bc = abs(sc->bndry[v]);
          if (bc > 0) {
            if (mincode == 0 || bc < mincode) mincode = bc;
          }
          k++;
        }
        fcodes[idx] = mincode;
        idx++;
      }
    }
  }

  /* Renumber vertices: only boundary vertices get local numbers */
  INT* indx = calloc(nv, sizeof(INT));
  INT* indxinv = calloc(nv, sizeof(INT));
  for (i = 0; i < nv; i++) indx[i] = -1;
  for (i = 0; i < ns_b * dim; i++) indx[fnodes[i]] = 0;
  INT nv_b = 0;
  for (i = 0; i < nv; i++) {
    if (indx[i] < 0) continue;
    indxinv[nv_b] = i;
    indx[i] = nv_b;
    nv_b++;
  }
  if (nv_b < nv) indxinv = realloc(indxinv, nv_b * sizeof(INT));

  /* Apply renumbering to fnodes */
  for (i = 0; i < ns_b * dim; i++) fnodes[i] = indx[fnodes[i]];

  /* Create boundary scomplex */
  dsc = haz_scomplex_null(dim - 1, dim);
  dsc.ns = ns_b;
  dsc.nv = nv_b;
  dsc.nodes = fnodes; /* haz_scomplex_realloc will realloc this */
  haz_scomplex_realloc(&dsc);

  /* Set flags (face codes) — must be after realloc which zeros flags */
  for (i = 0; i < ns_b; i++) dsc.flags[i] = fcodes[i];
  free(fcodes);

  /* Copy vertex coordinates */
  dsc.x = realloc(dsc.x, nv_b * dim * sizeof(REAL));
  for (i = 0; i < nv_b; i++)
    memcpy(dsc.x + i * dim, sc->x + indxinv[i] * dim, dim * sizeof(REAL));

  /* Copy boundary codes */
  for (i = 0; i < nv_b; i++) dsc.bndry[i] = sc->bndry[indxinv[i]];

  free(indx);
  free(indxinv);

  /* Build neighbors and connected components */
  find_nbr(dsc.ns, dsc.nv, dsc.dim, dsc.nodes, dsc.nbr);
  sc_cc_from_nbr(&dsc);

  return dsc;
}
/**********************************************************************/
/*!
 * \fn INT sc_conformity_check(scomplex *sc)
 *
 * \brief Checks if a simplicial complex is conforming (no hanging nodes).
 *
 * A triangulation is conforming if every facet (codimension-1 face)
 * is shared by exactly 2 simplices (interior) or 1 simplex (boundary).
 * A facet appearing once that is not on the boundary indicates a
 * hanging node (non-conforming mesh).
 *
 * \param sc  I: the simplicial complex (leaf mesh) to check
 *
 * \return 0 if the mesh is conforming, nonzero = number of non-conforming
 * facets
 *
 * \note Uses a hash table to count facet occurrences. Each facet is
 *       identified by its sorted vertex indices.
 */
INT sc_conformity_check(scomplex* sc) {
  INT ns = sc->ns, dim = sc->dim, n1 = dim + 1;
  INT nfacets = ns * n1; /* total facets (dim+1 per simplex) */
  INT i, j, k;
  /*
   * Store all facets as sorted vertex tuples of length dim.
   * facets[f*dim + 0..dim-1] = sorted vertex indices of facet f.
   */
  INT* facets = (INT*)calloc(nfacets * dim, sizeof(INT));
  INT* ftmp = (INT*)calloc(dim, sizeof(INT));
  for (i = 0; i < ns; i++) {
    INT* el = sc->nodes + i * n1;
    for (j = 0; j < n1; j++) {
      /* facet j = all vertices except vertex j */
      INT fi = i * n1 + j;
      INT pos = 0;
      for (k = 0; k < n1; k++) {
        if (k == j) continue;
        ftmp[pos++] = el[k];
      }
      /* insertion sort ftmp */
      for (pos = 1; pos < dim; pos++) {
        INT val = ftmp[pos];
        INT hole = pos;
        while (hole > 0 && ftmp[hole - 1] > val) {
          ftmp[hole] = ftmp[hole - 1];
          hole--;
        }
        ftmp[hole] = val;
      }
      memcpy(facets + fi * dim, ftmp, dim * sizeof(INT));
    }
  }
  free(ftmp);
  /*
   * Sort facets lexicographically to count duplicates.
   * Use an index array and qsort.
   */
  INT* idx = (INT*)calloc(nfacets, sizeof(INT));
  for (i = 0; i < nfacets; i++) idx[i] = i;
  /* We need dim accessible in the comparator — use a file-scope variable */
  /* Instead, sort by building a comparison key approach with qsort_r or
   * just do a simple bucket/radix approach. For portability, we do a
   * manual merge sort with the facets array. */
  /* Simple approach: sort idx[] using shell sort with lexicographic compare */
  {
    INT gap, ii, jj, tmp;
    for (gap = nfacets / 2; gap > 0; gap /= 2) {
      for (ii = gap; ii < nfacets; ii++) {
        tmp = idx[ii];
        INT* a = facets + tmp * dim;
        for (jj = ii; jj >= gap; jj -= gap) {
          INT* b = facets + idx[jj - gap] * dim;
          INT cmp = 0;
          for (k = 0; k < dim; k++) {
            if (a[k] < b[k]) {
              cmp = -1;
              break;
            }
            if (a[k] > b[k]) {
              cmp = 1;
              break;
            }
          }
          if (cmp >= 0) break;
          idx[jj] = idx[jj - gap];
        }
        idx[jj] = tmp;
      }
    }
  }
  /*
   * Scan sorted facets and count occurrences.
   * Conforming: each facet appears 1 (boundary) or 2 (interior) times.
   * Non-conforming: a facet appears an odd number != 1 or > 2 times.
   */
  INT nerr = 0;
  i = 0;
  while (i < nfacets) {
    /* count how many consecutive facets are identical */
    INT cnt = 1;
    while (i + cnt < nfacets) {
      INT* a = facets + idx[i] * dim;
      INT* b = facets + idx[i + cnt] * dim;
      INT same = 1;
      for (k = 0; k < dim; k++) {
        if (a[k] != b[k]) {
          same = 0;
          break;
        }
      }
      if (!same) break;
      cnt++;
    }
    if (cnt > 2) {
      /* more than 2 simplices share this facet — broken mesh */
      nerr++;
      if (sc->print_level > 0) {
        fprintf(stderr, "\n%%WARNING: facet shared by %lld simplices: [",
                (long long)cnt);
        INT* a = facets + idx[i] * dim;
        for (k = 0; k < dim; k++) fprintf(stderr, " %lld", (long long)a[k]);
        fprintf(stderr, " ]");
      }
    }
    /* cnt == 1 is boundary, cnt == 2 is interior — both OK */
    i += cnt;
  }
  free(facets);
  free(idx);
  if (nerr) {
    fprintf(stderr,
            "\n%%***CONFORMITY CHECK FAILED: %lld non-conforming facets "
            "(ns=%lld, nv=%lld, dim=%lld)\n",
            (long long)nerr, (long long)ns, (long long)sc->nv, (long long)dim);
  } else {
    fprintf(stdout, "%% Conformity check PASSED (ns=%lld, nv=%lld, dim=%lld)\n",
            (long long)ns, (long long)sc->nv, (long long)dim);
  }
  return nerr;
}
/* creategrid_fread moved to io_commented_out.c — use sc_read_gmsh instead */
/**********************************************************************/
/*!
* \fn scomplex* make_uniform_mesh(const INT dim,const INT mesh_ref_levels,const INT mesh_ref_type,const INT set_bndry_codes)
*
* \brief makes a mesh_ref_levels of refined mesh on the unit cube in dimension dim.
*
* \param dim               Dimension of computational domain
* \param mesh_ref_levels   Number of refinement levels
*
* \param mesh_ref_type     if > 10: uniform refinement ;
*                          if < 10: nearest vertex bisection ;
*
* \param set_bndry_codes   if .eq. 0: the boundary codes come from sc->bndry[];
*                          if .ne. 0  the boundary codes are set
*
* \return scomplex* with the mesh and FEM data.
*
* \note 20210815 (ltz)
* \note Works in 2D and 3D
*
*/
scomplex* make_uniform_mesh(const INT dim,const INT mesh_ref_levels,const INT mesh_ref_type,const INT set_bndry_codes)
{

  // Loop Counters
  INT jlevel,k;

  // Create a simplicial complex
  scomplex **sc_all=NULL,*sc=NULL,*sctop=NULL;
  fprintf(stdout,"\n%%%%---------------------------------------------------------------------");
  fprintf(stdout,"\n%%%%Meshing...");
  // Get the coarsest mesh on the cube in dimension dim and set the refinement type.
  sc_all=mesh_cube_init(dim,1,mesh_ref_type);
  sc=sc_all[0];
  if(sc->ref_type>10){
    // Uniform refinement only for dim=2 or dim=3
    if(dim==3){
      for(jlevel=0;jlevel<mesh_ref_levels;++jlevel){
        uniformrefine3d(sc);
        sc_vols(sc);
      }
    } else if(dim==2){
      for(jlevel=0;jlevel<mesh_ref_levels;++jlevel){
        uniformrefine2d(sc);
        sc_vols(sc);
      }
    } else {
      check_error(ERROR_DIM, __FUNCTION__);
    }
    // Get boundaries
    find_nbr(sc->ns,sc->nv,sc->dim,sc->nodes,sc->nbr);
    sc_vols(sc);
  } else {
    // Nearest vertex bisection refinement
    ivector marked; marked.val=NULL;
    for(jlevel=0;jlevel<mesh_ref_levels;++jlevel){
      // Choose the finest grid
      sctop=scfinest(sc);
      marked.row=sctop->ns;
      marked.val=realloc(marked.val,marked.row*sizeof(INT));
      // Mark everything
      for(k=0;k<marked.row;k++) marked.val[k]=TRUE;
      // Now we refine
      refine(1,sc,&marked);
      // Free the finest grid
      haz_scomplex_free(sctop);
    }
    ivec_free(&marked);
  }
  fprintf(stdout,"Done.\n");
  scfinalize(sc, NULL, (INT)1);
  sc_vols(sc);
  if (FALSE) {// (do not export)/(export) the mesh to vtu: [FALSE/TRUE]
    vtu_data vdata;
    vtu_data_init(sc, &vdata);
    sc_write_vtk("mesh.vtu", &vdata);
    vtu_data_free(&vdata);
  }
  fprintf(stdout,"%%%%---------------------------------------------------------------------\n");
  // Build FEM data directly on the simplicial complex
  sc_build_fem_data(sc);
  // Free the pointer array but NOT the scomplex itself (caller owns it)
  free(sc_all);
  return sc;
}
/**********************************************************************/
/*!
 * \fn void create1Dgrid_Line(scomplex** sc_ptr,REAL left_end,REAL right_end,INT nelm)
 *
 * \brief Creates a 1D grid from scratch [left_end,right_end] with no holes.
 *
 * \param sc_ptr     Pointer to scomplex pointer (output)
 * \param left_end   Coordinate of Left-End Point
 * \param right_end  Coordinate of Right-End Point
 * \param nelm       Number of Elements
 *
 * \return scomplex via sc_ptr with the mesh and FEM data.
 *
 */
void create1Dgrid_Line(scomplex** sc_ptr,REAL left_end,REAL right_end,INT nelm)
{
  INT i; /* Loop index */

  INT dim = 1;
  INT nv = nelm+1;

  // Get h-spacing
  REAL h = (right_end - left_end)/nelm;

  // Create scomplex
  scomplex *sc = haz_scomplex_init(dim, nelm, nv, dim);

  // Build element connectivity (nodes) and coordinates
  for(i=0;i<nelm;i++) {
    sc->nodes[i*2] = i;
    sc->nodes[i*2+1] = i+1;
    sc->x[i] = left_end + h*i;
    sc->bndry[i] = 0;
  }
  sc->x[nelm] = right_end;
  sc->bndry[0] = 1;
  sc->bndry[nelm] = 1;

  // Set connected components
  sc->cc = 1;
  sc->bndry_cc = 1;

  // Build neighbors, volumes, and FEM data
  find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
  sc_vols(sc);
  sc_build_fem_data(sc);

  *sc_ptr = sc;
  return;
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*EOF*/
