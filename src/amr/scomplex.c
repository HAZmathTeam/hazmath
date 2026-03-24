/*! \file src/amr/scomplex.c
 *
 *  Authors: James Adler, Xiaozhe Hu, and Ludmil Zikatanov
 *           HAZmath (https://hazmath.net)
 *           Created with the help of Claude (Anthropic)
 *
 *  Created 20170715.  Copyright 2017__HAZMATH__. All rights reserved.
 *
 *  \note Simplicial complex: initialization, geometry, volumes,
 *        FEM data, boundary, conformity check, neighbor finding.
 *        Refinement routines moved to amr_core.c (20260322).
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
  /* Free all iCSRmat pointers via fem. el_v->JA aliases sc->nodes
     (freed by haz_scomplex_free), so null it before icsr_free. */
  if (fem->el_v) { fem->el_v->JA = NULL; icsr_free(fem->el_v); free(fem->el_v); fem->el_v = NULL; }
  if (fem->el_ed) { icsr_free(fem->el_ed); free(fem->el_ed); fem->el_ed = NULL; }
  if (fem->ed_v) { icsr_free(fem->ed_v); free(fem->ed_v); fem->ed_v = NULL; }
  if (fem->el_f && fem->el_f != fem->el_ed) { icsr_free(fem->el_f); free(fem->el_f); }
  fem->el_f = NULL;
  if (fem->f_v && fem->f_v != fem->ed_v) { icsr_free(fem->f_v); free(fem->f_v); }
  fem->f_v = NULL;
  if (fem->f_ed) { icsr_free(fem->f_ed); free(fem->f_ed); fem->f_ed = NULL; }
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
  if (fem->coded_faces) free(fem->coded_faces);
  if (fem->coded_f_btype) free(fem->coded_f_btype);
  if (fem->dwork) free(fem->dwork);
  if (fem->iwork) free(fem->iwork);
  free(fem);
  sc->fem = NULL;
  /* Free inc pointer array only — entries were freed via fem above */
  if (sc->inc) { free(sc->inc); sc->inc = NULL; }
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
    /* Compute number of faces by counting from the neighbor array.
       Each interior face appears twice in nbr (once per element),
       each boundary face once. */
    INT nface = 0;
    if (dim == 2) {
      nface = nedge;
    } else {
      INT n_bndry_slots = 0;
      for (INT ii = 0; ii < ns_leaf * n1; ii++) {
        if (sc->nbr[ii] < 0) n_bndry_slots++;
      }
      nface = (ns_leaf * n1 - n_bndry_slots) / 2 + n_bndry_slots;
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
    /* Overwrite f_flag with macroface codes from bndry_f2v.
       For each boundary face in bndry_f2v, find the matching global face
       in f_v by vertex-set comparison, and set f_flag to the face code. */
    if (sc->bndry_f2v && sc->bndry_f2v->row > 0) {
      iCSRmat *fv = fem->f_v;
      for (INT bf = 0; bf < sc->bndry_f2v->row; ++bf) {
        INT ba = sc->bndry_f2v->IA[bf], bb = sc->bndry_f2v->IA[bf + 1];
        INT nvf = bb - ba;
        INT fcode = sc->bndry_f2v->val[ba];
        if (!fcode || nvf != dim) continue;
        /* Find matching face in f_v */
        for (INT gf = 0; gf < nface; ++gf) {
          INT ga = fv->IA[gf], gb = fv->IA[gf + 1];
          if (gb - ga != nvf) continue;
          INT match = 1;
          for (INT p = ba; p < bb && match; ++p) {
            INT found = 0;
            for (INT q = ga; q < gb; ++q)
              if (sc->bndry_f2v->JA[p] == fv->JA[q]) { found = 1; break; }
            if (!found) match = 0;
          }
          if (match) { f_flag[gf] = fcode; break; }
        }
      }
    }
    /* Build coded-faces list: faces with nonzero f_flag.
       Also determine boundary/interior using el_f transpose (face degree). */
    {
      iCSRmat f_el;
      icsr_trans(fem->el_f, &f_el);
      INT ncf = 0;
      for (INT f = 0; f < nface; ++f) if (f_flag[f]) ncf++;
      fem->n_coded_faces = ncf;
      fem->coded_faces = (INT *)malloc(ncf * sizeof(INT));
      fem->coded_f_btype = (INT *)malloc(ncf * sizeof(INT));
      INT ci = 0;
      for (INT f = 0; f < nface; ++f) {
        if (!f_flag[f]) continue;
        fem->coded_faces[ci] = f;
        /* boundary if face has only 1 adjacent element */
        INT deg = f_el.IA[f + 1] - f_el.IA[f];
        fem->coded_f_btype[ci] = (deg >= 2) ? 1 : 0;
        ci++;
      }
      icsr_free(&f_el);
    }
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
/* haz_add_simplex, set_color, dgs_initialize, v2s_build, v2s_free_data,
   dgs_type_N, dgs_lvl_N, dgs_get_bse, haz_bisect_new, haz_bisect_reuse,
   v2s_add, haz_refine_simplex, refine  -->  moved to amr_core.c */
/* make_uniform_mesh  -->  moved to amr_core.c */
/*********************************************************************/
/* REMOVED BLOCK: lines from haz_add_simplex body through refine() -- see amr_core.c */
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
  // end: computation of the number of boundary faces. now store the vertices
  // for every face, and look up face codes from bndry_f2v.
  INT* fnodes = calloc(ns_b * dim, sizeof(INT));
  INT* fcodes = calloc(ns_b, sizeof(INT)); /* face code per boundary face */
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
        /* Look up face code from bndry_f2v by matching vertex set */
        fcodes[ns_b1] = 0;
        if (sc->bndry_f2v) {
          INT *fv = fnodes + ns_b1 * dim;
          for (INT bf = 0; bf < sc->bndry_f2v->row; ++bf) {
            INT ba = sc->bndry_f2v->IA[bf], bb = sc->bndry_f2v->IA[bf + 1];
            if (bb - ba != dim) continue;
            /* check if all vertices match (order-independent) */
            INT match = 1;
            for (INT p = 0; p < dim && match; ++p) {
              INT found = 0;
              for (INT q = ba; q < bb; ++q) {
                if (fv[p] == sc->bndry_f2v->JA[q]) { found = 1; break; }
              }
              if (!found) match = 0;
            }
            if (match) { fcodes[ns_b1] = sc->bndry_f2v->val[ba]; break; }
          }
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
  /* Set boundary face codes as element flags */
  for (i = 0; i < ns_b; i++) dsc.flags[i] = fcodes[i];
  free(fcodes);
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
