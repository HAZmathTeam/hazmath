/*! \file examples/amr_grids/amr_grids_supporting.h
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program
 *
 * \note
 */
/************************************************************/
static char *str_add_dim(const INT dim, const char *prefix, const char *suffix)
{
  INT k=dim;
  if(k<1) k=1;
  else if(k>3) k=3;  
  char *dir0=calloc(1023,sizeof(char));
  snprintf(dir0,1023,"%s%d%s",prefix,k,suffix);
  trim_str(&dir0,1);
  return dir0;
}
static char *fname_set(const char *dir, const char *fname_in) {
  // combine names: fname_in[0]=dir/fname_in[0]
  size_t ldir0 = strlen(dir) + 1;
  size_t lfname = strlen(fname_in);
  char *fname = strndup(dir, ldir0);
  fname = realloc(fname, (lfname + ldir0 + 1) * sizeof(char));
  strncat(fname, fname_in, lfname);
  trim_str(&fname,1);
  return fname;
}
/************************************************************/
static ivector mark_around_pts(scomplex *sc, scomplex *scglobal, INT nstar,
                               REAL *xstar, iCSRmat *node_ins,
                               const INT max_nodes) {
  /* scglobal is the global sc that contains all refinements */
  iCSRmat ins_node;
  dCSRmat xins_node, xnode_ins;
  INT *ii, *jj;
  INT dim = scglobal->n, dim1 = dim + 1, n = sc->n, n1 = n + 1, ns = sc->ns;
  INT nzw = -10, jstar, jk, p, pj, j = -1, iwa, iwb;
  /* un-mark everything */
  ivector marked =
      ivec_create(sc->ns);  // this is only created on the top level....
  for (j = 0; j < ns; j++) marked.val[j] = FALSE;
  //  for(j=0;j<ns;j++) marked.val[j]=TRUE;
  //  bail out if no points
  if (!nstar || (xstar == NULL)) return marked;
  INT *scnjn = NULL;
  REAL *xstar0 = NULL;
  //  fprintf(stdout,"\nlevel=%d",scglobal->level);
  if (!(scglobal->level)) {
    // on the coarsest level we need to set the matrix as for a first time:
    node_ins->row = nstar;
    node_ins->col = scglobal->ns;
    node_ins->IA = calloc((nstar + 1), sizeof(INT));
    nzw = 0;
    node_ins->IA[0] = nzw;
    for (jstar = 0; jstar < nstar; ++jstar) {
      xstar0 = xstar + jstar * dim;
      for (j = 0; j < scglobal->ns; j++) {
        scnjn = scglobal->nodes + j * dim1;
        // check if the point is inside the simplex.
        if (!xins(dim, scnjn, scglobal->x, xstar0)) {
          nzw++;
        }
      }
      node_ins->IA[jstar + 1] = nzw;
    }
    node_ins->nnz = node_ins->IA[nstar];
    // now run again, this time filling in the points:
    node_ins->JA = calloc(node_ins->nnz, sizeof(INT));
    nzw = 0;
    for (jstar = 0; jstar < nstar; ++jstar) {
      xstar0 = xstar + jstar * dim;
      for (j = 0; j < scglobal->ns; j++) {
        scnjn = scglobal->nodes + j * dim1;
        // check if the point is inside the simplex.
        if (!xins(dim, scnjn, scglobal->x, xstar0)) {
          node_ins->JA[nzw] = j;
          marked.val[j] = TRUE;
          nzw++;
        }
      }
    }
    //
  } else {
    //    haz_scomplex_print(scglobal,0,__FUNCTION__);
    ins_node = icsr_create(node_ins->col, node_ins->row, node_ins->nnz);
    xins_node = dcsr_create(0, 0, 0);
    xins_node.row = ins_node.row;
    xins_node.col = ins_node.col;
    xins_node.nnz = ins_node.nnz;
    xins_node.IA = ins_node.IA;
    xins_node.JA = ins_node.JA;
    xins_node.val = NULL;
    // input:
    xnode_ins = dcsr_create(0, 0, 0);
    xnode_ins.row = node_ins->row;
    xnode_ins.col = node_ins->col;
    xnode_ins.nnz = node_ins->nnz;
    xnode_ins.IA = node_ins->IA;
    xnode_ins.JA = node_ins->JA;
    xnode_ins.val = NULL;
    // transpose
    dcsr_transz(&xnode_ins, NULL, &xins_node);
    //    fprintf(stdout,"\nNS=%d",scglobal->ns);fflush(stdout);
    //    icsr_write_icoo("ins_node.txt",&ins_node);
    /* INT nzw=0; */
    INT pjiter;
    nzw = 0;
    // first run compute nonzeroes:
    for (j = 0; j < scglobal->ns; j++) {
      // we do the finest level only!
      if (scglobal->child0[j] < 0 || scglobal->childn[j] < 0) {
        pj = scglobal->parent[j];
        p = pj;
        if (pj < 0) pj = j;
        //	fprintf(stdout,"\nj=%d,p=%d;pj=%d,ins_node.row=%d",j,p,pj,ins_node.row);fflush(stdout);
        pjiter = 0;
        while (pj >= ins_node.row) {
          pj = scglobal->parent[pj];
          pjiter++;
          if (pjiter > scglobal->level) {
            // stop w err if this is not working.
            fprintf(stderr,
                    "%%%%****ERROR (%s): ancestor not found in refinement tree "
                    "for element %lld on level=%lld\n%%%%****Exiting....",
                    __FUNCTION__, (long long int)j,
                    (long long int)scglobal->level);
            exit(16);
          }
        }
        /* if(pjiter>0){ */
        /*   fprintf(stdout,"\nNew:pj[%d]==%d,pjiter=%d",j,pj,pjiter);fflush(stdout);
         */
        /* } */
        iwa = ins_node.IA[pj];
        iwb = ins_node.IA[pj + 1];
        if ((iwb - iwa) <= 0) {
          // empty simplex do nothing
          continue;
        }
        scnjn = scglobal->nodes + j * n1;
        for (jk = iwa; jk < iwb; jk++) {
          jstar = ins_node.JA[jk];
          xstar0 = xstar + jstar * n;
          // check if the point is inside the simplex.
          if (!xins(n, scnjn, scglobal->x, xstar0)) {
            // this is in, so we count it.
            nzw++;
          } else {
            // fprintf(stdout, "\nchk_node=%d", jstar);
            // print_full_mat(1, n, xstar0, "xstar0");
          }
        }
      }
    }
    // reallocate as some nodes may belong to more than 1 simplex... otherwise
    // this can be simplified.
    if (nzw > ins_node.nnz) {
      ins_node.val = realloc(ins_node.val, nzw * sizeof(INT));
      ins_node.JA = realloc(ins_node.JA, nzw * sizeof(INT));
    }
    ins_node.nnz = nzw;
    if (nzw > node_ins->nnz) {
      node_ins->JA = realloc(node_ins->JA, nzw * sizeof(INT));
    }
    node_ins->nnz = nzw;
    for (j = 0; j < ins_node.nnz; ++j) {
      ins_node.val[j] = -1;
      node_ins->JA[j] = -1;
    }
    // now run again to get the JA in val.
    jj = node_ins->JA;
    ii = ins_node.val;
    nzw = 0;
    for (j = 0; j < scglobal->ns; j++) {
      /* // check if we have marked and do nothing if we are marked: */
      /* if(marked.val[j]) {continue;} */
      // we do the finest level only!
      if (scglobal->child0[j] < 0 || scglobal->childn[j] < 0) {
        pj = scglobal->parent[j];
        if (pj < 0) pj = j;
        pjiter = 0;
        while (pj >= ins_node.row) {
          pj = scglobal->parent[pj];
          pjiter++;
          if (pjiter > scglobal->level) {
            // stop w err if this is not working.
            fprintf(stderr,
                    "%%%%****ERROR (%s): ancestor not found in refinement tree "
                    "for element %lld on level=%lld\n%%%%****Exiting....",
                    __FUNCTION__, (long long int)j,
                    (long long int)scglobal->level);
            exit(16);
          }
        }
        iwa = ins_node.IA[pj];
        iwb = ins_node.IA[pj + 1];
        if ((iwb - iwa) <= 0) {
          // empty simplex do nothing
          continue;
        }
        scnjn = scglobal->nodes + j * n1;
        for (jk = iwa; jk < iwb; jk++) {
          jstar = ins_node.JA[jk];
          xstar0 = xstar + jstar * n;
          // check if the point is inside the simplex.
          if (!xins(n, scnjn, scglobal->x, xstar0)) {
            // this is in, so we add it.
            ii[nzw] = j;      // note: ii is alias for ins_node.val
            jj[nzw] = jstar;  // note: jj is alias for node_ins->JA
                              // (node_ins->val does not exist).
            nzw++;
          }else{
            // fprintf(stdout, "\nchk_node2=%d", jstar);
            // print_full_mat(1, n, xstar0, "xstar0_2");
          }
        }  // end-for points in j-th simplex
      }    // if on finest level
    }      // endfor over all simplices
    // now we need to convert and clean up (use dcoo_2_dcsr for our particular
    // case with tmp->val used as jj.
    int i, iind, jind;
    ins_node.row = scglobal->ns;
    ins_node.col = nstar;
    ins_node.nnz = nzw;
    INT *ind = calloc((ins_node.row + 1), sizeof(INT));
    ins_node.IA = realloc(ins_node.IA, (ins_node.row + 1) * sizeof(INT));
    // initins_node.IAlize
    memset(ind, 0, sizeof(INT) * (ins_node.row + 1));
    // count number of nonzeros in each row
    for (i = 0; i < nzw; ++i) ind[ii[i] + 1]++;
    // set row pointer
    ins_node.IA[0] = 0;
    for (i = 1; i <= ins_node.row; ++i) {
      ins_node.IA[i] = ins_node.IA[i - 1] + ind[i];
      ind[i] = ins_node.IA[i];
    }

    // set column index and values
    for (i = 0; i < nzw; ++i) {
      iind = ii[i];
      jind = ind[iind];
      ins_node.JA[jind] = jj[i];
      ind[iind] = ++jind;
    }
    free(ind);
    // now all should be ready we transpose and we are done:
    ins_node.nnz = nzw;
    xins_node.row = ins_node.row;
    xins_node.col = xins_node.col;
    xins_node.nnz = ins_node.nnz;
    node_ins->col = ins_node.row;  // transpose should also have this,
    node_ins->row = ins_node.col;
    node_ins->nnz = ins_node.nnz;
    // here we need to assign the right pointers:
    xins_node.IA = ins_node.IA;
    xins_node.JA = ins_node.JA;
    xnode_ins.IA = node_ins->IA;
    xnode_ins.JA = node_ins->JA;
    dcsr_transz(&xins_node, NULL, &xnode_ins);
    //    icsr_write_icoo("ins_node.txt",&ins_node);
    //    icsr_write_icoo("node_ins.txt",node_ins);
    for (j = 0; j < sc->ns; j++) {
      marked.val[j] = FALSE;
    }
    nzw = 0;  // working
    for (j = 0; j < scglobal->ns; j++) {
      if (scglobal->child0[j] < 0 || scglobal->childn[j] < 0) {
        nzw++;
        if ((ins_node.IA[j + 1] - ins_node.IA[j]) > max_nodes) {
          if (scglobal->level > 100) {
            fprintf(stdout, "\ntet=%d; num_nodes=%d;nodes=[", j,ins_node.IA[j + 1] - ins_node.IA[j]);
            for (i = ins_node.IA[j]; i < ins_node.IA[j + 1]; ++i) {
              fprintf(stdout, "%d", ins_node.JA[i]);
            }
            fprintf(stdout, "];\n");
          }
          p = abs((scglobal->child0[j] + 1));
          //	  fprintf(stdout,"\np=%d",p);fflush(stdout);
          marked.val[p] = TRUE;
        }
      }
    }
    icsr_free(&ins_node);
  }
  /* if(scglobal->level > 256){ */
  /*   fprintf(stdout,"\nnode_ins(row)=%d,node_ins(col)=%d,node_ins(nnz)=%d\n",node_ins->row,node_ins->col,node_ins->nnz);fflush(stdout);
   */
  /*   icsr_write_icoo("z123.txt",node_ins); */
  /* } */
  return marked;
}
static void data_1d_init(const INT dimbig,const char *idir, const char *odir,data_1d *g)
{
  const char *fnames[] = {
    "001_coordinates_of_vertices",
    "002_segments_definition",
    "003_number_of_points_by_segment",
    "004_coordinates_of_all_inner_points_:including_vertices:",
    "005_point_thickness",
    "009_thickness",
    "\0"};
  //  char
  //////////////////////////////////////////////////////////////////////////////
  g->odir=strdup(odir);
  g->idir=strdup(idir);
  //////////////////////////////////////////////////////////////////////////////
  g->dim = 1;
  g->dimbig = dimbig;
  if(g->dimbig<2) g->dimbig=2;
  else if(g->dimbig>3) g->dimbig=3;  
  g->nv = 0;
  g->nseg = 0;
  g->xv = NULL;
  g->seg_radius = NULL,
  g->pt_thickness = NULL;
  g->seg = NULL, g->divisions = NULL;
  // filenames:
  g->fv_coords = fname_set(idir, fnames[0]);
  g->fseg = fname_set(idir, fnames[1]);
  g->fdivisions = fname_set(g->idir, fnames[2]);
  g->fvtmp_coords = fname_set(g->idir, fnames[3]);
  g->fpt_thickness = fname_set(g->idir, fnames[4]);
  g->fseg_radius = fname_set(g->idir, fnames[5]);
  //
  g->fvtu_3d = str_add_dim(g->dimbig,g->odir,"d_grid.vtu");
  g->fvtu_1d = str_add_dim(g->dim,g->odir,"d_grid.vtu");
  ////////////////////////////////////////////////////////////////////
  fprintf(stdout,"\n%%%%INPUT FILENAMES=\n%s\n%s\n%s\n%s\n%s\n%s\n",\
	  g->fv_coords,g->fseg,g->fdivisions,g->fvtmp_coords,g->fpt_thickness,g->fseg_radius);
  fprintf(stdout,"\n%%%%OUTPUT FILENAMES:\n%s\n%s\n",g->fvtu_3d,g->fvtu_1d);
  //
}
/************************************************************/
/* void getdata_1d(data_1d *g)
 * read 1d data
 *
 *
 *
 */
/*************************************************************/
static void getdata_1d(data_1d *g) {
  FILE *fp;
  //
  REAL *xvtmp = NULL;  // working
  INT nvsize, nvaddsize, nsegsize;
  read_eof(g->fv_coords, (void **)&g->xv, &nvsize, 'R', 0);
  free(g->fv_coords);g->fv_coords=NULL;
  read_eof(g->fvtmp_coords, (void **)&xvtmp, &nvaddsize, 'R', 0);
  free(g->fvtmp_coords);g->fvtmp_coords=NULL;
  //
  g->nv = (INT)nvsize / g->dimbig;
  g->nvadd = (INT)nvaddsize / g->dimbig;
  //
  fprintf(stdout, "\nnvsize=%lld,nvaddsize=%lld\n", (long long)nvsize,
          (long long)nvaddsize);
  g->xv = realloc(g->xv, (nvsize + nvaddsize) * sizeof(REAL));
  memcpy((g->xv + nvsize), xvtmp, nvaddsize * sizeof(REAL));
  free(xvtmp);
  xvtmp = (g->xv + nvsize);
  /// integers
  read_eof(g->fseg, (void **)&g->seg, &nsegsize, 'I', 0);
  free(g->fseg);g->fseg=NULL;
  //
  g->nseg = nsegsize / (g->dim + 1);
  // for the rest we know the sizes so we alocate and then read:
  g->divisions = calloc(g->nseg, sizeof(INT));
  fp = fopen(g->fdivisions, "r");
  rveci_(fp, g->divisions, &g->nseg);
  fclose(fp);
  free(g->fdivisions);g->fdivisions=NULL;
  //
  g->seg_radius = calloc(g->nseg, sizeof(REAL));
  fp = fopen(g->fseg_radius, "r");
  rvecd_(fp, g->seg_radius, &g->nseg);
  fclose(fp);
  free(g->fseg_radius);g->fseg_radius=NULL;
  //
  g->pt_thickness = calloc(g->nvadd, sizeof(REAL));
  fp = fopen(g->fpt_thickness, "r");
  rvecd_(fp, g->pt_thickness, &g->nvadd);
  fclose(fp);
  free(g->fpt_thickness);g->fpt_thickness=NULL;
  //  fprintf(stdout, "\nnseg=%d\n", g->nseg);
  return;
}
/********************************************************************************/
static void data_1d_free(data_1d *g)
{  
  free(g->divisions);
  free(g->seg);
  free(g->pt_thickness);
  free(g->seg_radius);
  /// filenames if not freed:
  free(g->fv_coords);
  free(g->fvtmp_coords);
  free(g->fseg);
  free(g->fdivisions);
  free(g->fseg_radius);
  free(g->fpt_thickness);
  //
  free(g->fvtu_1d);
  free(g->fvtu_3d);
  return;
}
/*********************************************************************************************/
static INT init_pts(const INT dim, const INT npts, REAL *pts, scomplex *sc, const REAL scale)
{
  INT k = 0, i, j;
  cube2simp *c2s = cube2simplex(dim);  // now we have the vertices of the unit cube in bits
  REAL *vc = calloc(4 * dim * c2s->nvcube, sizeof(REAL));
  REAL *xmintmp = vc;                                 // maps to [0...0]
  REAL *xmaxtmp = xmintmp + dim * (c2s->nvcube - 1);  // last vertex
  REAL *xmin = xmaxtmp + dim * (c2s->nvcube - 1);     // last vertex
  REAL *xmax = xmin + dim * (c2s->nvcube - 1);        // last vertex
  INT kdimi;
  for (i = 0; i < dim; i++) {
    xmax[i] = pts[i];
    xmin[i] = xmax[i];
    kdimi = dim + i;
    for (k = 1; k < npts; k++) {
      if (pts[kdimi] > xmax[i]) {
        xmax[i] = pts[kdimi];
      }
      if (pts[kdimi] < xmin[i]) {
        xmin[i] = pts[kdimi];
      }
      kdimi += dim;
    }
  }
  for (i = 0; i < dim; i++) {
    xmaxtmp[i] = xmax[i] + (scale - 1e0) * (xmax[i] - xmin[i]);
    xmintmp[i] = xmin[i] - (scale - 1e0) * (xmax[i] - xmin[i]);
  }
  for (j = 1; j < c2s->nvcube - 1; j++) {
    for (i = 0; i < dim; i++) {
      vc[j * dim + i] =
          xmintmp[i] + (xmaxtmp[i] - xmintmp[i]) * (c2s->bits[dim * j + i]);
    }
  }
  mapit(sc, vc);
  free(vc);
  cube2simp_free(c2s);
  return 0;
}
/****************************************************************************************/
static void special_1d(scomplex *sc, data_1d *g, dvector *seg_r) {
  // construct 1-homogenous simplicial complex embedded in 2d or 3d
  INT i, j, k, l, m, dimbig;
  dimbig = g->dimbig;
  REAL *xvtmp = g->xv + g->nv * dimbig;
  REAL r = -1e20;
  INT i1, i2, ko, bego, endo, ptrn, ptre;
  bego = 0, ptrn = 0;
  ptre = 0;
  for (i = 0; i < g->nseg; ++i) {
    l = g->divisions[i] - 2;
    i1 = g->seg[2 * i];
    i2 = g->seg[2 * i + 1];
    r = g->seg_radius[i];
    if (!l) {
      endo = bego + 1;
      bego = endo + 1;
      sc->nodes[2 * ptre] = i1;
      sc->nodes[2 * ptre + 1] = i2;
      seg_r->val[ptre] = r;
      ptre++;
      continue;
    }
    endo = bego + l + 1;
    sc->nodes[2 * ptre] = i1;
    sc->nodes[2 * (ptre + l) + 1] = i2;
    seg_r->val[ptre + l] = r;
    k = ptrn;
    ko = bego + 1;
    for (j = ptre; j < (ptre + l); ++j) {
      sc->nodes[2 * j + 1] = k + g->nv;
      sc->nodes[2 * (j + 1)] = k + g->nv;
      // fprintf(stdout,"\nt123(%d,1:2)=[%d
      // %d];",j+1,sc->nodes[2*j]+1,sc->nodes[2*j+1]+1);
      seg_r->val[j] = r;
      for (m = 0; m < dimbig; ++m) {
        xvtmp[dimbig * k + m] = xvtmp[dimbig * ko + m];
      }
      k++;
      ko++;
    }
    // fprintf(stdout,"\nt123(%d,1:2)=[%d
    // %d];",ptre+l+1,sc->nodes[2*(ptre+l)]+1,sc->nodes[2*(ptre+l)+1]+1);
    bego = endo + 1;
    ptrn += l;
    ptre += (l + 1);
  }
  sc->nv = ptrn + g->nv;
  sc->ns = ptre;
  seg_r->val = realloc(seg_r->val, (sc->ns) * sizeof(REAL));
  sc->nodes = realloc(sc->nodes, (sc->n + 1) * (sc->ns) * sizeof(INT));
  free(sc->x);
  g->xv = realloc(g->xv, sc->nv * dimbig * sizeof(REAL));  //
  sc->x = g->xv;
  return;
}
/**************************************************************************************/
static void read_and_setup(const char *finput_solver,const char *dir_matrices, \
		    input_param *inparam, dCSRmat *A,   dvector *b, dvector *x)
{
  /* matrix and right hand side: reads block, returns monolitic */
  block_dCSRmat Ablk;

  INT i,j;
  fprintf(stdout,"\n===========================================================================\n");
  fprintf(stdout,"Reading the matrix, right hand side, and parameters\n");
  fprintf(stdout,"===========================================================================\n");
  /* set Parameters from Reading in Input File */
  param_input_init(inparam);
  param_input(finput_solver,inparam);
  /* read the matrix and right hand side (2 by 2) */
  INT brow = 2;
  INT bcol = 2;
  //bdcsr_alloc_minimal(brow, bcol, &Ablk);
  bdcsr_alloc(brow, bcol, &Ablk);
  unsigned char fmt='B'; // 'b' or 'B' is for binary format
  // Read the 00 block of the stiffness matrix
  /************************************************************/
  const char *fnames_mat[] = {"A.npy","B.npy","Bt.npy","C.npy","b0.npy","b1.npy","idofs.npy","\0"};
  //
  char *fmata  = fname_set(dir_matrices, fnames_mat[0]);
  char *fmatb  = fname_set(dir_matrices, fnames_mat[1]);
  char *fmatbt = fname_set(dir_matrices, fnames_mat[2]);
  char *fmatc  = fname_set(dir_matrices, fnames_mat[3]);
  char *fb0    = fname_set(dir_matrices, fnames_mat[4]);
  char *fb1    = fname_set(dir_matrices, fnames_mat[5]);
  char *fidofs = fname_set(dir_matrices, fnames_mat[6]);
  fprintf(stdout,"\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
	  fmata,				\
	  fmatb,				\
	  fmatbt,				\
	  fmatc,				\
	  fb0,					\
	  fb1,					\
	  fidofs);fflush(stdout);

// reading
  Ablk.blocks[0]=dcoo_read_eof_dcsr_p(fmata,NULL,fmt);
  Ablk.blocks[3]=dcoo_read_eof_dcsr_p(fmatc,NULL,fmt);
  INT size[2];
  size[0]=Ablk.blocks[0]->row;  size[1]=Ablk.blocks[3]->row;
  Ablk.blocks[1]=dcoo_read_eof_dcsr_p(fmatbt,size,fmt);
  size[1]=Ablk.blocks[0]->row;  size[0]=Ablk.blocks[3]->row;
  Ablk.blocks[2]=dcoo_read_eof_dcsr_p(fmatb,size,fmt);
  //
  dvector **b_blk=malloc(bcol*sizeof(dvector));
  b_blk[0]=dvector_read_eof_p(fb0,fmt);
  b_blk[1]=dvector_read_eof_p(fb1,fmt);
  // ivector *idofs=(ivector*)malloc(sizeof(ivector));
  //  idofs = ivector_read_eof_p(fidofs,fmt);
  //
  // clean filenames
  free(fmata);  free(fmatb);  free(fmatbt); free(fmatc);  free(fb0);    free(fb1); free(fidofs); 
  fmata=NULL; fmatb=NULL; fmatbt=NULL; fmatc=NULL; fb0=NULL; fb1=NULL;  fidofs=NULL; 
  b->row=b_blk[0]->row+b_blk[1]->row;
  b->val=calloc(b->row,sizeof(REAL));
  memcpy(b->val,b_blk[0]->val,b_blk[0]->row*sizeof(REAL));
  memcpy(&b->val[b_blk[0]->row],b_blk[1]->val,b_blk[1]->row*sizeof(REAL));
  for(i=0;i<brow;++i)
    free(b_blk[i]);// not needed any longer
  free(b_blk);
  /* set initial guess */
  dvec_alloc(b->row,x);
  dvec_set(b->row,x, 0e0);
  /*************** *************************************/
  A[0]=bdcsr_2_dcsr(&Ablk);
  /* Ablk is not needed */
  for(i=0;i<brow;++i){
    for(j=0;j<bcol;++j){
      /*these were allocated with _p, so they are just freed as linear vectors*/
      dcsr_free(Ablk.blocks[i*bcol+j]);
    }
  }
  bdcsr_free(&Ablk);
  /* write to check: do not do it for large matrices*/
  /* dcsr_write_dcoo("Ablk.ztxt",Ablk.blocks[0]); */
  /* dcsr_write_dcoo("C.ztxt",Ablk.blocks[3]); */
  /* dcsr_write_dcoo("Bt.ztxt",Ablk.blocks[1]); */
  /* dcsr_write_dcoo("B.ztxt",Ablk.blocks[2]); */
  /* dcsr_write_dcoo("A.ztxt",&A); */
  /* exit(33);   */
  return;
}
/*************************************************************************/
