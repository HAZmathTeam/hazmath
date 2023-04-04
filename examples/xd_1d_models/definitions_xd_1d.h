/*! \file examples/amr_grids/amr_grids_supporting.h
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov 2019/01/09.
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief This program
 *
 * \note
 */
/****************************************************************************/
typedef struct /* n-homogenous simplicial complex */
{
  INT dimbig;     /* the dimension of the space in which SC is embedded */
  INT dim;        /* the dimension of SC */
  INT nv;         /* number of 0-dimensional simplices */
  INT nvadd;      /* number of 0-dimensional simplices added on every segment */
  INT nseg;       /* number of segments */
  INT *seg;       /* nv boundary codes for vertices */
  INT *divisions; /*divisions per segment */
  REAL *xv; /*(nv times dimbig) array to hold the coordinates of vertices */
  REAL *pt_thickness; /* points attribute */
  REAL *seg_radius;   /* segments attribute */
  char *idir; /* directory with input files*/
  char *odir; /* directory for the output files */
  char *fv_coords;   /*input_file: coords of bifurcations*/
  char *fseg;       /*input_file: segments definition */
  char *fdivisions; /*input_file: divisions per segment*/
  char *fvtmp_coords; /*input_file: coordinates of all points (biffurcations or not */
  char *fpt_thickness;/*input_file: segment thickness*/
  char *fseg_radius;/*input_file: radius */
  char *fvtu_3d;  /*output_file: for the 3d grid in vtu*/
  char *fvtu_1d;  /*output_file: for the 1d grid on vtu */
} data_1d;
/*********************************************************************/
/************************************************************/
static char *str_add_dim(const INT dim, const char *prefix, const char *suffix);
/************************************************************/
static ivector mark_around_pts(scomplex *sc, scomplex *scglobal, INT nstar, \
			       REAL *xstar, iCSRmat *node_ins,		\
			       const INT max_nodes);
static void data_1d_init(const INT dimbig,const char *idir, const char *odir,data_1d *g);
static void getdata_1d(data_1d *g);
static INT init_pts(const INT dim, const INT npts, REAL *pts, scomplex *sc, const REAL scale);
/****************************************************************************************/
static void special_1d(scomplex *sc, data_1d *g, dvector *seg_r);
/**************************************************************************************/
