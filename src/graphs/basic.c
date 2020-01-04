/*! \file src/graphs/basic.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 01/03/2020.
 *  Copyright 2020__HAZMATH__. All rights reserved.
 *
 *  \note basic subroutines for graphs;
 *
 */

#include "hazmath.h"

/*!
 * \fn get_adjacency_from_transition(dCSRmat *P, dvector *weighted_degree)
 *
 * \brief get adjacency matrix from transition matrix (A = D*P)
 *
 * \param P                pointer to the transition matrix
 * \param weighted_degree  pointer to the weighted degree
 *
 * \return A               pointer to the adjacency matrix
 *
 */
dCSRmat get_adjacency_from_transition(dCSRmat *P,
                                      dvector *weighted_degree)
{
  // local variables
  dCSRmat A;

  // copy transition matrix
  dcsr_alloc(P->row, P->col, P->nnz, &A);
  dcsr_cp(P, &A);

  // A = D*P
  dcsr_row_scale(&A, weighted_degree);

  // return adjacency matrix
  return A;
}

/*!
 * \fn get_graphLaplacian_from_adj_wdeg(dCSRmat *A, dvector *weighted_degree)
 *
 * \brief get graph Laplacian from adjacency matrix and weighted degree matrix (L = D-A)
 *
 * \param A                pointer to the adjacency matrix
 * \param weighted_degree  pointer to the weighted degree
 *
 * \return L               pointer to the graph Laplacian
 *
 */
dCSRmat get_graphLaplacian_from_adjacency(dCSRmat *A,
                                          dvector *weighted_degree)
{
  // local variables
  dCSRmat L;

  // form diagonal weighted degree matrix
   dCSRmat D = dcsr_create_diagonal_matrix(weighted_degree);

  // L =  D-A
  dcsr_add(&D, 1.0, A, -1.0, &L);

  // clean
  dcsr_free(&D);

  // return adjacency matrix
  return L;
}

/*!
 * \fn get_normalizedgraphLaplacian_from_L_wdeg_inv(dCSRmat *L, dvector *weighted_degree_half_inv)
 *
 * \brief get normalized graph Laplacian from graph Laplacian and inverse half weighted degree matrix (N = D^-1/2*L*D^-1/2)
 *
 * \param L                pointer to the graph Laplacian
 * \param weighted_degree_half_inv  pointer to the inverse half weighted degree
 *
 * \return N               pointer to the normalized graph Laplacian
 *
 */
dCSRmat get_normalizedgraphLaplacian_from_L_wdeg_inv(dCSRmat *L,
                                                     dvector *weighted_degree_half_inv)
{
  // local variables
  dCSRmat N;

  // form diagonal weighted degree matrix
  dCSRmat D_half_inv = dcsr_create_diagonal_matrix(weighted_degree_half_inv);

  // N = D^-1/2*L*D^-1/2
  dcsr_rap(&D_half_inv, L, &D_half_inv, &N);

  // clean
  dcsr_free(&D_half_inv);

  // return adjacency matrix
  return N;
}
