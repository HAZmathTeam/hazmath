/*! \file src/mesh/mesh_create.c
*
* \brief Coordinate allocation utilities.
*
*  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 1/9/15.
*  Copyright 2015__HAZMATH__. All rights reserved.
*
*/

#include "hazmath.h"

/* creategrid_fread, make_uniform_mesh, create1Dgrid_Line moved to scomplex.c */

/*!
* \fn struct coordinates *allocatecoords(INT ndof,INT mydim)
*
* \brief Allocates memory and properties of coordinates struct
*
* \param ndof      Total number of coordinates
* \param mydim     Dimension of coordinates
*
* \return A        Coordinate struct
*
*/
struct coordinates *allocatecoords(INT ndof,INT mydim)
{
  // Flag for errors
  SHORT status;

  struct coordinates *A = malloc(sizeof(struct coordinates));
  assert(A != NULL);

  A->x = (REAL *) calloc(mydim*ndof,sizeof(REAL));
  switch (mydim)
  {
    case 1:
    A->y=NULL;
    A->z=NULL;
    break;
    case 2:
    A->y = A->x + ndof;
    A->z=NULL;
    break;
    case 3:
    A->y = A->x + ndof;
    A->z = A->y + ndof;
    break;
    default:
    status = ERROR_DIM;
    check_error(status, __FUNCTION__);
  }
  A->n = ndof;
  return A;
}

/*!
* \fn void free_coords(coordinates* A)
*
* \brief Frees memory of arrays of coordinate struct
*
* \param A         Pointer to coordinates struct to be freed
*
*/
void free_coords(coordinates* A)
{
  if (A==NULL) return;

  if(A->x) {
    free(A->x);
    A->x = NULL;
  }

  return;
}
