//
//  mesh.h
//
//
//  Created by James Adler, Xiaozhe Hi, and Ludmil Zikatanov 2015-01-09.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _mesh_h
#define _mesh_h

#include "sparse.h"
#include "vec.h"

/**
 * \struct coords
 * \brief Row-major coordinate storage: x[i*dim+j] = component j of node i
 */
typedef struct coords{

  //! coordinate values (ordered x(0) y(0) z(0) x(1) y(1) z(1) etc...))
  REAL* x;

  //! Size of arrays (number of nodes)
  INT n;

  //! dimension
  INT dim;

} coords;

/**
 * \struct coordinates
 * \brief Column-major coordinate storage for FE DOF locations and quadrature.
 *        x[i] = x-coord of node i, y[i] = y-coord, z[i] = z-coord.
 *        (y = x + n, z = y + n in memory)
 */
typedef struct coordinates{

  //! x values
  REAL* x;

  //! y values
  REAL* y;

  //! z values (if in 3D)
  REAL* z;

  //! Size of arrays (number of nodes)
  INT n;

} coordinates;


#endif
