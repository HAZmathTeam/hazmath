//
//  dense.h
//
//
//  Created by Hu, Xiaozhe on 01/04/20.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _dense_h
#define _dense_h

/**
 * \struct dDENSEmat
 * \brief Dense matrix of REAL type in DENSE format
 *
 * DENSE Format in REAL
 *
 */
typedef struct dDENSEmat{

  //! number of rows
  INT row;

  //! number of columns
  INT col;

  //! array of values
  REAL *val;

} dDENSEmat;

/**
 * \struct iDENSEmat
 * \brief Dense matrix of INT type in DENSE format
 *
 * DENSE Format in INT
 *
 */
typedef struct iDENSEmat{

  //! number of rows
  INT row;

  //! number of columns
  INT col;

  //! array of values
  INT *val;

} iDENSEmat;

#endif
