//
//  quad.h
//  
//
//  Created by Adler, James on 1/9/15.
//
//
//  Routines and Structs related to Quadrature

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _quad_h
#define _quad_h

/**
 * \struct qcoordinates
 * \brief Returns coordinates and weights of quadrature nodes
 */
typedef struct qcoordinates{

  //! x values
  dvector x;

  //! y values
  dvector y;

  //! z values (if in 3D)
  dvector z;

  //! weights 
  dvector w;

  //! Size of arrays (number of nodes)
  INT n;
	
} qcoordinates;


#endif
