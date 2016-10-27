/*! \file src/utilities/coefficient.c
 *
 *  \brief This code create generic functions for coefficients that are commonly used
 *        such as zero scalars and vectors and coefficients of 1.
 *  
 *  Created by James Adler and Xiaozhe Hu on 2/22/16.
 *  Copyright 2016__HAZMAT__. All rights reserved.
 *
 *  \note modified by Xiaozhe Hu 10/27/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/27/2016
 *
 */

#include "hazmat.h"

/****************************************************************************************/
void constant_coeff_scal(REAL *val,
                         REAL *x,
                         REAL constval)
{
   /*!
    * \fn void constant_coeff_scal(REAL *val, REAL *x, REAL constval)
    *
    * \brief Assigns constant value (uses "time" slot in assembly routines)
    *
    * \param val        Pointer to a scalar REAL number
    * \param x          Pointer to an REAL array
    * \param constval   constant value (scalar REAL number)
    *
    */

    *val = constval;
}
/****************************************************************************************/

/****************************************************************************************/
void constant_coeff_vec2D(REAL *val,
                          REAL* x,
                          REAL constval)
{
    /*!
     * \fn void constant_coeff_vec2D(REAL *val,REAL* x,REAL constval)
     *
     * \brief Assigns constant value (uses "time" slot in assembly routines) in 2D
     *
     * \param val        Pointer to an REAL array
     * \param x          Pointer to an REAL array
     * \param constval   constant value (scalar REAL number)
     *
     */

    val[0] = constval;
    val[1] = constval;
}
/****************************************************************************************/

/****************************************************************************************/
void constant_coeff_vec3D(REAL *val,
                          REAL* x,
                          REAL constval)
{
    /*!
     * \fn void constant_coeff_vec3D(REAL *val,REAL* x,REAL constval)
     *
     * \brief Assigns constant value (uses "time" slot in assembly routines) in 3D
     *
     * \param val        Pointer to an REAL array
     * \param x          Pointer to an REAL array
     * \param constval   constant value (scalar REAL number)
     *
     */

    val[0] = constval;
    val[1] = constval;
    val[2] = constval;
}
/****************************************************************************************/

/****************************************************************************************/
void zero_coeff_scal(REAL *val,
                     REAL* x,
                     REAL time)
{
    /*!
     * \fn void zero_coeff_scal(REAL *val,REAL* x,REAL time)
     *
     * \brief Assigns zero value
     *
     * \param val        Pointer to a scalar REAL number
     * \param x          Pointer to an REAL array
     * \param constval   constant value (scalar REAL number)
     *
     */

    *val = 0.0;
}
/****************************************************************************************/

/****************************************************************************************/
void zero_coeff_vec2D(REAL *val,
                      REAL* x,
                      REAL time)
{
    /*!
     * \fn void zero_coeff_vec2D(REAL *val,REAL* x,REAL time)
     *
     * \brief Assigns zero value in 2D
     *
     * \param val        Pointer to an REAL array
     * \param x          Pointer to an REAL array
     * \param constval   constant value (scalar REAL number)
     *
     */

    val[0] = 0.0;
    val[1] = 0.0;
}
/****************************************************************************************/

/****************************************************************************************/
void zero_coeff_vec3D(REAL *val,REAL* x,REAL time) 
{
    /*!
     * \fn void zero_coeff_vec3D(REAL *val,REAL* x,REAL time)
     *
     * \brief Assigns zero value in 3D
     *
     * \param val        Pointer to an REAL array
     * \param x          Pointer to an REAL array
     * \param constval   constant value (scalar REAL number)
     *
     */

    val[0] = 0.0;
    val[1] = 0.0;
    val[2] = 0.0;
}
/****************************************************************************************/

/****************************************************************************************/
void one_coeff_scal(REAL *val,
                    REAL* x,
                    REAL time)
{
    /*!
     * \fn void one_coeff_scal(REAL *val,REAL* x,REAL time)
     *
     * \brief Assigns value 1
     *
     * \param val        Pointer to a scalar REAL number
     * \param x          Pointer to an REAL array
     * \param time       REAL number
     *
     */

    *val = 1.0;
}
/****************************************************************************************/

/****************************************************************************************/
void one_coeff_vec2D(REAL *val,
                     REAL* x,
                     REAL time)
{
    /*!
     * \fn void one_coeff_vec2D(REAL *val,REAL* x,REAL time)
     *
     * \brief Assigns value 1 in 2D
     *
     * \param val        Pointer to a scalar REAL number
     * \param x          Pointer to an REAL array
     * \param time       REAL number
     *
     */

    val[0] = 1.0;
    val[1] = 1.0;
}
/****************************************************************************************/

/****************************************************************************************/
void one_coeff_vec3D(REAL *val,
                     REAL* x,
                     REAL time)
{
    /*!
     * \fn void one_coeff_vec3D(REAL *val,REAL* x,REAL time)
     *
     * \brief Assigns value 1 in 3D
     *
     * \param val        Pointer to a scalar REAL number
     * \param x          Pointer to an REAL array
     * \param time       REAL number
     *
     */

    val[0] = 1.0;
    val[1] = 1.0;
    val[2] = 1.0;
}
/*************************************  END  *******************************************/

