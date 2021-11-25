/*! \file src/utilities/timing.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 10/06/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/27/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/27/2016 & 08/28/2021
 */

#include "hazmath.h"

/*************************************************************************************/
/*!
 * \fn get_time (REAL *time)
 *
 * \brief Get system time
 *
 * \author Xiaozhe Hu
 * \date   10/06/2015
 *
 */
void get_time (REAL *time)
{
    if ( time != NULL ) {
        *time = (REAL) clock() / CLOCKS_PER_SEC;
    }
}

/******************************* END **************************************************/
