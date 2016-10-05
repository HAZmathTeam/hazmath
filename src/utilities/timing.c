/*! \file timing.c
 *
 *  Created by James Adler and Xiaozhe Hu on 10/06/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

void gettime (REAL *time)
{
    /**
     * \fn gettime (REAL *time)
     *
     * \brief Get system time
     *
     * \author Xiaozhe
     * \date   10/06/2015
     *
     */
    
    if ( time != NULL ) {
        *time = (REAL) clock() / CLOCKS_PER_SEC;
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
