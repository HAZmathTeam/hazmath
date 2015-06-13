/*
 *  input.c
 *
 *  Created by James Adler and Xiaozhe Hu on 13/6/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include <math.h>
#include <time.h>

// Our Includes
#include "macro.h"
#include "grid.h"
#include "sparse.h"
#include "vec.h"
#include "functs.h"
#include "fem.h"
#include "param.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
void param_input (const char *filenm,
                       input_param *inparam)
{
    /**
     * \fn void param_input (const char *filenm, input_param *inparam)
     *
     * \brief Read input parameters from disk file
     *
     * \param filenm    File name for input file
     * \param inparam   Input parameters
     *
     * \author Xiaozhe Hu
     * \date   06/13/2015 
     *
     */
    
    char     buffer[500]; // Note: max number of char for each line!
    int      val;
    SHORT    status = SUCCESS;
    
    // set default input parameters  -- to be added, Xiaozhe

    // if input file is not specified, use the default values
    if (filenm==NULL) return;
    
    FILE *fp = fopen(filenm,"r");
    if (fp==NULL) {
        printf("### ERROR: Could not open file %s...\n", filenm);
    }
    
    while ( status == SUCCESS ) {
        int     ibuff;
        double  dbuff;
        char    sbuff[500];
    
        val = fscanf(fp,"%s",buffer);
        if (val==EOF) break;
        if (val!=1){ status = ERROR_INPUT_PAR; break; }
        if (buffer[0]=='[' || buffer[0]=='%' || buffer[0]=='|') {
            fgets(buffer,500,fp); // skip rest of line
            continue;
        }
    
        // match keyword and scan for value
        if (strcmp(buffer,"workdir")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            strncpy(inparam->workdir,sbuff,128);
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"print_level")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->print_level = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"output_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->output_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"gridfile")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            strncpy(inparam->gridfile,sbuff,128);
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"dim")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->dim = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"nquad")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->nquad = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"FE_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->FE_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"FE_type_velocity")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->FE_type_velocity = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"FE_type_pressure")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->FE_type_pressure = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"time_step_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->time_step_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"time_steps")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->time_steps = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"time_step_size")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->time_step_size = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"nonlinear_itsolver_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->nonlinear_itsolver_maxit = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"nonlinear_itsolver_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->nonlinear_itsolver_tol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"linear_itsolver_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_itsolver_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }

        
        else if (strcmp(buffer,"linear_itsolver_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_itsolver_maxit = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"linear_itsolver_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_itsolver_tol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }

        else {
            printf("### WARNING: Unknown input keyword %s!\n", buffer);
            fgets(buffer,500,fp); // skip rest of line
        }
    }
    
    fclose(fp);
    
    // sanity checks -- to be added,  Xiaozhe
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
