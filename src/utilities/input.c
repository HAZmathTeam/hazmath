/*
 *  input.c
 *
 *  Created by James Adler and Xiaozhe Hu on 13/6/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 */

#include "hazmat.h"

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
        
        else if (strcmp(buffer,"linear_stop_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_stop_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"linear_restart")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_restart = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"linear_precond_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_precond_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_tol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_maxit = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"C")==0)||(strcmp(buffer,"c")==0))
                inparam->AMG_type = CLASSIC_AMG;
            else if ((strcmp(buffer,"SA")==0)||(strcmp(buffer,"sa")==0))
                inparam->AMG_type = SA_AMG;
            else if ((strcmp(buffer,"UA")==0)||(strcmp(buffer,"ua")==0))
                inparam->AMG_type = UA_AMG;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_levels")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_levels = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_cycle_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"V")==0)||(strcmp(buffer,"v")==0))
                inparam->AMG_cycle_type = V_CYCLE;
            else if ((strcmp(buffer,"W")==0)||(strcmp(buffer,"w")==0))
                inparam->AMG_cycle_type = W_CYCLE;
            else if ((strcmp(buffer,"A")==0)||(strcmp(buffer,"a")==0))
                inparam->AMG_cycle_type = AMLI_CYCLE;
            else if ((strcmp(buffer,"NA")==0)||(strcmp(buffer,"na")==0))
                inparam->AMG_cycle_type = NL_AMLI_CYCLE;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_smoother")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"JACOBI")==0)||(strcmp(buffer,"jacobi")==0))
                inparam->AMG_smoother = SMOOTHER_JACOBI;
            else if ((strcmp(buffer,"GS")==0)||(strcmp(buffer,"gs")==0))
                inparam->AMG_smoother = SMOOTHER_GS;
            else if ((strcmp(buffer,"SGS")==0)||(strcmp(buffer,"sgs")==0))
                inparam->AMG_smoother = SMOOTHER_SGS;
            else if ((strcmp(buffer,"CG")==0)||(strcmp(buffer,"cg")==0))
                inparam->AMG_smoother = SMOOTHER_CG;
            else if ((strcmp(buffer,"SOR")==0)||(strcmp(buffer,"sor")==0))
                inparam->AMG_smoother = SMOOTHER_SOR;
            else if ((strcmp(buffer,"SSOR")==0)||(strcmp(buffer,"ssor")==0))
                inparam->AMG_smoother = SMOOTHER_SSOR;
            else if ((strcmp(buffer,"GSOR")==0)||(strcmp(buffer,"gsor")==0))
                inparam->AMG_smoother = SMOOTHER_GSOR;
            else if ((strcmp(buffer,"SGSOR")==0)||(strcmp(buffer,"sgsor")==0))
                inparam->AMG_smoother = SMOOTHER_SGSOR;
            else if ((strcmp(buffer,"POLY")==0)||(strcmp(buffer,"poly")==0))
                inparam->AMG_smoother = SMOOTHER_POLY;
            else if ((strcmp(buffer,"L1DIAG")==0)||(strcmp(buffer,"l1diag")==0))
                inparam->AMG_smoother = SMOOTHER_L1DIAG;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_smooth_order")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"NO")==0)||(strcmp(buffer,"no")==0))
                inparam->AMG_smooth_order = NO_ORDER;
            else if ((strcmp(buffer,"CF")==0)||(strcmp(buffer,"cf")==0))
                inparam->AMG_smooth_order = CF_ORDER;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_relaxation")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_relaxation=dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_presmooth_iter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_presmooth_iter = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_postsmooth_iter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_postsmooth_iter = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarse_dof")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_coarse_dof = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarse_solver")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_coarse_solver = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_coarse_scaling")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
                inparam->AMG_coarse_scaling = ON;
            else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
                inparam->AMG_coarse_scaling = OFF;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggregation_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_aggregation_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_pair_number")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_pair_number = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_quality_bound")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_quality_bound = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_strong_coupled")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_strong_coupled = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_max_aggregation")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_max_aggregation = ibuff;
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
