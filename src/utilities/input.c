/*! \file src/utilities/input.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 13/6/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \note: modified by Xiaozhe Hu on 10/29/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/29/2016 & 08/28/2021
 *
 *  \todo: check errors at the end -- Xiaozhe Hu
 *
 */

#include "hazmath.h"

/***********************************************************************************************/
/*!
 * \fn input_param *param_input_p (const char *filenm)
 *
 * \brief Read input parameters from disk file
 *
 * \param filenm   Pointer to file name of the input file
 * \param inparam     Pointet to input_param structure (OUTPUT)
 *
 */
input_param *param_input_p(const char *filenm)
{
  input_param *inparam=malloc(1*sizeof(input_param));
  param_input(filenm,inparam);
  return inparam;
}

/***********************************************************************************************/
/*!
 * \fn void param_input (const char *filenm, input_param *inparam)
 *
 * \brief Read input parameters from disk file
 *
 * \param filenm   Pointer to file name of the input file
 * \param inparam     Pointet to input_param structure
 *
 */
void param_input (const char *filenm,		\
                  input_param *inparam)
{
  INT maxb=512; // max number of char for each line is maxb
  INT      val;
  SHORT    status = SUCCESS;
  char     *buffer=malloc(maxb*sizeof(char));
  if(!buffer){
    status= ERROR_ALLOC_MEM;
    check_error(status, __FUNCTION__);
  }

    // set default input parameters
  param_input_init(inparam);
    // if input file is not specified, use the default values
    if (filenm==NULL) return;
    // ltz: if we are here, the "inparam->inifile" must be the same as filenm:
    strcpy(inparam->inifile,filenm);
    // end ltz:
    FILE *fp = fopen(filenm,"r");
    if (fp==NULL) {
        status = ERROR_OPEN_FILE;
        check_error(status, __FUNCTION__);
    }
    // only read when successfully open the file
    while ( status == SUCCESS ) {
        int     ibuff;
        double  dbuff;
        char    sbuff[maxb];

        val = fscanf(fp,"%s",buffer);
        if (val==EOF) break;
        if (val!=1){ status = ERROR_INPUT_PAR; break; }
        if (buffer[0]=='[' || buffer[0]=='%' || buffer[0]=='|') {
            fgets(buffer,maxb,fp); // skip rest of line
            continue;
        }

        // match the keyword and read in the value

        // -----------
        // output
        // -----------
        if (strcmp(buffer,"print_level")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->print_level = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"output_dir")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            strncpy(inparam->output_dir,sbuff,128);
            fgets(buffer,maxb,fp); // skip rest of line
        }

        // ---------------
        // mesh stuff
        // ---------------
        else if (strcmp(buffer,"read_mesh_from_file")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->read_mesh_from_file = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"gridfile")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0){
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            strncpy(inparam->gridfile,sbuff,128);
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"spatial_dim")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->spatial_dim = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"refinement_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->refinement_type = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"refinement_levels")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->refinement_levels = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"boundary_codes")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->boundary_codes = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        // ---------------
        // finite element
        // ---------------
        else if (strcmp(buffer,"nquad")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->nquad = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"FE_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->FE_type = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"Mass_lump")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->Mass_lump = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        // --------------
        // time stepping
        // --------------
        else if (strcmp(buffer,"time_start")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->time_start = dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"time_step_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->time_step_type = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"time_steps")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->time_steps = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"time_step_size")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->time_step_size = dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"rhs_time_dep")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->rhs_time_dep = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        // -----------------
        // nonlinear solver
        // -----------------
        else if (strcmp(buffer,"nonlinear_itsolver_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->nonlinear_itsolver_type = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"nonlinear_itsolver_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->nonlinear_itsolver_maxit = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"nonlinear_itsolver_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->nonlinear_itsolver_tol = dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"nonlinear_itsolver_toltype")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->nonlinear_itsolver_toltype = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"fas_presmoothers")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->fas_presmoothers = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"fas_postsmoothers")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->fas_postsmoothers = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"fas_smooth_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->fas_smooth_tol = dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        // --------------
        // linear solver
        // --------------
        else if (strcmp(buffer,"linear_itsolver_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_itsolver_type = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"linear_itsolver_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_itsolver_maxit = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"linear_itsolver_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_itsolver_tol = dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"linear_stop_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_stop_type = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"linear_restart")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_restart = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"linear_precond_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->linear_precond_type = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        // -----------
        // AMG
        // -----------
        else if (strcmp(buffer,"AMG_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_tol = dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_maxit = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }

            if ((strcmp(buffer,"UA")==0)||(strcmp(buffer,"ua")==0))
                inparam->AMG_type = UA_AMG;
            else if ((strcmp(buffer,"SA")==0)||(strcmp(buffer,"sa")==0))
                inparam->AMG_type = SA_AMG;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_levels")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_levels = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
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
            else if ((strcmp(buffer,"ADD")==0)||(strcmp(buffer,"add")==0))
                inparam->AMG_cycle_type = ADD_CYCLE;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_amli_degree")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_amli_degree = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_nl_amli_krylov_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_nl_amli_krylov_type = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
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
            else if ((strcmp(buffer,"FJACOBI")==0)||(strcmp(buffer,"fjacobi")==0))
                inparam->AMG_smoother = SMOOTHER_FJACOBI;
            else if ((strcmp(buffer,"FGS")==0)||(strcmp(buffer,"fgs")==0))
                inparam->AMG_smoother = SMOOTHER_FGS;
            else if ((strcmp(buffer,"FSGS")==0)||(strcmp(buffer,"fsgs")==0))
                inparam->AMG_smoother = SMOOTHER_FSGS;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_relaxation")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_relaxation=dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_Schwarz_levels")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
              status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_Schwarz_levels = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_presmooth_iter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_presmooth_iter = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_postsmooth_iter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_postsmooth_iter = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_coarse_dof")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_coarse_dof = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_coarse_solver")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_coarse_solver = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
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
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_fpwr")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_fpwr=dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_aggregation_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_aggregation_type = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_strong_coupled")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_strong_coupled = dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_max_aggregation")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_max_aggregation = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        //-------------------
        // SA AMG
        //-------------------
        else if (strcmp(buffer,"AMG_tentative_smooth")==0) {
           val = fscanf(fp,"%s",buffer);
           if (val!=1 || strcmp(buffer,"=")!=0) {
               status = ERROR_INPUT_PAR; break;
           }
           val = fscanf(fp,"%lf",&dbuff);
           if (val!=1) { status = ERROR_INPUT_PAR; break; }
           inparam->AMG_tentative_smooth = dbuff;
           fgets(buffer,500,fp); // skip rest of line
       }

       else if (strcmp(buffer,"AMG_smooth_filter")==0) {
           val = fscanf(fp,"%s",buffer);
           if (val!=1 || strcmp(buffer,"=")!=0) {
               status = ERROR_INPUT_PAR; break;
           }
           val = fscanf(fp,"%s",buffer);
           if (val!=1) { status = ERROR_INPUT_PAR; break; }

           if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
               (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
               inparam->AMG_smooth_filter = ON;
           else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                    (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                    (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                    (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
               inparam->AMG_smooth_filter = OFF;
           else
               { status = ERROR_INPUT_PAR; break; }
           fgets(buffer,500,fp); // skip rest of line
       }

        // ------------------
        // Schwarz method
        // ------------------
        else if (strcmp(buffer,"Schwarz_mmsize")==0) {
          val = fscanf(fp,"%s",buffer);
          if (val!=1 || strcmp(buffer,"=")!=0) {
            status = ERROR_INPUT_PAR; break;
          }
          val = fscanf(fp,"%d",&ibuff);
          if (val!=1) { status = ERROR_INPUT_PAR; break; }
          inparam->Schwarz_mmsize = ibuff;
          fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"Schwarz_maxlvl")==0)
        {
          val = fscanf(fp,"%s",buffer);
          if (val!=1 || strcmp(buffer,"=")!=0) {
            status = ERROR_INPUT_PAR; break;
          }
          val = fscanf(fp,"%d",&ibuff);
          if (val!=1) {status = ERROR_INPUT_PAR; break; }
          inparam->Schwarz_maxlvl = ibuff;
          fgets(buffer,maxb,fp); // skip rest of line
        }

        else if (strcmp(buffer,"Schwarz_type")==0)
        {
          val = fscanf(fp,"%s",buffer);
          if (val!=1 || strcmp(buffer,"=")!=0) {
            status = ERROR_INPUT_PAR; break;
          }
          val = fscanf(fp,"%d",&ibuff);
          if (val!=1) { status = ERROR_INPUT_PAR; break; }
          inparam->Schwarz_type = ibuff;
          fgets(buffer,maxb,fp); // skip rest of line
        }
        else if (strcmp(buffer,"Schwarz_blksolver")==0)
        {
          val = fscanf(fp,"%s",buffer);
          if (val!=1 || strcmp(buffer,"=")!=0) {
            status = ERROR_INPUT_PAR; break;
          }
          val = fscanf(fp,"%d",&ibuff);
          if (val!=1) { status = ERROR_INPUT_PAR; break; }
          inparam->Schwarz_blksolver = ibuff;
          fgets(buffer,maxb,fp); // skip rest of line
        }

        // ------------------
        // HX-preconditioner
        // ------------------
        else if (strcmp(buffer,"HX_smooth_iter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->HX_smooth_iter = ibuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        // ------------------
        // BSR-preconditioner
        // ------------------
        else if (strcmp(buffer,"BSR_alpha")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->BSR_alpha = dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }
        else if (strcmp(buffer,"BSR_omega")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->BSR_omega = dbuff;
            fgets(buffer,maxb,fp); // skip rest of line
        }

        else {
            printf("### HAZMATH WARNING: Unknown input keyword %s!\n", buffer);
            fgets(buffer,maxb,fp); // skip rest of line
        }
    }

    fclose(fp);
    if(buffer) free(buffer);
    // check errors -- to be added,  Xiaozhe

}

/************************************ END ******************************************************/
