/*! \file src/utilities/wrapper.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 09/02/16.
 *  Copyright 2016__HAZMATH__. All rights reserved.
 *
 *  \note: modified by ltz1 on 08/29/2018
 *  \note: modified by Xiaozhe Hu on 10/27/2016
 *  \note: done cleanup for releasing -- Xiaozhe Hu 10/27/2016
 *
 */

#include "hazmath.h"
INT main(int *argc, char **argv) {
    dCOOmat         *Acoo;      // matrix
    dCSRmat         *Acsr;      // matrix
    dvector         *rhs, *sol; // right-hand-side, solution
    /***************************************************/
    AMG_param       amgparam; // parameters for AMG
    linear_itsolver_param  itparam;  // parameters for linear itsolver    
    input_param inparam;
    INT i;
    INT print_lvl=4;
    INT maxit=1000;
    REAL tol=1e-11;
    /***************************************************/
    Acoo=(dCOOmat *)malloc(sizeof(dCOOmat));
    FILE *fin = HAZ_fopen("LS/matrix3.ijv","r");
    i=fscanf(fin,"%i",&(Acoo->row));
    i=fscanf(fin,"%i",&(Acoo->col));
    i=fscanf(fin,"%i",&(Acoo->nnz));
    Acoo->rowind=calloc(Acoo->nnz,sizeof(INT));
    Acoo->colind=calloc(Acoo->nnz,sizeof(INT));
    Acoo->val=calloc(Acoo->nnz,sizeof(REAL));
    fprintf(stdout,"\nReading the matrix...");
    for(i=0;i<Acoo->nnz;i++){
      fscanf(fin,"%i %i %lg",(Acoo->rowind+i),(Acoo->colind+i),(Acoo->val+i));
      Acoo->rowind[i]--;Acoo->colind[i]--;
      //      fprintf(stdout,"\n%i: %i %i %23.16e",i,Acoo->rowind[i],Acoo->colind[i],Acoo->val[i]);
    }
    fprintf(stdout,"... %d nonzeroes: DONE.\n",Acoo->nnz);
    fclose(fin);
    fin = HAZ_fopen("LS/rhs3.dat","r");
    rhs=(dvector *)malloc(sizeof(dvector));
    rhs->row = Acoo->row; rhs->val = calloc(rhs->row,sizeof(REAL));
    fprintf(stdout,"\nReading the rhs...");
    rvecd_(fin,rhs->val,&(rhs->row));
    fprintf(stdout,"... %d rows: DONE.\n",rhs->row);
    fclose(fin);
    /***************************************************/
    Acsr=(dCSRmat *)malloc(sizeof(dCSRmat));
    Acsr->row=Acoo->row;
    Acsr->col=Acoo->col;
    dcoo_2_dcsr(Acoo,Acsr);
    //    FILE *fid = HAZ_fopen("matrix.nnn","w");
    //    csr_print_matlab(fid,Acsr);
    //    fclose(fid);
    //    exit(255);
    free(Acoo->val);free(Acoo->rowind);free(Acoo->colind);free(Acoo);
    sol=(dvector *)malloc(sizeof(dvector));
    sol->row = Acsr->col; sol->val = calloc(Acsr->row,sizeof(REAL));
    //
    param_input_init(&inparam);
    param_input("./input.dat", &inparam);
    
    // Set parameters for linear iterative methods
    param_linear_solver_init(&itparam);
    param_linear_solver_set(&itparam, &inparam);
    if (print_lvl > PRINT_MIN) param_linear_solver_print(&itparam);
    
    // Set parameters for algebriac multigrid methods
    param_amg_init(&amgparam);
    param_amg_set(&amgparam, &inparam);
    if (print_lvl > PRINT_MIN) param_amg_print(&amgparam);
        
    amgparam.print_level          = print_lvl;
    itparam.linear_tol            = tol;
    itparam.linear_print_level    = print_lvl;
    itparam.linear_maxit          = maxit;
    /************************************************************/    
    linear_solver_dcsr_krylov_amg(Acsr, rhs, sol, &itparam, &amgparam);
    //    directsolve_UMF(Acsr, rhs, sol, print_lvl);
    dvec_write("sol.dat",sol);
}
/***************************** END ***************************************************/

