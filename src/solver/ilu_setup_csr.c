/*
 *  ilu_setup_csr.c
 *
 *  Setup incomplete LU decomposition for dCSRmat matrices
 *
 *  Created by James Adler and Xiaozhe Hu on 12/26/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 *  Ref Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller
 *        Appendix P475 A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *        Academic Press Inc., San Diego, CA, 2001.
 */

#include "hazmat.h"

/* declarations for ilu.for */
#ifdef __cplusplus 
extern "C" {void iluk_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,REAL *alu,
                       INT *jlu,INT *iwk,INT *ierr,INT *nzlu);}
extern "C" {void ilut_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,
                       const REAL *droptol,REAL *alu,INT *jlu,INT *iwk,
                       INT *ierr,INT *nz);}
extern "C" {void ilutp_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,
                        const REAL *droptol,const REAL *permtol,const INT *mbloc,
                        REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nz);}
#else
extern void iluk_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,REAL *alu,
                  INT *jlu,INT *iwk,INT *ierr,INT *nzlu);
extern void ilut_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,const REAL *droptol,
                  REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nz);
extern void ilutp_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,const REAL *droptol,
                   const REAL *permtol,const INT *mbloc,REAL *alu,INT *jlu,INT *iwk,
                   INT *ierr,INT *nz);
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

SHORT ilu_dcsr_setup (dCSRmat *A,
                      ILU_data *iludata,
                      ILU_param *iluparam)
{

    /**
     * \fn SHORT ilu_dcsr_setup (dCSRmat *A, ILU_data *iludata, ILU_param *iluparam)
     *
     * \brief Get ILU decomposition of a CSR matrix A
     *
     * \param A         Pointer to dCSRmat matrix
     * \param iludata   Pointer to ILU_data
     * \param iluparam  Pointer to ILU_param
     *
     * \return  SUCCESS if successed; otherwise, error information.
     *
     * \author Shiquan Zhang, Xiaozhe Hu
     * \date   12/27/2009
     */

    const INT   type=iluparam->ILU_type, print_level=iluparam->print_level;
    const INT   n=A->col, nnz=A->nnz, mbloc=n;
    const REAL  ILU_droptol=iluparam->ILU_droptol;
    const REAL  permtol=iluparam->ILU_permtol;
    
    // local variable
    INT    lfil=iluparam->ILU_lfil, lfilt=iluparam->ILU_lfil;
    INT    ierr, iwk, nzlu, nwork, *ijlu;
    REAL  *luval;
    
    REAL   setup_start, setup_end, setup_duration;
    SHORT  status = SUCCESS;
    
    gettime(&setup_start);
    
    // Expected amount of memory for ILU needed and allocate memory 
    switch (type) {
        case ILUt:
            iwk=10*nnz;     // iwk is the maxim possible nnz for ILU
            lfilt=floor(n*0.5)+1;
            break;
        case ILUtp:
            iwk=10*nnz;     // iwk is the maxim possible nnz for ILU
            lfilt=floor(n*0.5)+1;
            break;
        default: // ILUk
            if (lfil == 0) iwk=nnz+500;
            else iwk=(lfil+5)*nnz;
            break;
    } 
    
    nwork  = 4*n;
    
    // setup ILU preconditioner
    iludata->row=iludata->col=n;    
    ilu_data_alloc(iwk, nwork, iludata);
    
    // ILU decomposition
    ijlu=iludata->ijlu;
    luval=iludata->luval;
    
    switch (type) {
        case ILUt:
            
            ilut_(&n,A->val,A->JA,A->IA,&lfilt,&ILU_droptol,luval,ijlu,&iwk,&ierr,&nzlu);
            break;
            
        case ILUtp:

            ilutp_(&n,A->val,A->JA,A->IA,&lfilt,&ILU_droptol,&permtol,
                   &mbloc,luval,ijlu,&iwk,&ierr,&nzlu);
            break;
            
        default: // ILUk
            
            iluk_(&n,A->val,A->JA,A->IA,&lfil,luval,ijlu,&iwk,&ierr,&nzlu);
            break;
    } 
    
    dcsr_shift(A, -1);
    
    iludata->nzlu=nzlu;
    iludata->nwork=nwork;
    
    if (ierr!=0) {
        printf("### ERROR: ILU setup failed (ierr=%d)!\n", ierr);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if (iwk<nzlu) {
        printf("### ERROR: Need more memory for ILU %d!\n", iwk-nzlu);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if (print_level>PRINT_NONE) {
        gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        
        switch (type) {
            case ILUt:
                printf("ILUt setup costs %f seconds.\n", setup_duration);    
                break;
            case ILUtp:
                printf("ILUtp setup costs %f seconds.\n", setup_duration);    
                break;
            default: // ILUk
                printf("ILUk setup costs %f seconds.\n", setup_duration);    
                break;
        }     
    }
    
FINISHED:     
    
    return status;
    
}

SHORT mem_iludata_check (ILU_data *iludata)
{
    /**
     * \fn SHORT mem_iludata_check (ILU_data *iludata)
     *
     * \brief Check wether a ILU_data has enough work space
     *
     * \param iludata    Pointer to be cheked
     *
     * \return           SUCCESS if success, else ERROR (negative value)
     *
     * \author Xiaozhe Hu, Chensong Zhang
     * \date   11/27/09
     */
    
    const INT memneed = 2*iludata->row; // estimated memory usage
    
    if ( iludata->nwork >= memneed ) {
        return SUCCESS;
    }
    else {
        printf("### ERROR: ILU needs %d RAM, only %d allocated!\n",
               memneed, iludata->nwork);
        return ERROR_ALLOC_MEM;
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
