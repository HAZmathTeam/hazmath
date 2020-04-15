/* two functions to go in io.c*/
/***********************************************************************************************/
/***********************************************************************************************/
/*!
 * \fn SHORT dcoo_2_dcsr (dCOOmat *A, dCSRmat *B)
 *
 * \brief Transform a dCOOmat matrix to a dCSRmat format.
 *
 * \param A   Pointer to dCOOmat matrix
 * \param B   Pointer to dCSRmat matrix
 *
 * \return    SUCCESS if successed; otherwise, error information.
 *
 */
dCSRmat *dcoo_2_dcsr_w (dCOOmat *A)
{
    // get size
    const INT m=A->row, n=A->col, nnz=A->nnz;
    // allocate 
    dCSRmat *B=dcsr_create_w(m,n,nnz);
    // local variable
    INT *ia = B->IA;
    INT *ja = B->JA;
    REAL *Bval = B->val;
    INT *row_idx = A->rowind;
    INT *col_idx = A->colind;
    REAL *Aval = A->val;
    INT i, iind, jind;

    INT *ind = (INT *) calloc(m+1,sizeof(INT));

    // initialize
    memset(ind, 0, sizeof(INT)*(m+1));

    // count number of nonzeros in each row
    for (i=0; i<nnz; ++i) ind[row_idx[i]+1]++;

    // set row pointer
    ia[0] = 0;
    for (i=1; i<=m; ++i) {
        ia[i] = ia[i-1]+ind[i];
        ind[i] = ia[i];
    }

    // set column index and values
    for (i=0; i<nnz; ++i) {
        iind = row_idx[i];
        jind = ind[iind];
        ja[jind] = col_idx[i];
        Bval[jind] = Aval[i];
        ind[iind] = ++jind;
    }

    if (ind) free(ind);

    return B;
}
/**************************************************************************************/
/**
 * \fn void dcoo_read_eof_dcsr(FILE *fp, dCSRmat *A, INT *size)
 *
 * \brief Read A from matrix disk file in I,J,V format and convert to CSR format
 *        first reads until the end of file to count the number of nonzeroes
 *        then reqinds the file and reads all nonzeroes.
 *        The number of rows is the max of all integers from I[] and size[0]
 *        The number of columns is the max of all integers from J[] and size[1]
 *        If size is NULL is the same as size[0]=size[1]=0
 *
 * \param filename  File name for matrix
 * \param A         Pointer to the CSR matrix
 * \param size      Pointer to an array with two integers or NULL.
 *
 * \note File format:
 *   - nrow ncol nnz     % number of rows, number of columns, and nnz
 *   - i  j  a_ij        % i, j a_ij in each line
 *
 */
dCSRmat *dcoo_read_eof_dcsr_w (FILE *fp,INT *size)
{
    int  i,j,k,m,n,nnz,ichk;
    REAL value;
    dCSRmat *A=NULL;
    if(size){
      m=size[0];
      n=size[1];
    } else {
      m=-1;
      n=-1;
    }
    nnz=0;
    while(!feof(fp)){
      ichk=fscanf(fp,"%d %d %lg", &i,&j, &value);
      if(ichk<0) break;
      if(i>m) m=i;
      if(j>n) n=j;
      nnz++;
    }
    fprintf(stdout,"%s FOUND %d records.\n",__FUNCTION__,nnz);
    // check if we read anything...
    if(!nnz) check_error(ERROR_WRONG_FILE, __FUNCTION__);
    //
    m++;n++; //bc indices strat from zero so row and coll are with one more.
    dCOOmat *Atmp=dcoo_create_w(m,n,nnz);
    //    fprintf(stdout,"\nCheck structure: %d %d %d\n\n",Atmp->row,Atmp->col,Atmp->nnz);fflush(stdout);
    rewind(fp);
    for ( k = 0; k < nnz; k++ ) {
        if ( fscanf(fp, "%d %d %lg", &i, &j, &value) != EOF ) {
            Atmp->rowind[k]=i; Atmp->colind[k]=j; Atmp->val[k] = value;
        } else {
          check_error(ERROR_WRONG_FILE, __FUNCTION__);
        }
    }
    A=dcoo_2_dcsr_w(Atmp);
    free(Atmp);
    return A;
}
/***********************************************************************************************/
/**
 * \fn dvector *dvector_read_eof(FILE *fp)
 *
 * \brief Read b from a disk file until EOF in array format
 *
 * \param filename  File name for vector b
 * \param b         Pointer to the dvector b (output)
 *
 * \note File Format:
 *   - nrow
 *   - val_j, j=0:nrow-1
 *
 */
dvector *dvector_read_eof_w(FILE *fp)
{
    int  i, n,ichk;
    REAL value;
    dvector *b=NULL;
    n=0;
    while(!feof(fp)){
      ichk=fscanf(fp,"%lg", &value);
      if(ichk<0) break;
      n++;
    }
    fprintf(stdout,"%s: FOUND %d records.\n",__FUNCTION__,n);
    if(!n)
      check_error(ERROR_WRONG_FILE, __FUNCTION__);
    b=dvec_create_w(n);
    rewind(fp);
    for ( i = 0; i < n; ++i ) {
        fscanf(fp, "%lg", &value);
        b->val[i] = value;
        if ( value > BIGREAL ) {
            printf("### ERROR: Wrong value = %lf\n", value);
            fclose(fp);
            if(b) free(b);
            exit(ERROR_INPUT_PAR);
        }
    }
    return b;
}
/*END OF SOLVERS_READ.h*/
