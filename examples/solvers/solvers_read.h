/* two functions to go in io.c*/
/***********************************************************************************************/
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
void dcoo_read_eof_dcsr (FILE *fp,dCSRmat *A, INT *size)
{
    int  i,j,k,m,n,nnz,ichk;
    REAL value;
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
    dCOOmat Atmp=dcoo_create(m,n,nnz);
    //    fprintf(stdout,"\nCheck structure: %d %d %d\n\n",Atmp->row,Atmp->col,Atmp->nnz);fflush(stdout);
    rewind(fp);
    for ( k = 0; k < nnz; k++ ) {
        if ( fscanf(fp, "%d %d %lg", &i, &j, &value) != EOF ) {
            Atmp.rowind[k]=i; Atmp.colind[k]=j; Atmp.val[k] = value;
        } else {
          check_error(ERROR_WRONG_FILE, __FUNCTION__);
        }
    }
    dcoo_2_dcsr(&Atmp,A);
    free((void *)Atmp.rowind);
}
/***********************************************************************************************/
/**
 * \fn void dvector_read_eof(FILE *fp, dvector *b)
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
void dvector_read_eof(FILE *fp,dvector *b)
{
    int  i, n,ichk;
    REAL value;
    n=0;
    while(!feof(fp)){
      ichk=fscanf(fp,"%lg", &value);
      if(ichk<0) break;
      n++;
    }
    fprintf(stdout,"%s: FOUND %d records.\n",__FUNCTION__,n);
    if(!n)
      check_error(ERROR_WRONG_FILE, __FUNCTION__);
    dvec_alloc(n,b);
    rewind(fp);
    for ( i = 0; i < n; ++i ) {
        fscanf(fp, "%lg", &value);
        b->val[i] = value;
        if ( value > BIGREAL ) {
            printf("### ERROR: Wrong value = %lf\n", value);
            dvec_free(b);
            fclose(fp);
            exit(ERROR_INPUT_PAR);
        }
    }
    return;
}
/*END OF SOLVERS_READ.h*/
