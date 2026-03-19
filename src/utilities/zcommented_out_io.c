/*! \file src/utilities/io_commented_out.c
 *
 * \brief Commented-out I/O functions preserved for reference.
 *        These are not compiled (wrapped in #if 0).
 *        Use sc_read_gmsh / sc_write_gmsh / sc_write_vtk instead.
 *
 */

#include "hazmath.h"

#if 0 /* hazw — old hazmath format writer, replaced by sc_write_gmsh */
/********************************************************************************/
/*!
 * \fn void hazw(char *nameout,scomplex *sc, const INT shift)
 *
 * \brief Write a simplicial complex to a file in a "hazmath" format.
 *
 * \param nameout   File name
 * \param sc        Pointer to a simplicial complex
 * \param shift     integer added to the elements of arrays such as sc->nodes.
 *
 * \note The data is organized as follows:
 *
 * 0. num_simplices,num_vertices,dimension,connected_components(bndry)-1;
 *
 * 1. As every simplex has (dim+1) vertices, the next
 *    (dim+1)*num_simplices integers are: 1st vertex for each simplex
 *    (num_simplices integers); 2nd vertex for every simplex, etc.
 *
 * 2. num_simplices integers are the flags (tags) associated with
 *    every element: could be something characterizing the element.
 *
 * 3. (dim)*num_vertices REALs (coordinates of the vertices) e.g.:
 *    num_vertices REALs with 1st coordinate, num_vertices REALs with
 *    2nd coordinate,...
 *
 * 4. num_vertices integers with tages (boundary codes) for every vertex.
 *
 */
/********************************************************************************/
void hazw(char *nameout,scomplex *sc, const INT shift)
{
  // WRITING in HAZMATH format.
  FILE *fmesh;
  INT n=sc->nv,ns=sc->ns, dim=sc->dim,ndl=sc->dim+1;
  INT *je = sc->nodes, *ib=sc->bndry;
  REAL *x = sc->x;
  INT k=-10,j=-10,kndl=-10;
  fmesh=HAZ_fopen(nameout,"w");
  fprintf(fmesh,"%lld %lld %lld %lld\n",(long long )ns,(long long )n,(long long )dim,(long long )(sc->bndry_cc-1));
  for (j=0;j<ndl;j++) {
    for (k=0;k<ns;k++){
      kndl=ndl*k+j;
      fprintf(fmesh," %lld ",(long long )(je[kndl]+shift));
    }
    fprintf(fmesh,"\n");
  }
  for (k=0;k<ns;k++){
    fprintf(fmesh," %lld ", (long long )sc->flags[k]);
  }
  fprintf(fmesh,"\n");
  for(j=0;j<dim;j++){
    for(k=0;k<n;k++){
      fprintf(fmesh," %23.16g ",x[k*dim+j]);
    }
    fprintf(fmesh,"\n");
  }
  for(k=0;k<n;k++){
    fprintf(fmesh," %lld ", (long long )ib[k]);
  }
  fprintf(fmesh,"\n");
  fprintf(stdout,"\n%%Output (hazmath) written on:%s\n",nameout);
  fclose(fmesh);
  return;
}
#endif /* hazw */
/********************************************************************************/

#if 0 /* haz_scomplex_read — old hazmath format reader, replaced by sc_read_gmsh */
/**********************************************************************/
/*!
 * \fn scomplex *haz_scomplex_read(FILE *fp, INT print_level)
 *
 * \brief Read a simplicial complex from a file in hazmath format.
 *
 */
scomplex* haz_scomplex_read(FILE* fp, INT print_level) {
  INT i, j, k;
  long long ns_, nv_, n_, nholes_;
  fscanf(fp, "%lld %lld %lld %lld", &ns_, &nv_, &n_, &nholes_);
  INT ns = (INT)ns_, nv = (INT)nv_, n = (INT)n_, nholes = (INT)nholes_;
  INT n1 = n + 1;
  INT nbig = n;
  scomplex* sc = (scomplex*)haz_scomplex_init(n, ns, nv, n);
  INT one_zero_flag = 1;
  long long readint;
  for (j = 0; j < n1; j++) {
    for (k = 0; k < ns; k++) {
      INT n1kj = n1 * k + j;
      fscanf(fp, " %lld ", &readint);
      sc->nodes[n1kj] = (INT)readint;
      if (sc->nodes[n1kj] == 0 && one_zero_flag == 1)
        one_zero_flag = 0;
    }
  }
  if (one_zero_flag == 1)
    for (i = 0; i < ns * n1; i++)
      sc->nodes[i] -= 1;
  for (k = 0; k < ns; k++) {
    fscanf(fp, " %lld ", (long long*)sc->flags + k);
  }
  for (j = 0; j < nbig; j++) {
    for (i = 0; i < nv; i++) {
      fscanf(fp, "%lg", sc->x + i * nbig + j);
    }
  }
  for (i = 0; i < nv; i++) {
    fscanf(fp, "%lld", (long long*)(sc->bndry + i));
  }
  sc->cc = 1;
  sc->bndry_cc = (nholes == 0) ? 1 : nholes + 1;
  sc->print_level = print_level;
  return sc;
}
#endif /* haz_scomplex_read */
/**********************************************************************/

#if 0 /* creategrid_fread — old wrapper, replaced by sc_read_gmsh + find_nbr + sc_vols + sc_build_fem_data */
/**********************************************************************/
/*!
* \fn scomplex* creategrid_fread(FILE *gfid,INT file_type)
*
* \brief Creates grid by reading in from file, returns an scomplex.
*
* \param gfid      Grid FILE ID
* \param file_type Type of File Input: 0 - haz format
*
* \return scomplex* with the mesh and FEM data.
*
*/
scomplex* creategrid_fread(FILE *gfid,INT file_type)
{
  if(file_type!=0) {
    fprintf(stderr,"Unknown mesh file type, %lld. Try using native (.haz) format. -Exiting\n",(long long )file_type);
    exit(255);
  }
  scomplex *sc = haz_scomplex_read(gfid, 0);
  fprintf(stdout,"reading complete...\n");fflush(stdout);
  find_nbr(sc->ns, sc->nv, sc->dim, sc->nodes, sc->nbr);
  sc_vols(sc);
  sc_build_fem_data(sc);
  return sc;
}
#endif /* creategrid_fread */
/**********************************************************************/
