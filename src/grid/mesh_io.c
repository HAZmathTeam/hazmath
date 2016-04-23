/*
 *  mesh_input.c   
 *  
 *  Created by James Adler and Xiaozhe Hu on 1/9/15.
 *  Copyright 2015__HAZMAT__. All rights reserved.
 *
 *  Obtains routines for reading in meshes via original format and vtk format.
 *
 */

#include "hazmat.h"

/******************************************************************************/
void read_grid_old(FILE *gfid,trimesh *mesh) 
{
  /* Reads in gridfile of the following form:
   *
   * First line:	nelms nnodes dim nholes
   * Next four lines:	Each element will have 3 vertices in 2D and 4 in 3D
   * 
   *	      line 1:	First node of all elements
   *	      line 2: 	Second node of all elements
   *	      line 3:   Third node of all elements
   *          line 4:	Dummy line in 2D or 4th node for 3D
   *			Columns are node-element map not sparse
   *
   * Next two (or 3 in 3D) lines:	x,y,z coordinates
   *			line 1:		x coordinates of all nodes
   *			line 2:		y coordinates of all nodes
   *			line 3:	        z coordinates of all nodes (only in 3D)
   *
   * Next line:         List of nodes and whether they land on boundary 
   *                    (binary) -1 indicates boundary of hole if any
   *
   * INPUT: 
   *         gfid    Grid File ID
   *
   * OUTPUT: 
   *         mesh    Properties of mesh:
   *                  dim, nelm, nvert, nbvert, v_per_elm, 
   *                  coordinates, and element to vertex map
   *
   */
	
  INT i,j,k; /* Loop indices */

  // Get Number of elements, nodes and boundary edges first
  INT* line1 = calloc(4,sizeof(INT));
  INT lenhead = 4;
  rveci_(gfid,line1,&lenhead);
  INT dim = line1[2];
  INT nelm = line1[0];
  INT nv = line1[1];
  INT nholes = line1[3];
  free(line1);

  // Get number of vertices per element
  INT v_per_elm = dim+1;
			
  // Allocate arrays to read in other data such as coordinate information
  INT* element_vertex = (INT *) calloc(nelm*v_per_elm,sizeof(INT));
  coordinates *cv = allocatecoords(nv,dim);
		
  // Get next 3-4 lines Element-Vertex Map
  INT* line2 = (INT *) calloc(nelm,sizeof(INT));
  for (i=0; i<v_per_elm; i++) {
    rveci_(gfid,line2,&nelm);
    for (j=0; j<nelm; j++) {
      element_vertex[j*v_per_elm+i] = line2[j];
    }
  }
  free(line2);
	
  // Get next 2-3 lines for coordinate map
  rvecd_(gfid,cv->x,&nv);
  rvecd_(gfid,cv->y,&nv);
  if(dim==3)
    rvecd_(gfid,cv->z,&nv);
	
  // Get next 1-2 lines for boundary nodes
  INT nbv = 0;
  INT* v_bdry = (INT *) calloc(nv,sizeof(INT));
  rveci_(gfid,v_bdry,&nv);
  for(i=0;i<nv;i++) {
    if(v_bdry[i]) {
      nbv++;
    }
  }
        
  /* Element Vertex Map */
  iCSRmat el_v = convert_elmnode(element_vertex,nelm,nv,v_per_elm);
  if(element_vertex) free(element_vertex);

  // This format only allows for 1 connected region and 
  // up to 2 connected boundaries (1 hole).
  INT nconn_reg = 0;
  INT nconn_bdry = 0;
  if(nholes==0) {
    nconn_reg = 1;
    nconn_bdry = 1;
  } else if(nholes==1) {
    nconn_reg = 1;
    nconn_bdry = 2;
  }

  // Update mesh with known quantities
  mesh->dim = dim;
  mesh->nelm = nelm;
  mesh->nv = nv;
  mesh->nbv = nbv;
  mesh->nconn_reg = nconn_reg;
  mesh->nconn_bdry = nconn_bdry;
  mesh->cv = cv;
  *(mesh->el_v) = el_v;
  mesh->v_per_elm = v_per_elm;
  mesh->v_bdry = v_bdry;

  return;
}
/******************************************************************************/

/******************************************************************************/
void read_grid_vtk(FILE *gfid,trimesh *mesh) 
{
  /* Reads in gridfile in vtk (really vtu) format.  Example follows:
   *
   * <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
   * <UnstructuredGrid>
   * <Piece NumberOfPoints="9" NumberOfCells="8">
   * <Points>
   * <DataArray type="Float64" NumberOfComponents="3" Format="ascii">                        
   * 0                        0                        0                                            
   * 0.5                      0                        0                      
   * 1                        0                        0                        
   * 0                      0.5                        0                      
   * 0.5                    0.5                        0                      
   * 1                      0.5                        0                         
   * 0                        1                        0                      
   * 0.5                      1                        0                      
   * 1                        1                        0 
   * </DataArray>
   * </Points>
   * <PointData Scalars="scalars">
   * <DataArray type="Int64" Name="v_bdry" Format="ascii"> 
   * 1  1  1  1  0  1  1  1  1  
   * </DataArray>
   * <DataArray type="Int64" Name="connectedcomponents" Format="ascii"> 
   * -1  -1  -1  -1  1  -1  -1  -1  -1  
   * </DataArray>
   * </PointData>
   * <Cells>
   * <DataArray type="Int64" Name="offsets" Format="ascii"> 
   * 3  6  9  12  15  18  21  24 
   * </DataArray>
   * <DataArray type="Int64" Name="connectivity" Format="ascii">
   * 0  1  3  3  1  4  1  2  4  4  2  5  3  4  6  6  4  7  4  5  7  7  5  8 
   * </DataArray>
   * <DataArray type="Int64" Name="types" Format="ascii"> 
   * 5  5  5  5  5  5  5  5 
   * </DataArray>
   * </Cells>
   * </Piece>
   * </UnstructuredGrid>
   * </VTKFile>

   *
   * INPUT: 
   *         gfid    Grid File ID
   *
   * OUTPUT: 
   *         mesh    Properties of mesh:
   *                  dim, nelm, nvert, nbvert, v_per_elm, 
   *                  coordinates, and element to vertex map
   *
   */
	
  fprintf(stderr,"I don't really want to process this vtk file just yet...Come back later...\n\n");
  exit(255);

  return;
}
/******************************************************************************/


/******************************************************************************/
void dump_mesh_vtk(char *namevtk,trimesh *mesh)
{

  /* Dumps mesh data to vtk format 
  *
  * Input:
  *  mesh:     Mesh struct to dump
  * 
  * Output:
  *  namevtk  File name of vtk file
  *
  */

  // Basic Quantities
  INT nv = mesh->nv;
  INT nelm = mesh->nelm;
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;

  // VTK needed Quantities
  INT tcell=-10;
  INT k=-10,j=-10,kndl=-10;
  char *tfloat="Float64", *tinto="Int64", *endian="LittleEndian";
   
  /* 
     What endian?:

     Intel x86; OS=MAC OS X: little-endian
     Intel x86; OS=Windows: little-endian
     Intel x86; OS=Linux: little-endian
     Intel x86; OS=Solaris: little-endian
     Dec Alpha; OS=Digital Unix: little-endian
     Dec Alpha; OS=VMS: little-endian
     Hewlett Packard PA-RISC; OS=HP-UX: big-endian
     IBM RS/6000; OS=AIX: big-endian
     Motorola PowerPC; OS=Mac OS X:  big-endian
     SGI R4000 and up; OS=IRIX: big-endian
     Sun SPARC; OS=Solaris: big-endian
  */

  /* 
     Types of cells for VTK 

     VTK_VERTEX (=1) 
     VTK_POLY_VERTEX (=2)
     VTK_LINE (=3)
     VTK_POLY_LINE (=4)
     VTK_TRIANGLE(=5)
     VTK_TRIANGLE_STRIP (=6)
     VTK_POLYGON (=7) 
     VTK_PIXEL (=8) 
     VTK_QUAD (=9)
     VTK_TETRA (=10)
     VTK_VOXEL (=11)
     VTK_HEXAHEDRON (=12)
     VTK_WEDGE (=13) 
     VTK_PYRAMID (=14)
  */
  const INT TRI=5;  
  const INT TET=10;
  
  if(dim==2) 
    tcell=TRI; /* triangle */
  else 
    tcell=TET; /* tet */

  // Open File for Writing
  FILE* fvtk = HAZ_fopen(namevtk,"w");  

  // Write Headers
  fprintf(fvtk, \
	  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n", \
	  endian);
  fprintf(fvtk,"<UnstructuredGrid>\n");
  fprintf(fvtk,"<Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n",nv,nelm);
  fprintf(fvtk,"<Points>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" NumberOfComponents=\"3\" Format=\"ascii\">", \
  	  tfloat);

  // Dump coordinates
  if(dim == 2) {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %24.16g %24.16g %24.16g ",mesh->cv->x[k],mesh->cv->y[k],0e0);
    }
  } else {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %24.16g %24.16g %24.16g ",mesh->cv->x[k],mesh->cv->y[k], \
  	      mesh->cv->z[k]);
    }
  }
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Points>\n");

  // Dump v_bdry Data to indicate if vertices are boundaries
  fprintf(fvtk,"<PointData Scalars=\"scalars\">\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"v_bdry\" Format=\"ascii\">",tinto);
  for(k=0;k<nv;k++) fprintf(fvtk," %i ",mesh->v_bdry[k]);
  fprintf(fvtk,"</DataArray>\n");

  // Dump information about connected components.
  // Positive integers indicate connected components of a domain
  // Negative integers indicate connected components of the boundaries
  // Example: A cube (1 connected domain and 1 connected boundary)
  //            would be 1 on the interior and -1 on points on the boundary
  //          A cube with a hole (1 connected domain and 2 connected boundaries)
  //          would have 1 on the points in the interior and
  //          -1 on points on the outer boundary and -2 on the inner boundary
  // If NULL, then one connected region and boundary.
  if(mesh->v_component) {
    fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectedcomponents\" Format=\"ascii\">",tinto);
    for(k=0;k<nv;k++) fprintf(fvtk," %i ",mesh->v_component[k]);
    fprintf(fvtk,"</DataArray>\n");
  }
  fprintf(fvtk,"</PointData>\n");

  // Dump el_v map
  fprintf(fvtk,"<Cells>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"offsets\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++) fprintf(fvtk," %i ",mesh->el_v->IA[k]-1);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  for(k=0;k<nelm;k++){
    kndl=k*v_per_elm;
    for(j=0;j<v_per_elm;j++) fprintf(fvtk," %i ",mesh->el_v->JA[kndl + j]-1);
  }
  fprintf(fvtk,"</DataArray>\n");

  // Dump Element Type
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"types\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++)
    fprintf(fvtk," %i ",tcell);

  // Put in remaining headers
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Cells>\n");
  fprintf(fvtk,"</Piece>\n");
  fprintf(fvtk,"</UnstructuredGrid>\n");
  fprintf(fvtk,"</VTKFile>\n");

  fclose(fvtk);

  return;
}
/******************************************************************************/

