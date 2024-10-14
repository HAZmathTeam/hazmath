/*! \file src/mesh/mesh_io.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 1/9/15.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \brief Obtains routines for reading in meshes via original format and vtk format.
 *
 *   \note updated by James Adler 07/25/2018
 *
 */

#include "hazmath.h"

/******************************************************************************/
/*!
 * \fn void read_grid_haz(FILE *gfid,mesh_struct *mesh)
 *
 * \brief Reads in gridfile of the following form:
 *
 *        First line:             nelms nnodes dim nholes
 *        Next two-four lines:	Each element will have 2 vertices in 1D, 3 in 2D, and 4 in 3D.
 *                                Columns are node-element map not sparse
 *
 *	      Line 1:	First node of all elements
 *	      Line 2: 	Second node of all elements
 *	      Line 3:   Third node of all elements (only 2D & 3D)
 *          Line 4:	4th node for 3D
 *
 *        Next one, two, or 3 lines:	x,y,z coordinates
 *
 *          Line 1:	x coordinates of all nodes
 *		  Line 2:	y coordinates of all nodes (only in 2 or 3D)
 *		  Line 3:	z coordinates of all nodes (only in 3D)
 *
 *        Next line:  List of nodes and whether they land on boundary
 *                    Binary or Flagged (i.e., -1 indicates boundary of hole)
 *
 * \param gfid             Grid FILE ID.
 * \param mesh             Pointer to mesh struct
 *
 * \return mesh.dim        Dimension of problem
 * \return mesh.nelm       Number of elements in mesh
 * \return mesh.nvert      Number of vertices in mesh
 * \return mesh.nbvert     Number of boundary vertices in mesh
 * \return mesh.v_per_elm  Number of vertices per element
 * \return mesh.cv         Coordinates of mesh vertices
 * \return mesh.el_v       Element to vertex map
 *
 * \note This code will check if the file starts counting nodes at 0 or 1.
 *       We assume that if there is a hole, the 0th or 1st node is not eliminated
 *
 */

void read_grid_haz(FILE *gfid,mesh_struct *mesh)
{

  // Loop indices
  INT i,j,k;

  // Flag to check if file started counting at 1 or 0.
  INT one_zero_flag = 1;

  // Get basic data
  long long nelm_,nv_,dim_,nholes_;
  //fscanf(gfid,"%lld %lld %lld %lld",  (long long *)&nelm,  (long long *)&nv,  (long long *)&dim,  (long long *)&nholes);
  fscanf(gfid,"%lld %lld %lld %lld", &nelm_,  &nv_,  &dim_,  &nholes_);
  INT nelm = (INT )nelm_,nv=(INT) nv_,dim=(INT) dim_,nholes=(INT) nholes_;
  printf("%d\t%d\t%d\n",nelm,nv,dim);
  // exit(666);
  // Get number of vertices per element
  INT v_per_elm = dim+1;

  long long readint;

  // Element-Vertex Map
  mesh->el_v=malloc(sizeof(iCSRmat));
  mesh->el_v->row=nelm;
  mesh->el_v->col=nv;
  mesh->el_v->nnz=nelm*v_per_elm;
  mesh->el_v->IA = (INT *)calloc(nelm+1, sizeof(INT));
  mesh->el_v->JA = (INT *)calloc(mesh->el_v->nnz, sizeof(INT));
  for(i=0;i<nelm+1;i++) {
    mesh->el_v->IA[i] = v_per_elm*i;
  }
  for (i=0;i<v_per_elm;i++) {
    for (j=0;j<nelm;j++){
      k=v_per_elm*j+i;
    //  fscanf(gfid,"%lld",   (long long *)(mesh->el_v->JA+k));
      fscanf(gfid,"%lld",   &readint);
      mesh->el_v->JA[k] = (INT ) readint;
      if(mesh->el_v->JA[k]==0 && one_zero_flag==1)
        one_zero_flag = 0;
    }
  }
  if(one_zero_flag==1)
    for(i=0;i<mesh->el_v->nnz;i++)
      mesh->el_v->JA[i]-=1;

  mesh->el_v->val=NULL;

  // Read in Flags for each element
  mesh->el_flag = (INT *) calloc(nelm,sizeof(INT));
  rveci_(gfid,mesh->el_flag,&nelm);

  // Get next 2-3 lines for coordinate map
  mesh->cv = allocatecoords(nv,dim);
  INT nvdim=nv*dim;
  rvecd_(gfid,mesh->cv->x,&nvdim); // This actually reads in all dimensions

  // Get next 1-2 lines for boundary flags
  INT nbv = 0;
  mesh->v_flag = (INT *) calloc(nv,sizeof(INT));
  rveci_(gfid,mesh->v_flag,&nv);
  for(i=0;i<nv;i++) {
    if(mesh->v_flag[i]>0) {
      nbv++;
    }
  }

  // Check for holes
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
  mesh->v_per_elm = v_per_elm;

  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void read_grid_vtk(FILE *gfid,mesh_struct *mesh)
 *
 * \brief Reads in gridfile in vtk (really vtu) format.  Example follows:
 *
 *        <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
 *        <UnstructuredGrid>
 *        <Piece NumberOfPoints="9" NumberOfCells="8">
 *        <Points>
 *        <DataArray type="Float64" NumberOfComponents="3" Format="ascii">
 *        0                        0                        0
 *        0.5                      0                        0
 *        1                        0                        0
 *        0                      0.5                        0
 *        0.5                    0.5                        0
 *        1                      0.5                        0
 *        0                        1                        0
 *        0.5                      1                        0
 *        1                        1                        0
 *        </DataArray>
 *        </Points>
 *        <PointData Scalars="scalars">
 *        <DataArray type="Int64" Name="v_flag" Format="ascii">
 *        1  1  1  1  0  1  1  1  1
 *        </DataArray>
 *        <DataArray type="Int64" Name="connectedcomponents" Format="ascii">
 *        -1  -1  -1  -1  1  -1  -1  -1  -1
 *        </DataArray>
 *        </PointData>
 *        <Cells>
 *        <DataArray type="Int64" Name="offsets" Format="ascii">
 *        3  6  9  12  15  18  21  24
 *        </DataArray>
 *        <DataArray type="Int64" Name="connectivity" Format="ascii">
 *        0  1  3  3  1  4  1  2  4  4  2  5  3  4  6  6  4  7  4  5  7  7  5  8
 *        </DataArray>
 *        <DataArray type="Int64" Name="types" Format="ascii">
 *        5  5  5  5  5  5  5  5
 *        </DataArray>
 *        </Cells>
 *        </Piece>
 *        </UnstructuredGrid>
 *        </VTKFile>
 *
 * \param gfid             Grid FILE ID.
 * \param mesh             Pointer to mesh struct
 *
 * \return mesh.dim        Dimension of problem
 * \return mesh.nelm       Number of elements in mesh
 * \return mesh.nvert      Number of vertices in mesh
 * \return mesh.nbvert     Number of boundary vertices in mesh
 * \return mesh.v_per_elm  Number of vertices per element
 * \return mesh.cv         Coordinates of mesh vertices
 * \return mesh.el_v       Element to vertex map
 *
 */
void read_grid_vtk(FILE *gfid,mesh_struct *mesh)
{
  fprintf(stderr,"I don't really want to process this vtk file just yet...Come back later...\n\n");
  exit(ERROR_OPEN_FILE);

  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void dump_mesh_haz(char *namehaz,mesh_struct *mesh)
 *
 * \brief Dumps mesh data to haz format
 *
 *        First line:             nelms nnodes dim nholes
 *        Next two-four lines:	Each element will have 2 vertices in 1D, 3 in 2D, and 4 in 3D.
 *                                Columns are node-element map not sparse
 *
 *	      Line 1:	First node of all elements
 *	      Line 2: 	Second node of all elements
 *	      Line 3:   Third node of all elements (only 2D & 3D)
 *          Line 4:	4th node for 3D
 *
 *        Next one, two, or 3 lines:	x,y,z coordinates
 *
 *          Line 1:	x coordinates of all nodes
 *		  Line 2:	y coordinates of all nodes (only in 2 or 3D)
 *		  Line 3:	z coordinates of all nodes (only in 3D)
 *
 *        Next line:  List of nodes and whether they land on boundary
 *                    Binary or Flagged (i.e., -1 indicates boundary of hole)
 *
 * \param namehaz          Output file name
 * \param mesh             Pointer to mesh struct
 *
 * \return namehaz.haz     File with mesh data.
 *
 * \note Dumps with numbering starting at 0.
 *
 */
void dump_mesh_haz(char *namehaz,mesh_struct *mesh)
{
  // Basic Quantities
  INT i,j,ja,jb,jcnt;
  INT nv = mesh->nv;
  INT nelm = mesh->nelm;
  INT dim = mesh->dim;
  INT v_per_elm = mesh->v_per_elm;

  // Open File for Writing
  FILE* fhaz = HAZ_fopen(namehaz,"w");

  // Write First Line
  fprintf(fhaz,"%lld %lld %lld %lld\n",				\
	  (long long )nelm,(long long )nv,(long long )dim,	\
	  (long long )mesh->nconn_bdry-(long long )mesh->nconn_reg);

  // Write Next 2-4 lines for element->vertex map
  INT* nodes = (INT *) calloc(nelm*v_per_elm,sizeof(INT));
  for(i=0;i<nelm;i++) {
    ja = mesh->el_v->IA[i];
    jb = mesh->el_v->IA[i+1];
    jcnt=0;
    for(j=ja;j<jb;j++) {
      nodes[i*v_per_elm+jcnt] = mesh->el_v->JA[j];
      jcnt++;
    }
    // Print first node of elements
    fprintf(fhaz," %lld ",  (long long )nodes[i*v_per_elm]);
  }
  fprintf(fhaz,"\n");

  // Print Second Node
  for(i=0;i<nelm;i++) {
    fprintf(fhaz," %lld ",  (long long )nodes[i*v_per_elm+1]);
  }
  fprintf(fhaz,"\n");
  // Print Third if in 2 or 3D
  if(dim==2 || dim==3) {
    for(i=0;i<nelm;i++) {
      fprintf(fhaz,"%lld ",(long long )nodes[i*v_per_elm+2]);
    }
    fprintf(fhaz,"\n");
  }
  // Print Fourth node if in 3D
  if(dim==3) {
    for(i=0;i<nelm;i++) {
      fprintf(fhaz," %lld ",  (long long )nodes[i*v_per_elm+3]);
    }
    fprintf(fhaz,"\n");
  }

  // Dump coordinates
  // x
  for(i=0;i<nv;i++) {
    fprintf(fhaz,"%23.16e ",mesh->cv->x[i]);
  }
  fprintf(fhaz,"\n");
  //y
  if(dim == 2 || dim==3) {
    for(i=0;i<nv;i++) {
      fprintf(fhaz,"%23.16e ",mesh->cv->y[i]);
    }
    fprintf(fhaz,"\n");
  }
  if(dim==3) {
    for(i=0;i<nv;i++) {
      fprintf(fhaz,"%23.16e ",mesh->cv->z[i]);
    }
    fprintf(fhaz,"\n");
  }

  free(nodes);

  // Dump v_flag Data to indicate if vertices are boundaries
  for(i=0;i<nv;i++) fprintf(fhaz,"%lld ",  (long long )mesh->v_flag[i]);
  fprintf(fhaz,"\n");

  fclose(fhaz);

  return;
}
/******************************************************************************/

/******************************************************************************/
/*!
 * \fn void dump_mesh_vtk(char *namevtk,mesh_struct *mesh)
 *
 * \brief Dumps mesh data to vtk format. Example follows:
 *
 *        <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
 *        <UnstructuredGrid>
 *        <Piece NumberOfPoints="9" NumberOfCells="8">
 *        <Points>
 *        <DataArray type="Float64" NumberOfComponents="3" Format="ascii">
 *        0                        0                        0
 *        0.5                      0                        0
 *        1                        0                        0
 *        0                      0.5                        0
 *        0.5                    0.5                        0
 *        1                      0.5                        0
 *        0                        1                        0
 *        0.5                      1                        0
 *        1                        1                        0
 *        </DataArray>
 *        </Points>
 *        <PointData Scalars="scalars">
 *        <DataArray type="Int64" Name="v_flag" Format="ascii">
 *        1  1  1  1  0  1  1  1  1
 *        </DataArray>
 *        <DataArray type="Int64" Name="connectedcomponents" Format="ascii">
 *        -1  -1  -1  -1  1  -1  -1  -1  -1
 *        </DataArray>
 *        </PointData>
 *        <Cells>
 *        <DataArray type="Int64" Name="offsets" Format="ascii">
 *        3  6  9  12  15  18  21  24
 *        </DataArray>
 *        <DataArray type="Int64" Name="connectivity" Format="ascii">
 *        0  1  3  3  1  4  1  2  4  4  2  5  3  4  6  6  4  7  4  5  7  7  5  8
 *        </DataArray>
 *        <DataArray type="Int64" Name="types" Format="ascii">
 *        5  5  5  5  5  5  5  5
 *        </DataArray>
 *        </Cells>
 *        </Piece>
 *        </UnstructuredGrid>
 *        </VTKFile>
 *
 * \param namevtk          Output file name
 * \param mesh             Pointer to mesh struct
 *
 * \return namevtk.vtu     File with mesh data.
 *
 */
void dump_mesh_vtk(char *namevtk,mesh_struct *mesh)
{
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
  const INT LINE=3;
  const INT TRI=5;
  const INT TET=10;

  if(dim==1) /* Line */
    tcell=LINE;
  else if(dim==2)
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
  fprintf(fvtk,"<Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",  (long long )nv,  (long long )nelm);
  fprintf(fvtk,"<Points>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" NumberOfComponents=\"3\" Format=\"ascii\">", \
          tfloat);

  // Dump coordinates
  if(dim == 1) {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],0e0,0e0);
    }
  } else if(dim == 2) {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],mesh->cv->y[k],0e0);
    }
  } else {
    for(k=0;k<nv;k++) {
      fprintf(fvtk," %23.16e %23.16e %23.16e ",mesh->cv->x[k],mesh->cv->y[k], \
              mesh->cv->z[k]);
    }
  }
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"</Points>\n");

  // Dump v_flag Data to indicate if vertices are boundaries
  fprintf(fvtk,"<PointData Scalars=\"scalars\">\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"v_flag\" Format=\"ascii\">",tinto);
  for(k=0;k<nv;k++) fprintf(fvtk," %lld ",  (long long )mesh->v_flag[k]);
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
    for(k=0;k<nv;k++) fprintf(fvtk," %lld ",  (long long )mesh->v_component[k]);
    fprintf(fvtk,"</DataArray>\n");
  }
  fprintf(fvtk,"</PointData>\n");

  // Dump el_v map
  fprintf(fvtk,"<Cells>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"offsets\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++) fprintf(fvtk," %lld ",  (long long )mesh->el_v->IA[k]);
  fprintf(fvtk,"</DataArray>\n");
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"connectivity\" Format=\"ascii\">\n",tinto);
  for(k=0;k<nelm;k++){
    kndl=k*v_per_elm;
    for(j=0;j<v_per_elm;j++) fprintf(fvtk," %lld ",  (long long )mesh->el_v->JA[kndl + j]);
  }
  fprintf(fvtk,"</DataArray>\n");

  // Dump Element Type
  fprintf(fvtk,"<DataArray type=\"%s\" Name=\"types\" Format=\"ascii\">",tinto);
  for(k=1;k<=nelm;k++)
    fprintf(fvtk," %lld ",  (long long )tcell);

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

// Create Mesh from Scratch Routines
/******************************************************************************/
/*!
 * \fn void create1Dgrid_Line(mesh_struct* mesh,REAL left_end,REAL right_end,INT nelm)
 *
 * \brief Creates a 1D grid from scratch [left_end,right_end] with no holes.
 *
 * \param left_end    Coordinate of Left-End Point
 * \param right_end   Coordinate of Right-End Point
 * \param nelm        Number of Elements
 *
 * \return mesh       Mesh struct and all its properties.
 *
 */
void create1Dgrid_Line(mesh_struct* mesh,REAL left_end,REAL right_end,INT nelm)
{
  // Initialize mesh for read in.
  initialize_mesh(mesh);
  mesh->el_v = malloc(sizeof(struct iCSRmat));

  INT i,j; /* Loop indices */

  INT dim = 1;
  INT nv = nelm+1;

  // Get number of vertices per element
  INT v_per_elm = 2;

  // Get h-spacing
  REAL h = 1.0/nelm;

  // Allocate arrays for element-vertex map,coordinates, and boundaries
  iCSRmat el_v = icsr_create(nelm,nv,v_per_elm*nelm);
  coordinates *cv = allocatecoords(nv,dim);
  INT* v_flag = (INT *) calloc(nv,sizeof(INT));

  // Build element-vertex map, coordinates, and v_flag
  j=0;
  for(i=0;i<nelm;i++) {
    el_v.IA[i]=i*v_per_elm;
    el_v.JA[j]=i;
    el_v.JA[j+1]=i+1;
    j=j+2;
    cv->x[i]=left_end+h*i;
    v_flag[i]=0;
  }
  el_v.IA[nelm] = nelm*v_per_elm;
  cv->x[nelm] = right_end;
  v_flag[0] = 1;
  v_flag[nelm] = 1;
  INT nbv = 2;

  // Assume no holes
  INT nconn_reg = 1;
  INT nconn_bdry = 1;

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
  mesh->v_flag = v_flag;

  // Build rest of mesh
  build_mesh_all(mesh);

  return;
}
/******************************************************************************/
