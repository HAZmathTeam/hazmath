/*
 *  creategrid.c
 *  
 *  Created by James Adler and Xiaozhe Hu on 1/9/15.
 *  Copyright 2015_JXLcode__. All rights reserved.
 *
 */

// Standard Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Our Includes
#include "grid.h"

/**************************************************************************************************************************************************/
/********* Read in arbritray grid and create FEM structure for the grid ***************************************************************************/
/**************************************************************************************************************************************************/
void creategrid(FILE *gfid,INT *element_node,REAL *xn,REAL *yn,INT *bdry_n,INT nelm, INT n, INT nbedge,INT element_order) 
{
	/* Reads in gridfile of the following form:
	 *
	 * First line:				nelms nnodes nboundaryedges nboundaryedges (This will be read outside this routine)
	 * Next four lines:		       	Each element will have 3 vertices in 2D and 4 in 3D
	 *			line 1:	       	First node of all elements
	 *			line 2:		Second node of all elements
	 *			line 3:	        Third node of all elements
	 *			line 4:		Dummy line (ignore) (or 4th node for 3D)
	 *					Columns are node-element map not sparse
	 * Next two (or 3 in 3D) lines:		x,y,z coordinates
	 *			line 1:					x coordinates of all nodes
	 *			line 2:					y coordinates of all nodes
	 *			(line 3:				z coordinates of all nodes (only in 3D))
	 * Next two lines:					List of nodes that land on boundary
	 *			line 1:					Corresponds to first part of edge on bondary
	 *			line 2:					Corresponds to second part of edge on boundary
	 *									Columns of this are edge_node map of boundary
	 * Rest of lines are dummy and not used
	 *
	 * Outputs: element to node map (not sparse)
	 *			node to coordinate map xn,yn,zn
	 *			Boundary Edge to Node Correspondence bdry_n(nbedge,2)
	 *
	 */
	
	INT i,j; /* Loop indices */
		
	// Get next 3 lines Element-Node Map
	INT* line2 = calloc(nelm,sizeof(INT));
	for (i=0; i<element_order; i++) {
		rveci_(gfid,line2,&nelm);
		for (j=0; j<nelm; j++) {
			element_node[j*element_order+i] = line2[j];
		}
	}
	// Next line is dummy
	rveci_(gfid,line2,&nelm);
	free(line2);
	
	// Get next 2 lines X and Y Node map
	rvecd_(gfid,xn,&n);
	rvecd_(gfid,yn,&n);
	
	// Get next 2 lines for boundary nodes
	INT* line3 = calloc(nbedge,sizeof(INT));
	INT* line4 = calloc(nbedge,sizeof(INT));
	rveci_(gfid,line3,&nbedge);
	rveci_(gfid,line4,&nbedge);
	for (i=0; i<nbedge; i++) {
		bdry_n[i*2+0] = line3[i];
		bdry_n[i*2+1] = line4[i];
	}
	free(line3);
	free(line4);
	
	return;
}
/**************************************************************************************************************************************************/

/**************************************************************************************************************************************************/
/********* Read in arbritray grid (3D Triangulation) **********************************************************************************************/
/**************************************************************************************************************************************************/
void readgrid3D(FILE *gfid,INT *element_node,REAL *xn,REAL *yn,REAL *zn,INT *n_bdry,INT nelm, INT n,INT element_order,INT *nbvert) 
{
	/* Reads in gridfile of the following form:
	 *
	 * First line:						nelms nnodes nboundaryedges nboundaryedges (This will be read outside this routine)
	 * Next four lines:					Each element will have 3 vertices in 2D and 4 in 3D
	 *			line 1:					First node of all elements
	 *			line 2:					Second node of all elements
	 *			line 3:					Third node of all elements
	 *			line 4:					Dummy line (ignore) (or 4th node for 3D)
	 *									Columns are node-element map not sparse
	 * Next two (or 3 in 3D) lines:		x,y,z coordinates
	 *			line 1:					x coordinates of all nodes
	 *			line 2:					y coordinates of all nodes
	 *			line 3:					z coordinates of all nodes (only in 3D)
	 * Next	line:						List of nodes and whether they land on boundary
	 *
	 * Rest of lines are dummy and not used
	 *
	 * Outputs: element to node map (not sparse)
	 *			node to coordinate map xn,yn,zn
	 *			Boundary Edge to Node Correspondence bdry_n(nbedge,2)
	 *			nbvert Number of Nodes on Boundaries
	 *
	 */
	
	INT i,j,cnt; /* Loop indices */
	
	// Get next 4 lines Element-Node Map
	INT* line2 = calloc(nelm,sizeof(INT));
	for (i=0; i<element_order; i++) {
		rveci_(gfid,line2,&nelm);
		for (j=0; j<nelm; j++) {
			element_node[j*element_order+i] = line2[j];
		}
	}
	free(line2);

	// Get next 3 lines X and Y and Z Node map
	rvecd_(gfid,xn,&n);
	rvecd_(gfid,yn,&n);
	rvecd_(gfid,zn,&n);
	
	// Get next line for boundary nodes
	INT* line3 = calloc(n,sizeof(INT));
	cnt = 0;
	rveci_(gfid,line3,&n);
	for (i=0; i<n; i++) {
		n_bdry[i] = line3[i];
		if(n_bdry[i]) {
			cnt++;
		}
	}
	*nbvert = cnt;
	free(line3);
	return;
}
/**************************************************************************************************************************************************/

/**************************************************************************************************************************************************/
/******** Make Default 2D Grid ***********************************************************************************************************************/
/**************************************************************************************************************************************************/
void std_tri_grid(INT *element_node,REAL *xn,REAL *yn,INT nex,INT ney,REAL xL,REAL xR,REAL yB,REAL yT,INT *bdry_n) 
{
	
	INT i,j; /* Loop indices */
	
	/* Creates the defualt triangular mesh using T3 elements and edges
	 
	 *    9----10----11----12
	 *    |\ 8  |\10  |\12  |
	 *    | \   | \   | \   |
	 *    |  \  |  \  |  \  |
	 *    |   \ |   \ |   \ |
	 *    |  7 \|  9 \| 11 \|
	 *    5-e11-6-e12-7-e13-8
	 *    |\ 2  |\ 4  |\ 6  |
	 *    | \   | \   | \   |
	 *    e4 e5 e6 e7 e8 e9 e10
	 *    |   \ |   \ |   \ |
	 *    |  1 \|  3 \|  5 \|
	 *    1-e1--2-e2--3-e3--4
	 *
	 *    Input:   nex,ney   Number of trigangle elements in x and y directions
	 *             xL,xR     Endpoints of domain in x
	 *             yB,yT     Endpoints of domain in y
	 *
	 *    Output:  element_node(elms,nodes)    Element to Node Map size (nelm,nodes)
	 *             xn,yn(nodes)					Coordinates of Nodes
	 *			   bdry_n(nbedge,2)			   Nodes on given edge that lies on boundary
	 */
	
	/******* Define grid *******************************************************/
	
	// Get number of nodes in each direction and total
	INT nx = nex/2 + 1;
	INT ny = ney/2 + 1;
	//INT n = nx*ny;
	
	// Get total number of elements
	//INT nelm = 2*(nex/2)*(ney/2);
	//INT nedge = 3*nelm/2 + nex/2 + ney/2;
	
	/********* Node Matrix *********/	
	for (j=0;j<ny;j++) {
		for (i=0;i<nx;i++) {
			xn[i+j*nx] = xL + i*(xR-xL)/(nx-1);
			yn[i+j*nx] = yB + j*(yT-yB)/(ny-1);
		}
	}
	
	/********** Element-Node Matrix, Element-Edge Matrix, and Edge-Node Matrix */
	INT k = 0;
	INT sw,se,nw,ne,south,west,mid,east,north;
	
	INT m=0;
	for (j=1;j<=ney/2;j++) {
		for (i=1;i<=nex/2;i++) {
			
			sw = i+(j-1)*nx;
			se = sw+1;
			nw = i+j*nx;
			ne = nw+1;
	
			south = i+(j-1)*(1.5*nex+1);
			west = south+(nex/2)+(i-1);
			mid = west+1;
			east = mid+1;
			north = i+j*(1.5*nex+1);
	
			element_node[k*3 + 0] = sw;
			element_node[k*3 + 1] = se;
			element_node[k*3 + 2] = nw;
			element_node[(k+1)*3 + 0] = ne;
			element_node[(k+1)*3 + 1] = nw;
			element_node[(k+1)*3 + 2] = se;
	
			//element_edge[k*3 + 0] = south;
			//element_edge[k*3 + 1] = mid;
			//element_edge[k*3 + 2] = west;
			//element_edge[(k+1)*3 + 0] = north;
			//element_edge[(k+1)*3 + 1] = mid;
			//element_edge[(k+1)*3 + 2] = east;
	
			k = k+2;
			
			// Boundaries
			
			if(j==1){
				bdry_n[m*2+0] = sw;
				bdry_n[m*2+1] = se;
				m++;
			}
			if (i==1) {
				bdry_n[m*2+0] = sw;
				bdry_n[m*2+1] = nw;
				m++;
			}
			if (j==ney/2) {
				bdry_n[m*2+0] = nw;
				bdry_n[m*2+1] = ne;
				m++;
			}
			if (i==nex/2) {
				bdry_n[m*2+0] = se;
				bdry_n[m*2+1] = ne;
				m++;
			}
		}
	}
	return;
}
/**************************************************************************************************************************************************/
	
