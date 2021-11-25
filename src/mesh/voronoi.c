/*! \file src/mesh/voronoi.c
 *
 *  Created by Casey Cavanaugh and James Adler on 2/18/20.
 *  Copyright 2015__HAZMATH__. All rights reserved.
 *
 *  \brief Obtains routines for generating voronoi mesh data given a
 *         Delaunay mesh.  Mostly used for mimetic finite-difference
 *         implementations of Maxwell's Equations.
 *
 *	NOTE: Only works in R3 for now!
 *
 */

#include "hazmath.h"

/******************************************************************************/

/*******************************************************************************/
/*!
* \fn void compute_Voronoi_nodes(mesh_struct* mesh, coordinates* cv_vor)
*
* \brief Computes the Voronoi nodes by finding circumcenters of Delaunay tets
*
* \param mesh      Delaunay triangulation mesh struct
*
* \return cv_vor   coord struct with Voronoi nodes (in ordering of Del tets)
*/

void compute_Voronoi_nodes(mesh_struct* mesh, coordinates* cv_vor)
{
  //mesh info of Delaunay mesh
  INT nelm = mesh->nelm;
  iCSRmat* el_v = mesh->el_v;
  
  //create coord struct for Voronoi nodes
  coordinates* cv_del = mesh->cv;
  
  //temp variables
  INT i,j,k;
  INT index[4];
  //coords of each element
  REAL elm_coords[12];
  //matrices we have to take determinants of
  REAL A[16];
  REAL Dx[16];
  REAL Dy[16];
  REAL Dz[16];
  //determinants
  REAL a = 666;
  REAL dx = -666;
  REAL dy = -666;
  REAL dz = -666;


  for(i = 0; i < nelm; i++){
	  
	//Delaunay node-element map
    get_incidence_row(i,el_v,index);	
	
	
    for(j=0; j<4; j++){							

      elm_coords[3*j] = cv_del->x[index[j]];
      elm_coords[3*j+1] = cv_del->y[index[j]];
      elm_coords[3*j+2] = cv_del->z[index[j]];
    }
	
	//fill in matrices we need to compute det of
    for(j=0; j<4; j++){							
      for(k=0; k<4; k++){
        if(k == 0){
          A[4*j+k] = elm_coords[3*j+k];
          Dx[4*j+k] = elm_coords[3*j]*elm_coords[3*j] + elm_coords[3*j+1]*elm_coords[3*j+1] + elm_coords[3*j+2]*elm_coords[3*j+2];
          Dy[4*j+k] = elm_coords[3*j]*elm_coords[3*j] + elm_coords[3*j+1]*elm_coords[3*j+1] + elm_coords[3*j+2]*elm_coords[3*j+2];
          Dz[4*j+k] = elm_coords[3*j]*elm_coords[3*j] + elm_coords[3*j+1]*elm_coords[3*j+1] + elm_coords[3*j+2]*elm_coords[3*j+2];
        }
        else if(k==1){
          A[4*j+k] = elm_coords[3*j+k];
          Dx[4*j+k] = elm_coords[3*j+k];
          Dy[4*j+k] = elm_coords[3*j];
          Dz[4*j+k] = elm_coords[3*j];
        }
        else if(k == 3){
          A[4*j+k]  = 1;
          Dx[4*j+k] = 1;
          Dy[4*j+k] = 1;
          Dz[4*j+k] = 1;
        }
        else{
          A[4*j+k] = elm_coords[3*j+k];
          Dx[4*j+k] = elm_coords[3*j+k];
          Dy[4*j+k] = elm_coords[3*j+k];
          Dz[4*j+k] = elm_coords[3*j+k-1];
        }
      }
    }

	//get determinants
    find_det_4(A, &a);
    find_det_4(Dx, &dx);
    find_det_4(Dy, &dy);
    find_det_4(Dz, &dz);

	//store voronoi pts
    cv_vor->x[i] = dx/(2*a);
    cv_vor->y[i] = -dy/(2*a);
    cv_vor->z[i] = dz/(2*a);
  }
	  

	  
  return;
}

/*!
* \fn  REAL compute_Voronoi_edges(mesh_struct * mesh, coordinates* cv_vor, dvector* vor_edge_length)
*
* \brief Computes the length of Voronoi edges using subtraction (interior edges)
		 or projections (boundary edges)
*
* \param mesh				Delaunay mesh struct
* \param cv_vor				coord struct with Voronoi nodes
*
* \return vor_edge_length		vector of voronoi edge length (in ordering of Del faces)
*/
void compute_Voronoi_edges(mesh_struct * mesh, coordinates* cv_vor, dvector* vor_edge_length)
{
  //Delaunay mesh info
  INT nface = mesh->nface;
  iCSRmat f_el;
  icsr_trans(mesh->el_f, &f_el);
  coordinates* cv_del = mesh->cv;
  
  //temp variables
  INT i;
  REAL x, y, z, vx, vy, vz;
  REAL nx, ny, nz;
  INT index[2];
  INT index_f[3];
  INT flag;
  REAL proj;

  for(i=0; i<nface; i++){
	
	//Delaunay elm-face maps = vor node-edge map
	get_incidence_row(i,&f_el,index); 
	
	//check if boundary del face = vor boundary edge
    flag = mesh->f_flag[i];		
	
	//if face is not on boundary, |edge| = |p1 - p2|
    if(flag == 0){	
								
      x = cv_vor->x[index[0]] - cv_vor->x[index[1]];
      y = cv_vor->y[index[0]] - cv_vor->y[index[1]];
      z = cv_vor->z[index[0]] - cv_vor->z[index[1]];
	  
      vor_edge_length->val[i] = sqrt(x*x + y*y + z*z);
    }
	//face is on boundary, so edge is chopped off
    //compute distance from point to the plane using projections
    else{

	  //look at face to vertex map, grab a point on the face, normal vec n
      get_incidence_row(i,mesh->f_v,index_f);

      vx = cv_vor->x[index[0]] - cv_del->x[index_f[0]];
      vy = cv_vor->y[index[0]] - cv_del->y[index_f[0]];
      vz = cv_vor->z[index[0]] - cv_del->z[index_f[0]];

      nx = mesh->f_norm[3*i];
      ny = mesh->f_norm[3*i+1];
      nz = mesh->f_norm[3*i+2];

	  proj = vx*nx+ vy*ny + vz*nz;

	  if(proj < 0){
		proj = -proj;
	  }

	  vor_edge_length->val[i] = proj;

    }
  }
  icsr_free(&f_el);

  return;
}



/*!
* \fn void neighbor_elm (mesh_struct * mesh, INT n, INT m, INT* ind, INT* f)
*
* \brief Determines whether two elements, n and m, are neighbors
*
* \param mesh 				Delaunay triangulation mesh struct
* \param n					index of one element
* \param m					index of another element
*
* \return ind				1 if they share a face, 0 if not
* \return f_el				index of the face they share (otherwise -666)
*
*/
void neighbor_elm (mesh_struct * mesh, INT n, INT m, INT* ind, INT* f)
{
  //if same element, return ind = 0, not neighbors
  if (n == m){	
    *ind =0;
    return;
  }
  //elm to face map
  iCSRmat* el_f = mesh->el_f;	
  //indices of 4 faces on n
  INT face_n[4];
  //indices of 4 faces on m  
  INT face_m[4];
  //temp variables  
  INT i = 0;
  INT j = 0;
  INT temp = 0;
  *f = -666;
  //get faces on n
  get_incidence_row(n,el_f,face_n);	
  //get faces on m
  get_incidence_row(m,el_f,face_m);	

  while (i<4 && temp == 0){
    j=0;
    while(j<4 && temp ==0){
	  //find if face i of elm n is in face set of m
      if (face_n[i] == face_m[j]){		
        temp = 1;
        *f = face_n[i];
      }
      j++;
    }
    i++;
  }
  *ind = temp;
  return;
}



/*!
* \fn REAL compute_Voronoi_faces(mesh_struct* mesh,coordinates* cv_vor, REAL* pt_on_face, dvec* vor_face_area)
*
* \brief Computes the area of Voronoi faces by partitioning the polygon into triangles
			and computing the area of each triangle using cross products
*
* \param mesh 				Delaunay triangulation mesh struct
* \param cv_vor				Voronoi coord struct
*
* \return vor_face_area		vector of Voronoi face area (in ordering of Del edges)
* \return pt_on_face		a point on each face--use for compute_Voronoi_volumes
*
* \note
*
*/
void compute_Voronoi_faces(mesh_struct* mesh,coordinates* cv_vor, REAL* pt_on_face, dvector* vor_face_area)
{
//Delaunay mesh info	
INT nedge = mesh->nedge;
// Get transposes of incident matrices
iCSRmat ed_el;
icsr_trans(mesh->el_ed,&ed_el);
iCSRmat ed_f;
icsr_trans(mesh->f_ed, &ed_f);
iCSRmat f_el;
icsr_trans(mesh->el_f, &f_el);


//temp variables
INT i,j;
REAL cross[3];
REAL u[3];
REAL v[3];
REAL mag = 0;
REAL area;
INT ind = 0;
INT f = -666;
INT a,k,l, elmi, nnz;
REAL x, y,z, nx, ny, nz;
INT num_nodes = 0;
REAL* coords = NULL;
INT* index = NULL;
INT* order = NULL;
INT* temp = NULL;
INT count;
REAL vx,vy,vz;

//get coord of an endpt of the edge
INT endpt[2];
//faces for boundary edges (only 2)
INT faces[2];
//temp vector for potential boundary faces
INT f_temp[10];


for(i=0; i<nedge; i++){
	
	//number of del elements incident to edge = num vor pts per face
	nnz	= ed_el.IA[i+1] - ed_el.IA[i];		
	//indices of del elements = vor pts
	index =(INT *)calloc(nnz,sizeof(INT));	
	get_incidence_row(i,&ed_el,index);
	
	//stores ordered nodes
	order = (INT *)calloc(nnz,sizeof(INT));	
	//marker for ordered or not
	temp =(INT *)calloc(nnz,sizeof(INT));	
	//choose p1 = first elm in index
	order[0] = index[0];	
	//mark it ordered	
	temp[0] = 1;						
	//count of points I've ordered
	count = 0;
	//do ordering of points I have
	while (count < nnz-1){			
		//grab last one I've ordered	
		a = order[count];		
		j = -1;
		ind = 0;	
		//loop over all elms in index		
		while (ind == 0 && j<nnz-1){
			j++;
			//if it hasn't been marked, see if it's a neighbor
			if (temp[j] == 0){	
			//check if a and ind[j] share a face			
			neighbor_elm(mesh, a, index[j], &ind, &f);	
			}
		}
		//when ind = 1, add index j to ordered points
		order[count+1] = index[j];	
		//label it as ordered
		temp[j] = 1;	
		//update counter		
		count++;						
	}
	
	//boundary edge--find intersection of face with boundary
	if(mesh->ed_flag[i] == 1){		

			
			//we need 3 extra points
			num_nodes = nnz+3;			
			//stores coords of points of face
			coords =(REAL *)calloc(3*(nnz+3),sizeof(REAL));		
			
			//grab a point on the del edge = pt on vor face
			get_incidence_row(i,mesh->ed_v,endpt);		
			x = mesh->cv->x[endpt[0]];
			y = mesh->cv->y[endpt[0]];
			z = mesh->cv->z[endpt[0]];
			
			//get the incident faces
			//don't know how many, so fill f_temp with -1's
			for(k=0; k<10; k++) {
				f_temp[k] = -1;
			}
			get_incidence_row(i,&ed_f,f_temp);
			//only store boundary faces in faces
			l=0;
			for(k=0; k<10; k++){
				if(f_temp[k] >-1){
					if(mesh->f_flag[f_temp[k]] == 1){
						faces[l] = f_temp[k];	
						l++;
					}
				}
			}
			if(l!=2){
				printf("HAZMATH DANGER: MORE BOUNDARY FACES THAN POSSIBLE!!\n");
				exit(0);
			}

		//compute points on faces -- add to beginning and end
		for(j=0; j<2; j++){
			
			f = faces[j];		

			nx = mesh->f_norm[3*f];
			ny = mesh->f_norm[3*f+1];
			nz = mesh->f_norm[3*f+2];
			
			//if boundary face, only belongs to one element
			INT elm[1];					
			get_incidence_row(f,&f_el,elm);
			vx = cv_vor->x[elm[0]];
			vy = cv_vor->y[elm[0]];
			vz = cv_vor->z[elm[0]];

			mag = nx* (vx-x) + ny*(vy-y) + nz*(vz-z);

			//determine whether we add the point to the beg or the end
			if(elm[0] == order[0]){    
				coords[3] = mag*nx + vx ;
				coords[4] = mag*ny + vy;
				coords[5] = mag*nz + vz;
			}
			else{
				coords[3*(nnz+3)-3]  = mag*nx + vx ;
				coords[3*(nnz+3)-2]  = mag*ny + vy;
				coords[3*(nnz+3)-1] = mag*nz + vz;
			}
		}
			//compute point on the boundary edge -- just the midpoint.
			coords[0] = mesh->ed_mid[3*i];	
			coords[1] = mesh->ed_mid[3*i +1];
			coords[2] = mesh->ed_mid[3*i+2];


		for(k=0; k< nnz; k++){			
			//get coords of pts we ordered		
			elmi = order[k];
			coords[3*k+6] = cv_vor->x[elmi];
			coords[3*k+6+1] = cv_vor->y[elmi];
			coords[3*k+6+2] = cv_vor->z[elmi];
			}
		
	}

	else{		
		//get coords of the ordered points
		num_nodes = nnz;
		coords =(REAL *)calloc(3*nnz,sizeof(REAL));
		for(k=0; k<nnz; k++){
			elmi = order[k];
			coords[3*k] = cv_vor->x[elmi];
			coords[3*k+1] = cv_vor->y[elmi];
			coords[3*k+2] = cv_vor->z[elmi];
		}
	}

	area = 0;
	for(j=0; j<num_nodes-2; j++){	
	
			//compute area of vor face using cross products
			u[0] = coords[3*(j+1)] - coords[0];
			u[1] = coords[3*(j+1)+1] - coords[1];
			u[2] = coords[3*(j+1)+2] - coords[2];
			v[0] = coords[3*(j+2)] - coords[0];
			v[1] = coords[3*(j+2)+1] - coords[1];
			v[2] = coords[3*(j+2)+2] - coords[2];
			cross_product(u, v, cross, &mag);
			area += .5*mag;

}

	vor_face_area->val[i] = area;	

	pt_on_face[3*i] = coords[0];
	pt_on_face[3*i+1] = coords[1];
	pt_on_face[3*i+2] = coords[2];


  // Free temporary arrays
	if(index) {
      free(index);
      index=NULL;
    }
	if(coords) {
      free(coords);
      coords=NULL;
    }
	if(order) {
      free(order);
      order=NULL;
    }
	if(temp) {
      free(temp);
      temp=NULL;
    }

}
	icsr_free(&ed_el);
	icsr_free(&f_el);
	icsr_free(&ed_f);
return;

}

/*!
* \fn void compute_Voronoi_volumes(mesh_struct* mesh,coordinates* cv_vor, dvec* vor_face_area, REAL* pt_on_face, dvector* vor_el_vol)
*
* \brief Computes the volume of the Voronoi polyhedra by partitioning Voronoi polyhedra into tetrahedra
*
* \param mesh 				Delaunay triangulation mesh struct
* \param cv_vor 			Voronoi coord struct
* \param vor_face_area 		vector of Voronoi face areas
* \param pt_on_face			point on each Voronoi face
*
* \return vor_el_vol			vector of Voronoi element volumes (in ordering of Del nodes)
*
*/
void compute_Voronoi_volumes(mesh_struct* mesh, coordinates* cv_vor, dvector* vor_face_area, REAL* pt_on_face, dvector* vor_el_vol)
{
   //Delaunay and Voronoi mesh info
  INT nv_del = mesh->nv;
  coordinates* cv_del = mesh->cv;
  iCSRmat v_ed;
  icsr_trans(mesh->ed_v,&v_ed);
  //temp variables
  INT i,j;
  INT* index = NULL; 		
  REAL tx, ty, tz, x, y, z;
  REAL volume;
  INT f, nnz;
  REAL height;

  for(i=0; i<nv_del; i++){
    nnz	= v_ed.IA[i+1] - v_ed.IA[i];
    index = (INT*)calloc(nnz, sizeof(INT));
	
	//get indices of incident edges
    get_incidence_row(i,&v_ed,index);	
    volume = 0;
	//loop over faces/poly = edges/del node
    for(j=0; j<nnz; j++){
		
	  //what face we're on	
      f = index[j];	
	  
	  // get tangent vector to Delaunay edge = normal toCoronoi face	  
      tx = mesh->ed_tau[3*f]; 	
      ty = mesh->ed_tau[3*f+1];
      tz = mesh->ed_tau[3*f+2];
	  
	  //find vector from face to del pt
      x = cv_del->x[i] - pt_on_face[3*f];	
      y = cv_del->y[i] - pt_on_face[3*f+1];
      z = cv_del->z[i] - pt_on_face[3*f+2];
	  
	  //project it in direction of normal to face (tan to Del edge)
	  height = x*tx + y*ty + z*tz;			
	  
	   if(height<0){
		  height = -height;
	  }
	  
	  volume+= height*vor_face_area->val[f]/3;
    }
	
    vor_el_vol->val[i] = volume;
	

    if(index) {
      free(index);
      index=NULL;
    }
  }
  icsr_free(&v_ed);
  return;
}
