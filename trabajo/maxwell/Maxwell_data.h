/*! \file Maxwell_data.h
*
*  Created by Adler, Hu, Zikatanov on 09/29/2017.
*  Edited by Casey Cavanaugh 5/20/19
*  Copyright 2016_HAZMATH__. All rights reserved.
*
* \brief This contains all the Data parameters and coefficients
*        for the Maxwell example.  This includes exact solutions,
*        RHS functions, coefficients, and boundary conditions.
*
*/

// PDE Coefficients
void permitivity(REAL *val,REAL* x,REAL time,void *param) {
  *val = 1.0;
}
void permeability(REAL *val,REAL* x,REAL time,void *param) {
  *val = 1.0;
}
void oneovermu(REAL *val,REAL* x,REAL time,void *param) {
  REAL mu = -6.66;
  permeability(&mu,x,time,param);
  *val = 1.0/mu;
}


// True Solution (if you have one)
void truesol(REAL *val,REAL* x,REAL time,void *param) {

  REAL a = 1.0;
  REAL myexp = exp(-a*time);

  //ordering (E1,E2,E3,B1,B2,B3,p)

  //Real test problem
  val[0] = -cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]) *myexp/M_PI;
  val[1] = sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]) *myexp/M_PI;
  val[2] = 0.0;
  val[3] = -sin(M_PI*x[0])*cos(M_PI*x[1])*cos(M_PI*x[2]) *myexp;
  val[4] = -cos(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*x[2]) *myexp;
  val[5] = 2*cos(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]) *myexp;
  val[6] = 0.0;

  /* val[0] = myexp*(x[0]*x[0] - x[0]*x[1] + 2*x[1]*x[2]);
  val[1] = myexp*(2*x[2]*x[2] - 2*x[0]*x[1]);
  val[2] = myexp*(x[1]*x[2] - 3*x[1]*x[1]);
  val[3] = myexp*(-6*x[1] - 3*x[2]);
  val[4] = myexp*2*x[1];
  val[5] = myexp*(x[0] - 2*x[1] - 2*x[2]);
  val[6] = 0.0;
  */

}
void Etrue(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
}
void Btrue(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  val[0] = myu[3];
  val[1] = myu[4];
  val[2] = myu[5];
}
void ptrue(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  *val = myu[6];
}

// Right-hand Side
//NEED TO INPUT -j as RHS
void current_density(REAL *val,REAL* x,REAL time,void *param) {
  REAL a = 1.0;
  REAL myexp = exp(-a*time);

  val[0] = myexp*(3*M_PI + 1/M_PI)*cos(M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[2]);
  val[1] = myexp*(3*M_PI + 1/M_PI)*-sin(M_PI*x[0])*cos(M_PI*x[1])*sin(M_PI*x[2]);
  val[2] = 0.0;

  /* val[0] = -myexp*(x[0]*x[0] - x[0]*x[1] + 2*x[1]*x[2] - 2);
  val[1] = -myexp*(2*x[2]*x[2] - 2*x[0]*x[1] - 4);
  val[2] = -myexp*(x[1]*x[2] - 3*x[1]*x[1] + 6);
  */
}


// Boundary Conditions
void bc(REAL *val,REAL* x,REAL time,void *param) {
  REAL myu[7];
  truesol(myu,x,time,param);
  val[0] = myu[0];
  val[1] = myu[1];
  val[2] = myu[2];
  val[3] = myu[3];
  val[4] = myu[4];
  val[5] = myu[5];
  val[6] = myu[6];
}
void bc_test(REAL *val,REAL* x,REAL time,void *param) {
  *val = 0.0;
}

/***********************************************************************************************/
/*!
* \fn void find_det_4( REAL* A, REAL deta)
*
* \brief find det of 4x4 matrix
*
* \param A            vectorized matrix
*
* \return deta        determinant
*
*/

void find_det_4( REAL* A, REAL* deta)
{

  REAL s0p, s0m, s1p, s1m, s2p, s2m, s3p, s3m;

  s0p = A[5]*A[10]*A[15] + A[6]*A[11]*A[13] + A[7]*A[9]*A[14];
  s0m = A[13]*A[10]*A[7] + A[14]*A[11]*A[5] + A[15]*A[9]*A[6];

  s1p = A[4]*A[10]*A[15] + A[6]*A[11]*A[12] + A[7]*A[8]*A[14];
  s1m = A[12]*A[10]*A[7] + A[14]*A[11]*A[4] + A[15]*A[8]*A[6];

  s2p = A[4]*A[9]*A[15] + A[5]*A[11]*A[12] + A[7]*A[8]*A[13];
  s2m = A[12]*A[9]*A[7] + A[13]*A[11]*A[4] + A[15]*A[8]*A[5];

  s3p = A[4]*A[9]*A[14] + A[5]*A[10]*A[12] + A[6]*A[8]*A[13];
  s3m = A[12]*A[9]*A[6] + A[13]*A[10]*A[4] + A[14]*A[8]*A[5];

  *deta = A[0]*(s0p-s0m) - A[1]*(s1p-s1m) + A[2]*(s2p-s2m) - A[3]*(s3p-s3m);

  return;

}


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
  INT nelm = mesh->nelm;
  INT nv = mesh->nv;
  INT dim = mesh->dim;
  coordinates* cv_del = mesh->cv;
  iCSRmat* el_v = mesh->el_v;



  INT i,j,k;
  REAL elm_coords[12];
  REAL* A = (REAL *)calloc(16,sizeof(REAL));
  REAL* Dx =(REAL *)calloc(16,sizeof(REAL));
  REAL* Dy =(REAL *)calloc(16,sizeof(REAL));
  REAL* Dz =(REAL *)calloc(16,sizeof(REAL));
  REAL a = 666;
  REAL dx = -666;
  REAL dy = -666;
  REAL dz = -666;
  INT p[4];
  REAL piv[4];
  INT* index = (INT *)calloc(4,sizeof(INT));


  for(i = 0; i < nelm; i++){
    get_incidence_row(i,el_v,index);	//find Del nodes on each elm
	
    for(j=0; j<4; j++){							//get Del coords on each elm
      elm_coords[3*j] = cv_del->x[index[j]];
      elm_coords[3*j+1] = cv_del->y[index[j]];
      elm_coords[3*j+2] = cv_del->z[index[j]];
    }
	
    for(j=0; j<4; j++){							//fill in matrices we need to compute det of
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
	
	//fill in voronoi pts
    cv_vor->x[i] = dx/(2*a);
    cv_vor->y[i] = -dy/(2*a);
    cv_vor->z[i] = dz/(2*a);
  }

  return;
}

/*!
* \fn  REAL compute_Voronoi_edges(mesh_struct * mesh, coordinates* cv_vor, REAL* vor_edge_length)
*
* \brief Computes the length of Voronoi edges
*
* \param mesh				Delaunay mesh struct
* \param cv_vor				coord struct with Voronoi nodes
*
* \return vor_edge_length		vector of voronoi edge length (in ordering of Del faces)
*/
void compute_Voronoi_edges(mesh_struct * mesh, coordinates* cv_vor, REAL* vor_edge_length)
{
  INT nface = mesh->nface;
  iCSRmat f_el;
  icsr_trans(mesh->el_f, &f_el);
  coordinates* cv_del = mesh->cv;
  INT i,j;
  REAL x, y, z, vx, vy, vz;
  REAL nx, ny, nz;
  INT index[2];
  INT index_f[3];
  INT flag;
  REAL proj;

  for(i=0; i<nface; i++){
	 
    get_incidence_row(i,&f_el,index); //incident del elm on face = vor nodes on edge
 
    flag = mesh->f_flag[i];		//check if boundary del face = vor boundary edge
	
    if(flag == 0){	//if face is not on boundary, |edge| = |p1 - p2|
								
      x = cv_vor->x[index[0]] - cv_vor->x[index[1]];
      y = cv_vor->y[index[0]] - cv_vor->y[index[1]];
      z = cv_vor->z[index[0]] - cv_vor->z[index[1]];
	  
      vor_edge_length[i] = sqrt(x*x + y*y + z*z);

    }
    else{
      //face is on boundary, so edge is chopped off
      //compute distance from point to the plane using projections
      //f_v face to vertex map, grab a point on the face, normal vec n

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

	  vor_edge_length[i] = proj;

    }
  }
  
  return;
}

/*!
* \fn void cross_product(REAL u, REAL v, REAL cross, REAL mag)
*
* \brief Computes the cross product and it's magnitude of two vectors R3
*
* \param u 			vector in R3
* \param v 			vector in R3
*
* \return cross		  u x v
* \return mag 		| u x v |
*
*/
void cross_product(REAL* u, REAL* v, REAL* cross, REAL* mag)
{
  cross[0] = u[1]*v[2] - u[2]*v[1];
  cross[1] = u[2]*v[0] - u[0]*v[2];
  cross[2] = u[0]*v[1] - u[1]*v[0];

  *mag = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

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
  if (n == m){	//if same element, return ind = 0, not neighbors
    *ind =0;
    return;
  }
  iCSRmat* el_f = mesh->el_f;	//elm to face map
  INT* face_n =(INT *)calloc(4,sizeof(INT));	//indices of 4 faces on n
  INT* face_m =(INT *)calloc(4,sizeof(INT));	//indices of 4 faces on m
  INT i = 0;
  INT j = 0;
  INT temp = 0;
  *f = -666;
  get_incidence_row(n,el_f,face_n);	//get faces on n
  get_incidence_row(m,el_f,face_m);	//get faces on m

  while (i<4 && temp == 0){
    j=0;
    while(j<4 && temp ==0){
      if (face_n[i] == face_m[j]){		//find if face i of elm n is in face set of m
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
* \fn REAL compute_Voronoi_faces(mesh_struct* mesh,coordinates* cv_vor, REAL* pt_on_face, REAL* vor_face_area)
*
* \brief Computes the area of Voronoi faces
*
* \param mesh 				Delaunay triangulation mesh struct
* \param cv_vor				Voronoi coord struct
*
* \return vor_face_area		vector of Voronoi face area (in ordering of Del edges)
* \return pt_on_face		a point on each face--use for compute_Voronoi_volumes
*
*/
void compute_Voronoi_faces(mesh_struct* mesh,coordinates* cv_vor, REAL* pt_on_face, REAL* vor_face_area)
{
	
INT nface = mesh->nface;
INT nedge = mesh->nedge;

iCSRmat ed_el;
icsr_trans(mesh->el_ed,&ed_el);

INT i,j;
REAL cross[3];
REAL* u = (REAL *)calloc(3,sizeof(REAL));
REAL* v =(REAL *)calloc(3,sizeof(REAL));
REAL mag = 0;
REAL area;
INT ind = 0;
INT f = -666;
INT a,k, elmi, nnz;
REAL x, y,z, nx, ny, nz;
INT num_nodes = 0;
REAL* coords = NULL;
INT* index = NULL;
INT* order = NULL;
INT* temp = NULL;
INT count;

iCSRmat ed_f;
icsr_trans(mesh->f_ed, &ed_f);

iCSRmat f_el;
icsr_trans(mesh->el_f, &f_el);

//get coord of an endpt of the edge
INT* endpt =(INT *)calloc(2,sizeof(INT));
INT* faces =(INT *)calloc(2,sizeof(INT));

for(i=0; i<nedge; i++){
	//printf("face1\n");
	
	nnz	= ed_el.IA[i+1] - ed_el.IA[i];		//number of del elements incident to edge = num vor pts per face
	index =(INT *)calloc(nnz,sizeof(INT));	//indices of del elements = vor pts
	get_incidence_row(i,&ed_el,index);
	//printf("face2\n");
	
	order = (INT *)calloc(nnz,sizeof(INT));	//stores ordered nodes
	temp =(INT *)calloc(nnz,sizeof(INT));	//marker for ordered or not

	order[0] = index[0];			//choose p1 = first elm in index
	temp[0] = 1;					//mark it ordered

	//do ordering of points I have
	count = 0;
	while (count < nnz-1){			//count of points I've ordered
		//printf("face3\n");	
		a = order[count];		//grab last one I've ordered
		j = -1;
		ind = 0;									//loop over all elms in index
		while (ind == 0 && j<nnz-1){
			j++;
			if (temp[j] == 0){		//if it hasn't been marked, see if it's a neighbor
			neighbor_elm(mesh, a, index[j], &ind, &f);	//check if a and ind[j] share a face
			}
		}
		order[count+1] = index[j];		//when ind = 1, add index j to ordered points
		temp[j] = 1;					//label it as ordered
		count++;						//update counter
	}
	
	
	if(mesh->ed_flag[i] == 1){			//boundary edge--find intersection of face with boundary
		//	printf("face4\n");
			num_nodes = nnz+3;			//we need 3 extra points

			coords =(REAL *)calloc(3*(nnz+3),sizeof(REAL));		//stores coords of points of face

			get_incidence_row(i,mesh->ed_v,endpt);		//grab a point on the del edge = pt on vor face
			x = mesh->cv->x[endpt[0]];
			y = mesh->cv->y[endpt[0]];
			z = mesh->cv->z[endpt[0]];
			
			//get the incident faces
			get_incidence_row(i,&ed_f,faces);
			
			//printf("face5\n");
		for(j=0; j<2; j++){		//compute points on faces -- add to beginning and end
			//printf("face6\n");
			REAL vx,vy,vz;
			
			f = faces[j];
			nx = mesh->f_norm[3*f];
			ny = mesh->f_norm[3*f+1];
			nz = mesh->f_norm[3*f+2];

			INT elm[1];					//if boundary face, only belongs to one element
			get_incidence_row(f,&f_el,elm);
			//printf("face7\n");
			vx = cv_vor->x[elm[0]];
			vy = cv_vor->y[elm[0]];
			vz = cv_vor->z[elm[0]];

			mag = nx* (vx-x) + ny*(vy-y) + nz*(vz-z);

			if(elm[0] == order[0]){    //determine whether we add the point to the beg or the end
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

			coords[0] = mesh->ed_mid[3*i];	//compute point on the boundary edge -- just the midpoint.
			coords[1] = mesh->ed_mid[3*i +1];
			coords[2] = mesh->ed_mid[3*i+2];

		for(k=0; k< nnz; k++){			//get coords of pts we ordered		
			elmi = order[k];
			coords[3*k+6] = cv_vor->x[elmi];
			coords[3*k+6+1] = cv_vor->y[elmi];
			coords[3*k+6+2] = cv_vor->z[elmi];
			}
	}

	else{		//get coords of the ordered points
		num_nodes = nnz;
		coords =(REAL *)calloc(3*nnz,sizeof(REAL));
		for(k=0; k<nnz; k++){
			elmi = order[k];
			coords[3*k] = cv_vor->x[elmi];
			coords[3*k+1] = cv_vor->y[elmi];
			coords[3*k+2] = cv_vor->z[elmi];
		}
		//printf("face8\n");
	}

	area = 0;
	for(j=0; j<num_nodes-2; j++){		//compute area of vor face using cross products
			//printf("face9\n");
			u[0] = coords[3*(j+1)] - coords[0];
			u[1] = coords[3*(j+1)+1] - coords[1];
			u[2] = coords[3*(j+1)+2] - coords[2];
		//	printf("u: %f  %f  %f\n", u[0], u[1], u[2]);
			v[0] = coords[3*(j+2)] - coords[0];
			v[1] = coords[3*(j+2)+1] - coords[1];
			v[2] = coords[3*(j+2)+2] - coords[2];
		//	printf("v: %f  %f  %f\n", v[0], v[1], v[2]);
			cross_product(u, v, cross, &mag);
		//	printf("mag u x v: %f\n", mag);
			area += .5*mag;

}

	vor_face_area[i] = area;	
	pt_on_face[3*i] = coords[0];
	pt_on_face[3*i+1] = coords[1];
	pt_on_face[3*i+2] = coords[2];


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
return;

}

/*!
* \fn void compute_Voronoi_volumes(mesh_struct* mesh,coordinates* cv_vor, REAL* vor_face_area, REAL* pt_on_face, REAL* vor_el_vol)
*
* \brief Computes the volume of the Voronoi polyhedra
*
* \param mesh 				Delaunay triangulation mesh struct
* \param cv_vor 			Voronoi coord struct
* \param vor_face_area 		vector of Voronoi face areas
* \param pt_on_face			point on each Voronoi face
*
* \return vor_el_vol			vector of Voronoi element volumes (in ordering of Del nodes)
*
*/
void compute_Voronoi_volumes(mesh_struct* mesh, coordinates* cv_vor, REAL* vor_face_area, REAL* pt_on_face, REAL* vor_el_vol)
{

  INT nv_del = mesh->nv;
  coordinates* cv_del = mesh->cv;
  iCSRmat* el_ed = mesh->el_ed;
  iCSRmat v_ed;
  icsr_trans(mesh->ed_v,&v_ed);

  INT i,j;
  INT* index = NULL; 		//preallocate?
  REAL tx, ty, tz, x, y, z;
  REAL volume;
  INT f, nnz;
  REAL height;
  
  for(i=0; i<nv_del; i++){
    nnz	= v_ed.IA[i+1] - v_ed.IA[i];
	//printf("%d \n", nnz);
    index = (INT*)calloc(nnz, sizeof(INT));
    get_incidence_row(i,&v_ed,index);	//get indices of incident edges
    volume = 0;
    for(j=0; j<nnz; j++){		//loop over faces/poly = edges/del node
      f = index[j];				//what face we're on
      tx = mesh->ed_tau[3*f]; 	// get tangent vector to edge = normal to face
      ty = mesh->ed_tau[3*f+1];
      tz = mesh->ed_tau[3*f+2];
      x = cv_del->x[i] - pt_on_face[3*f];	//find vector from face to del pt
      y = cv_del->y[i] - pt_on_face[3*f+1];
      z = cv_del->z[i] - pt_on_face[3*f+2];
	  
	  height = x*tx + y*ty + z*tz;			//project it in direction of normal to face (tan to del edge)
	  
	   if(height<0){
		  height = -height;
	  }
	  
	  volume+= height*vor_face_area[f]/3;
	  // printf("height = %f, base = %f, vol = %f  \n", height, vor_face_area[f], volume);
    }
	
    vor_el_vol[i] = volume;
	
    if(index) {
      free(index);
      index=NULL;
    }
  }
  return;
}

/***********************************************************************************************/
/*!
* \fn dCSRmat dcsr_create_diagonal_matrix (const INT m, const INT n, const INT index_start, REAL* vec)
*
* \brief Create a dCSRmat sparse matrix that is diagonal
*
* \param m             Number of rows
* \param index_start   Number from which memory is indexed (1 for fortran, 0 for C)
*
* \return A            the new dCSRmat matrix
*
*/
dCSRmat dcsr_create_diagonal_matrix(const INT m,
  const INT index_start, REAL* vec)
  {
    dCSRmat A;

    A.IA = (INT *)calloc(m+1, sizeof(INT));
    A.JA = (INT *)calloc(m, sizeof(INT));
    A.val = (REAL *)calloc(m, sizeof(REAL));

    A.row = m; A.col = m; A.nnz = m;

    INT i;
    for(i=0;i<m;i++){
      A.IA[i]=i + index_start;
      A.JA[i]=i + index_start;
      A.val[i] = vec[i];
    }
    A.IA[m] = m + index_start;

    return A;
  }
