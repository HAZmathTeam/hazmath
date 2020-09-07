/*********** HAZMATH FUNCTIONS and INCLUDES ***************************************/
#include "hazmath.h"
//#include "elasticity_data_2.h"
#include "elasticity_data.h"
//#include "elasticity_system.h"
/*********************************************************************************/

void main()
{
  // mesh read
  INT cycle=0;
  char filename_per_cycle[100]={'\0'};
  sprintf(filename_per_cycle, "2dSQ-%d.haz",cycle);
  FILE *gfid = HAZ_fopen(filename_per_cycle,"r");
  if(gfid == NULL){
    perror("Could not find and open the file !!!! ");
  }

  // Create the mesh (now we assume triangles in 2D or tetrahedra in 3D)
  // File types possible are 0 - old format; 1 - vtk format
  INT mesh_type = 0;
  mesh_struct *mesh=(mesh_struct *)malloc(1*sizeof(mesh_struct)); //
  printf(" --> loading grid from file: %s\n",filename_per_cycle);
  initialize_mesh(mesh);   // Q1. Why only here? 
  creategrid_fread(gfid,mesh_type,mesh);
  fclose(gfid);
  ///////////////////////////////////////////////////////
  // here is the action of interest
  iCSRmat *f_el=(iCSRmat *)malloc(1*sizeof(iCSRmat)); // face_to_element;
  icsr_trans(mesh->el_f,f_el); // f_el=transpose(el_f); 
  INT elm,nbr,jk,face,nbr0,pq; // current simplex and neighbor
  iCSRmat *el_f=mesh->el_f;//place holder
  ////////////////////////////////////////////////////////
  INT facein=0; //num faces interior
  INT facebd=0;// num faces on boundary
  for(elm=0;elm<mesh->nelm;elm++) {
    for(jk=el_f->IA[elm];jk<el_f->IA[elm+1];jk++){
      face=el_f->JA[jk];
      // now we use the transposed and the fact that we have at most
      // two elements sharing a face; we look at the two (or less)
      // elements attached to a face and the element with number not
      // equal to elm is the neighbor, or we are on the boundary.
      nbr=-1;
      for(pq=f_el->IA[face];pq<f_el->IA[face+1];pq++){	
	nbr0=f_el->JA[pq];
	if(nbr0==elm) continue;
	nbr=nbr0;
      }
      if(nbr<0) {
	// boundary
	facebd++;
	fprintf(stdout,"\nface %d in element %d is on the boundary",	\
		face,elm);
      } else {
	if(nbr>elm)facein++; // only count each face once, just for
			     // counting
	fprintf(stdout,"\nface %d in element %d is shared with neighbor %d", \
		face,elm,nbr);
      }
    }
  }
  // each of the faces is 
  fprintf(stdout,"\nfaces(interior)=%d; faces(boundary)=%d\n\n",	\
	  facein,facebd);
  return;
}
