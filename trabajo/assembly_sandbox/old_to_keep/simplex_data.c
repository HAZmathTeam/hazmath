#include "hazmath.h"
#include "biot_headers.h"
/********************************************************************/
INT calc_ndofs(INT fetype,INT dim,					\
	       const INT nv,const INT ne,const INT nf,const INT ns)
{
  // nv is num vertices, ne num edges, nf num faces, ns num simplices. 
    switch(fetype){
    case 1:
      return nv;
      break;
    case 20:
      return ne;
      break;
    case 30:
      return nf;
      break;
    case 60:
      return nv*dim;
      break;
    case 61:
      return nf;
      break;
    default:
      return ns;
      break;
    }
}
/*************************************************************************/    
simplex_data *simplex_data_init(INT dim,void *fe_in, const INT isit_block)
{
  //input is a (block) fespace; this function creates structure for
  //holding element data. This should only be called once before any
  //loops over elements. It does not change dynamically. The dynamic
  //changes are through the "simplex_data_update" function.
  // allocate the structure. 
  simplex_data *s=malloc(1*sizeof(simplex_data));
  //
  block_fespace *fe;
  INT i,j,k,l,fetype;
  s->dim=dim;
  if(isit_block){
    fe=(block_fespace *)fe_in;
  } else {
    fe = (block_fespace *)malloc(1*sizeof(block_fespace));
    fe->nspaces=1;
    fe->var_spaces[0]=(fespace *)fe_in;
  }   
  s->ns=1; // number of elements is 1.
  s->nv=dim+1; 
  s->ne=(INT )(dim*(dim+1))/2;
  s->nf=s->nv;// facets are opposite to vertices
  s->nspaces=fe->nspaces;// how many FE spaces;
  s->fetypes=calloc(2*s->nspaces,sizeof(INT)); // types of elements (could be only one type);
  s->ndofs=s->fetypes+s->nspaces; // array of number of local dofs per space.
  s->num_dofs=0;// counter for degrees of freedom. 
  for(j=0;j<s->nspaces;j++){
    s->fetypes[j]=fe->var_spaces[j]->FEtype;
    s->ndofs[j]=calc_ndofs(s->fetypes[j],dim,s->nv,s->ne,s->nf,1); // one simplex only
    
    s->num_dofs+=s->ndofs[j];      
  }
  // the ints:
  //allocate the intergers:
  INT memi = s->nv + s->nf + s->nf*dim + s->nf + 2*s->num_dofs; // vertices + faces + 2*(all dofs)
  /*************************************************************/
  INT *iwrk=calloc(memi,sizeof(INT)); // vertices in element
  s->v    = iwrk; // vertices in element
  s->f    = s->v + s->nv;// faces in an element;
  s->fv    = s->f + s->nf;// face to vertex map: s->nf faces * dim
			  // vertices per face.
  s->f_flag = s->fv + s->nf*dim;// flags on faces
  s->dofs   = s->f_flag + s->nf;// array with ALL degrees of freedom
				// space by space
  s->gdofs   = s->dofs + s->num_dofs;// array with ALL GLOBAL degrees
				     // of freedom
  INT *iwrkend   = s->gdofs + s->num_dofs;//where it all ends
  /*************************************************************/
  //allocate the reals:
  INT memr=(s->nv + s->nf + 1 + s->nf)*dim + s->nf ;//for the simplex
						    //(vcoords,barycenter,
						    //faces normals,
						    //faces
						    //barycenters)
  REAL *wrk=calloc(memr,sizeof(REAL)); 
  s->x=wrk; // vertex coords
  s->f_norm = s->x   + s->nv*dim;// normal vectors to the faces
  s->mid  = s->f_norm + s->nf*dim;// barycenter
  s->f_mid  = s->mid + dim;//barycenters for the faces;
  s->f_area = s->f_mid + s->nf*dim;// areas of faces.
  REAL *wrkend=s->f_area+s->nf; // end
  /*******basis*******************************************/
  memr = (1 + dim + dim*dim)*s->num_dofs; //values of basis;
  REAL *wrkp=calloc(memr,sizeof(REAL)); 
  s->p      = wrkp; // basis function values;
  s->dp     = s->p+s->num_dofs; // gradient of a basis;
  s->ddp    = s->dp+s->num_dofs*dim; // hessian of the basis;
  REAL *wrkpend=s->ddp+s->num_dofs*dim*dim; // end
  /*******basis*******************************************/
  //edge data is null in our application;
  s->ed=NULL;// edges in an element
  s->edv=NULL;
  s->fed=NULL;
  s->ed_len=NULL;
  s->ed_tau=NULL;
  s->ed_mid=NULL;
  // numdofs=(dim+1)+dim*(dim+1)+(dim+1)+1
  INT memwrk=6*(s->num_dofs + s->num_dofs*dim + s->num_dofs*dim*dim);
  // 5 should be enough....
  INT memwrk0=6*((dim+1)*dim+(dim+1)*dim*dim+ (dim+1)*dim*dim);
  if(memwrk<memwrk0) memwrk=memwrk0;
  s->wrk=(void *)calloc(memwrk,sizeof(REAL));
  // clean up
  if(!isit_block) free(fe);
  // return ; end
  return s;
}
/*************************************************************************/
void simplex_data_update(INT snum,simplex_data *s, mesh_struct *mesh,	\
			 void *fe_in,const INT isit_block)
{
  INT i,j,k,l;
  INT dim=s->dim;
  s->vol=mesh->el_vol[snum];// volume
  s->flag=mesh->el_flag[snum];
  get_incidence_row(snum,mesh->el_v,s->v);
  get_incidence_row(snum,mesh->el_f,s->f);
  for(i=0;i<s->nf;i++){
    k=s->f[i];
    get_incidence_row(k,mesh->f_v,(s->fv+i*dim));    
  }
  for(i=0;i<s->nv;i++){
    s->x[i*dim+0]=mesh->cv->x[s->v[i]];
    if(dim<2) continue;
    s->x[i*dim+1]=mesh->cv->y[s->v[i]];
    if(dim<3) continue;
    s->x[i*dim+2]=mesh->cv->z[s->v[i]];
  }    
  for(i=0;i<s->nf;i++){
    k=s->f[i]*dim;
    l=i*dim;
    s->f_area[i] = mesh->f_area[s->f[i]];
    s->f_flag[i] = mesh->f_flag[s->f[i]];
    for(j=0;j<dim;j++){
      s->f_norm[l+j] = mesh->f_norm[k+j];
      s->f_mid[l+j]  = mesh->f_mid[k+j];//barycenters for the faces;
    }
  }  
  for(j=0;j<dim;j++)
    s->mid[j] =mesh->el_mid[snum*dim+j];
  //edge data is null in our application;
  s->ed=NULL;// edges in an element
  s->edv=NULL;
  s->fed=NULL;
  s->ed_len=NULL;
  s->ed_tau=NULL;
  s->ed_mid=NULL;
  // FE spaces; this should be improved....
  // get all dofs. 
  block_fespace *fe; 
  if(isit_block){
    fe=(block_fespace *)fe_in;
  } else {
    fe = (block_fespace *)malloc(1*sizeof(block_fespace));
    fe->nspaces=1;
    fe->var_spaces[0]=(fespace *)fe_in;
  }
  INT *ias,*jas,rowa,rowb;
  i =0;
  l=0;
  for(k=0;k<fe->nspaces;k++) {
    ias = fe->var_spaces[k]->el_dof->IA;
    rowa = ias[snum];
    rowb = ias[snum+1];
    jas=fe->var_spaces[k]->el_dof->JA;
    for (j=rowa; j<rowb; j++) {
      //      fprintf(stdout,"\ni=%d:%d,jas=%d",i,s->num_dofs,jas[j]);
      s->dofs[i] = jas[j];
      s->gdofs[i] = l + jas[j];
      i++;
    }
    l+=fe->var_spaces[k]->ndof;
  }
  if(!isit_block)free(fe);
  return;
}
/*************************************************************************/
local_vec **dof_data_init(simplex_data *s)
{
  // local (element=snum) data for a vector. this is application 
  INT j,dim=s->dim,ntime_steps=2;
  REAL *ueptr=NULL,*ueall=NULL;
  local_vec **ue=malloc(ntime_steps*sizeof(local_vec *));
  for(j=0;j<(ntime_steps+1);j++)
    ue[j]=malloc(sizeof(local_vec));  
  local_vec *ue1=ue[0];
  local_vec *uet=ue[1];
  //  local_vec *ue0=ue[2];
  ueall=calloc(3*s->num_dofs,sizeof(REAL));
  //ue1 (current solution at time t)
  ue1->b = ueall;
  ue1->u = ue1->b   + s->ndofs[0];
  ue1->w = ue1->u   + s->ndofs[1]*dim;
  ue1->p = ue1->w   + s->ndofs[3];
  ueptr  = ue1->p   + s->ndofs[4];
  // time derivatives:
  //  fprintf(stdout,"\nsize=%ld vs %d\n",ueptr-ue1->b,s->num_dofs); fflush(stdout);
  uet->b = ueptr;
  uet->u = uet->b  + s->ndofs[0];
  uet->w = uet->u  + s->ndofs[1]*dim;
  uet->p = uet->w  + s->ndofs[3];
  ueptr  = uet->p + s->ndofs[4];
  //  fprintf(stdout,"\n2size=%ld vs %d\n",ueptr-uet->b,s->num_dofs); fflush(stdout);
  /* // just for computing time derivatives. */
  /* ue0->b  = ueptr; */
  /* ue0->u  = ue0->b  + s->ndofs[0]; */
  /* ue0->w  = ue0->u  + s->ndofs[1]*dim; */
  /* ue0->p  = ue0->w  + s->ndofs[3]; */
  /* ueptr   = ue0->p  + s->ndofs[4]; */
  /* fprintf(stdout,"\n3size=%ld vs %d\n",ueptr-uet->b,s->num_dofs); fflush(stdout); */
  return ue;
}
/*************************************************************************/
void dof_data_update(local_vec **ue,				\
		     INT snum, simplex_data *s,			\
		     dvector *uh0, dvector *uh1, const REAL dt)
{
  /* 
     local (element=snum) data for a vectors uh1 and uh0. if uh1 and uh0 are not
     null returns approximation of the difference on a simplex of
     their degrees of freedom. 
  */
  INT j;
  local_vec *ue1=ue[0];
  local_vec *uet=ue[1];
  /* local_vec *ue0=ue[2]; */
  if(uh1 != NULL){
    for(j=0;j<s->num_dofs;j++){
      ue1->b[j]=uh1->val[s->gdofs[j]];
      //      fprintf(stdout,"\nuh(%d)=%.15f",s->gdofs[j],uh1->val[s->gdofs[j]]);
    }
  } else { 
    for(j=0;j<s->num_dofs;j++){
      ue1->b[j]=0e0;
    }
  }
  /**/
  if(uh0 != NULL){
    for(j=0;j<s->num_dofs;j++){
      uet->b[j]=(ue1->b[j] - uh0->val[s->gdofs[j]])/dt;
    }
  } else {
    for(j=0;j<s->num_dofs;j++){
       uet->b[j]=ue1->b[j];
    }
  }
  return;
}
/*************************************************************************/
void simplex_data_print(INT snum, simplex_data *s, local_vec **ue)
{
  fprintf(stdout,
	  "\n================================================");
  fprintf(stdout,
	  "\nSimplex=%d; flag=%d; dim=%d; nf=%d; nv=%d;  vol=%e", \
	  snum,s->flag,s->dim,s->nf,s->nv,s->vol);
  print_full_mat_int(1,s->nv,s->v,"nop");
  print_full_mat_int(1,s->nf,s->f,"faces");
  print_full_mat_int(s->nf,s->dim,s->fv,"fv");
  print_full_mat_int(1,s->nf,s->f_flag,"f_flag");
  print_full_mat_int(1,s->num_dofs,s->dofs,"DOFs");
  print_full_mat_int(1,s->num_dofs,s->gdofs,"GLOB_DOFs");
  print_full_mat(s->nv,s->dim,s->x,"xv");
  print_full_mat(1,s->dim,s->mid,"bary");
  print_full_mat(s->nf,s->dim,s->f_norm,"normals");
  print_full_mat(1,s->nf,s->f_area,"area");
  print_full_mat(s->nf,s->dim,s->f_mid,"f_centers");
  if(ue!=NULL){
    print_full_mat(1,s->num_dofs,ue[0]->b,"ue1"); 
    print_full_mat(1,s->num_dofs,ue[1]->b,"uet");
  }
  fprintf(stdout,
	  "\n--------------------------------------------------");
  fflush(stdout);
  return;
}
/*************************************************************************/
void simplex_data_free(simplex_data *s)
{
  free(s->fetypes);  // fe
  free(s->v);  // integers
  free(s->x);  // real 
  free(s->p);  // basis and derivatives
  free(s->wrk);// working area
}
/*************************************************************************/
void dof_data_free(local_vec **ue)
{
  INT j;
  //  print_full_mat(1,13,ue[0]->b,"ueue");
  free(ue[0]->b);// these are all REAL pointers
  free(ue[0]);
  free(ue[1]);
  free(ue);
  return;
}
