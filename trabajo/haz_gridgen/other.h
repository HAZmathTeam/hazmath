/************************************************************/
static unsigned int cmp_simplex(INT n, INT sim1, INT sim2,	\
			 INT *sv1, INT *sv2, INT *stos1, INT *stos2)
{
  //sim1 and sim2 are pointers to the neighboring lists of two
  //simplices stos1 and stos2 are pointers to rows in
  //simplex-to-simplex incidence; sv1 and sv2 are pointers to rows in
  //simplex-vertex incidence for two neighboring simplices.
  //
  //compare two sets of length n. they should overlap at exactly n-1 members.  If they do
  //not, 0 is returned, otherwise 1. no check if the members of the
  //sets are distinct (they should be for the purpose of comparing two
  //simplices.
  unsigned int fnd;
  unsigned int nf=0;
  INT i1,k1=-10,i2,k2=-10;
  INT i,j;
  for (i=0;i<n;i++){
    fnd=0;
    i1=sv1[i];
    for (j = 0;j < n;j++){
      if(i1 == sv2[j]){
	fnd=1;
	break;
      }
    }
    if(fnd) {
      continue;
    } else {
      // not found
      nf++;
      if(nf > 1) return 0;
      k1=i;
    }
  }
  // same with sv1 and sv2 swapped.
  nf=0;
  for (i=0;i<n;i++){
    fnd=0;
    i2=sv2[i];
    for (j=0;j<n;j++){
      if(i2 == sv1[j]){
	fnd=1;
	break;
      }
    }
    if(fnd) {
      continue;
    } else {
      // not found
      nf++;
      if(nf > 1) return 0;
      k2=i;
    }
  }
  /* NOW put the neightbors at the right places ******* */
  if(k1<0||k2<0){
    fprintf(stderr,"\n***ERROR in %s; k1=%d,k2=%d\n\n",__FUNCTION__,k1,k2);
    exit(65);
  } else {
    stos1[k1]=sim2;
    stos2[k2]=sim1;
    return 1;
  }
}
void find_nbr(INT ns,INT nv,INT n,INT *sv,INT *stos)
{
  // find neighboring list
  INT* ivs=NULL, *jvs=NULL, *stosi=NULL, *stosk=NULL, *svi=NULL, *svk=NULL;
  INT i,j,k,jp,jia,jib,iabeg,iaend,ibbeg,ibend,kn1;
  INT n1 = n+1,nv1=nv+1,ns1=ns+1,nsv=ns*n1;
  /* allocate */
  ivs=(INT *) calloc(nv1,sizeof(INT));
  jvs=(INT *) calloc(nsv,sizeof(INT));
  /* transpose to get the vertex--simplex incidence */
  for (i = 0; i < nv1; ++i)
    ivs[i] = 0;
  //
  for (k = 0; k < ns; ++k) {
    kn1=k*n1;
    for (i = 0; i < n1; i++) {
      j = sv[kn1+i] + 2;
      if (j < nv1) ivs[j]++;
    }
  }
  ivs[0] = 0;
  ivs[1] = 0;
  if (nv != 1) {
    for (i= 2; i < nv1; ++i) {
      ivs[i] += ivs[i-1];
    }
  }
  for (i = 0; i < ns; ++i) {
    iabeg = i*n1;
    iaend = i*n1+n1;
    for (jp = iabeg; jp < iaend; ++jp) {
      j = sv[jp]+1;
      k = ivs[j];
      jvs[k] = i;
      ivs[j] = k+1;
    }
  }
  /* for (i = 0; i < nv; ++i) { */
  /*   iabeg = ivs[i]; */
  /*   iaend = ivs[i+1]; */
  /*   fprintf(stdout,"row %d: ",i+1); */
  /*   for (jia = iabeg; jia < iaend; ++jia) { */
  /*     fprintf(stdout,"%d ",jvs[jia]+1); */
  /*   } */
  /*   fprintf(stdout,"\n"); */
  /* } */
  /*
   */
  INT *icp=(INT *) calloc(ns,sizeof(INT));
  INT j123=-10;
  for (i = 0; i < ns; ++i) icp[i] = -1;
  for (i = 0; i < nsv; ++i) stos[i] = -1;
  for (i = 0; i < ns; ++i) {
      iabeg = i*n1;
      iaend = iabeg + n1;
      stosi=stos+iabeg; svi=sv+iabeg;
      for (jia = iabeg; jia < iaend; ++jia) {
	j = sv[jia]; // vertex number
	ibbeg = ivs[j];
	ibend = ivs[j+1];
	// loop over all simplices with this j as a vertex.
	for (jib = ibbeg; jib < ibend; ++jib) {
	  k = jvs[jib];
	  if(k<=i) continue;
	  if (icp[k] != i) {
	    icp[k] = i;
	    kn1=k*n1;
	    stosk=stos+kn1; svk=sv+kn1;
	    if(!cmp_simplex(n1,i,k,svi,svk,stosi,stosk)) continue;
	  }//if
	} //for
      } //for
    } //for (i...
    if (icp) free(icp);
    if (ivs) free(ivs);
    if (jvs) free(jvs);
    //  }
    return;
}

/**********************************************************************/
/*!
 * \fn void faces_attr(subscomplex *subsc)
 *
 * \brief
 *
 * \param
 *
 * \return
 *
 * \note
 *
 */
void faces_attr0(subscomplex *subsc)
{
  INT dim=subsc->nbig, dim1=dim+1;
  size_t nbits=dim*sizeof(REAL);
  scomplex *sc=subsc->parent;
  INT *nbr=sc->nbr, *elf=subsc->elf;
  INT is,it,di,i,j,k;
  INT *fflags=subsc->flags;
  INT l,node,bflag,j1,jf;
  //  find the minimal vertex bounday flag:
  INT minvflag=sc->bndry[0];
  for(i=1;i<sc->nv;i++)
    if(minvflag>sc->bndry[i])  minvflag=sc->bndry[i];
  for(jf=0;jf<subsc->ns;jf++){fflags[jf]=minvflag-1;}
  //  fprintf(stdout,"\nmin vertex flag=%d\n",minvflag);
  /* face_vertex map*/
  REAL *xf=(REAL *)calloc(dim*dim1,sizeof(REAL));
  REAL *snsn=(REAL *)calloc(dim*dim1,sizeof(REAL));
  REAL *areas=(REAL *)calloc(dim1,sizeof(REAL));
  // length of the work space needed : 3 matrices(dim*dim), mass
  // center(dim); pivots(dim),integer permutation(dim) ;
  INT nlength=3*dim*dim*sizeof(REAL)+2*dim*sizeof(REAL)+dim*sizeof(INT);
  void *wrk=(void *)calloc(nlength,sizeof(char));
  REAL fact=sc->factorial;
  REAL s;
  for(it=0;it<sc->ns;it++){
    di=dim1*it;
    //    fprintf(stdout,"\nel=%d, faces:",it);fflush(stdout);
    for(j=0;j<dim1;j++){
      is=nbr[di+j];
      jf=elf[di+j];
      //      fprintf(stdout,"%d(%d) ",jf,fflags[jf]);fflush(stdout);
      j1=jf*((subsc->n+1));
      // vertices in it which are opposite to jf
      if(fflags[jf]>=minvflag)
	{
	  fflags[jf]=0; //this is an interior face because we have
			//already been here.
	  continue;
	}
      l=0;
      bflag=minvflag;
      for(k=0;k<dim1;k++){
	node=sc->nodes[di+k];
	if(k==j) continue;
	subsc->nodes[j1+l]=node;
	if(abs(bflag)<abs(sc->bndry[node])) bflag=sc->bndry[node];
	l++;
      }
      fflags[jf]=bflag;
    }
    for(j=0;j<dim1;j++){
      node=sc->nodes[di+j];
      memcpy((xf+j*dim),(sc->x+node*dim),nbits);
    }
    /* xf[9]=0.; xf[10]=0.; xf[11]=0.; */
    /* xf[6]=1.; xf[7]=0.; xf[8]=0.; */
    /* xf[0]=0.; xf[1]=1.; xf[2]=0.; */
    /* xf[3]=0.; xf[4]=0.; xf[5]=1.; */
    area_face0(dim, fact, xf, snsn,		\
	      areas,(sc->vols+it),		\
	      wrk);
    // get the signs of the normals to agree is it>is;
    for(j=0;j<dim1;j++){
      is=nbr[di+j];
      //      jf=elf[di+j];
      s=chk_sign0(it,is);
      for(i=0;i<dim;i++){snsn[j*dim+i]*=s;}
      //      if(fflags[jf]==22) fprintf(stdout,"\n(%d,%d); s=%f",it,is,s);
    }
    //    print_full_mat(dim1,dim,snsn,"snsn");
    for(j=0;j<dim1;j++){
      jf=elf[di+j];
      memcpy(subsc->normals+jf*dim,(snsn+j*dim),nbits);
      subsc->areas[jf]=areas[j];
    }
    //    exit(11);
  }
  //  haz_scomplex_print(sc,0,"after attr");
  //  haz_subscomplex_print(subsc,0,__FUNCTION__);
  return;
}
/*ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ*/
// fio=fopen(nameio,"w");
// if(fio){
  //    fprintf(stdout,"I/O on file \"%s\"",nameio);
  // for (i=0;i<nv;i++){
  //   fprintf(stdout,"\n%10i: (",i+1);
  //   for(j=0;j<dim-1;j++){
  //     fprintf(stdout,"%16.8e,",x[i*dim+j]);
  //   }
  //   fprintf(stdout,"%16.8e)",x[i*dim+dim-1]);
  // }
  // for (i=0;i<nv;i++){
  //   fprintf(stdout,"\n%10i-->%10i: (",i,p[i]);
  //   for(j=0;j<dim-1;j++){
  //     fprintf(stdout,"%16.8e,",x[p[i]*dim+j]);
  //   }
  //   fprintf(stdout,"%16.8e)",x[p[i]*dim+dim-1]);
  // }
  // fclose(fio);
// } else {
//   fprintf(stderr,"****ERROR:Could not open file \"%s\" for I/O\n",nameio);
//   //return;
//   exit(127);
// }
