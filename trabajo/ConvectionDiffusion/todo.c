
/**************************************************/
// tangent vectors are the columns of bs[].  then we take i and j
// always in one and the same direction by doing i<j or things like
// that...
    iaa=ia[i]-1;
    iab=ia[i+1]-1;
    for (jk=iaa; jk<iab; jk++){
      j=ja[jk]-1;
      /*      fprintf(stdout,"aaaaa %i,%i\n",i,j); */
      if(i != j){
	xj=xyz->x[j];
	yj=xyz->y[j];
	if(dim>2) 
	  zj=xyz->z[j];
	/* 
	   compute the advection field at the middle of the edge and
	   then the bernoulli function
	*/
	te[0]  = xi - xj;
	te[1]  = yi - yj;
	xm[0] = (xi + xj)*0.5e+0;
	xm[1] = (yi + yj)*0.5e+0;
	if(dim>2){
	  te[2]  = zi - zj;
	  xm[2] = (zi + zj)*0.5e+0;
	}
	vector_val_ad(ad,xm,0.0,&edge_flag);
	bte=0.;
	for(jdim=0;jdim<dim;jdim++){
	  bte += ad[jdim]*te[jdim];
	}
	/* alpe=a(xmid)\approx harmonic_average=|e|/(int_e 1/a); */
	/* should be computed by quadrature in general for 1/a(x). */
	/* a_{ij}=B(beta.t_e/harmonic_average)*harmonic_average*omega_e; */
	/* 	  for (i,j):    B(beta.t_e/harmonic_average);    */
	/* 	  for (j,i): B(-beta.t_e/harmonic_average). */
	/* 
	   This is direction independent, because t_e=xi-xj; If places
	   of i and j are switched, then we are working on a_{ji} with
	   the t_e of opposite sign.
	 */	
	scalar_val_d(&alpe,xm,0.0,&edge_flag);
	a[jk] *= (alpe*(bernoulli(bte/alpe))); 
	//	fprintf(stdout,"\ndd=%22.14e",(alpe*(bernoulli(bte/alpe))));
	/* the diagonal is equal to the negative column sum;*/

	diag0.val[j]-=a[jk];
      }
    } 
  }
  /* Another loop to set up the diagonal equal to whatever it needs to
     equal. */
  for (i = 0; i < nv; i++){
    iaa=ia[i]-1;
    iab=ia[i+1]-1;
    for (jk=iaa; jk<iab; jk++){
      j = ja[jk]-1; 
      if(i == j) a[jk]=diag0.val[i];
    }
  }
  /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
