
/**************************************************/
for(i=0;i<dim1;i++){
  id=i*dim; // this is for xs, so only dim columns
  for(j=0;j<dim1;j++){
    if(i==j) continue;
    /* \beta\cdot \tau_e */
    ij=i*dim1+j; // this is for stiff local.
    jd=j*dim; // this is for xs, so only dim columns
    bte=0.; 
    if(nop(j)>nop(i))
      for(k=0;k<dim;k++) bte+=ad[k]*(xs[jd+k]-xs[id+k]);
    else if(nop(j)!=nop(i)) // reverse tau
      for(k=0;k<dim;k++) bte+=ad[k]*(xs[id+k]-xs[jd+k]);
      //old for the above... for(k=0;k<dim;k++) te[k]=xs[jd+k]-xs[id+k];
      //old for the above... for(k=0;k<dim;k++) te[k]=xs[id+k]-xs[jd+k];
    /*midpoints*/
    // mid points on the edges if needed!    for(k=0;k<dim;k++) xm[k]=0.5*(xs[jd+k]+xs[id+k]);
    aloc[ij] *= bernoulli(bte/alpe); 
    //	fprintf(stdout,"\ndd=%22.14e",(alpe*(bernoulli(bte/alpe))));
    /* the diagonal is equal to the negative column sum;*/    
    d[j]-=aloc[ij];
  }
 }
/* Another loop to set up the diagonal equal to whatever it needs to
   equal. */
for (i = 0; i < dim1; i++)
  aloc[i+dim1+i]=diag0.val[i];

//	vector_val_ad(ad,xm,0.0,&edge_flag);
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
/* scalar_val_d(&alpe,xm,0.0,&edge_flag); */
  /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
