void dcsr_print(dCSRmat *A){
  INT i=0, j=0, incre;
  for (i=0;i<A->row;i++){
    incre = A->IA[i+1]-A->IA[i]; 
    for (j=A->IA[i];j<(A->IA[i]+incre);j++){
      printf("(%d,%d)-elem is %.4f\t",i+1,A->JA[j]+1,A->val[j]);
    }
    printf("\n");
  }
  fprintf(stdout,"===========================\n");
}

void icsr_incidence_print(iCSRmat *A){
  INT i=0, j=0, incre;
  for (i=0;i<A->row;i++){
    incre = A->IA[i+1]-A->IA[i]; 
    for (j=A->IA[i];j<(A->IA[i]+incre);j++){
      printf("(%d,%d)-elem is %d\t",i+1,A->JA[j]+1,1);
    }
    printf("\n");
  }
  printf("===========================\n");
}

void dvec_print(dvector *v){
  INT i=0, j=0, n=v->row;
  for (i=0;i<n;i++){
      printf("The %dth-elem is %f\n",i+1,v->val[i]);
  }
  printf("===========================\n");
}

void ivec_print(ivector *v){
  INT i=0, j=0, n=v->row;
  for (i=0;i<n;i++){
      printf("The %dth-elem is %d\n",i+1,v->val[i]);
  }
  printf("===========================\n");
}

void squaremesh(dvector *node,ivector* elem,REAL x0, REAL x1, REAL y0, REAL y1, INT n) {
  REAL hx = (REAL)((x1-x0)/n), hy = (REAL)((y1-y0)/n);
  INT i, j, nv = (n+1)*(n+1), nt = 2*n*n;
  node->row = 2*nv; 
  elem->row = 3*nt;
  node->val = (REAL*)calloc(2*nv,sizeof(REAL));
  elem->val = (INT*)calloc(3*nt,sizeof(INT));
  for (i=0;i<n+1;i++){
    for (j=0;j<n+1;j++){
      node->val[i*(n+1)+j] = x0 + i*hx;
      node->val[nv+i*(n+1)+j] = y0 + j*hy;
      if ((i!=n)&&(j!=n)) {
        elem->val[2*n*i+2*j] = i*(n+1)+j;
        elem->val[nt+2*(n*i+j)] = (i+1)*(n+1)+j+1;
        elem->val[2*nt+2*(n*i+j)] = i*(n+1)+j+1;

        elem->val[2*(n*i+j)+1] = i*(n+1)+j;
        elem->val[nt+2*(n*i+j)+1] = (i+1)*(n+1)+j;
        elem->val[2*nt+2*(n*i+j)+1] = (i+1)*(n+1)+j+1;
      }
    }
  }
}

void icsr_squaremesh(dvector *node, iCSRmat *elem, REAL x0, REAL x1, REAL y0, REAL y1, INT n) {
  REAL hx = (REAL)((x1-x0)/n), hy = (REAL)((y1-y0)/n);
  INT i, j, nv = (n+1)*(n+1), nt = 2*n*n;
  node->row = 2*nv; 
  elem->row = nt;
  elem->col = nv;
  elem->nnz = 3*nt;
  node->val = (REAL*)calloc(2*nv,sizeof(REAL));
  elem->IA = (INT*)calloc(nt+1,sizeof(INT));
  elem->JA = (INT*)calloc(3*nt,sizeof(INT));
  
  for (i=0;i<nt;i++) {
    elem->IA[i+1] = elem->IA[i] + 3;
  }

  for (i=0;i<n+1;i++){
    for (j=0;j<n+1;j++){
      node->val[i*(n+1)+j] = x0 + i*hx;
      node->val[nv+i*(n+1)+j] = y0 + j*hy;
      if ((i!=n)&&(j!=n)) {
         elem->JA[6*(n*i+j)] = i*(n+1)+j;
         elem->JA[6*(n*i+j)+1] = (i+1)*(n+1)+j+1;
         elem->JA[6*(n*i+j)+2] = i*(n+1)+j+1;
    
         elem->JA[6*(n*i+j)+3] = i*(n+1)+j;
         elem->JA[6*(n*i+j)+4] = (i+1)*(n+1)+j;
         elem->JA[6*(n*i+j)+5] = (i+1)*(n+1)+j+1;
      }
    }
  }
}

void icsr_mesh(iCSRmat *elem, ivector* elem0, INT nv) {
  INT i, j, nt = elem0->row/3;
  elem->row = nt;
  elem->col = nv;
  elem->nnz = 3*nt;
  elem->IA = (INT*)calloc(nt+1,sizeof(INT));
  elem->JA = (INT*)calloc(3*nt,sizeof(INT));
  
  for (i=0;i<nt;i++){
    elem->IA[i+1] = elem->IA[i] + 3;
    elem->JA[i*3] = elem0->val[i];
    elem->JA[i*3+1] = elem0->val[i+nt];
    elem->JA[i*3+2] = elem0->val[i+2*nt];
    }
}
