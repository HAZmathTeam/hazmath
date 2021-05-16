// C program for array implementation of bfs.
void run_bfs(grapha *g,int *root,const int nroots,const int lvlmax)
{
  int iii,i,k,q,v,vi,qbeg,qend,lvl;  
  int *ia=g->ia,*ja=g->ja;
  int *mask=calloc(g->nv,sizeof(int));
  int *il=calloc((g->nv+1),sizeof(int));
  int *jl=calloc(g->nv,sizeof(int));
  int *anc=calloc(g->nv,sizeof(int));
  for(i=0;i<g->nv;++i) mask[i]=0;
  for(i=0;i<g->nv;++i) anc[i]=-1;
  lvl=0;
  il[lvl]=0;
  k=il[lvl];
  for(i=0;i<nroots;++i){
    mask[root[i]]=lvl+1;
    jl[k]=root[i];
    k++;
  }
  il[lvl+1]=k;
  // we need to repeat this
  while(1){
    qbeg=il[lvl];
    qend=il[lvl+1];
    for(q=qbeg;q<qend;++q){
      v = jl[q];
      for(vi=ia[v];vi<ia[v+1];++vi){
	i=ja[vi];
	if (!mask[i]){
	  mask[i]=lvl+1;
	  jl[k]=i;
	  k++;
	  anc[i]=v; // ancestor;
	}
	fprintf(stdout,"\nlvl=%d,v=%d; nbr=%d,mask=%d",lvl,v,i,mask[i]);fflush(stdout);
      }
    }
    lvl++;
    il[lvl+1]=k;
    if(k<=qend) break;
  }
  fprintf(stdout,"\nord (k=%d):",k);
  for(i=0;i<k;i++){
    v=jl[i];
    fprintf(stdout,"\nmask[%d]=%d",v,mask[v]);fflush(stdout);
  }
  for(i=0;i<(lvl+1);i++){
    fprintf(stdout,"\nil[%d]=%d",i,il[i]);
  }
  il=(int *)realloc(il,(lvl+1)*sizeof(int));
  fprintf(stdout,"\nAnother ord::");
  for(i=0;i<lvl;i++){
    qbeg=il[i];
    qend=il[i+1];
    fprintf(stdout,"\nlevel=%d (lvl=%d)::: ",i,lvl);
    for(q=qbeg;q<qend;++q){
      fprintf(stdout," v=%d(anc=%d;lvl=%d)",jl[q],anc[jl[q]],mask[jl[q]]);
    }
  }
  fprintf(stdout,"\n");
  free(anc);
  free(il);
  free(jl);
  free(mask);
  return;
}
int main()
{
  int nv=4;
  graphij *gij=graphij_init(nv);  
  edge_add(gij,0, 1,1.,1.);
  edge_add(gij,0, 2,1.,1.);
  edge_add(gij,1, 2,1.,1.);
  edge_add(gij,2, 0,1.,1.);
  edge_add(gij,2, 3,1.,1.);
  edge_add(gij,3, 3,1.,1.);
  print_graphij(gij);  
  grapha *ga=malloc(sizeof(grapha));
  coo2adj(gij,ga);
  print_grapha(ga);
  //
  int nroots=1;
  int *root=calloc(nroots,sizeof(int));
  root[0]=2;
  fprintf(stdout,"BFS:");
  run_bfs(ga,root,nroots,ga->nv+1);
  return 0;
}
