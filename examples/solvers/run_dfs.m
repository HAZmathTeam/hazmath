load vector_laplace_coarse.txt;
Ac=sparse(vector_laplace_coarse(:,1)+1,vector_laplace_coarse(:,2)+1,vector_laplace_coarse(:,3));
clear vector_laplace_coarse;
nc=size(Ac,1);
[iblkc,jblkc,nblkc,errc]=dfs0(Ac,nc);

load vector_laplace_fine.txt;
A=sparse(vector_laplace_fine(:,1)+1,vector_laplace_fine(:,2)+1,vector_laplace_fine(:,3));
clear vector_laplace_fine;
n=size(A,1);
[iblk,jblk,nblk,err0]=dfs0(A,n);
