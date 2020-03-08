function [AG,VT,n1,nT1] = tri_to_graph1(t);

[idummy,nel]=size(t);

iT=[[1:nel];[1:nel];[1:nel]];

iT=reshape(iT,3*nel,1);

t=reshape(t(1:3,1:nel),3*nel,1);

%% incidence vertex -- triangle 

VT=sparse(t,iT,ones(3*nel,1));

t=reshape(t,3,nel);t=transpose(t);

%% Get the stiffness matrix sparsity:

[n,nT]=size(VT)

n1=n;nT1=nT;
AG=sparse(n,n);
AG=sparse(VT*transpose(VT));
AG=AG-spdiags(diag(AG),0,n,n);
[i,j,s]=find(AG);s=-ones(size(s));
AG=sparse(i,j,s,n,n);
AG=AG-spdiags(AG*ones(n,1),0,n,n);



