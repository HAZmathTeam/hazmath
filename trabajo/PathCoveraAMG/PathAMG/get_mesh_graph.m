function [AG,xy,bpts,ipts]=get_mesh_graph(number_refinements)
%%
%% generates a graph corresponding to a triangulation of a domain. 
%% The triangulation is done with Matlab 
%%  decsg + initmesh + refinemesh.
%% of a domain given by "dgL" below. Right now the dgL is for L-shaped domain. 
%% It has 6 segments (see dgL(2)=number of segments). 
%% in the column array below from position 3 we have:  
%% first the x-coord of the segments and then the y coords. 
%---------------------------------------------------
dgL=[2; 
     6;
     0.;
     1.;
     1.;
     0.5;
     0.5;
     0.;
     0.;
     0.;
     0.5;
     0.5;
     1.;
     1.;
    ];

% to make the domain square
%dgL=[2; 
%     4;
%     0.;
%     1.;
%     1.;
%     0.;
%     0.;
%     0.;
%     1.;
%     1.;
%    ];

dlL=decsg(dgL);

%% Initialize mesh: hmax = inf, will result 
%% in triangulation with all verticies on the boundary.
hmax=0.4;
[xy,ebndry,t]=initmesh(dlL,'hmax',hmax);
xy=jigglemesh(xy,ebndry,t);

%figure(1)
%pdemesh(xy,ebndry,t); axis equal

nref=number_refinements;


for j=1:nref
 [xy,ebndry,t]=refinemesh(dlL,xy,ebndry,t,'longest');
 xy=jigglemesh(xy,ebndry,t);
 xy=jigglemesh(xy,ebndry,t);
end
xy=transpose(xy);
n=length(xy);
%% boundary points
bpts=unique([ebndry(1,:),ebndry(2,:)]); 
ipts=setdiff([1:n],bpts);
ipts=ipts'; bpts=bpts';
%% write T as a sparse matrix
%% 
[idummy,nel]=size(t);
iT=[[1:nel];[1:nel];[1:nel]];
iT=reshape(iT,3*nel,1);
t=reshape(t(1:3,1:nel),3*nel,1);
%% incidence vertex -- triangle 
VT=sparse(t,iT,ones(3*nel,1));
t=reshape(t,3,nel);t=transpose(t);
%% Get the stiffness matrix sparsity:
[n,nT]=size(VT);

AG=sparse(n,n);
AG=sparse(VT*transpose(VT)); %% adjacency matrix
%%generate the graph laplacian
d=diag(AG);
if size(d,1)==1
    d = d';
end
AG = AG-spdiags(d,0,n,n);
[i,j]=find(AG);
%AG = sparse(i,j,-ones(length(i),1),n,n);
AG = sparse(i,j,-1./rand(length(i),1),n,n);
AG = (AG+AG')*0.5;
d = sum(AG);
if size(d,1)==1
    d = d';
end
AG = AG-spdiags(d,0,n,n);

