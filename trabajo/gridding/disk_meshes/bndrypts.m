%%function [VT,n1,nT1] = tri_to_graph1(t);
R=1; level=0;
[t,x,y]=disk_mesh(R,level)
t=transpose(t); %% this is fake just to get the incidence matrices.
[idummy,nel]=size(t);
iT=[[1:nel];[1:nel];[1:nel]];
iT=reshape(iT,3*nel,1);
t=reshape(t(1:3,1:nel),3*nel,1);
%% incidence vertex -- triangle 
VT=sparse(t,iT,ones(3*nel,1));
t=reshape(t,3,nel);t=transpose(t);
[nv,nt]=size(VT)
%% incidence: edge--vertex
A=sparse(VT*transpose(VT));
[iw,jw]=find(triu(A,1));
ne=length(iw);
iev=[[1:ne]';[1:ne]'];jev=[iw;jw];
EV=sparse(iev,jev,ones(2*ne,1),ne,nv);
%% incidence: edge--triangle
[iw,jw]=find(EV*VT == 2);
ET=sparse(iw,jw,ones(length(iw),1),ne,nt);

%% we need to find the T,T matrix:
TT = transpose(VT)*VT;
[i,j,s] = find(TT==2);
TT=sparse(i,j,s,nt,nt);


