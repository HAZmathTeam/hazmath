function bcodes = bndrypts(t,c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% function bcodes = bndrypts(t,c);
%% finds the boundary nodes and sets their code to "c"
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=transpose(t); %% this is fake just to get the incidence matrices.
[dim1,nel]=size(t);
dim=dim1-1;
iT=[[1:nel];[1:nel];[1:nel]];
iT=reshape(iT,dim1*nel,1);
t=reshape(t(1:dim1,1:nel),dim1*nel,1);
%% incidence vertex -- triangle
VT=sparse(t,iT,ones(dim1*nel,1));
t=reshape(t,dim1,nel);
t=transpose(t);
[nv,nt]=size(VT);
%% incidence: edge--vertex
A=sparse(VT*transpose(VT));
[iw,jw]=find(triu(A,1));
ne=length(iw);
iev=[[1:ne]';[1:ne]'];jev=[iw;jw];
EV=sparse(iev,jev,ones(2*ne,1),ne,nv);
%% incidence: edge--triangle
[iw,jw]=find(EV*VT == dim);
ET=sparse(iw,jw,ones(length(iw),1),ne,nt);

ou=find((ET*ones(nt,1))<2);
[~,bp]=find(EV(ou,1:nv));
bcodes=zeros(nv,1);
bcodes(bp)=ones(length(bp),1);

%% no need to find the T,T matrix:
%TT = transpose(VT)*VT;
%[i,j,s] = find(TT==dim);
%TT=sparse(i,j,s,nt,nt);

return
