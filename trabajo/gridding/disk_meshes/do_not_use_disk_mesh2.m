function [t,x,y,bc]=disk_mesh(ntheta,r,R)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% function [t,x,y,bc]=annul_mesh(ntheta,r,R)
  %%
  %% creates a mesh on an annular domain.
  %% bc each has two columns,
  %% first column is vertex number, second column
  %% is the number along the circle (the polar angle).
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hth=(2*pi)/(ntheta-1);thmesh=transpose([0:hth:2*pi]);
delth=(2-hth)/(2+hth);
nr=fix((log(r)-log(R))/log(delth)-1e-05);
sth=1/delth;
rmesh=spdiags([-sth*ones(nr,1),ones(nr,1)],[-1,0],nr,nr)\[r;zeros(nr-1,1)];
nr=length(rmesh);
rmesh=[rmesh;(rmesh(nr)+R)*0.5;R];
nr=length(rmesh);
ith=[1:ntheta-1]'; jr=[1:nr-1]';
nq=length(ith)*length(jr);
np=(ntheta-1)*nr; %% vertices
[ith,jr]=meshgrid(ith,jr);
ith=reshape(ith,nq,1);
jr=reshape(jr,nq,1);

s=(ith-1)*nr+jr;

t = sparse(4*nq,3);

t(1:4:4*nq-3,1) = s; %(jj-1)*n+ii;
t(1:4:4*nq-3,2) = s+1; %(jj-1)*n+ii+1;

t(2:4:4*nq-2,1) = s+1; %(jj-1)*n+ii+1;
t(2:4:4*nq-2,2) = s+nr+1; %(jj)*n+ii+1;

t(3:4:4*nq-1,1) = s+nr+1; %(jj)*n+ii+1;
t(3:4:4*nq-1,2) = s+nr; %(jj)*n+ii;

t(4:4:4*nq  ,1) = s+nr; %(jj)*n+ii;
t(4:4:4*nq  ,2) = s; %(jj-1)*n+ii;

[i0,j0,s0]=find(t);s0=floor(s0);
jjs0=find((s0>(ntheta-1)*nr) & (s0<(ntheta*nr+1)));
s0(jjs0)=mod(s0(jjs0),(ntheta-1)*nr);
t = sparse(i0,j0,s0);clear i0;clear j0; clear s0;

%% vertices:
ith=[1:ntheta-1]'; jr=[1:nr]';
np=length(ith)*length(jr);
[ith,jr]=meshgrid(ith,jr);
ith=reshape(ith,np,1);jr=reshape(jr,np,1);

s1=(ith-1)*nr+jr;

%% boundary vertices, the second coordinate is the theta coordinate.
bc=[s1(find(jr==nr)),ith(find(jr==nr))]; %% outer circle

s2=np+[1:nq]';

t(1:4:4*nq-3,3) = s2; %center;
t(2:4:4*nq-2,3) = s2; %center;
t(3:4:4*nq-1,3) = s2; %center;
t(4:4:4*nq  ,3) = s2; %center;

t=full(t);
nt=size(t,1);

x(s1)=[rmesh(jr).*cos(thmesh(ith))];
y(s1)=[rmesh(jr).*sin(thmesh(ith))];

c1=sparse(t(:,3)-np,t(:,1),ones(nt,1),nq,np); d=c1*ones(np,1);
x(s2)=(c1*x')./d;
y(s2)=(c1*y')./d;
nv=size(x,1);
if(nv<2) x=x';y=y';nv=size(x,1);end
nv=nv+1;
x(nv)=y(nv)=0.
t=delaunay(x,y);
flag=print_mesh_hazmath(x,y,t);
flag
return
