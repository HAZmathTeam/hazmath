n0=7;
%%::old::%%n=20; x=rand(n,1); y=rand(n,1); tri=delaunay(x,y);

%% a uniform mesh on a triangle: 
% pick two random angles in pi/6, pi/2
iter=0;
while(1)
 a1=1/6+rand(1)*(0.5-1/6);
 a2=1/6+rand(1)*(0.5-1/6);
%%a1=70/180;a2=30/180;
a3=1-(a1+a2);
 iter=iter+1;
if((a3<=0.5 && a3>=1/6) || (iter>100)) break; end
end
a1=a1*pi;a2=a2*pi;a3=a3*pi;

[t,x,y]=do_tri(n0,a1,a2,a3);

[A,VT,nv,nt]=tri_to_graph1(t');

%% get the adjacency matrix;
[iw,jw]=find(triu(A,1));
ne=length(iw);
iev=[[1:ne]';[1:ne]'];jev=[iw;jw];
EV=sparse(iev,jev,ones(2*ne,1),ne,nv);
[iw,jw]=find(EV*VT == 2);
ET=sparse(iw,jw,ones(length(iw),1),ne,nt);

%% we need to find the T,T matrix:

At = transpose(VT)*VT;
[i,j,s] = find(At==2);
At=sparse(i,j,s,nt,nt);
xt=1/3*transpose(VT)*[x,y];
yt=xt(:,2);
xt=xt(:,1);
gplot(At,[xt,yt])


