n=32
new_graph=0;

fname_vertex='g_vertex.dat';
fname_simplex='g_simplex.dat';
fname_x='x_coord.dat';
if(new_graph==0)
    fid=fopen(fname_x,'r');
else
    fid=-1;
end
if(fid<0)    
    fid=fopen(fname_x,'w');
    x=randn(n,2);
    fprintf(fid,'%.16e %.16e\n',x'); 
else
    x=load(fname_x); 
end
tri = delaunay (x(:,1), x(:,2));
%%plot (x(:,1), x(:,2), "r*");
x1min=min(x(:,1));
x1max=max(x(:,1));
x2min=min(x(:,2));
x2max=max(x(:,2));
%%
hx=(x1max-x1min)/sqrt(n);
hy=(x2max-x2min)/sqrt(n);
dim=size(x,2);
[nvt,ntt,tri,vt,iwt]=clean_sc(tri,dim);
[nei,neb,tt,et,ev]=tt_ft_fv(tri,vt,nvt,ntt,dim);
vt=vt';
xt=(1/(dim+1))*(vt'*x);
vv=ev'*ev;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(fname_vertex,'w');
write_csr(vv,fid);
fclose(fid);
fid=fopen(fname_simplex,'w');
%%tt=tt+speye(size(tt,1));
write_csr(tt,fid);
fclose(fid);
% [vv,status]=read_csr(fname_vertex); vv=vv-diag(diag(vv));
% [tt,status]=read_csr(fname_simplex); tt=tt-diag(diag(tt));
commv=sprintf('./lexi %s > ordv.m\n',fname_vertex);
comms=sprintf('./lexi %s > ords.m\n',fname_simplex);
disp(commv)
%% this generates a ordv.m with the poles and residues. 
statusv=system(commv)
disp(comms)
%% this generates a ords.m with the poles and residues. 
statuss=system(comms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% defines p, pinv and level;
ordv
pv=p;pinvv=pinv;levelv=level;trv=tree;
%% defines p, pinv and level for the dual graph;
ords
%% defines p, pinv and level for the dual graph;
ps=p;pinvs=pinv;levels=level;trs=tree;
clear p; clear pinv; clear level; clear tree;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
levmax=max(levelv)
%plot_lexi(vv,x,pv,levelv,2,1)
%plot_lexi(tt,xt,ps,levels,2,2)
%%
%k=1;
%for k=1:length(p)
%    [p(k),tree(p(k))]
%end

