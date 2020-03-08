function [t,xy]=disk_mesh(R,level)
%% using distmesh: the mesh size is R*2^(-level-1)
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off','all');
addpath("distmesh/");
fd=@(p) sqrt(sum(p.^2,2))-R;
h=R*2^(-(level+1));
disp(['Disk domain: R=',num2str(R),'; step size h=',num2str(h)]);
[xy,t]=distmesh2d(fd,@huniform,h,[-R,-R;R,R],[]);
size(xy)
disp(['Elements=',int2str(size(t,1)),'; Verices=',int2str(size(xy,1))]);
return
