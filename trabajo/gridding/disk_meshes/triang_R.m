%%%%  READS a triangulation input by triang_W.m
  function [p,e,t]=traing_R(name0);
fid = fopen(name0,'r');
	     [ipar,count]=fscanf(fid,' %i %i %i %i\n',[1,4]);
n=ipar(2),nel=ipar(1),ned=ipar(3)
t=zeros(4,nel);
p=zeros(2,n);
e=zeros(7,ned);
for j = 1 : 4
[t(j:j,1:nel),cnt]=fscanf(fid,' %i ',[1,nel]);
end
[p(1,1:n),cnt]=fscanf(fid,' %e ',[1,n]);
[p(2,1:n),cnt]=fscanf(fid,' %e ',[1,n]);
for j = 1 : 2
[e(j,1:ned),cnt]=fscanf(fid,' %i ',[1,ned]);
end
for j = 3 : 4
[e(j,1:ned),cnt]=fscanf(fid,' %e ',[1,ned]);
end
for j = 5 : 7
[e(j,1:ned),cnt]=fscanf(fid,' %i ',[1,ned]);
end
statusi=fclose(fid)
