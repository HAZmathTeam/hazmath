clear all;
pack;
load fort.800;
X=fort(1,:)';
Y=fort(2,:)';
[temp,seq] = size(X)
kend = temp/3
X = reshape(X,[3,kend]);
Y = reshape(Y,[3,kend]);
Z=  reshape(fort(3,:)',[3,kend]);
clear fort
xmin = min(min(X))
ymin = min(min(Y))
xmax = max(max(X))
ymax = max(max(Y))
zmin = min(min(Z))
zmax = max(max(Z))
figure(1);
%axis([xmin xmax ymin ymax zmin zmax]);
%colormap(hsv)
colormap('default')
set(fill3(X,Y,Z,'w'),...
	'erasemode','background','edgecolor','k');
%patch(X,Y,'w')
for i = 1 : 3
for i = 1 : 2	
rot	
end
pause
end	
% Contour lines in regular grid case only. 
load fort.100
kk = sqrt(length(fort))
figure(2)
sol = reshape(fort,[kk,kk]);
contour(sol,30);
%clear all;
%pack;
