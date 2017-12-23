function plotgraf1(x,y,n);
xy = [x,y];
load fort.201
a = sparse(fort(:,1),fort(:,2),fort(:,3));
clear fort
gplot(a,xy);
if(n<300)
plot(x,y,'c.'); 
else
plot(x,y,'c.'); 
end
axis off
return