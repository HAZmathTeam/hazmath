function [ajac]=jacp(p);
dim=size(p,1);
ajac=speye(dim);
if (p(1)==0), return; end
ajac=sparse(dim,dim);
r=p(1);p(1)=1;
ajac(1:dim,1)=p2c(p)
pause
p(1)=r;
above_d=-p(1);
for k=2:dim
  angle=p(k);
  p(k)=0.5*pi-p(k);
  c=p2c(p);
  ajac(k:dim,k)=c(k:dim);
  p(k)=angle;
  above_d=above_d*sin(p(k));
  ajac(k-1,k)=above_d;
end
%if(dim==2)
%    rmanual=[cos(p(2)),-r*sin(p(2));...
%             sin(p(2)),r*cos(p(2))];
%else
%    if(dim==3)
%        rmanual=[cos(p(2)) , -r*sin(p(2)) ,    0;...
%                 sin(p(2))*cos(p(3)) , r*cos(p(2))*cos(p(3)) , -r*sin(p(2))*sin(p(3));...
%                 sin(p(2))*sin(p(3)) , r*cos(p(2))*sin(p(3)) , r*sin(p(2))*cos(p(3))];
%    end
%end
return;
