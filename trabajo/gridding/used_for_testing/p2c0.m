function [c]=p2c0(p);
dim=length(p);
c=zeros(dim,1);
if (p(1)==0), return; end
c(1)=p(1)*cos(p(2));
c(2)=p(1)*sin(p(2));
if(dim==3)
    c(3)=c(2)*sin(p(3));
    c(2)=c(2)*cos(p(3));
end
return;
