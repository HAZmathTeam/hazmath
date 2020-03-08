function [c]=p2c(p);
dim=length(p);
c=zeros(dim,1);
if (p(1)==0), return; end
rho=p(1);
cend=rho;
for i=1:(dim-1)
    c(i)=rho*cos(p(i+1));
    for j=1:(i-1)
        c(i)=c(i)*sin(p(j+1));
    end
    c
    cend=cend*sin(p(i+1));
end
c(dim)=cend;
return;
