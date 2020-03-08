x=[7.0710678118654757e-01 ; 6.1232339957367660e-17 ;-7.0710678118654746e-01];
p=[norm(x);pi/4;pi/4]
pi2=2*pi;
dim=size(p,1);
for iter=1:50
    f=x-p2c(p);
    e=norm(f)
    if(e<1e-14),break; end
    ajac=jacp(p);
    det0=det(ajac);
    full(ajac)
    %%if(abs(det0)<1e-8), ajac=ajac+speye(dim); end
    p=p+ajac\f;
    p(2:dim)=p(2:dim)-floor(p(2:dim)/pi2)*pi2;
    break
end