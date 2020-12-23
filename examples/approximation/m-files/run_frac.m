do_coeff=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 1.0;
s   = -0.5;
beta  = 0.0;
t   =  0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lev =  10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ea=0;
em=ea;
nh=2^lev-1;
x=[linspace(0,1,nh+2)]';
x=x(2:nh+1);
dim=1;
[A,M,f]=matrix_setup_mass(nh,do_coeff);

sm=dim*(dim+1)/min(diag(M)); sa=1/norm(A,inf);

if(nnz(A)/size(A,1)<3) 
    [U,d]=eig(sa*A,sm*M);
    check_nrm=norm(sa*A*U-sm*M*U*d)
    d=diag(d);
    Ds=spdiags(d.^s,[0],nh,nh);
    As=M*U*Ds*U'*M;
else
    %% load matrix
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnd0=0;
bnd1=sm*sa;
status0=system('make -C ..');
comm0=sprintf('../aaa.ex <<EOF_FRAC >../m-files/frac.m\n %.2Lf %.2Lf %.2Lf %.2Lf %.2f %.2f\nEOF_FRAC\n',s,t,alpha,beta,bnd0,bnd1);
%%      disp(comm0)
status1=system(comm0)
[res,pol,z,w,f,er]=frac();
m=length(res);
m1=m-1;
D=diag(d);
n=size(A,1);
I=speye(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs=randn(n,1);
%%
ff=M*rhs;
asf=res(m)*ff;
for j=1:m1
    asf=asf+res(j)*((sa*A-sm*pol(j)*M)\ff);
end
asf=M*asf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DD=res(m)*I;
% for j=1:m1
%     dd=res(j)./(d-pol(j)*ones(n,1));
%     DD=DD+spdiags(dd,[0],n,n);
% end
%for j=1:m1
%    dd=res(j)./(d-pol(j)*ones(n,1));
%    Dj=spdiags(dd,[0],n,n);
%    DD=DD+Dj;
%end
%AAs=M*U*DD*U'*M;
%%AAs=res(m)*I;
%%for j=1:m1
%%    AAs=AAs+res(j)*((sa*A-sm*pol(j)*M)\I);
%%end
%% AAs=M*AAs*M;
