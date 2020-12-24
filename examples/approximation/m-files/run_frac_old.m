do_coeff=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 1.0;
s   = -0.5;
beta  = 0.0;
t   =  0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ea=0;
em=ea;
dim=2;
%% load matrix
disp('Loading matrix A...')
load A.dat
nrowa=A(1,1); ncola=A(1,2); nnza=A(1,3); A=A(2:nnza+1,1:3);
A=sparse(A(:,1)+1,A(:,2)+1,A(:,3),nrowa,ncola);
disp('Loading matrix M...')
load M.dat
nrowm=M(1,1); ncolm=M(1,2); nnzm=M(1,3); M=M(2:nnzm+1,1:3);
M=sparse(M(:,1)+1,M(:,2)+1,M(:,3),nrowm,ncolm);
disp('Loading matrix b and u...')
load b.dat ; b=b(2:end);
load u.dat ; u=u(2:end);

sm=dim*(dim+1)/min(diag(M)); sa=1/norm(A,inf);
sm=1;
sa=1;
[U,d]=eig(sa*A,sm*M);
check_nrm=norm(sa*A*U-sm*M*U*d)
d=diag(d);
n=nrowa;
Ds=spdiags(d.^s,[0],n,n);

As=M*U*Ds*U'*M;

D=diag(d);
I=speye(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnd0=0;
bnd1=max(d);
status0=system('make -C ..');
comm0=sprintf('../aaa.ex <<EOF_FRAC >../m-files/frac.m\n %.2Lf %.2Lf %.2Lf %.2Lf %.2f %.2f\nEOF_FRAC\n',s,t,alpha,beta,bnd0,bnd1);
%%      disp(comm0)
status1=system(comm0)
[res,pol,z,w,f,er]=frac();
m=length(res);
m1=m-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%
% ff=M*b;
% asf=res(m)*ff;
% for j=1:m1
%     asf=asf+res(j)*((sa*A-sm*pol(j)*M)\ff);
% end
% asf=M*asf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DD=res(m)*I;
x00=min(d);
d01=res(m);
for j=1:m1
    d01=d01+res(j)/(x00-pol(j));
    dd=res(j)./(d-pol(j)*ones(n,1));
    Dj=spdiags(dd,[0],n,n);
    DD=DD+Dj;
end
AAs=res(m)*I;
for j=1:m1
    AAs=AAs+res(j)*((sa*A-sm*pol(j)*M)\I);
end
AAs=M*AAs*M;
