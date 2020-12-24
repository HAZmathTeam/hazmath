function [ur,As,Asinv_a]=approx_frac_debug(A,M,b,s_in,dim_in)
    %% DEBUGGING: more outputs 
    %% [ur,As,Asinv_a,Ds,Dsinv_a]=approx_frac_debug(A,M,b,s_in,dim_in)
    %%
    %% if s_in is specified, s=s_in;
    %% if s_in is not specified, s=-0.5
    %%
    %% D and U below satisfy A*U=M*U*D, U'*M*U=U*U'*M=I
    %% OUTPUT:
    %% As = (M*U*D^s*U*M);    
    %% ur APPROXIMATES the solution to As*x = b, i.e. (M*U*D^s*U*M)*x=b,
    %% Asinv_a rational function approx to inv(As);    
    %% Ds is D^s 
    %% Dsinv_a approximates inv(Ds) with rational function.    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha = 1.0; s=-0.5; beta  = 0.0; t   =  0.5;
    if(nargin<3)
        ur=NaN;
        return
    end
    if(nargin>3)
        s=s_in;
    end
    dim=2;
    if(nargin>4)
        dim=dim_in; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U,d]=eig(A,M);
    check_nrm=norm(A*U-M*U*d)
    d=diag(d);
    n=size(A,1);
    Ds=spdiags(d.^s,[0],n,n); %% exact Ds(matlab)
    As=M*U*Ds*U'*M; %% exact As(matlab)
    As_inv=U*(spdiags(d.^(-s),[0],n,n))*U';
    D=diag(d);
    I=speye(n);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OLD: bnd0=0; bnd1=3*max(d);
    sm=dim*(dim+1)/min(diag(M)); sa=1/norm(A,inf);
    bnd0=0;bnd1=sm/sa; %%norm(A,inf)*dim*(dim+1)/min(diag(M))
    status0=system('make -C .. clean ; make -C ..');
    %% watch s is (-s) to approximate inverses;
    comm0=sprintf('../aaa.ex <<EOF_FRAC >../m-files/frac.m \n %.2f %.2f %.2f %.2f %.2f %.2f\nEOF_FRAC\n',-s,t,alpha,beta,bnd0,bnd1);
    %%      disp(comm0)
    status1=system(comm0)
    [res,pol,z,w,f,er]=frac();
    m=length(res);
    m1=m-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dsinv_a=res(1)*I; %% approximation of D^(-s)=inv(Ds)
    % for j=1:m1
    %     dd=res(j+1)./(d-pol(j)*ones(n,1));
    %     Dj=spdiags(dd,[0],n,n);
    %     Dsinv_a=Dsinv_a+Dj;
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Asinv_a=res(1)*(M\I); %% approx (A^(-s))
    for j=1:m1
        Asinv_a=Asinv_a+res(j+1)*((A-pol(j)*M)\I);
    end
    ur=res(1)*(M\b);
    for j=1:m1
        ur=ur+res(j+1)*((A-pol(j)*M)\b);
    end
    return
end

