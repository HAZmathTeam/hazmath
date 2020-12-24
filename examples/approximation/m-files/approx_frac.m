function ur=approx_frac(A,M,b,s_in,dim_in)
    %% ur=approx_frac(A,M,s_in,dim_in)
    %% APPROXIMATES the solution to (M*U*D^s*U*M)*x=b,
    %% if s_in is specified, s=s_in;
    %% if s_in is not specified, s=-0.5 
    %% where A*U=M*U*D. 
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sa=norm(A,inf);
    sm=dim*(dim+1)/min(diag(M));
    bnd0=0;
    bnd1=sa*sm;
    status0=system('make -C ..');
    %% here we input (-s) as a power. 
    comm0=sprintf('../aaa.ex <<EOF_FRAC >../m-files/frac.m\n %.2Lf %.2Lf %.2Lf %.2Lf %.2f %.2f\nEOF_FRAC\n',-s,t,alpha,beta,bnd0,bnd1);
    disp(comm0)
    %% this generates a file m-files/frac.m with the poles and residues. 
    status1=system(comm0)
    [res,pol,z,w,f,er]=frac();
    m=length(res);
    m1=m-1;
    ur=res(1)*(M\b);
    for j=1:m1
        ur=ur+res(j+1)*((A-pol(j)*M)\b);
    end
    return
end

