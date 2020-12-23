function [A,M,f]=matrix_setup_mass(nh,do_coeff)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %% defines a mass and stiffness tridiagonal matrix corresponding to
  %% the discretization of the eigenproblem Ax = lambda Mx
  %%
  %% -(p(x)u')'=f, u(0)=u(1)=0
  %%
  %% do_coeff=0: p==1  random coefficient
  %% do_coeff ~= 0: p=random pw. constant
  %% 
  %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h=1/(nh+1);
  o1=ones(nh+2,1);
  if(do_coeff)
    o=rand(nh+2,1);
  else
    o=ones(nh+2,1);
  end
  A = (1/h)*spdiags(-o, [1], nh+2, nh+2);
  A = 0.5*(A+A');
  A=A+spdiags(-[sum(A)]',[0],nh+2,nh+2);
  A=A(2:nh+1,2:nh+1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  M = spdiags([o1,4*o1,o1], [-1:1], nh+2, nh+2);
  %%  M=M+spdiags([2*sum(A)]',[0],nh+2,nh+2);
  M=h/6*M;
  M=M(2:nh+1,2:nh+1);
  f=M*o1(2:nh+1); %% constant right hansd side;
return

end

