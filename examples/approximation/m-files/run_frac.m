%% APPROXIMATES inv(M*U*D^s*U*M), where A*U=M*U*D.
%% KEEP alpha beta and t as they are below or change 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 1.0; beta  = 0.0; t   =  0.5; s=-0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load matrices and vectors all symmetric:
A=get_matrix_haz('A',1);
M=get_matrix_haz('M',1);
%%%%%%%%%%%%%%%%%%%%%%%%%Af=get_matrix_haz('Af',1);
b=get_vector_haz('b');
u=get_vector_haz('u');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE APPROXIMATION:
dim=2;
ur=approx_frac(A,M,b,s,dim);
fprintf(1,'==========================\n |u-ur|=%.12e\n==========================\n',norm(ur-u,inf))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%