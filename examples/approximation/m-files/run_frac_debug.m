%% DEBUG for rational approx; 
%% KEEP alpha beta and t as they are below or change 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 1.0; beta  = 0.0; t   =  0.0; s=0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load matrices and vectors all symmetric:
A=get_matrix_haz('A',1);
M=get_matrix_haz('M',1);
%%%%%%%%%%%%%%%%%Af=get_matrix_haz('Af',1);
b=get_vector_haz('bm');
u=get_vector_haz('um');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE APPROXIMATION:
dim=2;
[ur,As,Asinv_a]=approx_frac_debug(A,M,b,s,dim);
I=speye(size(A,1),size(A,1));
%fprintf(1,'==========================\n |u-ur|=%.12e\n==========================\n',norm(ur-u,inf))
fprintf(1,'==================================\n |As*Asinv_a-I|=%.12e\n==================================\n',norm(As*Asinv_a-I,inf))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
