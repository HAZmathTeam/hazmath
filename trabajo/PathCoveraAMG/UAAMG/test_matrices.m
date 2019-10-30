function test_matrices( filename,solver_type, prec_type, nsize, epsilon)

% test AMG
% 
% @ Xiaozhe Hu, Tufts University

% load AMG package
addpath(genpath('~/01042018/UAAMG'));
addpath(genpath('~/01042018/PathAMG'));
addpath(genpath('~/DSD/20170515'));
addpath(genpath('~/UAAMG/newWSGraphs'));

%---------------------
% load matrices
%---------------------
load(filename)
display(filename)

%nsize
%epsilon
% %A = assembleLaplace(20);
%A = assembleGraphLaplace(nsize);
%A = assembleAnisotropicGraphLaplace(nsize,epsilon);


%-------------------------pre-process------------------%
% [ L_cell,b_cell ] = split_conncomp( 'L_and_b_with_2674_alternatives_13849_voters_of_function_1.mat' );
% A = L_cell{1};
% b = b_cell{1};
% load('../data/L1_largest_comp.mat');
%  A = L;
% 
% load('pmf_lhs')
% load('pmf_rhs')
% A = A(13201:end,13201:end);
% b = b(13201:end);
% 
%% get the largest connected component
G = graph(A);
S = conncomp(G);
i = mode(S);
largest_conn = find(S==i);
A = A(largest_conn,:);
A = A(:,largest_conn);

% get a undirected (symmetric) positive weighted graph
A_sym = (A + A')/2;
A = abs(A_sym);

% take out diagonals of A (if any)
n = size(A,1);
d = diag(A);
A = (A-spdiags(d,0,n,n));
d = sum(A);
if size(d,1)==1
  d = d';
end


L = spdiags(d,0,n,n)-A;

% load('A_c.mat')
% L = A_c;

% get normalized graph laplacian
% D_half = spdiags(1./sqrt(d),0,n,n);
% L = D_half * L * D_half;

% calculate condition number
%lambda_n = eigs(L,1,'lm');
%lambda_sm = eigs(L,10,'sm');
%ind = find(lambda_sm>1e-6);
%lambda_2 = lambda_sm(ind(1));
%cond = lambda_n/lambda_2
%return;

%L = A;

N = size(L,1);

%---------------------
% right hand side
%---------------------
%b = zeros(N,1);

% x_2
%B = sparse((eye(size(L,1))-(tril(L)+diag(diag(L)))*L)*(eye(size(L,1))-(triu(L)+diag(diag(L)))*L));
% B = sparse(eye(size(L,1))-diag(1./diag(L))*L);
% [V,~] = eigs(B,2,'lm');
% x_2 = V(:,2);
 %exact_x = rand(N,1);
 exact_x = zeros(N,1);
% %exact_x = x_2;
 b = L*exact_x;
%  load('r_c.mat')
%  b = r_c;

%---------------------
% AMG parameters
%---------------------
amgParam.print_level = 1; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
% setup phase parameters
amgParam.max_level = 20; % maximal number of level in AMG
amgParam.coarsest_size = 100; % size of the coarest level

amgParam.strong_connection = 0.08; 
amgParam.agg_type = 'HEC';  % 'HEC': heavy edge coarsening | 'MIS': maximial indepent set coarsening 
                                      % | 'path_matching': path cover coarsening
                                      % | 'path_matching_aff': path cover coarsening using affinity
                                      % | 'path_matching_adapt': adaptive path cover coarsening using affinity
                                      % | 'regular_matching': regular matching coarsening
amgParam.agg_radius = 1; % radius of the aggregates, only used when amgParam.agg_type = 'MIS'

% solve pahse parameters
amgParam.cycle_type = 'V'; % 'V': V-cycle | 'W': W-cycle | 'nV': n-fold V cycle | 'K': K-cycle (nonlinear AMLI-cycle)

amgParam.coarse_it = 2; % Number of iterations on coarse levels. Used when amgParam.cycle_type = 'nV' or 'K'
amgParam.coarse_krylov_type = 'GCG'; % Type of Krylove method when amgParam.cycle_type = 'K'
                                     % 'GCG': generalized CG | 'FGRMES': flexible GMRes 
                                     
amgParam.n_presmooth = 1; % number of presmoothing      
amgParam.n_postsmooth = 1; % number of postsmoothing

amgParam.prec_max_it = 1; % when AMG is used as preconditioner, how cycles will be applied
amgParam.max_it = 2500;  % when AMG is used as standalone solver, maximal number of iterations that is allowed
amgParam.tol = 1e-10;    % when AMG is used as standalone solver, tolerance for the reletive residual

amgParam.ILU_level = 0;  % how many levels use ILU smoother (from fine to coarse), 0 means no ILU smoother

% ILU parameters (only used when amgParam.ILU_level > 0)
amgParam.ILU_type = 'nofill';  % nofill, crout, ilutp
amgParam.droptol = 0.01; % drop tolerance for ILU when ILU_type = 'crout' or 'ilutp'
% do NOT touch the following three
amgParam.milu = 'off';  % row, col, off
amgParam.udiag = 1;
amgParam.thresh = 1;

%---------------------
% Iterative parameters
%---------------------
iterParam.solver_type = 'AMG_Adapt'; % 'SL': GS as standalone solver | 'AMG': AMG as standalone solver 
                              % | 'CG': conjugate gradiant | 'FGMRES':
                              % flexible GMRes | AMG_path: path cover AMG
                              % solver | AMG_Adapt: adaptive AMG solver
                              % 'AMG_path': path cover solver
iterParam.prec_type = 'NULL'; % 'NULL': no preconditioner | 'AMG': AMG precondtioner
iterParam.print_level = 0; % how much information to print for iterative solver
                           % 0: print nothing | positive number print information 

iterParam.max_it = 100; % maximal number of iterations that is allowed
iterParam.tol = 1e-8;  % Tolerance for relative residual 
iterParam.restart = 100; % restart number for flexible GMRes

% initial guess
% x = construct_x1_1(nsize);  
% x = construct_x1_2(20);
% x = construct_x1_3(nsize);
% x = construct_x1_4(20);
rng(5)
x = rand(N,1);
%x = zeros(N,1);

%---------------------
% setup phase
%---------------------
if ( strcmp(iterParam.solver_type, 'AMG') || strcmp(iterParam.prec_type, 'AMG') )
%if  strcmp(iterParam.prec_type, 'AMG') 
    amgData = AMG_Setup(L, amgParam);
end

if  strcmp(amgParam.agg_type, 'path_matching_adapt' )
    [amgData,x] = AMG_pathSetup_adapt(L, b, x, amgParam);
end
%---------------------
% solve phase
%---------------------

switch iterParam.solver_type
    
    case 'SL', % use simple linear iterative solver
        x = Simple_Solve(L, b, x, iterParam);

    case 'AMG', % use AMG directly
        x = AMG_Solve(amgData, b, x, amgParam);
    
    case 'AMG_path' % use AMG-path solver
        x = AMG_pathSolve(L, b, x, amgParam);
        
    case 'AMG_Adapt' % use adaptive AMG solver
        x = AMG_Solve_Adapt(L, b, x, amgParam, exact_x);

    otherwise, % use preconditioned Krylov methods
        x = Krylov_Solve(L, b, x, iterParam, amgParam, amgData);
        
end
