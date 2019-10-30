
% get matrix
nsize = 16;
A = assembleLaplace(nsize);

% get size
N = size(A,1);

% get adjancy matrix
D = spdiags(diag(A), 0, N, N);
Adj = D - A;

% get smooth error
iterParam.print_level = 0; % how much information to print for iterative solver
iterParam.max_it = 100; % maximal number of iterations that is allowed
iterParam.tol = 1e-8;  % Tolerance for relative residual 

b = zeros(N,1);
rng(1);
x = rand(N,1);
x = Simple_Solve(A, b, x, iterParam);
smooth_error = x/norm(x);

err0 = sqrt(x'*A*x)

% find level set
[~,max_idx] = max(smooth_error);
[~,min_idx] = min(smooth_error);
s = [max_idx, min_idx];
%s = max_idx;

%A2 = A*A;
%A2adj = A2 - spdiags(diag(A2), 0, N, N);
%[ level_set ] = BFS_level_set( A2adj, s );
[ level_set ] = BFS_level_set( Adj, s );

% get rid of edges between level_sets
A2 = A*A;
A2adj = A2 - spdiags(diag(A2), 0, N, N);
[ Adj_level_set ] = cut_edge_between_level_set( A2adj, level_set );

% generate paths
%[ path ] = level_set_to_path( level_set );
[ path ] = generate_BFS_path( Adj_level_set, level_set );

% get aggregation from path
num_agg = length(path);
aggregation = zeros(N,1);
for i=1:num_agg
    aggregation(path{i}) = i;
end

% form prolongation
%aggregation = level_set;
%num_agg = max(level_set);
[ P ] = generate_unsmoothed_P( aggregation, num_agg );

Ac = P'*A*P;

% test error
DL = tril(A);
DU = triu(A);

%sqrt(x'*A*x)
x = forward_gs(A, b, x, DL, 1);
%sqrt(x'*A*x)
r = b - A*x;
rc = P'*r;
ec = Ac\rc;
x = x + P*ec;
%sqrt(x1'*A*x1)
x = backward_gs(A, b, x, DU, 1);
err1 = sqrt(x'*A*x)


err1/err0
