function [ A ] = compute_affinity( X_k, A_ori, threshold )
% Calculate affinity score based on section 3.3.2 in OREN E. LIVNE and ACHI
% BRANDT's paper LEAN ALGEBRAIC MULTIGRID (LAMG): FAST GRAPH LAPLACIAN LINEAR SOLVER
%
%   Joanne Lin @ Tufts University Math Dept.

% input: X_k - a matrix of current estimation
%        [x_k(i) = f(x,y) = x | x_k(i) = f(x,y) = y | x_k(i) = f(x,y) = x + y | x_k(i) = f(x,y) = x - y], 
%        where f(x,y) is of dimention sqrt(length(x_k)) * sqrt(length(x_k)).
%        x = floor(i/sqrt(length(x_k))), y = mod(i, sqrt(length(x_k)))
%        We reshape each initial guess into a long vector of length n. so X_k is of size n * 4.

% output: A - Adjacency matrix of the new mesh after adding some edges has small 
% affinity c_uv between vertex u and vetex v (the smaller c_uv is, the closer are u 
% and v) 

n = size(X_k,1);
%A = zeros(n,n); % Adjacency matrix of the new mesh
%C = zeros(n,n); % Affinity scores

% find nonzeros of A^2
[row_idx, col_idx, ~] = find(triu(A_ori,1));

% compute affinity
% C = 1 - ((  sum(X_k(row_idx,:).*X_k(col_idx,:),2).^2 )./ ...
%     (sum(X_k(row_idx,:).*X_k(row_idx,:),2).*sum(X_k(col_idx,:).*X_k(col_idx,:),2) ));
% 
% % check affinity
% flag =  (C<threshold);

% construct A
% val = -1./(max(abs(X_k(row_idx(flag),:) - X_k(col_idx(flag),:)),[],2) + 0.001);
% A = sparse(row_idx(flag), col_idx(flag), val, n, n);

val = -1./(max(abs(X_k(row_idx,:) - X_k(col_idx,:)),[],2) + 1e-6);
A = sparse(row_idx, col_idx, val, n, n);
A = A + A';
%A = sparse(A);

end

