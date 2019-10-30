function [ A ] = compute_similarity_full( coord  )
% similarity (full matrix)
%
% @ Xiaozhe Hu, Tufts University

% get size
[n, ~] = size(coord);

% compute distance
S = dist(coord');

% find nonzeros
[row_idx, col_idx, val] = find(triu(S,1));

% construct A
val = -1./val;
A = sparse(row_idx, col_idx, val, n, n);
A = A + A';

end

