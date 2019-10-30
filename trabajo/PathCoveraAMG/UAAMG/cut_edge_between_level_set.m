function [ Adj_level_set ] = cut_edge_between_level_set( Adj, level_set )
% cut edges between level sets

% @ Xiaozhe Hu

% get size
N = size(Adj, 1);

% get nonzeros
[row_idx, col_idx, val] = find(Adj);

% get index for edges within level sets
idx = ( level_set(row_idx) == level_set(col_idx) );

% get adjancy matrix for level sets
Adj_level_set = sparse(row_idx(idx), col_idx(idx), val(idx), N, N);

end

