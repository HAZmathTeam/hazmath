function [ path ] = dfs_path1( A, startnode )
% find path given adjacency matrix using dfs
%   Inputs: A - adjacency matrix
%           startnode - where the path starts
%   Output: path - the path in A that starts with startnode
%   Joanne (Junyuan) Lin @Dept. of Math, Tufts University


% initialization
path = startnode;
ngbr = find(A(startnode,:));
% main loop
while ngbr ~= 0
    path(end+1) = ngbr;
    edge_endpoints = find(A(ngbr,:));
    test = sum(edge_endpoints);
    temp = ngbr;
    ngbr = sum(edge_endpoints) - startnode;
    startnode = temp;
end




end

