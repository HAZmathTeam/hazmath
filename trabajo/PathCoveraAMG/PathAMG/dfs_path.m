function [ path ] = dfs_path( vertex_ngbr, startnode,path_len )
% find path given adjacency matrix using dfs
%   Inputs: vertex_ngbr - an n by 1 cell with each vertex's
%           neighbors in the cell
%           startnode - where the path starts
%           path_len - length of path
%   Output: path - the path in A that starts with startnode
%   Joanne (Junyuan) Lin @Dept. of Math, Tufts University


% initialization
path = zeros(path_len,1);
path(1) = startnode;
%path = startnode;
ind = 2;
ngbr = vertex_ngbr{startnode};
% main loop
while ngbr ~= 0
    path(ind) = ngbr;
    %path(end + 1) = ngbr;
    edge_endpoints = vertex_ngbr{ngbr};
    temp = ngbr;
    ngbr = sum(edge_endpoints) - startnode;
    startnode = temp;
    ind = ind + 1;
end




end

