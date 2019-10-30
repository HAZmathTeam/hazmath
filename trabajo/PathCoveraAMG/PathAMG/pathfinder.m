function [ path,path_fwd,path_bkwd ] = pathfinder( AG )
% The algorithm is intended for covering undirected graphs based on simply
% growing the longest paths favoring heaviest edges (see test for video demo)
% Input: An undirected graph G = (V, E) and a weight function WG:E->Z+
% Output: path, a path of G containing u and v.

% Joanne Lin @Tufts University

%%set-up
n = size(AG,1);
d=diag(AG);
L = AG-spdiags(d,0,n,n);

%%first round path-finding
[M,row_ind] = max(-(L)); %start with the heaviest edge 
[w,col_ind] = max(M);
edge = [row_ind(col_ind) col_ind];

%favors less negative edge weights when dealing with negative weights
% if isempty(find(M))
%     [row,col,val] = find(L);
%     [~,I] = min(val);
%     edge = [row(I) col(I)];
% end

u = edge(1);
v = edge(2);

[ path_fwd ] = path_oneside( AG,u,v );

%%preparation
% get rid of the used nodes
del_nodes = setdiff(path_fwd,v);
for i = 1:length(del_nodes)
    L(del_nodes(i),:) = zeros(1,n);
    L(:,del_nodes(i)) = zeros(n,1);
end
% add back e(u,v)
L(u,v) = -w;
L(v,u) = -w;
% % regenerate the graph
% G = graph(L,'omitself');
% G.Edges.Weight = -G.Edges.Weight;

%%finding the path from the other side
[ path_bkwd ] = path_oneside( L,v,u );

%%join the two path and output a cover
add = flipud(path_bkwd);
add = add(3:end);
path = union(add,path_fwd,'stable');

    


end

