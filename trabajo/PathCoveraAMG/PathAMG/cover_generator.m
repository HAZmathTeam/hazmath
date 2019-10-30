function [ cover,isolated,L_cover,AG_ori ] = cover_generator( AG )
% The algorithm is intended for covering undirected graphs. 
% Input: A graph laplacian AG of undirected graph G = (V, E) and a weight function WG:E->Z+
% Output: cover, a cover of G.

% Joanne Lin @Tufts University

% %translate graph laplacian AG to graph G
% G_ori = graph(AG,'omitself');
% G_ori.Edges.Weight = -G_ori.Edges.Weight;
% 
% %construct the laplacian of the original graph
% B_ori = incidence(G_ori);
% L_ori = B_ori*diag(G_ori.Edges.Weight)*B_ori';
% n = numnodes(G_ori);


%get rid of the diagonals
n = size(AG,1);
d=diag(AG);
AG_ori = AG-spdiags(d,0,n,n);
AG = AG_ori;
L_cover = zeros(size(AG_ori));

%place holder for the path cover
cover = {};

%set-up
set = [];

while ~isempty(find(AG))
    %find a path of G
    [ path ] = pathfinder( AG );

    %construct the laplacian of the graph of cover
    for i = 1:(length(path)-1)
        L_cover(path(i),path(i+1)) = AG_ori(path(i),path(i+1));
        L_cover(path(i+1),path(i)) = AG_ori(path(i+1),path(i));
    end

    cover{end+1} = path; 

    set = [set;path];

    %get rid of the used nodes on the path
    del_nodes = path;
    for i = 1:length(del_nodes)
        AG(del_nodes(i),:) = zeros(1,n);
        AG(:,del_nodes(i)) = zeros(n,1);
    end
%     % regenerate the graph
%     G = graph(AG,'omitself');
%     G.Edges.Weight = -G.Edges.Weight;

end

isolated = setdiff(1:n,set);


% G = graph(L_cover);
% G.Edges.Weight = -G.Edges.Weight;
    


end

