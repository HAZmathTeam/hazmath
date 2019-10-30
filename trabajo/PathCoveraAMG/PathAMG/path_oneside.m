function [ path ] = path_oneside( AG,u,v )
% The algorithm is intended for covering undirected graphs. From Moran,
% Newman and Wolfstahl's paper of Approximation Algorithms for Covering a
% Graph by Vertex-Disjoint Paths of
% Maximum Total Weight
% Input: laplacian AG of an undirected graph G = ( V, E) and a weight function WG:E->Z+
% Output: A path from one side.

% Joanne Lin @Tufts University

% initialization
E_new = [];
n = size(AG,1);
d=diag(AG);
AG = AG-spdiags(d,0,n,n);

p = 1:n;
o = 1:n;
% If v in V is an endpoint of a path in the cover, then p ( v ) is the
% path covering v, and o(v) is the other endpoint of this path.
% Edges = G.Edges.EndNodes;
% Weight = G.Edges.Weight;
% [~,I] = ismember([u v],Edges,'rows');
% if I==0
%     [~,I] = ismember([v u],Edges,'rows');
% end
% edge = Edges(I,:);
% Edges(I,:) = zeros(1,2);
% Weight(I) = 0;

AG(u,v) = 0;
AG(v,u) = 0;

temp_w = 1;
while norm(temp_w)~=0
    % Choose an edge e = ( u , v ) in E, such that W(e) is maximum.
    if p(u)~=p(v)
%         E_new = [E_new;edge];
%         Edges(I,:) = zeros(1,2);
%         Weight(I) = 0;

        AG(u,v) = 0;
        AG(v,u) = 0;
%         if o(u)~=u
%             for i = 1:size(Edges,1)
%                 if Edges(i,1) == u 
%                     Edges(i,:) = zeros(1,2);
%                     Weight(i) = 0;
%                 elseif Edges(i,2) == u 
%                     Edges(i,:) = zeros(1,2);
%                     Weight(i) = 0;
%                 end
%             end
%         end
% since o(u)=u always
% This algorithm only enlongates the path from vertex v's end. Maybe change
% the algorithm to enlongates from both ends later.
        p(o(u))=p(v);
        t = o(u);
        o(v) = o(o(u));
        o(o(v))=t;
        if o(v)~=v
%            for i = 1:size(Edges,1)
%                 if Edges(i,1) == v
%                     Edges(i,:) = zeros(1,2);
%                     Weight(i) = 0;
%                 elseif Edges(i,2) == v
%                     Edges(i,:) = zeros(1,2);
%                     Weight(i) = 0;
%                 end

                AG(:,v) = 0;
                AG(v,:) = 0;
%            end
        end

%       temp_w = zeros(length(Weight),1);
%        for i = 1:size(Edges,1)
%             if Edges(i,1) == u 
%                 temp_w(i) = Weight(i);
%             elseif Edges(i,2) == u 
%                 temp_w(i) = Weight(i);
%             end
%        end
       [~,y,temp_w] = find(AG(u,:));
       if norm(temp_w)==0
           break;
       end
    end
    
    [~,I] = max(-(temp_w));
    
    edge = [u y(I)];
    v = u;
    u = setdiff(edge,v);
    
end
        
% for finding a path
isolated_pts = [];
for i = 1:length(p)
    if p(i)==o(i)
        isolated_pts = [isolated_pts;i];
    elseif p(i)==i
        start_pt = i;
    elseif o(i)==i
        end_pt = i;
    end
end
path = start_pt;
node = start_pt;
while node~=end_pt
    path = [path; o(node)];
    node = o(node);
end
end



