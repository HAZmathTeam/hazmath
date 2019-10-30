function [ cover, W_cover, p, o, AG_new ] = genCover( AG )
% The algorithm is intended for covering undirected graphs. 
% From Moran, Newman and Wolfstahl's paper of Approximation Algorithms for Covering a
% Graph by Vertex-Disjoint Paths of
% Maximum Total Weight
% Input: AG, graph laplacian of undirected graph G = (V, E) and a weight function WG:E->Z+
% Output: cover, a cover of G which has cell stucture and each cell contains a path.
%         W_cover, the weight of cover S on graph G 

% Joanne Lin @Tufts University

% Created 6/12/2017

%%%%%%%%%%%%%%%%%%%% Video Effect Module %%%%%%%%%%%%%%%%%%%%
%translate graph laplacian AG to graph G
% G_ori = graph(AG,'omitself');
% pl = plot(G_ori, 'EdgeLabel', -G_ori.Edges.Weight);
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialization

n = size(AG, 1);
D = diag(AG);
AG = abs(AG-spdiags(D,0,n,n)); %take diagonals out and use positive weights
V_new = [];
p = 1:n; % contains nodes where the paths start with
o = 1:n; % contains nodes where the paths end with
% for example path 1->2->3, p = [1,1,1], o = [3,3,3]

AG_new = zeros(n,n);


while ~isempty(find(AG))
    % update the new graph
        for i = 1:n
            if p(i)~=i && o(i)~=i
                AG(i,:) = 0; %delete redundent edges
                AG(:,i) = 0; %delete redundent edges
            end
        end
        
    % take the heaviest edge
    [w_col, id_row] = max(triu(AG));
    [w, v] = max(w_col);
    u = id_row(v);
    
    % get rid of the heaviest edge
    AG(u,v) = 0;
    AG(v,u) = 0;
    
    % main loop
    if p(u) ~= p(v)
        %highlight(pl,[u v],'EdgeColor','r','LineWidth',5);
        %pause(1.5)
        if isempty(intersect(u, V_new)) % u is the newly added node, e(u,v) is added either to the head or to the end
            if v == p(v) % if v is the old starting point
                old_start = p(v);
                for x = 1:n
                    if p(x) == old_start % take all the edges/paths incident e(u,v) on node v
                        p(x) = p(u); % change the edges' beginning nodes to u
                    end
                end
                o(u) = o(v);
            else % if v is the old ending point
                old_end = o(v);
                for x = 1:n
                    if o(x) == old_end % take all the edges/paths incident e(u,v) on node v
                        o(x) = o(u); % change the edges' beginning nodes to u
                    end
                end
                p(u) = p(v);
            end
%             AG(v,:) = 0;
%             AG(:,v) = 0;
        elseif isempty(intersect(v, V_new))% if v is the newly added node
            if u == p(u) % if v is the old starting point
                old_start = p(u);
                for x = 1:n
                    if p(x) == old_start; % take all the edges/paths incident e(u,v) on node v
                        p(x) = p(v); % change the edges' beginning nodes to u
                    end
                end
                o(v) = o(u);
            else
                old_end = o(u);
                for x = 1:n
                    if o(x) == old_end % take all the edges/paths incident e(u,v) on node v
                        o(x) = o(v); % change the edges' beginning nodes to u
                    end
                end
                p(v) = p(u);
            end
%             AG(u,:) = 0;
%             AG(:,u) = 0;
        else % both u and v are in the set
            if u == p(u)
                if v == p(v)
                    if o(u) > o(v)
                        old_end_sm = o(v);
                        old_end_lg = o(u);
                        for x = 1:n
                            if o(x) == old_end_sm || o(x) == old_end_lg% take all the edges/paths incident e(u,v) on node v
                                o(x) = old_end_lg; % change the edges' beginning nodes to u
                                p(x) = old_end_sm;
                            end
                        end
                    else
                        old_end_sm = o(u);
                        old_end_lg = o(v);
                        for x = 1:n
                            if o(x) == old_end_sm || o(x) == old_end_lg% take all the edges/paths incident e(u,v) on node v
                                o(x) = old_end_lg; % change the edges' beginning nodes to u
                                p(x) = old_end_sm;
                            end
                        end
                    end
                elseif v == o(v)
                    if o(u) > p(v)
                        old_end_sm = p(v);
                        old_end_lg = o(u);
                        for x = 1:n
                            if p(x) == old_end_sm || o(x) == old_end_lg% take all the edges/paths incident e(u,v) on node v
                                o(x) = old_end_lg; % change the edges' beginning nodes to u
                                p(x) = old_end_sm;
                            end
                        end
                    else
                        old_end_sm = o(u);
                        old_end_lg = p(v);
                        for x = 1:n
                            if o(x) == old_end_sm || p(x) == old_end_lg% take all the edges/paths incident e(u,v) on node v
                                o(x) = old_end_lg; % change the edges' beginning nodes to u
                                p(x) = old_end_sm;
                            end
                        end
                    end
                end
            else
                if v == p(v)
                    if p(u) > o(v)
                        old_end_sm = o(v);
                        old_end_lg = p(u);
                        for x = 1:n
                            if o(x) == old_end_sm || p(x) == old_end_lg% take all the edges/paths incident e(u,v) on node v
                                o(x) = old_end_lg; % change the edges' beginning nodes to u
                                p(x) = old_end_sm;
                            end
                        end
                    else
                        old_end_sm = p(u);
                        old_end_lg = o(v);
                        for x = 1:n
                            if p(x) == old_end_sm || o(x) == old_end_lg% take all the edges/paths incident e(u,v) on node v
                                o(x) = old_end_lg; % change the edges' beginning nodes to u
                                p(x) = old_end_sm;
                            end
                        end
                    end
                elseif v == o(v)
                    if p(u) > p(v)
                        old_end_sm = p(v);
                        old_end_lg = p(u);
                        for x = 1:n
                            if p(x) == old_end_sm || p(x) == old_end_lg% take all the edges/paths incident e(u,v) on node v
                                o(x) = old_end_lg; % change the edges' beginning nodes to u
                                p(x) = old_end_sm;
                            end
                        end
                    else
                        old_end_sm = p(u);
                        old_end_lg = p(v);
                        for x = 1:n
                            if p(x) == old_end_sm || p(x) == old_end_lg% take all the edges/paths incident e(u,v) on node v
                                o(x) = old_end_lg; % change the edges' beginning nodes to u
                                p(x) = old_end_sm;
                            end
                        end
                    end
                end
            end
                
        end
        
        
        AG_new(u,v) = w;
        AG_new(v,u) = w;
        V_new = union(V_new, u);
        V_new = union(V_new, v);
    end
end

W_cover = sum(sum(AG_new))/2;
G = graph(AG_new,'omitself');
starting_points = union(p,p);
cover = cell(length(starting_points),1);
for i = 1:length(starting_points)
    cover{i} = dfsearch(G,starting_points(i));
end


end

