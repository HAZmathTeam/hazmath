function [ cover, numUniquePath ] = genCover1( AG )
% The algorithm is intended for covering undirected graphs. 
% From Moran, Newman and Wolfstahl's paper of Approximation Algorithms for Covering a
% Graph by Vertex-Disjoint Paths of
% Maximum Total Weight
% Input: AG, graph laplacian of undirected graph G = (V, E) and a weight function WG:E->Z+
% Output: cover, a cover of G which has cell stucture and each cell contains a path.
%         W_cover, the weight of cover S on graph G 
%         AG_new, the graph laplacian of the cover

% Joanne Lin @Tufts University

% Created 6/25/2017

%%%%%%%%%%%%%%%%%%%% Video Effect Module %%%%%%%%%%%%%%%%%%%%
%translate graph laplacian AG to graph G
figure
G_ori = graph(AG,'omitself');
n = sqrt(numnodes(G_ori));
x = 0:1/(n-1):1;
x = repmat(x,n,1);
y = (0:1/(n-1):1)';
y = repmat(y, 1,n);
pl = plot(G_ori, 'XData', x(:), 'YData', y(:));
%pl = plot(G_ori);
%pl = plot(G_ori, 'EdgeLabel', -G_ori.Edges.Weight);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialization
%m = length(find(AG))/2;
n = size(AG, 1);
D = diag(AG);
AG = abs(AG-spdiags(D,0,n,n)); %take diagonals out and use positive weights

%AG_new = zeros(n,n);
vertex_ngbr = cell(n,1);
path_flag = zeros(n,1); % vertices on the same path have same value. 0 means 
                       % it hasn't been visited
path_ref = cell(nnz(AG)/2,2); % for ith row, first cell is the all the vertices in the path,
                      % the second cell is the length of the current path length
path_ref{1,1} = zeros(1,n);
path_ref{1,2} = 0;
numPath = 0;

%sort the edges in descending order
AG = triu(AG); %if the graph is undirected
[row,col,val] = find(AG); 
[edges,ind] = sort(val,'descend');
    
for i = 1:length(edges)
    u = row(ind(i));
    v = col(ind(i));
    
    if length(vertex_ngbr{u}) < 2 && length(vertex_ngbr{v}) < 2
    %if length(find(AG_new(u,:))) < 2 && length(find(AG_new(v,:))) < 2 %if u and v are both the endpoints
        %if ~(path_ref(path_flag(u),3) == path_ref(path_flag(v),3) && sum(path_ref(path_flag(u),1:2))~=2) % when u and v are not in the same nonzero path
        if ~(path_flag(u) == path_flag(v) && path_flag(u) ~= 0) % when u and v are not in the same nonzero path
%             AG_new(u,v) = val(ind(i));
%             AG_new(v,u) = val(ind(i));
%             AG_new(u,v) = 1;
%             AG_new(v,u) = 1;
            vertex_ngbr{u} = [vertex_ngbr{u} v];
            vertex_ngbr{v} = [vertex_ngbr{v} u];

            %%%%%%%%%%%%%%%%%%%% Video Effect Module %%%%%%%%%%%%%%%%%%%%
            highlight(pl,[u v],'EdgeColor','r','LineWidth',5);
            %pause
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if path_flag(u) == 0 && path_flag(v) ~= 0 % v is already in a path, but u is just visited
                path_flag(u) = path_flag(v);
                %path_ref{path_flag(v),1} = [path_ref{path_flag(v),1} u];
                path_ref{path_flag(v),2} = path_ref{path_flag(v),2} + 1;
                path_ref{path_flag(v),1}(path_ref{path_flag(v),2}) = u;
            elseif path_flag(v) == 0 && path_flag(u) ~= 0 % u is already in a path, but v is just visited
                path_flag(v) = path_flag(u);
                %path_ref{path_flag(u),1} = [path_ref{path_flag(u),1} v];
                path_ref{path_flag(u),2} = path_ref{path_flag(u),2} + 1;
                path_ref{path_flag(u),1}(path_ref{path_flag(u),2}) = v;
            elseif path_flag(v) == 0 && path_flag(u) == 0 % both u and v are not in the path
                numPath = numPath + 1;
                path_flag(u) = numPath;
                path_flag(v) = numPath;
                %path_ref{numPath,1} = [u v];
                %path_ref{numPath,1} = sparse(1,n);
                path_ref{numPath,1}(1) = u;
                path_ref{numPath,1}(2) = v;
                path_ref{numPath,2} = 2;
            else % u and v are in different nonzero paths
%                 path_endpoints(vertex_flag(u),:) = [sum(path_endpoints(vertex_flag(u),:))-u, sum(path_endpoints(vertex_flag(v),:))-v];
%                 path_endpoints(vertex_flag(v),:) = path_endpoints(vertex_flag(u),:);
%                 [path_id, row_id] = min([path_ref(path_flag(u),3), path_ref(path_flag(v),3)]);
%                 if row_id == 1    
%                     path_ref(path_flag(v),:) = [sum(path_ref(path_flag(v),1:2))-v, u, path_id];
%                 else
%                     path_ref(path_flag(u),:) = [sum(path_ref(path_flag(u),1:2))-u, v, path_id];
%                 end
                if path_ref{path_flag(u),2} > path_ref{path_flag(v),2}
                    %path_ref{path_flag(u),1} = [path_ref{path_flag(u),1} path_ref{path_flag(v),1}];
                    %path_ref{path_flag(u),2} = length(path_ref{path_flag(u),1});
                    path_ref{path_flag(u),1}((path_ref{path_flag(u),2} + 1) :(path_ref{path_flag(u),2} + path_ref{path_flag(v),2})) = path_ref{path_flag(v),1}(1:path_ref{path_flag(v),2});
                    path_ref{path_flag(u),2} = path_ref{path_flag(u),2} + path_ref{path_flag(v),2};
                    %path_flag(path_ref{path_flag(v),1}) = path_flag(u); % change all the vertices' path_flag previously in path containds v
                    path_flag(path_ref{path_flag(v),1}(1:path_ref{path_flag(v),2})) = path_flag(u);
                else
                    %path_ref{path_flag(v),1} = [path_ref{path_flag(u),1} path_ref{path_flag(v),1}];
                    %path_ref{path_flag(v),2} = length(path_ref{path_flag(v),1});
                    path_ref{path_flag(v),1}((path_ref{path_flag(v),2} + 1) :(path_ref{path_flag(v),2} + path_ref{path_flag(u),2})) = path_ref{path_flag(u),1}(1:path_ref{path_flag(u),2});
                    path_ref{path_flag(v),2} = path_ref{path_flag(v),2} + path_ref{path_flag(u),2};
                    path_flag(path_ref{path_flag(u),1}(1:path_ref{path_flag(u),2})) = path_flag(v); % change all the vertices' path_flag previously in path containds u
                end
                %numUniquePath = numUniquePath - 1;
            end
%             elseif vertex_flag(u) ~= vertex_flag(v) % u and v are in different nonzero paths
%                 ind2 = find(vertex_flag == vertex_flag(v));
%                 vertex_flag(ind2) = vertex_flag(u);
%             end
        end
    end
end


pause
%W_cover = sum(sum(AG_new))/2;
isolated_points = find(path_flag == 0);
unique_pathid = unique(path_flag);
if ~isempty(isolated_points)
    unique_pathid(1) = []; % delete 0 if there are any non-visited nodes
end
numUniquePath = length(unique_pathid);
all_numPath = numUniquePath + length(isolated_points);

%cover = cell(all_numPath,1);
cover = cell(numUniquePath,1);

for i = 1:numUniquePath
%     ind = find(vertex_flag==path_id(i));
%     found_start = 0;
%     for j = 1:length(ind)
%         %numNGBR = find(AG_new(ind(j),:));
%         numNGBR = vertex_ngbr{ind(j)};
%         if length(numNGBR) == 1 && found_start == 0
%             startnode = ind(j);
%             found_start = 1;
%         end
%     end
    %vertex_list = path_ref{unique_pathid(i),1};
    vertex_list = path_ref{unique_pathid(i),1}(1:path_ref{unique_pathid(i),2});
    flag = 0;
    for j = 1:path_ref{unique_pathid(i),2}
        if flag==0 && length(vertex_ngbr{vertex_list(j)})==1
            start_pt = vertex_list(j);
            flag=1;
        end
    end
    
    cover{i} = dfs_path( vertex_ngbr, start_pt,path_ref{unique_pathid(i),2});
end
        
        
% old code
% G_new = graph(AG_new);
% deg = degree(G_new);
% ind = find(deg==1);
% bin = conncomp(G_new);
% numPath = length(union(bin,bin));
% % path_start = union(vertex_flag,vertex_flag);
% % numPath = length(path_start);
% cover = cell(numPath,1);
% for i = 1:numPath
%     starting_point = find(bin(ind)==i);
%     if length(starting_point)==2
%         cover{i} = dfsearch(G_new,ind(starting_point(1)));
%         %cover{i} = find(vertex_flag==path_start(i));
%     end
% end

end

