function [ p ] = agg_plot( G_ori,aggregation,num_agg, iso_edges)
%plot of G that different vertices are in different colors
%p = plot(G_ori,'EdgeLabel',G_ori.Edges.Weight);
%p = plot(G_ori);
n = sqrt(numnodes(G_ori));
x = 0:1/(n-1):1;
x = repmat(x,n,1);
y = (0:1/(n-1):1)';
y = repmat(y, 1,n);
p = plot(G_ori, 'XData', x(:), 'YData', y(:));
G_ori.Nodes.NodeColors =aggregation;
p.NodeCData = G_ori.Nodes.NodeColors;

% get the isolated points and plot the isolated edges
n = length(iso_edges);
pt = [];
% for k = 1:n
%     pt = [pt; iso_edges{k}(1)]; %since all the isolated points are stored in first element
%     highlight(p,iso_edges{k},'EdgeColor','r','LineWidth',5)
% end

for i = 1:num_agg
    ind = find(aggregation==i);
    
    if length(ind) == 2
        highlight(p,ind,'EdgeColor','r','LineWidth',5)
    else
        ditch = intersect(ind,pt);
        remain = setdiff(ind,ditch);
        highlight(p,remain,'EdgeColor','r','LineWidth',5)
    end
    
end


end

