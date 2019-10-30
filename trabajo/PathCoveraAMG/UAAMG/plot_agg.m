function [ p ] = plot_agg( G_ori, aggregation,num_agg)
%plot of G that different vertices are in different colors
%p = plot(G_ori,'EdgeLabel',G_ori.Edges.Weight);
%p = plot(G_ori);
n = sqrt(numnodes(G_ori));
x = 1:n;
x = repmat(x,n,1);
y = (1:n)';
y = repmat(y, 1,n);
p = plot(G_ori, 'XData', x(:), 'YData', y(:));
G_ori.Nodes.NodeColors =aggregation;
p.NodeCData = G_ori.Nodes.NodeColors;


for i = 1:num_agg
    ind = find(aggregation==i);
    
    %if length(ind) == 2
        highlight(p,ind,'EdgeColor','r','LineWidth',5)
%     else
%         ditch = intersect(ind,pt);
%         remain = setdiff(ind,ditch);
%         highlight(p,remain,'EdgeColor','r','LineWidth',5)
    %end
    
end


end