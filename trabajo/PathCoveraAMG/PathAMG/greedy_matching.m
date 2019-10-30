function [ aggregation, num_agg, iso_edges ]=greedy_matching( L )
% Coarsening based on heavy edge coarsening

n = length(L(:,1));
%get rid of the diagonals
d=diag(L);
L = L-spdiags(d,0,n,n);


num_agg = 0;
iso_edges = {};

aggregation=zeros(n,1);

[row,col,val] = find(L);
[edges,ind] = sort(val,'ascend');

% main loop
for i = 1:length(edges)
    %check if the two endpoints are aggregated
    if aggregation(row(ind(i)))==0 && aggregation(col(ind(i)))==0
        num_agg = num_agg+1;
        aggregation(row(ind(i))) = num_agg;
        aggregation(col(ind(i))) = num_agg;
    end
end



% while ~isempty(find(L))
%     
%     %find the heaviest edge in L
%     [m,row_ind] = min(L);
%     [~,col_ind] = min(m);
%     
%     num_agg = num_agg+1;
%     aggregation(row_ind(col_ind)) = num_agg;
%     aggregation(col_ind) = num_agg;
% 
%     %get rid of the used nodes since they can't be paired up with
%     %other nodes
%     L(row_ind(col_ind),:) = 0;
%     L(:,row_ind(col_ind)) = 0;
%     L(col_ind,:) = 0;
%     L(:,col_ind) = 0;
%      
%     
% end

%for the isolated points
while ~isempty(find(aggregation == 0))
    left = find(aggregation == 0);
    %aggregate the left-over points with the heaviest edge
    for i = 1:length(left)
        [~,j] = max(-L(left(i),:));
        aggregation(left(i)) = aggregation(j);
        iso_edges{end+1} = [left(i) j];
    end
end

end



     
        
