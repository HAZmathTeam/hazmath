function [ aggregation,num_agg,iso_edges ] = pair_match_XH( cover,AG_ori )
%Greedy Algorithm for Maximal Matching
% Input: cover(cell)-each cell has a path
% Output: a maximal matching M

%initialization
m = size(AG_ori,1);
n = length(cover);
aggregation = zeros(m,1);
agg_size = [];
num_agg = 0;
iso_edges = {};
d=diag(AG_ori);
AG_ori = AG_ori-spdiags(d,0,m,m);

for i = 1:n
    cover_len = length(cover{i});
    
    for j = 1:floor(cover_len/2)
        num_agg = num_agg+1;
        aggregation(cover{i}(2*j-1))=num_agg;
        aggregation(cover{i}(2*j))=num_agg;
        agg_size(num_agg,1) = 2; 
    end
    
%     for j = 1:floor(cover_len/2)-1  % start from the second one
%         num_agg = num_agg+1;
%         aggregation(cover{i}(2*j+1))=num_agg;
%         aggregation(cover{i}(2*j))=num_agg;
%         
%     end
    
end
        
isolated = find(aggregation==0);

for i = 1:length(isolated)
% adding the isolated points based on the heaviest edge
    if (aggregation(isolated(i)) == 0) % make sure isolated(i) is actually isolated -- added by Xiaozhe
        [~,I] = min(AG_ori(:,isolated(i)));
        if aggregation(I)==0 % the heavist neighbor is also isolated
            num_agg = num_agg+1;
            aggregation(isolated(i)) = num_agg;
            aggregation(I) = num_agg;
            iso_edges{end+1} = [isolated(i) I];
            agg_size(num_agg,1) = 2;
        elseif agg_size(aggregation(I),1) >= 3
        %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
        % Add this to avoid large diameter (large aggregation size)
            agg_size(aggregation(I),1) = agg_size(aggregation(I),1) - 1;
            num_agg = num_agg+1;
            aggregation(I) = num_agg;
            aggregation(isolated(i)) = num_agg;
            agg_size(num_agg,1) = 2;
        else % the heavist neighbor is in an aggregation of size 2
            aggregation(isolated(i))=aggregation(I);
            agg_size(aggregation(I),1) = agg_size(aggregation(I),1) + 1;
            iso_edges{end+1} = [isolated(i) I];
        end
    end
%========================================================================%
% aggregate isolated point based on the heaviest edge that connects to the aggregated neighbor
%     temp_arr = full(AG_ori(isolated(i),:));
%     [~,I] = sort(abs(temp_arr),'descend');
%     nonzeros_ind = find(aggregation(I)); % in case the heaviest edge is agg 0
%     aggregation(isolated(i))=aggregation(I(nonzeros_ind(1)));
end

%========================================================================%
% % just leave the aggregations alone
% num_agg_new = num_agg+length(isolated);
% aggregation(isolated) = num_agg+1: num_agg_new;
% num_agg = num_agg_new;



end

