function [ aggregation, num_agg]=coarsening_agg_MWM(L, flag)
% Coarsening based on maximal weighted matching

n = length(L(:,1));
%get rid of the diagonals
d=diag(L);
L = L-spdiags(d,0,n,n);

num_agg = 0;
%iso_edges = {};

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

%for the isolated points
% put diagonal back
L = L+spdiags(d,0,n,n);

% find leftovers
left = find(aggregation == 0);

switch flag
    %aggregate the left-over points with the heaviest edge
    case 1,
        for i = 1:length(left)
            %[~,j] = max(-L(left(i),:));  % need to deal with positive off-diagonal entries 
            
            if aggregation(left(i))==0
                
                 [~,j] = min(L(:,left(i))); % only works for symmetric matrix!!!
                 
                  if aggregation(j)==0
            
                      num_agg = num_agg+1;
                      aggregation(j) = num_agg;
                      aggregation(left(i)) = num_agg;
                
                  else
                      
                      aggregation(left(i)) = aggregation(j);
                      
                  end
            
            end 
           
        end
        
     otherwise,
         
         aggregation(left) = (num_agg+1:num_agg+length(left))';
         num_agg = num_agg+length(left);
         
end      
             
%for the isolated points
% while ~isempty(find(aggregation == 0, 1))
%     
%     left = find(aggregation == 0);
%     
%     switch flag
%         %aggregate the left-over points with the heaviest edge
%         case 1,
%             for i = 1:length(left)
%                 %[~,j] = max(-L(left(i),:));  % need to deal with positive off-diagonal entries 
%                  [~,j] = min(L(:,left(i))); % only works for symmetric matrix!!!
%                 aggregation(left(i)) = aggregation(j);
%                % iso_edges{end+1} = [left(i) j];
%             end
%         otherwise,
%             aggregation(left) = (num_agg+1:num_agg+length(left))';
%             num_agg = num_agg+length(left);
%     end
%     
% end

end



     
        
