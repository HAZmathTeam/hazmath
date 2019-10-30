function [ aggregation, num_agg, coarse_path ] = path_matching( path, Adj )
% matching on the paths 
%
% @ Xiaozhe Hu

% get number of path
num_path = length(path);

% get size
N = size(Adj,1);

% initialize
aggregation = zeros(N,1);
num_agg = 0;

coarse_path = cell(num_path,1);

% main loop
for i=1:num_path
    
    % nodes on path
    node_path = path{i};
    
    % loca adjancy matrix
    A_path = Adj(node_path, node_path);
    
    % matching
    [ aggregation_path,num_agg_path ] = coarsening_MWM(A_path);
    
    % put local matching to global
    aggregation(node_path) = num_agg + aggregation_path;
    coarse_path{i} = num_agg+1:num_agg + num_agg_path;
    num_agg = num_agg + num_agg_path;
             
end

isolated = find(aggregation==0);
num_iso = length(isolated);

for i=1:num_iso
    
    if (aggregation(isolated(i)) == 0)
        
        [~,j] = min(Adj(:,isolated(i)));
        
        if aggregation(j)==0
            num_agg = num_agg + 1;
            aggregation(isolated(i)) = num_agg;
            aggregation(I) = num_agg;
        else
            aggregation(isolated(i))=aggregation(j);  
        end
            
    end
    
    
end

end

