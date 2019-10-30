function [ aggregation, num_agg, iso_edges ]=matching( L )
% Coarsening based on heavy edge coarsening

n = length(L(:,1));
%get rid of the diagonals
d=diag(L);
L = L-spdiags(d,0,n,n);
L_ori = L;

num_agg = 0;
iso_edges = {};

%p=symrcm(L);
p = amd(L);

aggregation=zeros(n,1);
used_nodes = [];

% main loop
for i=1:n
    
    if aggregation(p(i))==0
    
        [~,m]=min(L(:,p(i)));
        
        if aggregation(m)==0
        
            num_agg = num_agg+1;
            aggregation(m) = num_agg;
            aggregation(p(i)) = num_agg;
            
            %get rid of the used nodes since they can't be paired up with
            %other nodes
            L(p(i),:) = 0;
            L(:,p(i)) = 0;
            L(m,:) = 0;
            L(:,m) = 0;
        

         end
        
    end
    
end

%for the isolated points
while ~isempty(find(aggregation == 0))
    left = find(aggregation == 0);
    %aggregate the left-over points with the heaviest edge
    for i = 1:length(left)
        [~,j] = max(-L_ori(left(i),:));
        aggregation(left(i)) = aggregation(j);
        iso_edges{end+1} = [left(i) j];
    end
end


     
        
