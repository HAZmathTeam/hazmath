function [ path ] = generate_BFS_path( Adj, level_set, smooth_error)
% generate path cover from BFS level set 
%
% Note: Adj should be the adjancy matrix within level set only
%
% @ Xiaozhe Hu

% get size
num_level_set = max(level_set);

% initialize path
num_path = 0;

% main loop
for i = 1:num_level_set
    
    % nodes on this level
    node_level = find(level_set == i);
    
    % adjancy matrix on this level
    A_level = Adj(node_level, node_level);
    
    % find connected component
    [C_level, num_C_level, ~] = find_CC( A_level );
    
    % put into paths
    for j = 1:num_C_level
        
        % increase paths
        num_path = num_path+1;
        
        % get path
        path{num_path,1} = node_level(C_level{j});
        
    end
    
end


end

