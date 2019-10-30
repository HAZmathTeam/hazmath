function [ path ] = level_set_to_path( level_set )

% get number of paths
num_path = max(level_set);

% initilize path
path = cell(num_path,1);

% main loop
for i = 1:num_path
    
    path{i} = find(level_set == i);
    
end

end

