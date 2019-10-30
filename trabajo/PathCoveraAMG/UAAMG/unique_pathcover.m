function [ cover,  num_cover ] = unique_pathcover( cover )
% compress the paths 
%
% @ Xiaozhe Hu, Tufts University

% number of paths (might be empty)
total_cover = length(cover);

% local variable
num_keep = 0;      % number of paths that will be kept (nonempty paths)
num_delete = 0;    % number of paths that will be delete (empty paths)      

while (num_keep + num_delete < total_cover)
    
    if ( isempty(cover{num_keep+1}) )
        cover(num_keep+1,:) = [];
        num_delete = num_delete + 1;
    else
        num_keep = num_keep + 1;
    end
    
end

num_cover = num_keep;

end

