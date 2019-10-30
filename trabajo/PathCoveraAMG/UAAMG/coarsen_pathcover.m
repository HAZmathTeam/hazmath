function [ coarse_cover ] = coarsen_pathcover( cover, aggregation )
% coarsen the path covers based on the aggregations
%
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)

% initialize coarse_cover
num_cover = length(cover);
coarse_cover = cell(num_cover,1);

% main loop
for i=1:num_cover
    
    %!!!!!!!!!!!!!!!%
    % Assume the the matching on the each path starts from the first node
    % does not work if the matching starts from second node
    %!!!!!!!!!!!!!!!%
    n_cover = length(cover{i});
   
    if (mod(n_cover, 2)==0) % the length is even
        temp = aggregation(cover{i});
        coarse_cover{i} = temp(1:2:end);
    else 
        if (n_cover == 1) % path already is just one node
            coarse_cover{i} = aggregation(cover{i});
        else % the length is odd >1
            temp = aggregation(cover{i});
            coarse_cover{i} = temp(1:2:end-1);
        %[~,jj] = unique(temp);
        %jj = sort(jj);
        %coarse_cover{i} = temp(jj);
        end
    end
    
    %temp = aggregation(cover{i});
    %[~,jj] = unique(temp);
    %jj = sort(jj);
    %coarse_cover{i} = temp(jj);
    %coarse_cover{i} = temp(jj);
    
end


end

