function [ cover_flag ] = convert_cover_flag( cover, N )
% Conver list of path covers to a flag
%
% @ Xiaozhe Hu, Tufts University

cover_flag = zeros(N,1);
num_cover = length(cover);

for i = 1:num_cover
    
    cover_flag(cover{i}) = i;
    
end


end

