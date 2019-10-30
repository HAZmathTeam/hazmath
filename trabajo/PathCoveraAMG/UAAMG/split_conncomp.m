function [ L_cell,b_cell ] = split_conncomp( filename )
%preprocess the data
%   take L and b and split into cells containing connected components
% Junyuan Lin @Tufts Math Dept.

load(filename);

G = graph(L,'omitself');
C = conncomp(G);
count_C = union(C,C);
L_cell = cell(length(count_C),1);
b_cell = cell(length(count_C),1);

for i = 1:length(count_C)
    ind = find(C==count_C(i));
    L_cell{i} = L(ind,ind);
    b_cell{i} = b(ind);
end
    
end

