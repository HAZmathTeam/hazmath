function [ A ] = compute_similarity( coord, A, sparsity_level )
% similarity
%
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)

[n, ~] = size(coord);

S0 = spones(A);
S = S0;
for i=2:sparsity_level
    S = S*S0;
end

% find nonzeros of A^2
[row_idx, col_idx, ~] = find(triu(S,1));

% construct A
%val = -1./(max(abs(coord(row_idx,:) - coord(col_idx,:)),[],2) + 1e-6 );
val = -1./(max(abs(coord(row_idx,:) - coord(col_idx,:)),[],2) );
A = sparse(row_idx, col_idx, val, n, n);
A = A + A';

end

