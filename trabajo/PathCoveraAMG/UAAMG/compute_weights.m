function [ A ] = compute_weights(smooth_error, A, sparsity_level )
% commpute the weights based on XZ identity:  
% \| e - scaling P R e \|_{D}/ \| e \|_A 
%
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)

% get size
[n, n_smooth] = size(smooth_error);

% find nonzeros of A^sparsity_level
S0 = spones(A);
S = S0;
for i=2:sparsity_level
    S = S*S0;
end
[row_idx, col_idx, ~] = find(triu(S,1));
n_nonzeros = length(row_idx);

% compute weights: 0.5*|v_i - v_j| sqrt(a_ii^{-1} + a_jj^{-1})
val = zeros(n_nonzeros,1);

%invD = 1./diag(A);  %get diagonal inverse
%invD_weight = invD(row_idx) + invD(col_idx);  % a_ii^{-1} + a_jj^{-1}
D = diag(A);
D_weight = D(row_idx) + D(col_idx);  % a_ii + a_jj

err = (smooth_error(row_idx,:) - smooth_error(col_idx,:)).^2; % (v_i - v_j)^2

for i=1:n_smooth
    %val = val + 0.25*err(:,i).*invD_weight;
    val = val + 0.25*err(:,i).*D_weight;
end
val = sqrt(val)/sqrt(smooth_error'*A*smooth_error);
val = -1./val;

% construct A
A = sparse(row_idx, col_idx, val, n, n);
A = A + A';

end

