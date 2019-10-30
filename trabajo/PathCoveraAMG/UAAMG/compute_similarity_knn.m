function [ A ] = compute_similarity_knn( coord, k  )
% similarity (k nearest neighbor)
%
% assume just one coordinates !! (knn in 1D)
%
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)

% get size
[n, ~] = size(coord);

% order coordinates
[~, order] = sort(coord);

% local variable
temp = zeros(n, 2*k);

if (k==1)
    temp(1,:) = order(2:3)';
    temp(n,:) = order(n-2:n-1)';
    temp(2:n-1,1) = order(1:n-2);
    temp(2:n-1,2) = order(3:n);
    
    row = repmat(order,2,1);
    col = temp(:);
    clear temp;
    
    Atemp = sparse(row, col, ones(2*n,1), n, n);
    Atemp = Atemp + Atemp';
    clear row col;
    
    [row_idx, col_idx, ~] = find(triu(Atemp,1));
    clear temp;
    val = -1./abs(coord(row_idx) - coord(col_idx));
    A = sparse(row_idx, col_idx, val, n, n);
    A = A + A';  
end
    

end

