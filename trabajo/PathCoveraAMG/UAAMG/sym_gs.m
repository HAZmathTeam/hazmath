function [x] = sym_gs(A, b, x, DL, DU, nsmooth)
% Symmetric Gauss-Seidel smoother
%
% @ Junyuan Lin, Tufts University

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % GS iteration (forward)
    x = x + DL\(b - A*x);   
    % GS iteration (backward)
    x = x + DU\(b - A*x);   

end
    
end