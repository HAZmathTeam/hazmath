function [x] = backward_gs(A, b, x, DU, nsmooth)
% Backward Gauss-Seidel smoother
%
% @ Xiaozhe Hu, Tufts University

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % GS iteration
    x = x + DU\(b - A*x);   
     
end