function [x, err_hist] = jacobi(A, b, x, D, nsmooth)
% Forword Jacobi smoother
%
% @ Junyuan Lin, Tufts University

err_hist = zeros(1,nsmooth);

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    x_prev = x;
    
    % Jacobi iteration
    x = x + (b - A*x)./D;   
    
    % Check convergence based on
    % current and previous unknown values
    err_hist(i) = norm(x - x_prev) / (1 + max(norm(x), norm(x_prev)));
    if err_hist(i) < 1e-8
        fprintf('Convergence at iteration %d \n', i);
        break;
    end
end
    
end