function [x, err_hist] = forward_gs_revise(A, b, x, DL, nsmooth)
% Forword Gauss-Seidel smoother
%
% @ Xiaozhe Hu and Junyuan Lin, Tufts University

err_hist = zeros(1,nsmooth*length(b));
count = 0;
%----------------------------
% Step 1: Main loop
%----------------------------

B = inv(DL);

for t = 1:nsmooth
    
    %x_prev = x;
    
    for i = 1: length(b)
        x_prev = x;
        
        x(i) = x(i) + B(i,:)*(b - A*x);   

        % Check convergence based on
        % current and previous unknown values
        count = count + 1;
        err_hist(count) = norm(x - x_prev) / (1 + max(norm(x), norm(x_prev)));
        if err_hist(count) < 1e-8
            fprintf('Convergence at iteration %d \n', count);
            break;
        end
    end
    
end
    
end