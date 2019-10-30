function [ x, k, err ] = AMG_Solve(amgData, b, x, amgParam)
% Solve phase for AMG method
%
% @ Xiaozhe Hu and Junyuan Lin, Tufts University

% parameters
print_level = amgParam.print_level;
max_it = amgParam.max_it;
tol = amgParam.tol;

% prepare solve 
level = 1; 
err = zeros(max_it+1,1);

r = b - amgData(1).A*x;
err(1) = norm(r);  

% print
fprintf('----------------------------------------------------\n')
fprintf('              Calling AMG solver    \n');
fprintf('----------------------------------------------------\n')

if print_level > 0
    fprintf('----------------------------------------------------\n')
    display(sprintf(' # It |  ||r||/||r0||  |     ||r||      |   Rate.  |'));
    fprintf('----------------------------------------------------\n');
    display(sprintf(' %4d |  %e  |  %e  | %f |', 0, 1.0, err(1), 0.0));
end

% main loop
solve_start = tic;

for k = 1:max_it
    
    % call multigrid
    x = AMG_Cycle(amgData, b, x, level, amgParam);
    
    % compute residual
    r = b - amgData(level).A*x;
    %err_temp = x + amgData(level).DL / r;
%     err_temp = x;
%     err_temp = reshape(err_temp,[20,20]);
%     mesh(1:20,1:20,err_temp);
%     pause
    
    % compute error
    err(k+1) = norm(r);
    
    % display
    if print_level > 0
        display(sprintf(' %4d |  %e  |  %e  | %f |', k, err(k+1)/err(1), err(k+1), err(k+1)/err(k)));
    end
    
    if ((err(k+1)/err(1)) < tol)
        break;
    end
    
end

solve_duration = toc(solve_start);

% cut err
err = err(1:k+1);

% print
fprintf('----------------------------------------------------\n');
if k == max_it
    fprintf('        AMG reached maximal number of iterations \n');
else
    fprintf('        AMG converged \n');
    fprintf('        Number of iterations = %d \n', k);
end
fprintf('        Relative residual = %e \n', err(k+1)/err(1))
fprintf('----------------------------------------------------\n');
fprintf('        AMG solve costs %f seconds\n', solve_duration);
fprintf('----------------------------------------------------\n');





end
