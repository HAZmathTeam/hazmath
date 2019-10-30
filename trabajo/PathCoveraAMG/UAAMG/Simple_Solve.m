function [ x, k, err ] = Simple_Solve(A, b, x, iterParam)
% Solve phase for simple linear iterative solvers
%
% @ Xiaozhe Hu, Tufts University

% parameters
print_level =iterParam.print_level;
max_it = iterParam.max_it;
tol = iterParam.tol;

% prepare solve 
err = zeros(max_it+1,1);

r = b - A*x;
err(1) = norm(r);  

% print
fprintf('----------------------------------------------------\n')
fprintf('              Calling Simple Linear Iterative solver    \n');
fprintf('----------------------------------------------------\n')

if print_level > 0
    fprintf('----------------------------------------------------\n')
    display(sprintf(' # It |  ||r||/||r0||  |     ||r||      |   Rate.  |'));
    fprintf('----------------------------------------------------\n');
    display(sprintf(' %4d |  %e  |  %e  | %f |', 0, 1.0, err(1), 0.0));
end

% data
DL = tril(A);

% main loop
solve_start = tic;

for k = 1:max_it
    
    % call simple linear iterative method
    x = forward_gs(A, b, x, DL, 1);
    %x = jacobi(A, b, x, D, 1);
    
    %plot_solution(x);
    %pause
    
    % compute residual
    r = b - A*x;
    
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
    fprintf('        Simple linear iterative solver reached maximal number of iterations \n');
else
    fprintf('        Simple linear iterative solver converged \n');
    fprintf('        Number of iterations = %d \n', k);
end
fprintf('        Relative residual = %e \n', err(k+1)/err(1))
fprintf('----------------------------------------------------\n');
fprintf('        Simple linear iterative solve costs %f seconds\n', solve_duration);
fprintf('----------------------------------------------------\n');





end
