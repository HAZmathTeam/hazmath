function [ x, k, err ] = AMG_Solve_Adapt(A, b, x, amgParam, exact_x)
% Adaptive solve phase 
% 
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)
% 

% parameters
print_level = amgParam.print_level;
max_it = amgParam.max_it;
tol = amgParam.tol;

% local iterative method parameters
% iterParam.solver_type = 'CG'; % 'SL': GS as standalone solver | 'AMG': AMG as standalone solver 
%                               % | 'CG': conjugate gradiant | 'FGMRES': flexible GMRes
%                               % 'AMG_path': path cover solver
% iterParam.prec_type = 'AMG'; % 'NULL': no preconditioner | 'AMG': AMG precondtioner
% iterParam.print_level = 1; % how much information to print for iterative solver
%                            % 0: print nothing | positive number print information 
% 
% iterParam.max_it = 1; % maximal number of iterations that is allowed
% iterParam.tol = 1e-6;  % Tolerance for relative residual 
% iterParam.restart = 1; % restart number for flexible GMRes

% prepare solve 
level = 1; 
err = zeros(max_it+1,1);

r = b - A*x;
err(1) = norm(r);  

%setup amgData
amgParam.agg_type = 'MWM'; % this forces the aggregation on the superficial level is always HEC
[ amgData ] = AMG_Setup( A, amgParam ); %the printout of the setup phase is only for the initial AMG setup

if print_level > 0
    fprintf('----------------------------------------------------\n')
    fprintf('              Calling Adaptive AMG solver    \n');
    fprintf('----------------------------------------------------\n')
end

if print_level > 0
    fprintf('----------------------------------------------------\n')
    display(sprintf(' # It |  ||r||/||r0||  |     ||r||      |   Rate.  |'));
    fprintf('----------------------------------------------------\n');
    display(sprintf(' %4d |  %e  |  %e  | %f |', 0, 1.0, err(1), 0.0));
end

% main loop
solve_start = tic;

%resetup = 0;

for k = 1:max_it
    
   % keep the previous solution 
   x_old = x;

    % call multigrid
    x = AMG_Cycle(amgData, b, x, level, amgParam);
    
    %orthorgonalization
    x = x - x' * ones(amgData(1).N,1) / amgData(1).N;
    
    % get the smooth error
    %smooth_error = x-x_old;
    %smooth_error = x/norm(x);
    smooth_error = (exact_x- x)/norm(exact_x-x);
    
    approx_smooth_error = (x-x_old)/norm(x-x_old);
    
    % plot
%     if resetup == 1
%         figure(2)
%         subplot(1,2,1);
%         plot_solution(exact_x- x);
%         subplot(1,2,2);
%         plot_solution(x- x_old);
%         pause
%         resetup = 0;
%     end
        
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
 
    % re-aggregating only when the convergence becomes slower
    if (err(k+1)/err(k) > 0.5) % 1-0.5^#oflayers? 
    %if (abs(1-smooth_error'*approx_smooth_error) < 1e-2 )
    %if ( k>1 && abs(err(k+1)/err(k) - err(k)/err(k-1)) < 1e-3  )
    
         fprintf('-- resetup at iteration %3d --\n', k);
        
        %amgParam.print_level = 0;
        %if (k>=10 )
        %    amgParam.max_level = 10;
        %end  
        %if (k>=16 )
        %    amgParam.max_level = 20;
        %end
        
%         figure(1)
%         subplot(1,2,1);
%         plot_solution(smooth_error);
%         subplot(1,2,2);
%         plot_solution(approx_smooth_error);
%         
%         smooth_error'*approx_smooth_error
%         
%         pause;

        %[ amgData ] = AMG_Setup_Adapt( A, smooth_error, amgParam );
        %[ amgData ] = AMG_Setup_Adapt( A, approx_smooth_error, amgParam );
        
        [ amgData ] = AMG_Setup_pathcover( A, smooth_error, amgParam );
        %[ amgData ] = AMG_Setup_pathcover( A, approx_smooth_error, amgParam );

        %[ amgData ] = AMG_Setup_BFS_path( A, smooth_error, amgParam );
        %[ amgData ] = AMG_Setup_BFS_path( A, approx_smooth_error, amgParam );


        %resetup = 1;
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
