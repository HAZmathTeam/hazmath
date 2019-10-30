function [ x, k, err ] = AMG_pathSolve(A, b, x, amgParam)
% Solve phase for AMG method
% @ Joanne Lin based on codes of Xiaozhe Hu, Tufts University
% 

% %generate a random x
% N = size(amgData(1).A,1);
% x = rand(N,1);

% parameters
print_level = amgParam.print_level;
max_it = amgParam.max_it;
tol = amgParam.tol;

% local iterative method parameters
iterParam.solver_type = 'AMG'; % 'SL': GS as standalone solver | 'AMG': AMG as standalone solver 
                              % | 'CG': conjugate gradiant | 'FGMRES': flexible GMRes
                              % 'AMG_path': path cover solver
iterParam.prec_type = 'AMG'; % 'NULL': no preconditioner | 'AMG': AMG precondtioner
iterParam.print_level = 1; % how much information to print for iterative solver
                           % 0: print nothing | positive number print information 

iterParam.max_it = 1; % maximal number of iterations that is allowed
iterParam.tol = 1e-6;  % Tolerance for relative residual 
iterParam.restart = 1; % restart number for flexible GMRes

% prepare solve 
level = 1; 
err = zeros(max_it+1,1);

r = b - A*x;
err(1) = norm(r);  

%setup amgData
amgParam.agg_type = 'HEC'; % this forces the aggregation on the superficial level is always HEC
[ amgData ] = AMG_Setup( A, amgParam ); %the printout of the setup phase is only for the initial AMG setup

% % print
% fprintf('----------------------------------------------------\n')
% fprintf('              Calling AMG solver    \n');
% fprintf('----------------------------------------------------\n')

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
    x_old = x;
    x = AMG_Cycle(amgData, b, x, level, amgParam);
    %x = Krylov_Solve(A, b, x, iterParam, amgParam, amgData);
    %smooth_error = x-x_old;
    smooth_error = x;
        
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
    if (err(k+1)/err(k) > 0.6 || k==100) % 1-0.5^#oflayers?
        
         fprintf('-- resetup at iteration %3d --\n', k);
        
        % use path cover to construct multilevel
        %[ amgData ] = AMG_pathSetup_aff( A,x-x_old,amgParam );%nothing is printing from this subroutine
        %[ amgData ] = AMG_pathSetup_aff( A,v,amgParam );%nothing is printing from this subroutine
        % use HEC + path cover to construct multilevel
        %norm(x)
        %pause
        %[ amgData ] = AMG_pathSetup_adapt1(A, e_smoothed, amgParam); %path
        %cover + HEC
        %x = e_smoothed;
        %[ amgData,x ] = AMG_pathSetup_adapt(A, b, x, amgParam); %path cover only
        
        %amgParam.print_level = 0;
        %if (k>=10 )
        %    amgParam.max_level = 2;
        %end  
        %if (k>=16 )
        %    amgParam.max_level = 20;
        %end
        [ amgData ] = AMG_Setup_pathcover( A, smooth_error, amgParam );
        %x = x + (smooth_error'*r)/(smooth_error'*A*smooth_error)*smooth_error; 
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
