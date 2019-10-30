function [ A, amgParam, amgData ] = AMG_pathSetup_onelevel1(A, e_smoothed, amgParam)
% Adaptive setup by solving (HEC on multilevel then one layer down by pathcover)
% @ Joanne Lin based on codes of Xiaozhe Hu, Tufts University
% 

% %generate a random x
% N = size(amgData(1).A,1);
% x = rand(N,1);

% parameters
% print_level = amgParam.print_level;
% max_it = amgParam.max_it;
% tol = amgParam.tol;

% prepare solve 
%level = 1; 
% err = zeros(max_it+1,1);
% 
% r = b - A*x;
% err(1) = norm(r);

%x_smoothed = x;

%setup amgData
%amgParam.agg_type = 'HEC'; % this forces the aggregation on the superficial level is always HEC
%[ amgData ] = AMG_Setup( A, amgParam ); %the printout of the setup phase is only for the initial AMG setup
% % print
% fprintf('----------------------------------------------------\n')
% fprintf('              Calling AMG solver    \n');
% fprintf('----------------------------------------------------\n')

% if print_level > 0
%     fprintf('----------------------------------------------------\n')
%     display(sprintf(' # It |  ||r||/||r0||  |     ||r||      |   Rate.  |'));
%     fprintf('----------------------------------------------------\n');
%     display(sprintf(' %4d |  %e  |  %e  | %f |', 0, 1.0, err(1), 0.0));
% end

% main loop
%solve_start = tic;

%for k = 1:max_it
    
    % call multigrid
    %x_old = x;
    %x = AMG_Cycle(amgData, b, x, level, amgParam);
    
    %v = AMG_Cycle(amgData, zeros(length(b),1), x, level, amgParam);
    
%     % compute residual
%     r = b - amgData(level).A*x;
%     
%     % compute error
%     err(k+1) = norm(r);
    
%     % display
%     if print_level > 0
%         display(sprintf(' %4d |  %e  |  %e  | %f |', k, err(k+1)/err(1), err(k+1), err(k+1)/err(k)));
%     end
%         
%     if ((err(k+1)/err(1)) < tol)
%         break;
%     end
 
    % re-aggregating only when the convergence becomes slower
    %if (err(k+1)/err(1) > 0.68*err(k)/err(1)) % 1-0.5^#oflayers?
    %if k > 1
        % save HEC smoothed x
        %x_smoothed = x - x_old;
        
%         if length(x_smoothed)==400
%             err_temp = x_smoothed;
%             err_temp = reshape(err_temp,[20,20]);
%             mesh(1:20,1:20,err_temp);
%             pause
%         end
        
        %x_smoothed = x;
        old_amgParam = amgParam;
        amgParam.max_level = 2; % forces to setup only one coarse layer
        [ amgData ] = AMG_pathSetup_aff( A,e_smoothed,amgParam );%nothing is printing from this subroutine
        A = amgData(2).A;
%         x = amgData(1).R * x_smoothed;
%         b = amgData(1).R * b;
        amgParam = old_amgParam;
        amgParam.max_level = amgParam.max_level - 1;
        %break;
        %[ amgData ] = AMG_pathSetup_aff( A,v,amgParam );%nothing is printing from this subroutine
        
    %end
    
%end



% solve_duration = toc(solve_start);

% cut err
%err = err(1:k+1);

% % print
% fprintf('----------------------------------------------------\n');
% if k == max_it
%     fprintf('        AMG reached maximal number of iterations \n');
% else
%     fprintf('        AMG converged \n');
%     fprintf('        Number of iterations = %d \n', k);
% end
% fprintf('        Relative residual = %e \n', err(k+1)/err(1))
% fprintf('----------------------------------------------------\n');
% fprintf('        AMG solve costs %f seconds\n', solve_duration);
% fprintf('----------------------------------------------------\n');





end