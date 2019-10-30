function [ x, iter, residual ] = Krylov_Solve( A, b, x, iterParam, amgParam, amgData)
% Solve phase using Krylov method
%
% @ Xiaozhe Hu, Tufts University

max_it = iterParam.max_it;
tol = iterParam.tol;
print_level = iterParam.print_level;
restart = iterParam.restart;

fprintf('----------------------------------------------------\n');

solve_start = tic;

switch iterParam.solver_type
    
    case 'CG',
        
        fprintf('        Calling Preconditioned CG solver    \n');
        
        switch iterParam.prec_type
            
            case 'AMG'
                
                [x, iter, residual] = Prec_CG(A, b, x, @(r)AMG_prec(r, 1, amgData, amgParam), max_it, tol, print_level);

            otherwise
                
                display('        No preconditioner sepecified, run CG')
                %[x, iter, residual] = Prec_CG(A, b, x, [], max_it, tol, print_level);
                [x,~,~,iter,residual] = pcg(A,b,tol,max_it);
        end
        
    case 'FGMRES',
        
         fprintf('       Calling Preconditioned GMRES solver    \n');
        
        switch iterParam.prec_type
            
            case 'AMG'
                
                [x, iter, residual] = Prec_FGMRES(A, b, x, [], @(r)AMG_prec(r, 1, amgData, amgParam), max_it, restart, tol, print_level);
                
            otherwise
                
                display('       No preconditioner sepecified, run FGMRES')
                [x, iter, residual] = Prec_FGMRES(A, b, x, [], [], max_it, restart, tol, print_level);
        end   
        
    otherwise

        display('Wrong solver type!!')
end

solve_duration = toc(solve_start);
    
fprintf('----------------------------------------------------\n');

fprintf('----------------------------------------------------\n');
if iter == max_it
    fprintf('        Krylov method reached maximal number of iterations \n');
else
    fprintf('        Krylov method converged \n');
    fprintf('        Number of iterations = %d \n', iter);
end
fprintf('        Relative residual = %e \n', residual(iter)/residual(1))
fprintf('----------------------------------------------------\n');
fprintf('        Krylov method costs %f seconds\n', solve_duration);
fprintf('----------------------------------------------------\n');


end

