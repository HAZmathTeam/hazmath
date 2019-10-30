function [u, iter, residual] = Prec_FGMRES(A, f, u, M_l, M_r, maxit, restart, tol, print_level)
% Preconditioned Flexible General Minimal Residual Method (Left and Right Preconditioner)
%
% BRIEF: 
%   Solving Au=f by general minimal residual using flexible preconditioner
%   M (right preconditioning)
%
%
% INPUT: 
%   A:              Matrix (can be a function handle)
%   f:              Right hand side
%   u:              Initial guess
%   maxit:          Maximial number of itrations allowed
%   restart:        Restart number
%   tol:            Tolerance
%   M:              Precodnitoner (can be a function handle)
%   print_level:    How much information to print (=0: no output; >0 output)
%
% OUTPUT:
%   u:              Solution
%   iter:           Number of iterations
%   residual:          History of l^2 norm of residual
%
% USAGE:
%   [u, iter, residual] = PFGMRES(A, f, u, maxit, restart, tol, [],
%       print_level): PFGMRES without preconditoner
%   [u, iter, residual] = PFGMRES(A, f, u, maxit, restart, tol, M,
%       print_level): PFGMRES using preconditoner M
%
% COPYRIGHT:
% 
%   @ Xiaozhe Hu 02/03/2011, Penn State University

% TODO:
%   1. Check false convergence

%-------------------
% Preparation 
%-------------------
% size of the problem
N = size(f,1);

% restart number
if restart ~= 0
    restart = min(restart, N);
    restart = min(restart,maxit);
else
    restart = min(maxit,N);
end

% initalize memory
r = zeros(N,1);
V = zeros(N, restart+1);
Z = zeros(N, restart);
H = zeros(restart+1, restart);
b = zeros(restart+1,1);
R = zeros(restart,restart);
c = zeros(restart,1);
s = zeros(restart,1);
y = zeros(restart,1);
residual = zeros(maxit+1,1);

%local variables
iter = 1;
%norm_r = 0.0;

% computer the residual
if isa(A, 'double')
    r = f - A*u;
elseif isa(A, 'function_handle')
    r = f - A(u);
else
    error('A is neither a matrix or a function handle!!!');
end % end if

normr = norm(r);

% Left Preconditioning: r = M_l\r 
if isempty(M_l)
    
elseif isa(M_l, 'double')
    r = M_l\r;
elseif isa(M_l, 'function_handle')
    r = M_l(r);
else
    error('Preconditoner M_l is invalid!!!');
end % end if

% store the residual
residual(1) = norm(r);

% output if needed
if (print_level > 0)
    fprintf('----------------------------------------------------\n');
    display(sprintf(' # It |  ||r||/||r0||  |     ||r||      |   Rate.  |'));
    fprintf('----------------------------------------------------\n');
    display(sprintf(' %4d |  %e  |  %e  | %f |', 0, 1.0, residual(1), 0.0));
end

%-------------------
% Main loop 
%-------------------
while (iter < maxit)
    
    % reset converge  
    converge = 0;

    % first orthonormal basis v1
    norm_r = norm(r);
    V(:,1) = r/norm_r;
    
    % form right hand side b for the hessenberg system
    b(1) = norm_r; 
    
    % loop for restart
    for i = 1:restart
        
        % Right Preconditioning: z = M_r\v
        if isempty(M_r)
            Z(:,i) = V(:,i);
        elseif isa(M_r, 'double')
            Z(:,i) = M_r\V(:,i);
        elseif isa(M_r, 'function_handle')
            Z(:,i) = M_r(V(:,i));
        else
            error('Preconditoner M_r is invalid!!!');
        end % end if
        
        % w = Az
        if isa(A, 'double')
            V(:,i+1) = A*Z(:,i);
        elseif isa(A, 'function_handle')
            V(:,i+1) = A(Z(:,i));
        else
            error('A is neither a matrix or a function handle!!!');
        end % end if
        
        % Left Preconditioning: v = M_l\v
        if isempty(M_l)
    
        elseif isa(M_l, 'double')
            V(:,i+1) = M_l\V(:,i+1);
        elseif isa(M_l, 'function_handle')
            V(:,i+1) = M_l(V(:,i+1));
        else
            error('Preconditoner M_l is invalid!!!');
        end % end if
        
        %--------------------------------------
        % modified Gram-Schmidt 
        %--------------------------------------
        for k = 1:i
            
          H(k,i) = V(:,k)'*V(:,i+1);
          V(:,i+1) = V(:,i+1) - H(k,i)*V(:,k);
            
        end % end for k
        
        % new orthonormal basis 
        H(i+1,i) = norm(V(:,i+1));
        V(:,i+1) = V(:,i+1)/H(i+1,i); % becareful small H(i+1,i)
        
        %--------------------------------------
        % Use Givens transformation to get upper triangular system R 
        %--------------------------------------
        R(1,i) = H(1,i);
        
        % apply the previous Givens transformations
        if (i~=1)
            
            for k = 2:i
               
                temp = c(k-1)*R(k-1,i) + s(k-1)*H(k,i);
                R(k,i) = -s(k-1)*R(k-1,i) + c(k-1)*H(k,i);
                R(k-1,i) = temp;
                
            end % end for k
            
        end % end if (i~=1)
        
        % new Givens transformation
        delta = sqrt(R(i,i)^2 + H(i+1,i)^2);
        c(i) = R(i,i)/delta;
        s(i) = H(i+1,i)/delta;
        
        R(i,i) = c(i)*R(i,i) + s(i)*H(i+1,i);
        
        % apply Givens transformation to Right hand side b
        b(i+1) = -s(i)*b(i);
        b(i) = c(i)*b(i);
        
        % count iterations
        iter = iter + 1;
        
        % check convergence b(i+1) = || f-Au_k ||  
        residual(iter) = abs(b(i+1));
        
        % output
        if (print_level > 0)
            fprintf(' %4d |  %e  |  %e  | %f |\n', iter-1, residual(iter)/residual(1), residual(iter), residual(iter)/residual(iter-1));
        end
        
        if ((residual(iter)/residual(1)) < tol || iter > maxit)
            converge = 1;
            break;
        end
        
        %--------------------------------------
        
    end % end for i
    
    % solve the upper trangular matrix
    y(1:i) = R(1:i, 1:i)\b(1:i);
    
    % solution
    u = u + Z(:,1:i)*y(1:i);
    
    % update residual and restart
    if isa(A, 'double')
        r = f - A*u;
    elseif isa(A, 'function_handle')
        r = f - A(u);
    else
        error('A is neither a matrix or a function handle!!!');
    end % end if
    
    % check convergence
    if (converge && ((norm(r)/normr) < tol))
        break;
    else
        % Left Preconditioning: r = M_l\r 
        if isempty(M_l)
    
        elseif isa(M_l, 'double')
            r = M_l\r;
        elseif isa(M_l, 'function_handle')
            r = M_l(r);
        else
            error('Preconditoner M_l is invalid!!!');
        end % end if
    end
    
end % end while iter

% cut residual
iter = iter - 1;
residual = residual(1:iter+1);

end

