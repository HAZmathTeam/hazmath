function [ u, k, residual ] = Prec_CG(A, f, u, M, maxit, tol, print_level)
% Preconditioned Conjugate Gradient Method 
%
% @ Xiaozhe Hu, Tufts University

%-------------------
% Preparation 
%-------------------
% size of the problem
N = size(f,1);

r = zeros(N,1);
%p = zeros(N,1);
z = zeros(N,1);
Ap = zeros(N,1);
residual = zeros(maxit+1,1);

%alpha = 0.0;
%beta = 0.0;
%rz = 0.0;

if isa(A, 'double')
    r = f - A*u;
elseif isa(A, 'function_handle')
    r = f - A(u);
else
    error('A is neither a matrix or a function handle!!!');
end % end if
normr = norm(r);
residual(1) = normr;

if (print_level > 0)
    fprintf('----------------------------------------------------\n');
    display(sprintf(' # It |  ||r||/||r0||  |     ||r||      |   Rate.  |'));
    fprintf('----------------------------------------------------\n');
    display(sprintf(' %4d |  %e  |  %e  | %f |', 0, 1.0, residual(1), 0.0));
end

% Preconditioning: z = M\r
if isempty(M)
    z = r;
elseif isa(M, 'double')
    z = M\r;
elseif isa(M, 'function_handle')
    z = M(r);
else
    error('Preconditoner M is invalid!!!');
end % end if

p = z;

rz = r'*z;

%-------------------
% Main loop 
%-------------------
for k = 1:maxit
    
    % Ap
    if isa(A, 'double')
        Ap = A*p;
    elseif isa(A, 'function_handle')
        Ap = A(p);
    else
        error('A is neither a matrix or a function handle!!!');
    end % end if
    
    % alpha = (r,z)/(Ap,p)
    alpha = rz/(Ap'*p);
    
    % u = u + alpha*p
    u = u + alpha*p;
    
    % r = r - alpha*Ap
    r = r - alpha*Ap;
    residual(k+1) = norm(r);
 
    % display
    if (print_level > 0)
        display(sprintf(' %4d |  %e  |  %e  | %f |', k, residual(k+1)/residual(1), residual(k+1), residual(k+1)/residual(k)));
    end
        
    if ((residual(k+1)/residual(1)) < tol)
        break;
    end
    
    % Preconditioning: z = M\r
    if isempty(M)
        z = r;
    elseif isa(M, 'double')
        z = M\r;
    elseif isa(M, 'function_handle')
        z = M(r);
    else
        error('Preconditoner M is invalid!!!');
    end % end if
    
    % (r_new, z_new)
    rz_new = r'*z;
    
    % beta = (r_new, z_new)/(r,z)
    beta = rz_new/rz;
    rz = rz_new;
    
    % p = z + beta*p
    p = z + beta*p;
    
end %for

% cut residual
residual = residual(1:k+1);

end

