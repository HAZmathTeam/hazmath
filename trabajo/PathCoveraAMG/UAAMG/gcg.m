function [u] = gcg(A, b, tol, maxit, prec, print_level)
% Generilzed Conjugate Gradient Method 
%
% @ Xiaozhe Hu, Tufts University

% ---- setup ---- % 
n =size(A,1);
u0 = zeros(n,1);
u = u0;
r = b;

p = zeros(n,maxit);
beta = zeros(1, maxit);
% --------------- %

% ---- First iteration (steepest descent) ---- %
Br = feval(prec,r); 
p(:,1) = Br;
alpha = (p(:,1)'*r)/(p(:,1)'*A*p(:,1));
u = u + alpha*p(:,1);
r = r - alpha*A*p(:,1);
% -------------------------------------------- %

if (print_level > 0)
    error = zeros(maxit+1,1);
    error(1) = norm(r);
    display(sprintf('#It|  ||r||/||r0||  |  ||r||   |    Conv. Factor |'));
    display(sprintf(' %d |  %e  |  %e  | %f |', 0, 1.0, error(1), 0.0));
end
    
% ---- main loop ---- %
for k=2:maxit
    
    % --- form p --- %
    Br = feval(prec,r);
    p(:,k) = Br;
    
    for j = 1:k-1
       
        beta(j) = - (p(:,j)'*A*Br)/(p(:,j)'*A*p(:,j));
        
        p(:,k) = p(:,k) + beta(j)*p(:,j);
    
    end
    
    % ---- compute next iteration ---- %
    alpha = (p(:,k)'*r)/(p(:,k)'*A*p(:,k));
    u = u + alpha*p(:,k);
    
    % ---- compute residual ---- %
    r = r - alpha*A*p(:,k);
    
    if (print_level > 0)
        error(k) = norm(r);
        display(sprintf(' %d |  %e  |  %e  | %f |', k, error(k)/norm(b), error(k), error(k)/error(k-1)));
    end
        
    if ((norm(r)/norm(b)) < tol)
        %display(sprintf('gcg converged at iteration %d to a solution with relative residual %d', k, norm(r)/norm(b)));
        break;
    end
    
end

%if (k==maxit) display(sprintf('gcg does not converged within %d iterations', maxit));

end

