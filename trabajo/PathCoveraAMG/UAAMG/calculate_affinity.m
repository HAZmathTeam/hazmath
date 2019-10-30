function [ A,C ] = calculate_affinity( X_k, A_ori, threshold )
% Calculate affinity score based on section 3.3.2 in OREN E. LIVNE and ACHI
% BRANDT's paper LEAN ALGEBRAIC MULTIGRID (LAMG): FAST GRAPH LAPLACIAN LINEAR SOLVER
%
%   Joanne Lin @ Tufts University Math Dept.

% input: X_k - a matrix of current estimation
%        [x_k(i) = f(x,y) = x | x_k(i) = f(x,y) = y | x_k(i) = f(x,y) = x + y | x_k(i) = f(x,y) = x - y], 
%        where f(x,y) is of dimention sqrt(length(x_k)) * sqrt(length(x_k)).
%        x = floor(i/sqrt(length(x_k))), y = mod(i, sqrt(length(x_k)))
%        We reshape each initial guess into a long vector of length n. so X_k is of size n * 4.

% output: A - Adjacency matrix of the new mesh after adding some edges has small 
% affinity c_uv between vertex u and vetex v (the smaller c_uv is, the closer are u 
% and v) 

n = size(X_k,1);
A = zeros(n,n); % Adjacency matrix of the new mesh
C = zeros(n,n); % Affinity scores

A_square = A_ori * A_ori;

for i = 1:n
    for j = i+1:n
        C(i,j) = 1 - (X_k(i,:)*X_k(j,:)')^2/norm(X_k(i,:))^2/norm(X_k(j,:))^2;
        if C(i,j) < threshold
            
             xi = ceil(i/sqrt(n));
             yi = mod(i,sqrt(n));
             %yi = mod(i,n);
             if yi==0 
                 yi = sqrt(n);
             end
             xj = ceil(j/sqrt(n));
             %yj = mod(j,n);
             yj = mod(j,sqrt(n));
             if yj==0 
                 yj = sqrt(n);
             end
             d = norm([xi yi]-[xj yj]);
             if d < 2 || d==2 %look at graph A^2? tried but not that effective
           %if A_square(i,j)~=0
                A(i,j) = -1/(max(abs(X_k(i,:)-X_k(j,:)))+0.001); % should I just take the min from the four results?
            end
        end
    end
end

A = A + A';
A = sparse(A);

end

