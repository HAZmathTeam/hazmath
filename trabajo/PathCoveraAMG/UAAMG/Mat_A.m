function [A] = Mat_A(x_grid,DR)
%the right node is Dirichlet (DR = 1=> u(1) = 0) or if it is Neumann 
%(DR=0=>u'(1)=0)

%initialization
n = length(x_grid);
A = zeros(n);
A(1,1) = 1;

%If the right-node is Dirichlet also apply this to the last row u(x_n) = 0
if (DR)
    A(n,n) = 1;
end

% Loop over all the elements (intervals)
for elm=1:n-1
    leftnode = elm ; % Left endpoint of the element (x_{i-1})
    rightnode = elm+1; % Right endpoint of the element (x_{i})
    hi = x_grid(rightnode)-x_grid(leftnode) ; % Element length
    
    % Build the local stiffness matrix from class
    % | 1/ h_i  -1/ h_i |
    % | -1/ h_i  1/ h_i | on element x_{i-1}-x_i
    ALoc = ( 1/hi ).*ones(2);
    ALoc(1,2) = -ALoc(1,2);
    ALoc(2,1) = -ALoc(2,1);
    
    %Global matrix
    if (leftnode==1 &&(DR==1 && rightnode==n))
        % Do nothing. You have 1 element and two nodes and both values are
        % known
    elseif (leftnode==1)
        A(rightnode,rightnode) = A(rightnode,rightnode)+ALoc(2,2);
    elseif (DR==1 && rightnode==n)
        A(leftnode,leftnode) = A(leftnode,leftnode)+ALoc(1,1);
    else
        A(leftnode,leftnode) = A(leftnode,leftnode)+ALoc(1,1);
        A(leftnode,rightnode) = A(leftnode,rightnode)+ALoc(1,2);
        A(rightnode,leftnode) = A(rightnode,leftnode)+ALoc(2,1);
        A(rightnode,rightnode) = A(rightnode,rightnode)+ALoc(2,2);
    end
end
end
