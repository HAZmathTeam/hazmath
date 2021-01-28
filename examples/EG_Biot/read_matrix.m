format long

D= readtable('matrix.txt');

D = D{:,:};

n_dof_for_u = 290;
n_dof_for_p = 209;
timestep = 0.0005;

row_size_of_the_matrix = n_dof_for_u + n_dof_for_p;
col_size_of_the_matrix = n_dof_for_u + n_dof_for_p;

A = zeros(row_size_of_the_matrix,col_size_of_the_matrix);

for n=1:size(D)
    n
    row = D(n,1);
    col = D(n,2);
    value = D(n,3);
    A(row,col) = value;
end


A12 = A(1:n_dof_for_u,n_dof_for_u+1:end); 
A12 = transpose(A12);
A21 = A(n_dof_for_u+1:end,1:n_dof_for_u);
A21 = A21*(-timestep);

C = A12-A21;


[rows columns] = size(C);

for i=1:rows
    for j=1:columns
        if C(i, j)<=1e-10
            C(i, j) = 0;
        end
    end
end
%B=[1 2 3 4; 5 6 7 8 ; 9 10 11 12; 13 14 15 16];
%C=B(1:2,3:4)
%D=B(3:4,1:2)