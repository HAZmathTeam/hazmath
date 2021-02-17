format long

D= readtable('matrix.txt');

D = D{:,:};

n_dof_for_u = 1090;
n_dof_for_p = 512;


row_size_of_the_matrix = n_dof_for_u + n_dof_for_p;
col_size_of_the_matrix = n_dof_for_u + n_dof_for_p;

A = zeros(row_size_of_the_matrix,col_size_of_the_matrix);

for n=1:size(D)
    row = D(n,1);
    col = D(n,2);
    value = D(n,3);
    A(row,col) = value;
end


A12 = A(1:n_dof_for_u,n_dof_for_u+1:end); 
A12 = transpose(A12);
A21 = A(n_dof_for_u+1:end,1:n_dof_for_u);

A22 = A(n_dof_for_u+1:end,n_dof_for_u+1:end);

C = A12-A21;


[rows columns] = size(C);

for i=1:rows
    for j=1:columns
        if C(i, j) ~= 0
            print("C(%i,%i) = %f", i,j,C(i,j));
            error(0);
        end
    end
end

