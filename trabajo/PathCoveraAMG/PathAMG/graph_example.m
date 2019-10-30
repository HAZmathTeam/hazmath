% weighted graph example
edges = [1,2,3;
         1,4,2; 
         2,3,4;
         3,6,5;
         4,2,1;
         4,3,2;
         4,5,1;
         5,3,2;
         5,6,6;
         6,7,2;
         7,3,1];
     
 G = graph(edges(:,1)',edges(:,2)',edges(:,3)');
 
 % adjacency matrix
 A = [0,2,0,3,0,0,0; 2,0,4,1,0,0,0; 0,4,0,2,2,5,1; 3,1,2,0,1,0,0;...
     0,0,2,1,0,6,0; 0,0,5,0,6,0,2; 0,0,1,0,0,2,0];
 
 n = size(A, 1);
 D = diag(A);
 % graph laplacian
 AG = -(A-spdiags(D,0,n,n)); %take diagonals out and use positive weights

 
 clear edges;