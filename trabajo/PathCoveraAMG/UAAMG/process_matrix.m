function [Matrix] = process_matrix(Matrix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function takes in a matrix (sparse or dense) and returns a matrix
%that has one connected component, is symmetric, and has positive edge
%weights. 
%Input:
%   Matrix - 
%       an n x n sparse or dense matrix that will represent the adjacency
%       matrix of a graph. 
%Returns:
%   Matrix - 
%       the same matrix in memory, that will be altered to the
%       specifications mentioned in the commment
%Dependencies:
%   find_largest_conn_comp.m -
%       a function to find and return the adjacency matrix of the largest
%       connected component in the graph. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make symmetric
Matrix  = Matrix + Matrix';

%find largest connected component
Matrix  = find_largest_conn_comp(Matrix);


%make edgeweights positive
Matrix = abs(Matrix);

    
end