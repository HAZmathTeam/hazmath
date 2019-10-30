function [ Adjacency_Matrix ] = find_largest_conn_comp(Adjacency_Matrix)
%This function will use MATLAB graph objects and programs to find the
%largest connected component in the graph, and return the block matrix
%component corresponding to the subgraph from the original adjacency
%matrix. 

Graph = graph(Adjacency_Matrix);
connected_components = conncomp(Graph);

component_numbers = unique(connected_components);
number_of_components = max(size(component_numbers));

max_component_size = -inf;
max_component_index = -1;
for i = 1:number_of_components
    component_size = sum(connected_components == i);
    if component_size > max_component_size
        max_component_index = i;
        max_component_size = component_size;
    end
end

%find the indices corresponding to the largest connected component
largest_component = find(connected_components == max_component_index);

%replace matrix in memory with the submatrix of the largest connected
%component
Adjacency_Matrix = Adjacency_Matrix(largest_component,largest_component);


