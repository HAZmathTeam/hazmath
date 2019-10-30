function [ aggregation,num_agg,iso_edges ] = aggregate_matching( AG,M )
%find the aggregation based on matching

% input: AG - original graph laplacian
%        M - matching on graph of path cover, where M(i) = j and M(j) = i, means i and j get matched
%        together
% output: agg - aggregation, where agg(i) 

%initialization
m = size(AG,1);
d=diag(AG);
AG = AG-spdiags(d,0,m,m);
n = length(M);
aggregation = zeros(n,1);
iso_edges = {};
num_agg = 0;

%aggregate the matched nodes together
for i = 1:n
    neighbor = M(i);
    if neighbor~=0
        if M(neighbor)==i
            if aggregation(neighbor)==0
                num_agg = num_agg+1;
                aggregation(i) = num_agg;
                aggregation(neighbor) = num_agg;
            end
        end
    end
end

%aggregate left-over points
while ~isempty(find(aggregation == 0, 1))
    left = find(aggregation == 0);
    %aggregate the left-over points with the heaviest edge
    for i = 1:length(left)
        [~,j] = max(-AG(left(i),:));
        aggregation(left(i)) = aggregation(j);
        iso_edges{end+1} = [left(i) j];
    end
end


% check
if ~isempty(find(aggregation == 0, 1))
    fprintf('Error: some nodes are not aggregated');
end

end

