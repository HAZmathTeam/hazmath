function [ aggregation, num_agg ] = aggregate_MIS( A, isM )
% Generate aggregation based on given maximal indepdent set
%
% @ Xiaozhe Hu, Tufts University 

% prepare
N = size(A,1);
index = (1:N)';
M = index(isM);
num_agg = length(M);
aggregation = zeros(N,1);

% step 1: aggregate those connect to M
[i,j] = find(A(M,:));
aggregation(j) = i;

% step 2: aggregate left-overs
while ~isempty(find(aggregation == 0, 1))

    left = find(aggregation == 0);
    [i,j] = find(A(left,:));
    idx = ~(aggregation(j) == 0);
    i = i(idx);
    j = j(idx);
    aggregation(left(i)) = aggregation(j);

end

% check
if ~isempty(find(aggregation == 0, 1))
    fprintf('Error: some nodes are not aggregated');
end

end

