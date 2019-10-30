function [ P ] = generate_unsmoothed_P( aggregation, num_agg )
% Construct unsmoothed prolongation P
%
% @ Xiaozhe Hu, Tufts University

N = length(aggregation);
P = sparse((1:N)', aggregation,ones(N,1),N,num_agg);
%P = sparse((1:N)', aggregation,smooth_error,N,num_agg);

end

