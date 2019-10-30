function [ P ] = generate_unsmoothed_P_smooth_error( aggregation, num_agg, smooth_error)
% Construct unsmoothed prolongation P using smooth_error
%
% @ Xiaozhe Hu, Tufts University

N = length(aggregation);
P = sparse((1:N)', aggregation,smooth_error,N,num_agg);

end

