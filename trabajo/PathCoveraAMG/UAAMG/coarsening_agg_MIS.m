function [ aggregation, num_agg ] = coarsening_agg_MIS( A, agg_radius )
% Coarsening based on aggregations obtained from MIS
%
% @ Xiaozhe Hu, Tufts University 

%-----------------
% form power of A
%-----------------
A_radius = A;
for i = 1:agg_radius-1
    A_radius = A_radius*A;
end

A_diam = A_radius;
for i = 1:agg_radius
    A_diam = A_diam*A;
end

%-----------------
% find MIS of A_power
%-----------------
isM = MIS(A_diam);

%-----------------
% form aggregates
%-----------------
[ aggregation, num_agg ] = aggregate_MIS(A_radius, isM);

end

