function [ z ]=AMG_prec(r, level, amgData, amgParam)
% AMG preconditioner
%
% @ Xiaozhe Hu, Tufts University

%-------------------
% Preparation 
%-------------------
max_it = amgParam.prec_max_it;

N = size(r,1);
z = zeros(N,1); 

%-------------------
% Main loop 
%-------------------
for k = 1:max_it
    
    % call multigrid
    z = AMG_Cycle(amgData, r, z, level, amgParam);

end