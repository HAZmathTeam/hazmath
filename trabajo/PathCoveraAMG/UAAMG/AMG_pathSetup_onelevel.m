function [ A, b, x, amgParam, amgData, x_smoothed ] = AMG_pathSetup_onelevel(A, b, x, amgParam)
	
level = 1;

%amgParam.agg_type = 'HEC'; % this forces the aggregation on the superficial level is always HEC
%[ amgData ] = AMG_Setup( A, amgParam );

% call multigrid
%x_old = x;
%x_smoothed = AMG_Cycle(amgData, b, x, level, amgParam);
x_smoothed = x;

%v = AMG_Cycle(amgData, zeros(length(b),1), x, level, amgParam);

% compute residual
%r = b - amgData(level).A*x_smoothed;

% compute error
%err(k+1) = norm(r);

% re-aggregating after one HEC
% save HEC smoothed x
%e_smoothed = x_smoothed - x_old;

old_amgParam = amgParam;
amgParam.max_level = 2; % forces to setup only one coarse layer
[ amgData ] = AMG_pathSetup_aff( A, x, amgParam ); %nothing is printing from this subroutine
A = amgData(2).A;
%x = amgData(1).R * x_smoothed;
x = amgData(1).R * x;
b = amgData(1).R * b;
amgParam = old_amgParam;
amgParam.max_level = amgParam.max_level - 1;




end