function [ err, err_norm, err_norm_relative ] = check_approx( vector, pro, res )
% check the approximation property of given prolongationa and restriction
%
% only for unsmoothed aggregation!!
%
% @ Xiaozhe Hu, Tufts University

%Nc = size(pro,2);

temp = res * pro;
diag_inv = 1./sqrt(diag(temp));
     
err = vector - pro * (diag_inv.*(res*vector));
err_norm = norm(err);
err_norm_relative = norm(err)/norm(vector);

end

