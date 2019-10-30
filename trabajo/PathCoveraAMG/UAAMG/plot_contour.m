function plot_contour( solu )
% plot solution countor on a N times N mesh
%
% @ Xiaozhe Hu, Tufts University

N = sqrt(length(solu));
mat_solu = reshape(solu, [N,N]);

solu_max = max(solu);
solu_min = min(solu);
solu_inc = (solu_max - solu_min)/30;
solu_lev = solu_min:solu_inc:solu_max;

%contour(mat_solu, solu_lev, 'ShowText', 'On');
contour(mat_solu, solu_lev, 'LineWidth', 1);


end

