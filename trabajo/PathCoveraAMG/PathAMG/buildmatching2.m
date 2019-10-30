function [ aggregation, num_agg ] = buildmatching2( N )
%build matching for f(x,y) = y
% Input: number of rows/columns of the current level
% Output: a maximal matching M

aggregation = zeros(N^2,1);
num_agg = 1;

for i = 1:floor(N/2)
    for j = 1:N
        if aggregation(N*(j-1)+2*i) == 0 && aggregation(N*(j-1)+2*i-1)== 0
            aggregation(N*(j-1)+2*i) = num_agg;
            aggregation(N*(j-1)+2*i-1) = num_agg; % this is for y only
            num_agg = num_agg+1;
        end
    end
end

isolated_pts = find(aggregation == 0);
num_agg_new = num_agg-1+length(isolated_pts);
aggregation(isolated_pts) = num_agg: num_agg_new;
num_agg = num_agg_new;
end