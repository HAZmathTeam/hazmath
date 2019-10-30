function [ aggregation, num_agg ] = buildmatching1( N )
%build matching for f(x,y) = x
% Input: number of rows/columns of the current level
% Output: a maximal matching M

aggregation = zeros(N^2,1);
num_agg = 1;

for j = 1:floor(N/2)
    for i = 1:N
        if aggregation(N*2*(j-1)+i) == 0 && aggregation(N*2*(j-1)+i + N)== 0
            aggregation(N*2*(j-1)+i) = num_agg;
            aggregation(N*2*(j-1)+i + N) = num_agg; % this is for x only
            num_agg = num_agg+1;
        end
    end
end

isolated_pts = find(aggregation == 0);
num_agg_new = num_agg-1+length(isolated_pts);
aggregation(isolated_pts) = num_agg: num_agg_new;
num_agg = num_agg_new;
end