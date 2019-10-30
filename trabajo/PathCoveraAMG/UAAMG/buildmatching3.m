function [ aggregation, num_agg ] = buildmatching3( N )
%build matching for f(x,y) = x+y
% Input: number of rows/columns of the current level
% Output: a maximal matching M

aggregation = zeros(N^2,1);
num_agg = 1;

for j = 1:N-1
    for i = 2:N
        if aggregation(N*(j-1)+i) == 0 && aggregation(N*(j-1)+i + N-1)== 0
            aggregation(N*(j-1)+i) = num_agg;
            aggregation(N*(j-1)+i + N-1) = num_agg; % this is for x+y only
            num_agg = num_agg+1;
        end
    end
end

isolated_pts = find(aggregation == 0);
num_agg_new = num_agg-1+length(isolated_pts);
aggregation(isolated_pts) = num_agg: num_agg_new;
num_agg = num_agg_new;
end

