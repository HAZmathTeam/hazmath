function [ aggregation,num_agg ] = coarsening_MWM(A )
%Greedy Algorithm for Maximal Matching
% Input: matrix A (negative offdiagonals)
% Output: a maximal matching M

% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)
 
% get size
n = size(A,1);

% initilization
aggregation = zeros(n,1);
num_agg = 0;

% local variables
agg_size = [];

% modify A
D=diag(A);
A = A-spdiags(D,0,n,n);
clear D;

% find nearest neighbor
[~,NN] = min(A,[],1);

% main loop
for i = 1:n
    
    if aggregation(i) == 0  
        
        % adding the isolated points based on the heaviest edge
        %[~,j] = min(A(:,i));
        j = NN(i);
        if aggregation(j)==0 % the heavist neighbor is also isolated
            num_agg = num_agg+1;
            aggregation(i) = num_agg;
            aggregation(j) = num_agg;
            agg_size(num_agg) = 2;
        elseif (agg_size(aggregation(j)) >= 3)
            %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
            % Add this to avoid large diameter (large aggregation size)
            agg_size(aggregation(j)) = agg_size(aggregation(j)) - 1;
            num_agg = num_agg+1;
            aggregation(j) = num_agg;
            aggregation(i) = num_agg;
            agg_size(num_agg) = 2;
        else % the heavist neighbor is in an aggregation of size 2
            aggregation(i)=aggregation(j);
            agg_size(aggregation(j)) = agg_size(aggregation(j)) + 1;
        end
        
    end

end

end

