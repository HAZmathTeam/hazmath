function [ x1,f ] = construct_x1_2( n )
%CONSTRUCT_X1 construct initial guess x1 such that the level sets of each
%element f_i(x,y) = val(y)
%   input: n, the length away from 0, which make the whole matrix n * n
%   output: f_i(x,y) 
%           x1, vector of length n^2 
% Joanne Lin @Tufts University 8/08/2017

%side = -(n-1)/2 : (n-1)/2;
val = linspace(1,5,n);

% initialization
%x1 = zeros((2*n)^2, 1);
f = zeros(n, n);

for i = 1:n
    f(:,i) = val(i);
end

%mesh(side, side, f)

x1 = f(:);
end