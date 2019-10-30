function [ x1,f ] = construct_x1_4( n )
%CONSTRUCT_X1 construct initial guess x1 such that the level sets of each
%element f_i(x,y) = val(x) - val(y)
%   input: n, the length away from 0, which make the whole matrix n * n
%   output: f_i(x,y) 
%           x1, vector of length n^2 
% Joanne Lin @Tufts University 8/08/2017


val = linspace(1,5,n);

% initialization
%x1 = zeros((2*n)^2, 1);
f = zeros(n, n);

for i = 1:n
    for j = 1:n
        f(i,j) = (val(i) - val(j))+0.05; %shift from 0
    end
end

%mesh(side, side, f)

x1 = f(:);
end