function [ x1,f ] = construct_x1( n )
%CONSTRUCT_X1 construct initial guess x1 such that the level sets of each
%element f_i(x,y) are squares (i.e. f_i(x,y) = c, for (x, y) is a point on the square)
%   input: n, the length away from 0, which make the whole matrix n * n
%   output: f_i(x,y) are squares (i.e. f_i(x,y) = c, for (x, y) is a point on the square)
%           x1, vector of length n^2 
% Joanne Lin @Tufts University 6/28/2017

side = -(n-1)/2 : (n-1)/2;

% initialization
%x1 = zeros((2*n)^2, 1);
f = zeros(n, n);

if mod(n,2) == 1  % n is odd
    for i = 1 : (n-1)/2
        for j = 1 : i
            f((n-1)/2 + 1 + i, (n-1)/2 + 1 - j) = (n-1)/2 + 1 - i;
            f((n-1)/2 + 1 + i, (n-1)/2 + j) = (n-1)/2 + 1 - i;
        end
    end
    f = f + rot90(f, 1) + rot90(f, 2) + rot90(f, 3);
    f((n-1)/2 + 1, (n-1)/2 + 1) =  (n-1)/2 + 1;
else  % n is even
    for i = 1 : ceil((n-1)/2)
        for j = 1 : i
            f(ceil((n-1)/2) + i, ceil((n-1)/2) + 1 - j) = ceil((n-1)/2) + 1 - i;
            f(ceil((n-1)/2) + i, floor((n-1)/2) + j) = ceil((n-1)/2) + 1 - i;
        end
    end
    f = f + rot90(f, 1) + rot90(f, 2) + rot90(f, 3);
    %f(floor((n-1)/2) + 1, floor((n-1)/2) + 1) =  floor((n-1)/2) + 1;
end
mesh(side, side, f)

x1 = f(:);
end

