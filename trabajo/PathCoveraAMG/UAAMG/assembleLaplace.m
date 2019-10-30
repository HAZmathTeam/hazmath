function [ L ] = assembleLaplace(N)
%
% Copyright (C)  Xiaozhe Hu.

e = ones(N,1);

L1d = spdiags([-1*e 2*e -1*e], -1:1, N, N);
I = speye(N,N);

L =  (kron(L1d, I) + kron(I, L1d));


end

