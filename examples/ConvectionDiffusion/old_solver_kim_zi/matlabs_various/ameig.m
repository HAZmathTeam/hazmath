load fort.201
A = sparse(fort(1,:),fort(2,:),fort(3,:));
load fort.202
M = sparse(fort(1,:),fort(2,:),fort(3,:));
clear fort
e = eig(full(A));
ee = eig(full(A),full(M));
eesort = sort(ee);


% Exact eigvalues of Au=\lambda Mu.
%2 - 22.8657759367719
%3 - 20.5055448977079
%4 - 19.9297898422163
%5 - 
