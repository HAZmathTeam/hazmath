function [ As ] = construct_strong_connection( A, amgParam )
% Construct strongly connected matrix
%
% @ Xiaozhe Hu, Tufts University

% paramters
theta = amgParam.strong_connection;

% size of A 
N = size(A,1);

% Normalize A  (diagonal entries have to be positive)
Dinv = spdiags(1./sqrt(diag(A)),0,N,N);
Am = Dinv*A*Dinv;  

% delete weakly connect off-diagonal and diagonal
[im,jm,sm] = find(Am); 
idx = (-sm > theta); 

% matrix for strong connectness
As = sparse(im(idx),jm(idx),sm(idx),N,N); 
As = As + speye(N);    % add diagonal

end

