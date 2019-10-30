function [ L ] = assembleJumpLaplace(N, epsilon, jump_type, bd_type)
%
% Copyright (C)  Xiaozhe Hu.

h = 1/(N);
L = sparse(N^2,N^2);


% get Jump coefficients for off-diagonals
switch (jump_type)
    case 1
        jump = isolateIsland(N, epsilon);
    case 2
        jump = checkerboard(N, epsilon);
    case 3 
        jump = isolateIsland_rand(N, epsilon);
    case 4 
        jump = checkerboard_rand(N, epsilon);
    otherwise
        error('wrong jump type');
end
        
% generate off-diagonal (offset = -1)
sw = (-2.*jump(1:N,2:N+1).*jump(2:N+1,2:N+1))./(jump(1:N,2:N+1) + jump(2:N+1,2:N+1));
L = L + spdiags(sw(:),1,N^2,N^2);

% generate off-diagonal (offset = 1)
se = (-2.*jump(3:N+2,2:N+1).*jump(2:N+1,2:N+1))./(jump(3:N+2,2:N+1) + jump(2:N+1,2:N+1));
L = L + spdiags(se(:),-1,N^2,N^2);

% generate off-diagonal (offset = -N)
ss = (-2.*jump(2:N+1,1:N).*jump(2:N+1,2:N+1))./(jump(2:N+1,1:N) + jump(2:N+1,2:N+1));
L = L + spdiags(ss(:),N,N^2,N^2);

% generate off-diagonal (offset = N)
sn = (-2.*jump(2:N+1,3:N+2).*jump(2:N+1,2:N+1))./(jump(2:N+1,3:N+2) + jump(2:N+1,2:N+1));
L = L + spdiags(sn(:),-N,N^2,N^2);

% generate diagonal

if (bd_type == 1)
   
    % ghost
    jump(1,:) = jump(2,:);
    jump(N+2,:) = jump(N+1,:);
    jump(:,1) = jump(:,2);
    jump(:,N+2) = jump(:,N+1);
    
    sw = (-2.*jump(1:N,2:N+1).*jump(2:N+1,2:N+1))./(jump(1:N,2:N+1) + jump(2:N+1,2:N+1));
    sw(1,:) = 2*sw(1,:);
    se = (-2.*jump(3:N+2,2:N+1).*jump(2:N+1,2:N+1))./(jump(3:N+2,2:N+1) + jump(2:N+1,2:N+1));
    se(N,:) = 2*se(N,:);
    ss = (-2.*jump(2:N+1,1:N).*jump(2:N+1,2:N+1))./(jump(2:N+1,1:N) + jump(2:N+1,2:N+1));
    ss(:,1) = 2*ss(:,1);
    sn = (-2.*jump(2:N+1,3:N+2).*jump(2:N+1,2:N+1))./(jump(2:N+1,3:N+2) + jump(2:N+1,2:N+1));
    sn(:,N) = 2*sn(:,N);
    
end

sc = -1.*(sw + se + ss + sn);
L = L + spdiags(sc(:),0,N^2,N^2);

% get matrix
L = L'/(h^2);

end

%--------------------------
% generate Jump Coefficient
%--------------------------
function jump = isolateIsland(N,epsilon)
    
    jump = zeros(N+2, N+2);
    jump(2:N+1, 2:N+1) = 1;

    i = 2:2:N+1; j = 2:2:N+1;
    jump(i,j) = epsilon; 
    
end

function jump = checkerboard(N,epsilon)
    
    jump = zeros(N+2, N+2);
    jump(2:N+1, 2:N+1) = 1;

    i = 2:2:N+1; j = 2:2:N+1;
    jump(i,j) = epsilon; 

    theta = floor(log(epsilon)/log(10));
    
    if (N+1>3)
        i = 3:2:N+1; j = 3:2:N+1;
        jump(i,j) = 10^(-theta);
    end

end

function jump = isolateIsland_rand(N, epsilon)
    
    jump = zeros(N+2, N+2);
    jump(2:N+1, 2:N+1) = 1;
    
    theta = floor(log(epsilon)/log(10));

    i = 2:2:N+1; j = 2:2:N+1;
    jump(i,j) = 10.^(randi([theta, 0], length(i), length(j))); 

end

function jump = checkerboard_rand(N, epsilon)
    
    jump = zeros(N+2, N+2);
    jump(2:N+1, 2:N+1) = 1;
    
    theta = floor(log(epsilon)/log(10));

    i = 2:2:N+1; j = 2:2:N+1;
    jump(i,j) = 10.^(randi([theta, 0], length(i), length(j))); 

    if (N+1>3)
        i = 3:2:N+1; j = 3:2:N+1;
        jump(i,j) = 10.^(randi([theta, 0], length(i), length(j)));
    end
    
end

% function jump = total_rand(N, epsilon)
% 
%     jump = zeros(N+2, N+2);
%     
%     theta = floor(log(epsilon)/log(10));
% 
%     jump(2:N+1, 2:N+1) = 10.^(randi([theta, -theta], N, N)); 
% 
% end

