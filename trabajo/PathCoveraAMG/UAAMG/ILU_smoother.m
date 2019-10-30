function [ x ] = ILU_smoother( A, b, x, IL, IU, nsmooth )
% ILU smoother
%
% @ Xiaozhe Hu, Tufts University

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % ILU iteration
    x = x + IU\(IL\(b - A*x));   

end

end

