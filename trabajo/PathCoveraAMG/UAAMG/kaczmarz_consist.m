function [x, err_hist] = kaczmarz_consist(A, b, x, nsmooth, randomized)
% Kaczmarz algorithm for consistent system
% randomized takes in 1 and 0 to switch on/off the randomization
% @ Junyuan Lin, Tufts University

%------------------------------------------
% Step 1: Calculate probaility mass vector
%------------------------------------------
% 
tmp = A.^2;
pm = sum(tmp, 2) / sum(tmp(:));
% Calculate discrete CDF
F = cumsum(pm);

err_hist = zeros(1,nsmooth*length(x));
count = 0;

%Main loop

for t = 1:nsmooth
    %x_prev = x;

    for i = 1:length(x)
        % Choose row of A at random, with probability 
        % proportional to norm(A(rp, :))^2
        if randomized
            rm = find(F > rand, 1, 'first');
            %rm = randi([1,length(x)],1);
        else
            rm = i;
        end
        a = A(rm, :);
        % Update unknowns 
        x_prev = x;
        
        x = x + (b(rm) - a * x) * a' / (a * a');
        
        % Check convergence based on
        % current and previous unknown values
        count = count+1;
        err_hist(count) = norm(x - x_prev) / (1 + max(norm(x), norm(x_prev)));
        if err_hist(count) < 1e-8
            fprintf('Convergence at iteration %d \n', count);
            break;
        end

    end

    

end
%plot(1:t,err_hist)
end