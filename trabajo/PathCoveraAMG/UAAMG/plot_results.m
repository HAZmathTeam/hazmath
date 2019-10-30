function [ err_fgs,  err_k, err_rk ] = plot_results( A, b, x, nsmooth, consistent )
%Compare the results of forward gs, jacobi and kaczmarz in terms of relative error
%   solve AA't = b
%   nsmooth - smoothing steps
%   randomized - 0/1 whether use randomized kaczmarz or not
if consistent
    A_new = A*A';
    DL = tril(A_new);
    
    % forward Gauss Seidel
    fgs_time = tic;
    [~, err_fgs] = forward_gs_revise(A_new, b, x, DL, nsmooth);
    fgs_duration = toc(fgs_time)

    % kaczmarz method
    k_time = tic;
    [~, err_k] = kaczmarz_consist(A, b, x, nsmooth, 0);
    k_duration = toc(k_time)
    rk_time = tic;
    [~, err_rk] = kaczmarz_consist(A, b, x, nsmooth, 1);
    rk_duration = toc(rk_time)
else
    A_new = A' * A;
    DL = tril(A_new);
    b_new = A' * b;

    % forward Gauss Seidel
    fgs_time = tic;
    [~, err_fgs] = forward_gs_revise(A_new, b_new, x, DL, nsmooth);
    fgs_duration = toc(fgs_time)

    % kaczmarz method
    k_time = tic;
    [~, err_k] = kaczmarz_inconsist(A, b, x, nsmooth, 0);
    k_duration = toc(k_time)
    rk_time = tic;
    [~, err_rk] = kaczmarz_inconsist(A, b, x, nsmooth, 1);
    rk_duration = toc(rk_time)
end


% plot
figure
plot(1:nsmooth*size(A,1), err_fgs,'bo')
hold on
plot(1:nsmooth*size(A,1), err_k,'kd')
hold on
plot(1:nsmooth*size(A,1), err_rk,'g<')

legend('Forward GS', 'Kaczmarz', 'Randomized Kaczmarz')

end

