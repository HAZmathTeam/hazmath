function [ p ] = plot_different_methods( size_vec, methods )
%make a plot to compare the performances of differen methods
% inputs: size_vec - a vector containing size of testing matrix (# of row or col)
%         methods - each row represents the CPU time one method used to solve the
%         testing matrix (CG, AMG, Jacobi, LSQR)

% output: p - the plot that has all the methods

for i = 1:size(methods,1)
    p = loglog(size_vec, methods(i,:));
    %p = loglog(size_vec, methods(i,:));
    hold on
end

%p = plot(log(size_vec),2*log(size_vec)-7);
%p = loglog(size_vec,2*size_vec-7);
%hold on
%p = plot(log(size_vec),log(4)*log(size_vec)-4.5);
%p = plot(size_vec,size_vec);
%hold on
%p = plot(size_vec,log(4)*size_vec);

xlabel('\epsilon')
ylabel('CPU time')
legend('AMG w/ HEC','Adaptive AMG','Location','northwest')
%legend('AMG w/ HEC','Adaptive AMG','slope = 1','slope = log(4)','Location','northwest')

end

