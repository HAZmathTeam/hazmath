function [isM] = MIS(A)
% Find maximal indepent set of a given Graph (represented by matrix A)
%
% @ Xiaozhe Hu, Tufts University 

% local variables
N = size(A,1);
isM = false(N,1);       % independent set
isS = false(N,1);        % S: selected set
stop = 1;

% put isolated points in maximal independent set
deg = sum(spones(A)); 
deg = full(deg') - 1;
isS(deg == 0) = true;
isM(deg == 0) = true;

while stop
    
    % prepare for one sweep
    S = find(~isS); 
    n_left = length(S);
    As = A(S,S);
    deg = sum(spones(As)); 
    deg = full(deg') - 1;
    %maxdeg = max(deg);
    measure = deg + rand(n_left,1);
    count = zeros(n_left,1);
    
    % Find marked nodes with local maximum degree
    [i,j] = find(As);    % i,j: index for S
    idx = ( measure(i) > measure(j) );     % compare measure
    
    ii = i(idx);
    [repeat_time, ordered_node] = hist(ii, unique(ii));
    count(ordered_node) = repeat_time;
    
    tempM = S(count == deg);
    
    if sum(tempM) == 0
        stop = 0;
    end
    
    isM(tempM) = true;
    
    % Mark neighboring nodes as selected
    [~,j] = find(A(tempM,:)); 
    isS(j) = true;       
    
end

%fprintf('Number of Maximal Independent Points: %6.0u\n',sum(isM));