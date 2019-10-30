function [ amgData ] = AMG_pathSetup_aff( Af,x_temp,amgParam )
%Setup for AMG path cover method. 
%   Use affinity scores to change the weights on the graph
% 
% Setup phase for AMG method
%
% @ Xiaozhe Hu, Joanne Lin, Tufts University

%----------------
% local variable
%----------------
print_level = amgParam.print_level;

max_level = amgParam.max_level;
coarsest_size = amgParam.coarsest_size;

ILU_level = amgParam.ILU_level;
if ILU_level > 0 
    ILU_setup.type = amgParam.ILU_type;  % nofill, crout, ilutp
    ILU_setup.droptol = amgParam.droptol;
    ILU_setup.milu = amgParam.milu;  % row, col, off
    ILU_setup.udiag = amgParam.udiag;
    ILU_setup.thresh = amgParam.thresh;
end


agg_type = amgParam.agg_type;
agg_radius = amgParam.agg_radius;

level = 1;

%----------------
% AMG information
%----------------
P = cell(max_level,1);
R = cell(max_level,1);
A = cell(max_level,1);
N = cell(max_level,1);
DL = cell(max_level,1);
DU = cell(max_level,1);
D = cell(max_level,1);
IL = cell(max_level, 1);
IU = cell(max_level, 1);

%----------------
% finest level
%----------------
A{1} = Af;
N{1} = size(Af,1);
DL{1} = tril(Af);
DU{1} = triu(Af);
D{1} = diag(Af);

%----------------
% main loop
%----------------
% fprintf('----------------------------------------------------\n');
% fprintf('              Calling AMG setup    \n');
% fprintf('----------------------------------------------------\n');

setup_start = tic;

while ((level < max_level) && (N{level} > coarsest_size) )
    
    %----------------
    % form As
    %----------------
    %x1 = construct_x1_1(sqrt(N{level}));  
     %x2 = construct_x1_2(sqrt(N{level}));
     %x3 = construct_x1_3(sqrt(N{level}));
     %x4 = construct_x1_4(sqrt(N{level}));
     %x = [x1 x2 x3 x4];  %initial guess
     %x = [x2 x3];  %initial guess
     %x_temp = x4;
%    figure
%    temp = reshape(x,[20,20]);
%    mesh(1:20,1:20,temp);
%    pause
             %b = zeros(N{level}, size(x,2));

             %solve for several times
             %USE HEC instead of richardson
%              numTimes = 20;
%              %x_temp = x;
%              x_temp = rand(N{level},size(x,2));
%              for i = 1:numTimes;
%                 r = b - A{level}*x_temp;
%                 x_temp = x_temp + 0.6*r./D{level};
%                 x_temp = x_temp - x_temp' * ones(N{level},size(x,2)) / N{level};
% %                 err_temp = reshape(x_temp,[4,4]);
% %                 mesh(1:4,1:4,err_temp);
% %                 pause
%              end
%              err_temp = reshape(x_temp,[20,20]);
%              mesh(1:20,1:20,err_temp);
%              pause
             
             %use affinity scores to add edges and define weights
             %A_square = A{level} * A{level};
             %n = size(A_square,1);
             %d = diag(A_square);
             %A_square = (A_square-spdiags(d,0,n,n));
%              ind = find(A_square);
%              A_square(ind) = -1; % make A^2 unweighted 
             
             % scale to make everything above zero
             %add = max(abs(x_temp));
             %x_temp = x_temp+add+1;
%              err_temp = reshape(x_temp,[20,20]);
%              mesh(1:20,1:20,err_temp);
%              
             % For initial guesses  x_1 and x_2
             [ A_temp ] = compute_affinity( x_temp, A{level}, 0.01 );
             % For initial guesses  x_3 and x_4 and random initial guess
             %A_temp = compute_affinity( x_temp, A_square, 0.01 );
             
             %As = A{level} + A_temp;
             
             %As = A{level};
             As = A_temp;
             %As = A_square + A_temp;
             %As = A_temp;
          %   x_2 = forward_gs(A{level}, zeros(N{level},1), rand(N{level},1), DL{level}, 10);
%              [idx, jdx, ~] = find(A{level});
%              offdiag_idx = (~(idx == jdx));
            
            %=======using both u2 and u3 to build preconditioner=====%
%             As = sparse(idx(offdiag_idx), jdx(offdiag_idx), ...
%                  -1./(abs( x_2(idx(offdiag_idx)) - x_2(jdx(offdiag_idx)) )+abs( x_3(idx(offdiag_idx)) - x_3(jdx(offdiag_idx)) )), N{level}, N{level});

            %=======using only u2 to build preconditioner=====%
%             As = sparse(idx(offdiag_idx), jdx(offdiag_idx), ...
%                  -1./abs( x_2(idx(offdiag_idx)) - x_2(jdx(offdiag_idx)) ), N{level}, N{level});
            %As = -(As + As');
           
            %=======using x_1 (which the level set is square) to build preconditioner=====%
%              As = sparse(idx(offdiag_idx), jdx(offdiag_idx), ...
%                   -1./abs( x_1(idx(offdiag_idx)) - x_1(jdx(offdiag_idx)) + 0.001) , N{level}, N{level});

              As = -(As + As')/2;
              maximum = max(max(As));
%               ind = find( As );
%               threshold = prctile(As(ind),75);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               ind = find( As < maximum*0.75 );
%               As(ind) = 0;
              [row,col,val] = find( As > maximum*0.75 );
	          As = sparse(row, col, val, N{level}, N{level});
              %As = sparse(row, col, val, n, n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %----------------
    % form aggregation
    %----------------
    
    [ cover ] = genCover1( As );
    [ aggregation,num_agg,iso_edges ] = pair_match( cover, A_temp );
    %[ aggregation,num_agg ] = aggregate_matching( As,M );
%     figure
%     G_ori = graph(A_square,'omitself');
%     p = agg_plot(G_ori,aggregation,num_agg,iso_edges);
%      %[ p ] = agg_plot2( G_ori, aggregation,num_agg);
%     pause
%     scale_vec = ones(length(aggregation),1);
%     for i = 1:num_agg
%         ind = find(aggregation==num_agg);
%         size_agg = length(ind);
%         scale_vec(ind) = 1/sqrt(size_agg);
%     end
   
    clear As;
    
    %----------------
    % generate prolongation
    %----------------
    [ P{level} ] = generate_unsmoothed_P(aggregation, num_agg);
%     current_P = P{level};
%     for i = 1 : size(current_P,2)
%         ind = find(current_P(:,i));
%         size_agg = length(ind);
%         current_P(ind,i) = 1/sqrt(size_agg);
%     end


%+++++++++open this for checking scaled error on each level+++++++++++++++
    temp = P{level}' * P{level};
    diagonal_val = 1./sqrt(diag(temp));
    scale_factor = spdiags(diagonal_val,0,size(P{level},2),size(P{level},2));
    current_P = P{level} * scale_factor;
    P{level} = current_P;
    
    %----------------
    % generate restriction
    %----------------
    R{level} = P{level}';
    
    %----------------
    % generate x_temp on the coarse level
    %----------------
    x_temp = x_temp - x_temp' * ones(N{level},1) / N{level};
    err_temp = norm(x_temp - P{level} * R{level} * x_temp)
    x_temp = R{level}*x_temp;
    %pause
    %----------------
    % compute coarse grid matrix
    %----------------
    A{level+1} = R{level}*A{level}*P{level};
    N{level+1} = size(A{level+1},1);
    
    %----------------
    % extra information
    %----------------
    if level <= ILU_level
        [IL{level}, IU{level}] = ilu(A{level}, ILU_setup); 
    end
        
    DL{level+1} = tril(A{level+1});
    DU{level+1} = triu(A{level+1});
    D{level+1} = diag(A{level+1});
    
    %---------------
    % update
    %---------------
    if N{level+1}/N{level} >= 0.9
        amgParam.strong_connection = amgParam.strong_connection/8;
    elseif N{level}/N{level} >= 0.5 && N{level+1}/N{level} < 0.9
        amgParam.strong_connectionn = amgParam.strong_connection/4;
    end
    
    % smooth the error on the coarse level
    %numTimes = 10*level;
%     numTimes = 10;
%     for i = 1:numTimes;
%         r = zeros(N{level+1}, 1) - A{level+1}*x_temp;
%         x_temp = x_temp + 0.6*r./D{level+1};
%         x_temp = x_temp - x_temp' * ones(N{level+1},1) / N{level+1};
%     end
    
    level = level+1;
        
end

setup_duration = toc(setup_start);

% construct the data structure
amgData = struct('P',P,'R',R,'A',A,'DL',DL,'DU',DU,'D',D,...
    'IL',IL,'IU',IU,'N',N,'max_level',level);

amgData = amgData(1:level);



end

