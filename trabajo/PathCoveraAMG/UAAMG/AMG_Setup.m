function [ amgData ] = AMG_Setup( Af, amgParam )
% Setup phase for AMG method
%
% @ Xiaozhe Hu, Tufts University

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
    % strong connection
    %----------------
    %As = construct_strong_connection( A{level}, amgParam );
    As= A{level};
    
    %----------------
    % form aggregation
    %----------------
    switch agg_type
       
        case 'HEC',
            %[ aggregation, num_agg ] = coarsening_agg_HEM(A{level});
            [ aggregation, num_agg ] = coarsening_agg_HEC(As);
        case 'MIS',
            %[ aggregation, num_agg ] = coarsening_agg_MIS( A{level}, agg_radius );
            [ aggregation, num_agg ] = coarsening_agg_MIS( As, agg_radius );
        case 'MWM',
            %[ aggregation, num_agg ] = coarsening_agg_HEM(A{level});
            [ aggregation, num_agg ] = coarsening_MWM(As);
        case 'path_matching',
            
            %change weight only on the finest level, make it orthogonal
            %to constant vector
            
            %Use the second eigen vector of E to build the aggregation
             temp_size = size(As,1);
             [V,~] = eigs((eye(temp_size)-DL{level}\As)*(eye(temp_size)-DU{level}\As),10,'lm');
%             [V,~] = eigs((eye(temp_size)-DL{level}\As),10,'lm');
%             %[V,~] = eigs((eye(temp_size)-diag(D{level})\As),10,'lm');
%             x_3 = V(:,3);


%                [V,~] = eigs(As,5,'sm');
                x_2 = V(:,3); % get the e-vector corresponds to the secend smallest e-val
                x_2 = reshape(x_2,[40,40]);
                mesh(1:40,1:40,x_2);
               pause
             x_1 = construct_x1_3(40);
              
              
          %   x_2 = forward_gs(A{level}, zeros(N{level},1), rand(N{level},1), DL{level}, 10);
             A_square = A{level} * A{level};
            [idx, jdx, ~] = find(A_square);
             offdiag_idx = (~(idx == jdx));
            
            %=======using both u2 and u3 to build preconditioner=====%
%             As = sparse(idx(offdiag_idx), jdx(offdiag_idx), ...
%                  -1./(abs( x_2(idx(offdiag_idx)) - x_2(jdx(offdiag_idx)) )+abs( x_3(idx(offdiag_idx)) - x_3(jdx(offdiag_idx)) )), N{level}, N{level});

            %=======using only u2 to build preconditioner=====%
%             As = sparse(idx(offdiag_idx), jdx(offdiag_idx), ...
%                  -1./abs( x_2(idx(offdiag_idx)) - x_2(jdx(offdiag_idx)) ), N{level}, N{level});
            %As = -(As + As');
           
            %=======using x_1 (which the level set is square) to build preconditioner=====%
             As = sparse(idx(offdiag_idx), jdx(offdiag_idx), ...
                  -1./abs( x_1(idx(offdiag_idx)) - x_1(jdx(offdiag_idx)) + 0.001) , N{level}, N{level});
             As = -(As + As')/2;
%              ind = find( As < 500);
%              As(ind) = 0;
            %x = forward_gs(A{level}, zeros(N{level},1), rand(N{level},1), DL{level}, 100);
%             [idx, jdx, ~] = find(A{level});
%             offdiag_idx = (~(idx == jdx));
%             As = sparse(idx(offdiag_idx), jdx(offdiag_idx), ...
%                 -1./abs( x(idx(offdiag_idx)) - x(jdx(offdiag_idx)) ), N{level}, N{level});
%              
            %[ aggregation, num_agg ] = greedy_matching( As );
             [ cover ] = genCover1( As );
             [ aggregation,num_agg ] = pair_match( cover,As );
             
             
             
             G_ori = graph(As,'omitself');
             p = agg_plot(G_ori,aggregation,num_agg,[]);
             pause
             figure
             
%             
%             temp_P = zeros(size(As,1),num_agg);
%             for i = 1:num_agg
%                 [pos,~] = find(aggregation==i);
%                 temp_P(pos,i) = 1;
%             end
%             temp_D = sum(temp_P);
%             D_c = diag(1./temp_D);
%             norm1 = norm((temp_P*D_c*temp_P'-eye(size(temp_P,1)))*x_2,'fro')
          
        case 'path_matching_aff',
             %x1 = construct_x1_1(sqrt(N{level}));  
             %x2 = construct_x1_2(sqrt(N{level}));
             %x3 = construct_x1_3(sqrt(N{level}));
             %x4 = construct_x1_4(sqrt(N{level}));
             %x = [x1 x2 x3 x4];  %initial guess
             %x = [x2 x3];  %initial guess
             %x = x4;
%              figure
%              temp = reshape(x,[20,20]);
%              mesh(1:20,1:20,temp);
%              pause
             %b = zeros(N{level}, 1);

             %solve for several times
             %USE HEC instead of richardson
             %numTimes = 20;
             %x_temp = x;
             x_temp = rand(N{level},1);
%              for i = 1:numTimes;
%                 r = b - A{level}*x_temp;
%                 x_temp = x_temp + 0.6*r./D{level};
%                 x_temp = x_temp - x_temp' * ones(N{level},1) / N{level};
% %                 err_temp = reshape(x_temp,[4,4]);
% %                 mesh(1:4,1:4,err_temp);
% %                 pause
%              end
%              err_temp = reshape(x_temp,[20,20]);
%              mesh(1:20,1:20,err_temp);
%              pause
             
             %use affinity scores to add edges and define weights
             A_square = A{level} * A{level};
             n = size(A_square,1);
             d = diag(A_square);
             A_square = (A_square-spdiags(d,0,n,n));
%              ind = find(A_square);
%              A_square(ind) = -1; % make A^2 unweighted 
             
             % scale to make everything above zero
             add = max(abs(x_temp));
             x_temp = x_temp+add+1;
%              err_temp = reshape(x_temp,[20,20]);
%              mesh(1:20,1:20,err_temp);
             
             % For initial guesses  x_1 and x_2
             %[ A_temp ] = compute_affinity( x_temp, A{level}, 0.01 );
             % For initial guesses  x_3 and x_4
             A_temp = compute_affinity( x_temp, A_square, 0.01 );
             
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
              As = sparse(row, col, val, n, n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %x = forward_gs(A{level}, zeros(N{level},1), rand(N{level},1), DL{level}, 100);
%             [idx, jdx, ~] = find(A{level});
%             offdiag_idx = (~(idx == jdx));
%             As = sparse(idx(offdiag_idx), jdx(offdiag_idx), ...
%                 -1./abs( x(idx(offdiag_idx)) - x(jdx(offdiag_idx)) ), N{level}, N{level});
%              
            %[ aggregation, num_agg ] = greedy_matching( As );
              [ cover ] = genCover1( As ); % found cover on the cutted graph

              [ aggregation,num_agg,iso_edges ] = pair_match( cover, A_temp); %A is for adding back isolated points
              %[ aggregation, num_agg ] = buildmatching4( sqrt(N{level}) );
             
%              figure
%              G_ori = graph(A_square,'omitself');
%              p = agg_plot(G_ori,aggregation,num_agg,iso_edges);
%              %[ p ] = agg_plot2( G_ori, aggregation,num_agg);
%              pause
%             
%             temp_P = zeros(size(As,1),num_agg);
%             for i = 1:num_agg        
%                 [pos,~] = find(aggregation==i);
%                 temp_P(pos,i) = 1;
%             end
%             temp_D = sum(temp_P);
%             D_c = diag(1./temp_D);
%             norm1 = norm((temp_P*D_c*temp_P'-eye(size(temp_P,1)))*x_2,'fro')
          
                  
            
        case 'regular_matching',
            
            
%             temp_size = size(As,1);
            [ aggregation, num_agg ] = greedy_matching( As );
            
            %G_ori = graph(As,'omitself');
            %p = agg_plot(G_ori,aggregation,num_agg,[]);
            
%             [V,~] = eigs((eye(temp_size)-DL{level}\As),3,'lm');
%             x_2 = V(:,2); 
%             x_3 = V(:,3);
%             temp_P = zeros(size(As,1),num_agg);
%             for i = 1:num_agg
%                 [pos,~] = find(aggregation==i);
%                 temp_P(pos,i) = 1;
%             end
%             temp_D = sum(temp_P);
%             D_c = diag(1./temp_D);
%             norm2 = norm((temp_P*D_c*temp_P'-eye(size(temp_P,1)))*x_2,'fro')
            
            
        otherwise,
            display('Wrong aggregation type!');
            display('Use default type.')
            %[ aggregation, num_agg ] = coarsening_agg_MIS( A{level}, agg_radius );
            [ aggregation, num_agg ] = coarsening_agg_MIS( As, agg_radius );
       
    end
    
    clear As;
    
    %----------------
    % generate prolongation
    %----------------
    [ P{level} ] = generate_unsmoothed_P(aggregation, num_agg);
    
    %----------------
    % generate restriction
    %----------------
    R{level} = P{level}';
    
    %check error
%     new_x_1 = ((R{level} * x_1)' * R{level})';
%     x_2 = reshape(new_x_1,[20,20]);
%     mesh(1:20,1:20,x_2);
%     err = norm(x_1 - new_x_1 / 2, 2) %2 is the aggregation size (from matching is two) 
%     pause
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
    level = level+1;
        
end

setup_duration = toc(setup_start);

% construct the data structure
amgData = struct('P',P,'R',R,'A',A,'DL',DL,'DU',DU,'D',D,...
    'IL',IL,'IU',IU,'N',N,'max_level',level);

amgData = amgData(1:level);

% print information
if print_level > 0
   
    total_N = 0;
    total_NNZ = 0;
    
    fprintf('----------------------------------------------------\n');
    fprintf(' # Level |   # Row   |   # Nonzero  | Avg. NNZ/Row |\n');
    fprintf('----------------------------------------------------\n');
    
    for i = 1:level
            
            total_N = total_N + N{i};
            total_NNZ = total_NNZ + nnz(A{i});
            
            fprintf('   %2d    |%9d  | %10d   |   %7.3f    |\n', i, N{i}, nnz(A{i}), nnz(A{i})/N{i});
    end
    
    fprintf('----------------------------------------------------\n');
    
    fprintf(' Grid complexity: %0.3f | Operator complexity: %0.3f \n', total_N/N{1}, total_NNZ/nnz(A{1}));
    
    fprintf('----------------------------------------------------\n');
    
end

% print cputime
% fprintf('----------------------------------------------------\n');
% fprintf('        AMG setup costs %f seconds\n', setup_duration);
% fprintf('----------------------------------------------------\n');


end

