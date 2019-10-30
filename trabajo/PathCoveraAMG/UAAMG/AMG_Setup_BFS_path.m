function [ amgData ] = AMG_Setup_BFS_path( Af, smooth_error, amgParam )
% Setup phase for AMG method based on BFS paths
%
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)

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

if print_level > 0
    fprintf('----------------------------------------------------\n');
    fprintf('     Calling AMG setup based on path covering    \n');
    fprintf('----------------------------------------------------\n');
end
    
setup_start = tic;

%----------------
% get level set from smooth_error
%----------------
% get adjancy matrix
Adj = Af - spdiags(D{1}, 0, N{1}, N{1});

% get level set
[~,max_idx] = max(smooth_error);
[~,min_idx] = min(smooth_error);
s = [max_idx, min_idx];
[ level_set ] = BFS_level_set( Adj, s );

% compute strength matrix based on smooth_error
[ As ] = compute_similarity( smooth_error, Af, 2 );
%[ As ] = compute_weights(smooth_error, Af, 1 );

% get rid of edges between level_sets
%A2adj = As - spdiags(diag(As), 0, N{1}, N{1});
[ Adj_level_set ] = cut_edge_between_level_set( As, level_set );

% generate path
[ path ] = generate_BFS_path( Adj_level_set, level_set, smooth_error );
num_path = length(path);

% plot BFS path cover
% figure(1);
% subplot(1,2,1);
% plot_solution(smooth_error);
% 
% subplot(1,2,2);
% plot_contour(smooth_error);
% 
% figure(2);
% G = graph(A_affinity, 'OmitSelf');
% n = sqrt(N{1});
% xcoord = 1:n; xcoord = repmat(xcoord,n,1);
% ycoord = (1:n)'; ycoord = repmat(ycoord, 1,n);
% plot_cover( G, xcoord(:), ycoord(:), path );
% axis([1 n 1 n]);
% 
% pause

%----------------
% main loop
%----------------
while ((level < max_level) && (N{level} > max(coarsest_size, num_path)) )
    
    %----------------
    % form aggregation along the paths
    %----------------
    %[ aggregation,num_agg, ~ ] = pair_match_XH( path, A{level} );
    [ aggregation, num_agg, path ] = path_matching( path, As );

    %----------------
    % generate prolongation
    %----------------
    if  (level == 1)
        [ P{level} ] = generate_unsmoothed_P_smooth_error( aggregation, num_agg, smooth_error);
    else
        [ P{level} ] = generate_unsmoothed_P(aggregation, num_agg);
    end
        
    %----------------
    % generate restriction
    %----------------
    R{level} = P{level}';
    
%     %----------------
%     %+++++++++open this for checking scaled error on each level+++++++++++++++
%     temp = P{level}' * P{level};
%     diagonal_val = 1./sqrt(diag(temp));
%     scale_factor = spdiags(diagonal_val,0,size(P{level},2),size(P{level},2));
%     P{level} = P{level} * scale_factor;
%     R{level} = P{level}';
%     
%     %smooth_error = smooth_error - smooth_error' * ones(N{level},1) / N{level};
%     err_temp = norm(smooth_error - P{level} * R{level} * smooth_error)
%     %----------------
    
    %----------------
    % compute coarse grid matrix
    %----------------
    A{level+1} = R{level}*A{level}*P{level};
    N{level+1} = size(A{level+1},1);
    
    %----------------
    % update cover and affinity
    %----------------
    %[ path ] = coarsen_pathcover(path, aggregation);
    As = R{level}*As*P{level};
    
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

if print_level > 0
    fprintf('----------------------------------------------------\n');
    fprintf('    path cover AMG setup costs %f seconds\n', setup_duration);
    fprintf('----------------------------------------------------\n');
end

end

