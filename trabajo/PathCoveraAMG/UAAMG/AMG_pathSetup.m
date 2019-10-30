function [ amgData ] = AMG_pathSetup( Af,v,amgParam )
%Setup for AMG path cover method. 
%   Use v to change the weights on the graph
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
    % strong connection
    %----------------
    if (level == 1)
        [idx, jdx, ~] = find(A{level});
        offdiag_idx = (~(idx == jdx));
        As = sparse(idx(offdiag_idx), jdx(offdiag_idx), ...
           -1./abs( v(idx(offdiag_idx)) - v(jdx(offdiag_idx)) ), N{level}, N{level});             
    else
        As = construct_strong_connection( A{level}, amgParam );
    end
    
    %----------------
    % form aggregation
    %----------------
    
    [ cover ] = genCover1( A{level} );
    [ aggregation,num_agg ] = pair_match( cover,A{level} );
    %[ aggregation,num_agg ] = aggregate_matching( As,M );
              
    clear As;
    
    %----------------
    % generate prolongation
    %----------------
    [ P{level} ] = generate_unsmoothed_P(aggregation, num_agg);
    
    %----------------
    % generate restriction
    %----------------
    R{level} = P{level}';
    
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



end

