function [ amgData_adapt] = AMG_pathSetup_adapt1(Af, e_smoothed, amgParam)
% call AMG_pathSetup_onelevel to build one layer and HEC for the rest
% @ Joanne Lin, Tufts University
% 

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


setup_start = tic;
%while ((level < max_level) && (N{level} > coarsest_size) )
    [ Af, amgParam, amgData ] = AMG_pathSetup_onelevel1(Af, e_smoothed, amgParam);
%     max_level = amgParam.max_level; %amgParam.max_level is reduced by 1
%     each time when AMG_pathSetup_onelevel is called
%     if level == 1
%         x_HEC_smoothed = x_smoothed;
%     end
    % update amgData
    P{1} = amgData(1).P;
    R{1} = amgData(1).R;
    A{2} = amgData(2).A;
    N{2} = amgData(2).N;
    DL{2} = amgData(2).DL;
    DU{2} = amgData(2).DU;
    D{2} = amgData(2).D;
    
    %level = level+1;
%end
amgParam.agg_type = 'HEC';
[ amgData_HEC ] = AMG_Setup( Af, amgParam );
level_HEC = amgData_HEC(1).max_level;

for level = 2 : level_HEC
    P{level} = amgData_HEC(level-1).P;
    R{level} = amgData_HEC(level-1).R;
    A{level+1} = amgData_HEC(level).A;
    N{level+1} = amgData_HEC(level).N;
    DL{level+1} = amgData_HEC(level).DL;
    DU{level+1} = amgData_HEC(level).DU;
    D{level+1} = amgData_HEC(level).D;
end

setup_duration = toc(setup_start);

% % print cputime
% fprintf('----------------------------------------------------\n');
% fprintf('        AMG setup costs %f seconds\n', setup_duration);
% fprintf('----------------------------------------------------\n');

% construct the data structure
amgData_adapt = struct('P',P,'R',R,'A',A,'DL',DL,'DU',DU,'D',D,...
    'IL',IL,'IU',IU,'N',N,'max_level',level);

amgData_adapt = amgData_adapt(1:level_HEC+1);
end