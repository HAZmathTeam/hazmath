function [ amgData_adapt, x_HEC_smoothed ] = AMG_pathSetup_adapt(Af, b, x, amgParam)
% call AMG_pathSetup_onelevel and use path cover to build multilevel
% structure
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
x_HEC_smoothed = 0;

setup_start = tic;
while ((level < max_level) && (N{level} > coarsest_size) )
	[ Af, b, x, amgParam, amgData, x_smoothed ] = AMG_pathSetup_onelevel(Af, b, x, amgParam);
	if level == 1
	   x_HEC_smoothed = x_smoothed;
	end
	
	% update amgData
	P{level} = amgData(1).P;
	R{level} = amgData(1).R;
	A{level+1} = amgData(2).A;
	N{level+1} = amgData(2).N;
	DL{level+1} = amgData(2).DL;
	DU{level+1} = amgData(2).DU;
	D{level+1} = amgData(2).D;
	level = level+1;
end
setup_duration = toc(setup_start);

% % print cputime
% fprintf('----------------------------------------------------\n');
% fprintf('        AMG setup costs %f seconds\n', setup_duration);
% fprintf('----------------------------------------------------\n');

% construct the data structure
amgData_adapt = struct('P',P,'R',R,'A',A,'DL',DL,'DU',DU,'D',D,...
    'IL',IL,'IU',IU,'N',N,'max_level',level);

amgData_adapt = amgData_adapt(1:level);
end