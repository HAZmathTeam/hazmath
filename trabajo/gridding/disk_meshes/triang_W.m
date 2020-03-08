  function stat0=triang_W(t,p,fname0);
%%%%  Writes a triangulation t,p constructed with distmesh; 
n = length(p(1,:));
nel = length(t(1,:))
fid = fopen(fname0,'w');
if(fid > 0)
     disp(['Writing the mesh to: ',fname0]);
     %% write number of elements, number of nodes and twice number
     %of boundary edges
fprintf(fid,' %i %i %i %i\n',nel,n,ned,ned);
%%% Write T row by row;
for j = 1 : 4
 fprintf(fid,' %i ',t(j,1:nel));
 fprintf(fid,'\n');
end
%%% Write coordinates of the nodes;
fprintf(fid,' %23.17e ',p(1,1:n));
fprintf(fid,'\n');
fprintf(fid,' %23.17e ',p(2,1:n));
fprintf(fid,'\n');
for j = 1 : 2
%%% When j=1: Write the node numbers of the first end point of the boundary edges;
%%% When j=2: Write the node numbers of the SECOND end point of the boundary edges;
 fprintf(fid,' %i ',e(j,1:ned));
 fprintf(fid,'\n');
end
%% ignore these. 
for j = 3 : 4
fprintf(fid,' %23.17e ',e(j,1:ned));
 fprintf(fid,'\n');
end
for j = 5 : 7
 fprintf(fid,' %i ',e(j,1:ned));
 fprintf(fid,'\n');
end
stat0=fclose(fid);
else
disp(['error opening file **',fnam0,' **'])
end
