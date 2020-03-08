function stat0=m_hazw(t,p,fname0);
nssss=size(p)

%%%%  Writes a triangulation t,p constructed with distmesh;
nv = size(p,1)
nt = size(t,1)
dim1=size(t,2);
dim=dim1-1
fid = fopen(fname0,'w');
if(fid > 0)
     disp(['Writing the mesh to: ',fname0]);
     %% write number of elements, number of nodes and twice number
     %of boundary edges
     fprintf(fid,' %i %i %i %i\n',nt,nv,dim,0);
     %%% Write T row by row;
     for j = 1 : dim1
     fprintf(fid,' %i ',transpose(t(1:nt,j)));
     fprintf(fid,'\n');
     end
     %%% Write coordinates of the nodes;
     for j = 1 : dim
     fprintf(fid,' %23.17g ',transpose(p(1:nv,j)));
     fprintf(fid,'\n');
     end
     stat0=fclose(fid);
else
  disp(['error opening file **',fnam0,' **'])
end
