function stat0=m_hazw(fname0,t,p,bcodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% function stat0=m_hazw(fname0,t,p,bcodes);
%%
%% writes a hazmath mesh file "fname0". The input mesh is
%% given by:
%% simplices=t(nt,dim+1);
%% vertex coords=p(dim,nv);
%% boundary codes=bcodes(1:nv,1);
%% element codes=tcodes(1:nt,1) (set to zero not input)
%% values of a function evaluated at every vertex=f0(1:nv,1)
%% also we set f0=0 (no input).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Writes a triangulation t,p constructed with distmesh;
nv = size(p,1);
nt = size(t,1);
dim1=size(t,2);
dim=dim1-1;
fid = fopen(fname0,'w');
if(fid > 0)
  disp(['Writing the mesh to: ',fname0]);
  %% write number of elements, number of nodes dimension
  % and co-homology
     fprintf(fid,' %i %i %i %i\n',nt,nv,dim,0);
     %%% Write T row by row;
     for j = 1 : dim1
      fprintf(fid,' %i ',transpose(t(1:nt,j)));
      fprintf(fid,'\n');
     end
     disp(['...Wrote element-vertex table ']);
     %%% write element codes
     %%%     fprintf(fid,' %i ',transpose(tcodes));
     %%% now set to 0
     fprintf(fid,' %i ',zeros(1,nt));
     fprintf(fid,'\n');
     disp(['...Wrote element codes ']);
     %%% Write coordinates of the nodes;
     for j = 1 : dim
       fprintf(fid,' %23.17g ',transpose(p(1:nv,j)));
       fprintf(fid,'\n');
     end
     disp(['...Wrote coordinates of the vertices ']);
     %%% write codes for the vertices (non-zero is bndry node)
     fprintf(fid,' %i ',transpose(bcodes));
     fprintf(fid,'\n');
     %%% Write function values
     %%%     fprintf(fid,' %g ',transpose(f0));
     %%% now these are just zeroes;
     fprintf(fid,' %g ',zeros(1,nv));
     fprintf(fid,'\n');
      disp(['...Wrote function values ']);
     stat0=fclose(fid);
else
  disp(['error opening file **',fnam0,' **'])
  stat0=-1
end
return
