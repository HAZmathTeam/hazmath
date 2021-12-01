function stat0=write_coo_matrix(fname0,A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% function stat0=m_hazw(fname0,A);
%%
%% writes a sparse matrix A in csr format for hazmath
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Writes a triangulation t,p constructed with distmesh;
nrowa = size(A,1);
ncola = size(A,2);
nnza=nnz(A);
fid = fopen(fname0,'w');
if(fid > 0)
    disp(['Writing the matrix to: ',fname0]);
    fprintf(fid,' %i %i %i\n',nrowa,ncola,nnza);
    [ir,ic,val]=find(A);
    %%% Write A and use C numbering from 0;
    fprintf(fid,' %i %i %.16e\n',[ir-1,ic-1,val]');
    fclose(fid);
else
  disp(['error opening file **',fnam0,' **'])
  stat0=-1
end
end

