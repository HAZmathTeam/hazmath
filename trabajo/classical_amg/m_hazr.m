function [t,x,ib,el_code,status]=m_hazr(fname0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% function status=m_hazw(fname0,t,p,bcodes);
%%
%% reads a hazmath mesh file "fname0" and returns the triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Reads a triangulation t,p constructed with distmesh;
fid = fopen(fname0,'r');
if(fid > 0)
    disp(['READING mesh from: ',fname0]);
    %% read number of elements, number of nodes dimension
    % and co-homology
    nums=fscanf(fid,"%d",[1 4]);
    nt=nums(1);nv=nums(2);dim=nums(3);nholes=nums(4);
    dim1=dim+1;
    disp(['...Reading element-vertex table ']);
        t=fscanf(fid,"%d",[1 dim1*nt]);
    t=reshape(t,nt,dim1);
    disp(['...Reading element codes ']);
        el_code=fscanf(fid,"%d",[1 nt]);
    el_code=el_code';
    disp(['...Reading point coordinates ']);
        x=fscanf(fid,"%g",[1 nv*dim]);  
    x=reshape(x,nv,dim);
    disp(['...Reading boundary flags ']);
        ib=fscanf(fid,"%d",[1 nv]);
    ib=ib';
    status=fclose(fid);
else
    disp(['error opening file **',fnam0,' **'])
    status=-1
end
return
