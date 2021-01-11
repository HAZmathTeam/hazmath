function [status]=write_csr(A,fid)
    %% writes a matrix in a csr format 
    [ir,ic,val]=find(A);
    [ir,ip]=sort(ir);
    ic=ic(ip);val=val(ip);
    m=size(A,1);n=size(A,2);nz=length(ir);
    fprintf(fid,'%d ',[m,n,nz]);fprintf(fid,'\n');
    ia(1)=0;
    for k=1:m
        ia(k+1)=ia(k)+length(find((ir==k)));
    end
    fprintf(fid,'%d ',ia); fprintf(fid,'\n');
    fprintf(fid,'%d ',ic-1);  fprintf(fid,'\n');
    fprintf(fid,'%d ',val);  fprintf(fid,'\n');
    status=0;
end

