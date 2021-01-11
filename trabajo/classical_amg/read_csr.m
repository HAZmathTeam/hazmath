function [A,status]=read_csr(file_csr)
    fid=fopen(file_csr); %% this is in csr format. 
    if(fid)
        nums=fscanf(fid,"%d",[1 3]);
        n=nums(1);m=nums(2);nz=nums(3);n1=n+1;
        ia=fscanf(fid,"%d",[1 n1]);
        if(ia(1)==0) ish=1; end
        ia=ia+ish; 
        ii=zeros(1,nz);
        for i=1:n
            iaa=ia(i);iab=ia(i+1)-1;
            ii(iaa:iab)=i*ones(1,iab-iaa+1);
        end
        clear ia;
        jj=fscanf(fid,"%d",[1 nz]);
        jj=jj+ish;
        s=fscanf(fid,"%e",[1 nz+n]);
        if(!length(s))
            s=-ones(nz,1);
        else
            s=s(1:nz);%% the rest of the reals are rhs.
        end
        A=sparse(ii,jj,s,n,n);
        %%    clear ii;clear jj;clear s;
        ib=fscanf(fid,"%d",[1 n]);
        ib=setdiff([1:n],find(ib));
        %%remove boundary conditions and make graph laplacian;
        A=A(ib,ib);n=size(A,1);m=size(A,2);nz=nnz(A);
        A=A-spdiags(diag(A),[0],n,m);
        A=A-spdiags([sum(A)]',[0],n,m);
        status=fclose(fid);
        %%%%%%%%%%%%%%%%%%%%%%%%%%end read and conversion
    else
        disp(['could not open CSR file: ',file_csr,...
              ' continuing with 1d graph'])
        n=127;nz=2*n-2+n;
        o=ones(n,1);
        A=spdiags([-o,2*o,-o],[-1:1],n,n);
        A(1,1)=1; A(n,n)=1;
        status=1;
    end

end

