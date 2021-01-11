function [nfi,nfb,tt,ft,fv]=tt_ft_fv(t,tv,nv,nt,dim)
    %% in a simplicial mesh find
    %% the incidence matrices vertex_simplex;
    %% simplex_simptex (tt), face-vertex, edge-vertex    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% INPUT:  tv=incidence(dim,dim)
    %%         t is tv as a table
    %% OUTPUT: tt=incidence(dim,dim)
    %%         ft=incidence(dim-1,dim)
    %%         fv=incidence(dim-1,0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dim1=dim+1;
    %% incidence simplex -- simplex
    tt=sparse(tv*tv');
    tt=tt.*(tt==dim);
    [i,j,s]=find(tril(tt,-1));
    %k2=find(i<j);ne=length(k2);
    nfi=length(i);
    nfb=(dim1*nt) - 2*nfi;
    nf=nfi + nfb;
    ift=[[1:nfi]';[1:nfi]'];
    ft=sparse(ift,[i;j],[ones(nfi,1);-ones(nfi,1)],nf,nt);
    fv=abs(ft)*tv;
    %% this is dim independent, every vertex
    %% that was counted in both elements sharing the face.
    %% it finds all interior faces. 
    fv=fv.*(fv==2);     fv=(fv>0);
    it=[1:nt];
    for k=1:(dim-1), it=[it;[1:nt]]; end
    it=reshape(it,dim*nt,1);
    nb=nfi+1;
    for k=1:dim1
        indx0=setdiff([1:dim1],k);
        ttmp=t(:,indx0);
        ttmp=reshape(transpose(ttmp(1:nt,1:dim)),dim*nt,1);
        tv1=sparse(it,ttmp,ones(dim*nt,1),nt,nv);
        %%%%%%%%%%%%%%%%%%        if(dim==2),pause, end
        fb=tv1*[fv(1:nfi,:)]';
        fb=fb.*(fb==dim);
        indx0=find((fb*ones(nfi,1))==0);
        n0=length(indx0);
        if(n0==0) , continue, end
        %        fv(nb:nb+n0,:)=fb(indx0,:)
        fv(nb:nb+n0-1,:)=tv1(indx0,:);
        ft(nb:nb+n0-1,indx0)=speye(n0);
        nb = nb + n0;
    end
end