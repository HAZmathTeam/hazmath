function [nv,nt,tnew,tv,iw]=clean_sc(t,dim)
    %%    cleans up the table element-vertex by ?numbering vertices consecutively.
    dim1=dim+1;
    nv=max(max(t));
    [nt,idummy]=size(t);
    it=[1:nt];
    for k=1:dim, it=[it;[1:nt]]; end
    it=reshape(it,dim1*nt,1);
    t=reshape(transpose(t(1:nt,1:dim1)),dim1*nt,1);
    tv=sparse(it,t,ones(dim1*nt,1));
    iw=find(ones(1,nt)*tv);
    tv=tv(1:nt,iw);    
    [it,t]=find(tv); %% boundary faces
    [~,indx0]=sort(it);
    tnew=transpose(reshape(t(indx0),dim1,nt));
    [nt,nv]=size(tv);
end