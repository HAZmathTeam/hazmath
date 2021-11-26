function [iblk,jblk,nblk,errflag]=dfs0(A,n)
%
%    C...  This subroutine performs an algorithm due R.E.Tarjan 
%      for finding the
%    C...  strongly connected c omponents of a diggraph.  This is a
%    C...  modification of an implementation due  FRED GUSTAVSON  of the
%    C...  algorithm.  There is no warranty that such version would work
%    C...  properly.  The new numbering is in CSR format the pointers are in
%    C...  IBLK, and the actual numbering is in JBLK.
%    ---------------------------------------------------------------------*/
%  /*C...  Initialization.*/
    iblk=[];
    jblk=[];
    errflag=0;
    n1 = n + 1;
    n2 = n1 + 1;
    nblk = 0;
    nb = 0;
    count  = nb+1;
    k1 = count;
    vp = n1;
    sp = n1;
    for i = 1:n
        iedge(i) = 1;
        numb(i)  = 0;
        lowlink(i) = 0;
    end
    numb(n1) = 0;
    myflag=10;
    while (1)  %%10 continue
        %%    fprintf(stdout," myflag=%i %i, %i %i\n",myflag,nblk,count,n);
        if(myflag==10)%
                      %    C...  Get out when the  renumbering is done;
                      %    C...  Otherwise begin at vertex K1 
            if(count==n1) %
                iblk(nblk+1) = n1;
                iblk=iblk(1:nblk+1);
                return;
            end
            for ii=k1:n %
                i=ii;
                if(numb(ii)==0) %
                    myflag=30;
                    break; %%go to 30
                end
            end
            if(myflag ~=30)%
                disp(['There is an error in DEPTH FIRST SEARCH: i=  ', ...
                      int2str(i),'  n=',int2str(n)])
                errflag=254;
                return
            else %
                v = i;
                k1 = v + 1;
                myflag=50;
            end
        end
        if(myflag==40) 
            vp=vp-1;
            iblk(vp) = v;
            v=w;
            myflag=50;
        end
        if(myflag==50) 
            nb=nb + 1;
            numb(v) = nb;
            lowlink(v) = numb(v);
            sp=sp-1;
            jblk(sp) = v;
            v1 = v + 1;
            myflag=60;
        end
        %    /*...  */
        while(60) %% 60   continue;
            if(myflag == 60) %
                wp = iedge(v);
                %%% grab the indices corresponding to this row:
                ixxx=find(A(v,:));
                w = ixxx(wp);
                iedge(v) = wp+1;
                if(numb(w) >= numb(v)) 
                    myflag=70;
                elseif (numb(w)==0) 
                    myflag=40; 
                    break;
                else 
                    if(lowlink(v)>numb(w))
                        lowlink(v)=numb(w);
                    end
                    myflag=70;
                end
            end
            ixxx=find(A(v,:));
            lxxx=length(ixxx)+1;
            if(iedge(v) < lxxx) 
                myflag=60;
                continue;
            end
            %      /*...*/  
            if(lowlink(v) >= numb(v)) %%do not go to 90
                nblk = nblk + 1;
                iblk(nblk) = count;
                %	/**/
                while(80) %     %%    80 continue;
                    w = jblk(sp);
                    numb(w) = n1;
                    sp = sp + 1;
                    jblk(count) = w;
                    %%   count = count + 1
                    count=count+1;
                    if(v == w) break; end %%do not go to 80
                end
                if(sp==n1) 
                    myflag=10; 
                    break;
                end
            end 
            myflag=70;
            w = v;
            v = iblk(vp);
            vp = vp + 1;
            v1 = v + 1;
            if(lowlink(v)>lowlink(w)) lowlink(v)=lowlink(w); end
        end
    end
end
