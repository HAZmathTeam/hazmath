c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine mgilu(ja,a,lvl,ka)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd), dimension(10,*) :: ka
            real(kind=rknd), dimension(*) :: a
            real(kind=rknd), dimension(2) :: coeff
            integer(kind=iknd), allocatable, dimension(:) :: vtype
            real(kind=rknd), allocatable, dimension(:) :: dv,dw
cy
        n=ka(1,1)
        ispd=ka(1,lvl+1)
        allocate(dv(n),dw(n),vtype(n))
c
        a(n+1)=canorm(n,ispd,ja,a)
c
c       compute coarse graph matrices
c
        do level=lvl,1,-1
            call getptr(level,lvl,nf,nptr,japtr,iaptr,
     +          juptr,iuptr,jvptr,ivptr,iqptr,ibptr,nc,ncptr,ka)
            if(iuptr/=iaptr) then
                call snfilu(nf,ja(japtr),a(iaptr),ja(juptr),
     +              a(iuptr),ispd)
            endif
            if(level>1) then
                call getptr(level-1,lvl,nc,ncptr,jacptr,iacptr,
     +              jucptr,iucptr,jvcptr,ivcptr,iqcptr,ibcptr,
     1              ncc,nccptr,ka)
                call ceig(nf,ispd,ja(japtr),a(iaptr),ja(juptr),
     +              a(iuptr),dv,dw,0_iknd,vtype,coeff)
                call cwt(nf,nc,ispd,ja(jvptr),a(ivptr),vtype,
     +              ja(japtr),a(iaptr),ja(iqcptr),dv,dw)
                call a2ac(nf,ispd,ja(japtr),a(iaptr),nc,
     +              ja(jacptr),a(iacptr),ja(jvptr),a(ivptr))
            endif
        enddo
c
        deallocate(dv,dw,vtype)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine mginit(n,ispd,nblock,ib,maxja,ja,maxa,a,ncfact,
     +      maxlvl,maxfil,ka,lvl,dtol,method,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ib
            integer(kind=iknd), dimension(10,*) :: ka
            real(kind=rknd), dimension(*) :: a
            real(kind=rknd), dimension(2) :: coeff
            integer(kind=iknd), dimension(n) :: vtype,qc,p,q
            integer(kind=iknd), dimension(maxja) :: jc
            real(kind=rknd), dimension(n) :: dv,dw
cy
        iflag=0
c
        call chkja(n,ja,iflag)
        if(iflag/=0) return
        call chkib(n,nblock,ib,iflag)
        if(iflag/=0) return
c
        maxlvl=max(1,maxlvl)
        if(dtol<=0.0e0_rknd) maxlvl=1
        dtol=abs(dtol)
        ncfact=max(ncfact,2)
        if(method<0.or.method>2) method=0
c
        lenja=ja(n+1)-1
        if(ispd==1) then
            lena=lenja
        else
            lena=2*lenja-(n+1)
        endif
        minfil=((lenja-(n+1))/n)+3
c
        nf=n
        nptr=1
        japtr=1
        iaptr=1
        ibptr=japtr+lenja
        iqptr=ibptr+nf
        lvl=1
c
        ka(1,lvl)=nf
        ka(2,lvl)=nptr
        ka(3,lvl)=japtr
        ka(4,lvl)=iaptr
        ka(5,lvl)=0
        ka(6,lvl)=0
        ka(7,lvl)=0
        ka(8,lvl)=0
        ka(9,lvl)=iqptr
        ka(10,lvl)=ibptr
        a(n+1)=canorm(n,ispd,ja,a)
        do i=1,nblock
           do j=ib(i),ib(i+1)-1
              ja(lenja+j)=i
           enddo
        enddo
c
c
c
   10   nf=ka(1,lvl)
        nptr=ka(2,lvl)
        japtr=ka(3,lvl)
        iaptr=ka(4,lvl)
        iqptr=ka(9,lvl)
        ibptr=ka(10,lvl)
c
        if(iqptr+nf>maxja) then
            lvl=lvl-1
            if(lvl<1) iflag=20
            go to 20
        endif
c
c       compute ordering vector q, reorder ja and a
c
        rdtol=dtol
cc      if(lvl==maxlvl.and.maxlvl>1) rdtol=0.0e0_rknd
        call ja2jcb(nf,ispd,ja(japtr),a(iaptr),rdtol,maxja,jc,
     +      vtype,lena0,iflag)
        if(iflag/=0) then
            lvl=lvl-1
            if(lvl>=1) iflag=0
            go to 20
        endif
        lenja=ja(japtr+nf)-1
        call md(nf,lenja,ja(japtr),ja(iqptr),lenu0)
        call ja2jc1(nf,vtype,ja(iqptr))
        call ja2ja(nf,lenja,ja(japtr),ja(iqptr),ispd,a(iaptr),
     +      ja(ibptr))
        if(lvl>1) call vf2vf(ka(1,lvl-1),nf,ja(jvptr),ja(iqptr),qc)
c
        if(method==0) then
            juptr=iqptr+nf
            iuptr=iaptr+lena
            ka(5,lvl)=juptr
            ka(6,lvl)=iuptr
            lenju=maxja-juptr+1
            lenu=maxa-iuptr+1
            rdtol=dtol
cc          if(lvl==maxlvl.and.maxlvl>1) rdtol=0.0e0_rknd
            call sfilu(nf,ja(japtr),a(iaptr),lenju,ja(juptr),
     +          lenu,a(iuptr),ispd,rdtol,maxfil,iflag)
            if(iflag/=0) then
                lvl=lvl-1
                if(lvl>=1) iflag=0
                go to 20
            endif
c
            ka(7,lvl)=juptr+lenju
            ka(8,lvl)=iuptr+lenu
c
        else if(method==1) then
            juptr=japtr
            iuptr=iaptr+lena
            ka(5,lvl)=juptr
            ka(6,lvl)=iuptr
            ka(7,lvl)=iqptr+nf
            ka(8,lvl)=iuptr+lena
            if(iuptr+lena>maxa) then
                lvl=lvl-1
                if(lvl>=1) iflag=0
                go to 20
            endif
            call snfilu(nf,ja(japtr),a(iaptr),ja(juptr),a(iuptr),ispd)
c
        else if(method==2) then
            juptr=japtr
            iuptr=iaptr
            ka(5,lvl)=juptr
            ka(6,lvl)=iuptr
            ka(7,lvl)=iqptr+nf
            ka(8,lvl)=iuptr+lena
        endif
c
        ka(1,lvl+1)=0
        ka(2,lvl+1)=ka(2,lvl)+nf
        ka(3,lvl+1)=0
        ka(4,lvl+1)=0
        ka(5,lvl+1)=0
        ka(6,lvl+1)=0
        ka(7,lvl+1)=0
        ka(8,lvl+1)=0
        ka(9,lvl+1)=0
        ka(10,lvl+1)=0
        if(lvl>=maxlvl.or.nf<=1.or.lenu0<=lena0) go to 20
c
        jvptr=ka(7,lvl)
        ivptr=ka(8,lvl)
        lenjv=maxja-jvptr+1
        lenwt=maxa-ivptr+1
cc      call ja2jc(nf,ja(japtr),jc)
        rdtol=min(1.e-3_rknd,dtol)
        call ja2jf(nf,ispd,ja(japtr),a(iaptr),rdtol,jc,ja(ibptr))
        call crsncm(nf,nc,p,q,jc,maxja,vtype,qc,ncfact,iflag)
        call crsncr(nf,nc,p,q,jc,vtype,qc,ncfact,ispd,ja(japtr),
     +      a(iaptr),ja(juptr),a(iuptr),nblock,ja(ibptr),iflag)
        if(iflag/=0) then
            if(iflag==1) iflag=0
            go to 20
        endif
        call cvf(nf,nc,ispd,jc,lenjv,ja(jvptr),lenwt,
     +      qc,iflag)
        if(iflag/=0) then
            if(iflag==20) iflag=0
            go to 20
        endif
        if(nc<=0) go to 20
        call ceig(nf,ispd,ja(japtr),a(iaptr),ja(juptr),a(iuptr),
     +      dv,dw,0_iknd,vtype,coeff)
        call cwt(nf,nc,ispd,ja(jvptr),a(ivptr),vtype,
     +      ja(japtr),a(iaptr),qc,dv,dw)
        ka(1,lvl+1)=nc
        jacptr=ka(7,lvl)+lenjv
        iacptr=ka(8,lvl)+lenwt
        ka(3,lvl+1)=jacptr
        ka(4,lvl+1)=iacptr
c
        lenja=maxja-jacptr+1
        call ja2jac(nf,ja(japtr),jc,nc,ja(jacptr),ja(jvptr),
     +      lenja,iflag)
        if(iflag/=0) then
            if(iflag==20) iflag=0
            go to 20
        endif
        if(ispd==1) then
             lena=lenja
        else
             lena=2*lenja-(nc+1)
        endif
        if(iacptr+lena-1>maxa) then
            go to 20
        endif
        call a2ac(nf,ispd,ja(japtr),a(iaptr),nc,ja(jacptr),
     +      a(iacptr),ja(jvptr),a(ivptr))
        rdtol=min(1.e-3_rknd,dtol)
        call sfac(nc,ja(jacptr),a(iacptr),lenja,lena,ispd,
     +      rdtol,maxfil,minfil)
        ibcptr=jacptr+lenja
        iqcptr=ibcptr+nc
        ka(9,lvl+1)=iqcptr
        ka(10,lvl+1)=ibcptr
        if(iqcptr+nc>maxja) go to 20
        call mkib(nc,qc,ja(ibptr),ja(ibcptr))
c
        lvl=lvl+1
        go to 10
c
   20   ka(1,lvl+1)=ispd
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        function canorm(n,ispd,ja,a)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            real(kind=rknd), dimension(*) :: a
            real(kind=rknd) :: canorm
cy
c       compute anorm
c
        canorm=0.0e0_rknd
        eps=epsilon(1.0e0_rknd)
c
        do i=1,n
            canorm=max(canorm,abs(a(i)))
        enddo
        canorm=canorm*eps
        if(canorm>0.0e0_rknd) return
c
c       if diag is zero, try off diagonals
c
        nnz=ja(n+1)-ja(1)
        if(ispd/=1) nnz=2*nnz
        do i=1,nnz
            canorm=max(canorm,abs(a(ja(1)+i-1)))
        enddo
        canorm=canorm*eps
        if(canorm>0.0e0_rknd) return
c
c       if the matrix is zero
c
        canorm=eps
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ja2ja(n,lenja,ja,q,ispd,a,bindx)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,q,bindx
            integer(kind=iknd), dimension(lenja) :: link
            integer(kind=iknd) :: amtx
            real(kind=rknd), dimension(*) :: a
cy
c       compute linked list of new column indices
c
        if(ispd/=1) then
            amtx=ja(n+1)-ja(1)
        else
            amtx=0
        endif
c
        do i=1,n
            link(i)=0
        enddo
        do ii=1,n
            i=q(ii)
            do jj=ja(ii),ja(ii+1)-1
                j=q(ja(jj))
                if(i>j) then
                    irow=j
                    icol=i
                    aa=a(jj)
                    a(jj)=a(jj+amtx)
                    a(jj+amtx)=aa
                else
                    irow=i
                    icol=j
                endif
                ja(jj)=icol
                last=irow
   10           next=link(last)
                if(next==0) then
                    link(last)=jj
                    link(jj)=0
                else
                    if(icol<ja(next)) then
                        link(last)=jj
                        link(jj)=next
                    else
                        last=next
                        go to 10
                    endif
                endif
            enddo
        enddo
c
        ja(1)=n+2
        do i=1,n
            len=ja(i)
            last=i
            next=link(last)
   20       if(next>0) then
                last=next
                next=link(last)
                link(last)=len
                len=len+1
                go to 20
            endif
            ja(i+1)=len
            link(i)=q(i)
        enddo
c
c       reorder upper triangle
c
        do i=ja(1),ja(n+1)-1
   30       if(link(i)/=i) then
                ii=link(i)
                link(i)=link(ii)
                link(ii)=ii
                jj=ja(i)
                ja(i)=ja(ii)
                ja(ii)=jj
                a1=a(i)
                a2=a(i+amtx)
                a(i)=a(ii)
                a(i+amtx)=a(ii+amtx)
                a(ii)=a1
                a(ii+amtx)=a2
                go to 30
            endif
        enddo
c
c       diagonal of a
c
        do i=1,n
   40       if(link(i)/=i) then
                ii=link(i)
                link(i)=link(ii)
                link(ii)=ii
                kk=bindx(i)
                bindx(i)=bindx(ii)
                bindx(ii)=kk
                aa=a(i)
                a(i)=a(ii)
                a(ii)=aa
                go to 40
            endif
        enddo
c
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ihp(list,len)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: list
cy
c       reorder entries in list small to large
c
        if(len<=1) return
        n=len/2
        do m=n,1,-1
            k=m
            do
                kson=2*k
                if(kson>len) exit
                if(kson<len) then
                    if(list(kson)<list(kson+1)) kson=kson+1
                endif
                if(list(k)>=list(kson)) exit
                itemp=list(k)
                list(k)=list(kson)
                list(kson)=itemp
                k=kson
            enddo
        enddo
c
c
        do n=len,2,-1
            itemp=list(1)
            list(1)=list(n)
            list(n)=itemp
            k=1
            do
                kson=2*k
                if(kson>n-1) exit
                if(kson<n-1) then
                    if(list(kson)<list(kson+1)) kson=kson+1
                endif
                if(list(k)>=list(kson)) exit
                itemp=list(k)
                list(k)=list(kson)
                list(kson)=itemp
                k=kson
            enddo
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ja2jc(n,ja,jc)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,jc
cy
c       make jc data structure from ja data structure
c
        do i=1,n
            jc(i+1)=ja(i+1)-ja(i)
        enddo
c
c       compute new lengths
c
        do i=ja(1),ja(n+1)-1
            k=ja(i)+1
            jc(k)=jc(k)+1
        enddo
c
        jc(1)=n+2
        do i=2,n+1
            jc(i)=jc(i)+jc(i-1)
        enddo
c
        do i=1,n
            do jj=ja(i),ja(i+1)-1
                j=ja(jj)
                jc(jc(i))=j
                jc(i)=jc(i)+1
                jc(jc(j))=i
                jc(j)=jc(j)+1
            enddo
        enddo
c
        do i=n+1,2,-1
            jc(i)=jc(i-1)
        enddo
        jc(1)=n+2
c
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ja2jcb(n,ispd,ja,a,dtol,maxjc,jc,mark,lenjc,iflag)
c
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,jc,mark
            integer(kind=iknd), dimension(n) :: list
            real(kind=rknd), dimension(*) :: a
cy
c       make jc data structure from ja data structure
c
        iflag=0
        lenjc=0
        if(ispd==1) then
            ishift=0
        else
            ishift=ja(n+1)-ja(1)
        endif
        do i=1,n
            list(i)=0
            mark(i)=0
            jc(i+1)=0
        enddo
        anorm=a(n+1)
c
c       compute new lengths
c
        do i=1,n
            do 30 jj=ja(i),ja(i+1)-1
                j=ja(jj)
                ee=eltest(a(i),a(jj),a(jj+ishift),a(j))
                if(ee<dtol) go to 30
                jc(i+1)=jc(i+1)+1
                jc(j+1)=jc(j+1)+1
   30       enddo
        enddo
c
        jc(1)=n+2
        do i=2,n+1
            jc(i)=jc(i)+jc(i-1)
        enddo
c
        do i=1,n
            do 20 jj=ja(i),ja(i+1)-1
                j=ja(jj)
                ee=eltest(a(i),a(jj),a(jj+ishift),a(j))
                if(ee<dtol) go to 20
                jc(jc(i))=j
                jc(i)=jc(i)+1
                jc(jc(j))=i
                jc(j)=jc(j)+1
   20       enddo
        enddo
c
        do i=n+1,2,-1
            jc(i)=jc(i-1)
        enddo
        jc(1)=n+2
c
        num=0
        do 10 i=1,n
c
            if(mark(i)/=0) go to 10
            if(abs(a(i))>=anorm) go to 10
            if(jc(i)==jc(i+1)) go to 10
c
c       scan for zero diag entries, increase their degree relative to
c       largest connected entry with non-zero diag
c
            jx=jc(jc(i))
            ax=0.0e0_rknd
            do k=jc(i),jc(i+1)-1
               j=jc(k)
               if(a(j)/=0.0e0_rknd.and.mark(j)==0) then
                   call jamap(i,j,ij,ji,ja,ishift)
                   aa=abs(a(ij)*a(ji)/a(j))
                   if(aa>ax) then
                       jx=j
                       ax=aa
                   endif
               endif
            enddo
            if(ax>abs(a(i))) then
                mark(i)=jx
                mark(jx)=-i
                num=num+1
            endif
   10   continue
        lenjc=(jc(n+1)-jc(1))/2+n+1
        if(num==0) return
c
c       merge rows
c
        next=maxjc+1
        do i=n,1,-1
            i1=jc(i)
            i2=jc(i+1)-1
            jc(i+1)=next
c
c       shift row i
c
            list(i)=i
            len=0
            do k=i1,i2
                j=jc(k)
                if(list(j)==0) then
                    len=len+1
                    list(j)=list(i)
                    list(i)=j
                endif
                if(mark(j)<0) then
                    mj=-mark(j)
                    if(list(mj)==0) then
                        len=len+1
                        list(mj)=list(i)
                        list(i)=mj
                    endif
                endif
            enddo
            if(mark(i)>0) then
                mi=mark(i)
                do k=jc(mi),jc(mi+1)-1
                    j=jc(k)
                    if(list(j)==0) then
                        len=len+1
                        list(j)=list(i)
                        list(i)=j
                    endif
                    if(mark(j)<0) then
                        mj=-mark(j)
                        if(list(mj)==0) then
                            len=len+1
                            list(mj)=list(i)
                            list(i)=mj
                        endif
                    endif
                enddo
            endif
            if(next-len<i1) then
                iflag=20
                return
            endif
            do j=1,len
                next=next-1
                jc(next)=list(i)
                list(i)=list(jc(next))
                list(jc(next))=0
            enddo
            list(i)=0
        enddo
        jc(1)=next
        len=jc(n+1)-jc(1)
        ishift=next-(n+2)
        do i=1,len
            jc(n+1+i)=jc(next+i-1)
        enddo
        do i=1,n+1
            jc(i)=jc(i)-next+(n+2)
        enddo
c
        lenjc=(jc(n+1)-jc(1))/2+n+1
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ja2jf(n,ispd,ja,a,dtol,jf,bindx)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,jf,bindx
            real(kind=rknd), dimension(*) :: a
cy
c       make blockdiagonal jf data structure from ja data structure
c
        if(ispd==1) then
            ishift=0
        else
            ishift=ja(n+1)-ja(1)
        endif
        do i=1,n
            jf(i+1)=0
        enddo
c
c       compute new lengths
c
        do i=1,n
            do 10 jj=ja(i),ja(i+1)-1
                j=ja(jj)
                if(bindx(i)/=bindx(j)) go to 10
                ee=eltest(a(i),a(jj),a(jj+ishift),a(j))
                if(ee<dtol) go to 10
                jf(i+1)=jf(i+1)+1
                jf(j+1)=jf(j+1)+1
   10       enddo
        enddo
c
        jf(1)=n+2
        do i=2,n+1
            jf(i)=jf(i)+jf(i-1)
        enddo
c
        do i=1,n
            do 20 jj=ja(i),ja(i+1)-1
                j=ja(jj)
                if(bindx(i)/=bindx(j)) go to 20
                ee=eltest(a(i),a(jj),a(jj+ishift),a(j))
                if(ee<dtol) go to 20
                jf(jf(i))=j
                jf(i)=jf(i)+1
                jf(jf(j))=i
                jf(j)=jf(j)+1
   20       enddo
        enddo
c
        do i=n+1,2,-1
            jf(i)=jf(i-1)
        enddo
        jf(1)=n+2
c
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ja2jc1(n,mark,q)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: mark,q
cy
        do i=1,n
            if(mark(i)>0) then
                j=mark(i)
                if(q(i)<q(j)) then
                    k=q(i)
                    q(i)=q(j)
                    q(j)=k
                endif
            endif
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine vf2vf(nf,nc,vf,q,qc)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: vf,q,qc
            integer(kind=iknd), dimension(nc) :: list
cy
c
        do i=1,nf
            do j=vf(i),vf(i+1)-1
                vf(j)=q(vf(j))
            enddo
        enddo
c
c       save the fine to coarse mapping in q (md ordering not needed)
c
        do i=1,nc
            list(q(i))=qc(i)
        enddo
        do i=1,nc
            q(i)=list(i)
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine vf2vfc(n,nc,vf,vfc)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: vf,vfc
cy
c       make linked list of entries
c
        do i=1,nc+1
           vfc(i)=0
        enddo
        do i=1,n
            do jj=vf(i),vf(i+1)-1
                j=vf(jj)
                vfc(j+1)=vfc(j+1)+1
            enddo
        enddo
        vfc(1)=nc+2
        do i=2,nc+1
            vfc(i)=vfc(i)+vfc(i-1)
        enddo
c
        do i=1,n
            do jj=vf(i),vf(i+1)-1
                j=vf(jj)
                k=vfc(j)
                vfc(j)=k+1
                vfc(k)=i
            enddo
        enddo
        do i=nc+1,2,-1
            vfc(i)=vfc(i-1)
        enddo
        vfc(1)=nc+2
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ja2jac(n,ja,jc,nc,jac,vf,maxjac,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,jac,vf,jc
            integer(kind=iknd), dimension(n) :: list,mark
            integer(kind=iknd), dimension(maxjac) :: vfc
cy
c       make linked list of entries
c
        if(nc+1>maxjac) then
            iflag=20
            return
        endif
        iflag=0
        call ja2jc(n,ja,jc)
        call vf2vfc(n,nc,vf,vfc)
        do i=1,n
            mark(i)=0
        enddo
        jac(1)=nc+2
        do i=1,nc
            len=0
            do jj=vfc(i),vfc(i+1)-1
                j=vfc(jj)
                if(mark(j)/=-i) then
                     len=len+1
                     mark(j)=-i
                     list(len)=j
                endif
                do kk=jc(j),jc(j+1)-1
                    k=jc(kk)
                    if(mark(k)/=-i) then
                         len=len+1
                         mark(k)=-i
                         list(len)=k
                    endif
                enddo
            enddo
            next=jac(i)
            do jj=1,len
                j=list(jj)
                do kk=vf(j),vf(j+1)-1
                    k=vf(kk)
                    if(mark(k)/=i) then
                        if(next>maxjac) then
                            iflag=20
                            return
                        endif
                        jac(next)=k
                        next=next+1
                        mark(k)=i
                    endif
                enddo
            enddo
            jac(i+1)=next
            len=jac(i+1)-jac(i)
            if(len>1) call ihp(jac(jac(i)),len)
        enddo
        maxjac=jac(nc+1)-1
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine mkib(nc,qc,bindx,bindxc)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: qc,bindx,bindxc
cy
c       update ib data structure for coarse graph
c
        do i=1,nc
           bindxc(i)=bindx(qc(i))
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine a2ac(n,ispd,ja,a,nc,jac,ac,vf,wt)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,jac,vf
            integer(kind=iknd) :: amtx,acmtx,vmtx,wmtx
            real(kind=rknd), dimension(*) :: a,ac,wt
cy
c
        if(ispd==1) then
            wmtx=1-vf(1)
            vmtx=wmtx
            amtx=0
            acmtx=0
        else
            wmtx=1-vf(1)
            vmtx=wmtx+vf(n+1)-vf(1)
            amtx=ja(n+1)-ja(1)
            acmtx=jac(nc+1)-jac(1)
        endif
c
c       initialize
c
        do i=1,jac(nc+1)-1+acmtx
            ac(i)=0.0e0_rknd
        enddo
c
c       the main loop
c
        do i=1,n
c
c       diagonal entry
c
            aii=a(i)
            do kk=vf(i),vf(i+1)-1
                k=vf(kk)
                wtik=wt(kk+wmtx)
                wtki=wt(kk+vmtx)
                ac(k)=ac(k)+wtki*aii*wtik
                do mm=kk+1,vf(i+1)-1
                    m=vf(mm)
                    wtim=wt(mm+wmtx)
                    wtmi=wt(mm+vmtx)
                    call jacmap(m,k,mk,km,jac,acmtx)
                    if(mk>0) then
                        aa=ac(mk)+wtmi*aii*wtik
                        ac(km)=ac(km)+wtki*aii*wtim
                        ac(mk)=aa
                    endif
                enddo
             enddo
c
c       off diagonal entries
c
            do j=ja(i),ja(i+1)-1
                aij=a(j)
                aji=a(j+amtx)
                do kk=vf(ja(j)),vf(ja(j)+1)-1
                    k=vf(kk)
                    wtjk=wt(kk+wmtx)
                    wtkj=wt(kk+vmtx)
                    do mm=vf(i),vf(i+1)-1
                        m=vf(mm)
                        wtim=wt(mm+wmtx)
                        wtmi=wt(mm+vmtx)
                        if(k==m) then
                            ac(k)=ac(k)+wtkj*aji*wtim+wtmi*aij*wtjk
                        else
                            call jacmap(m,k,mk,km,jac,acmtx)
                            if(mk>0) then
                                aa=ac(mk)+wtmi*aij*wtjk
                                ac(km)=ac(km)+wtkj*aji*wtim
                                ac(mk)=aa
                            endif
                        endif
                    enddo
                 enddo
            enddo
        enddo
c
c       compute anorm
c
        ac(nc+1)=canorm(nc,ispd,jac,ac)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine jacmap(i,j,ij,ji,ja,amtx)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd) :: amtx
cy
c       compute location of a(i,j) and a(j,i)
c
        if(i<j) then
            imin=ja(i)
            imax=ja(i+1)-1
   10       imid=(imin+imax)/2
            if(ja(imid)==j) then
                ij=imid
                ji=ij+amtx
                return
            else if(imid==imax) then
                ij=0
                ji=0
                return
            else if(ja(imid)<j) then
                if(imid==imin) imid=imax
                imin=imid
                go to 10
            else
                imax=imid
                go to 10
            endif
c
        else
            jmin=ja(j)
            jmax=ja(j+1)-1
   20       jmid=(jmin+jmax)/2
            if(ja(jmid)==i) then
                ji=jmid
                ij=ji+amtx
                return
            else if(jmid==jmax) then
                ij=0
                ji=0
                return
            else if(ja(jmid)<i) then
                if(jmid==jmin) jmid=jmax
                jmin=jmid
                go to 20
            else
                jmax=jmid
                go to 20
            endif
        endif
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine sfac(n,ja,a,lenja,lena,ispd,dtol,maxfil,minfil)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd), dimension(400) :: ibin
            integer(kind=iknd) :: amtx
            real(kind=rknd), dimension(*) :: a
cy
c       sparse numeric factorization
c
        nbin=400
        if(ispd/=1) then
            amtx=ja(n+1)-ja(1)
        else
            amtx=0
        endif
        fact=1.0e4**(1.0e0_rknd/real(nbin))
        alf=log(fact)
        eps=epsilon(1.0e0_rknd)
        rtol=max(eps,dtol)
        lenja0=ja(n+1)-1
        qn=real(max(0,maxfil))*real(n)+real(n+1)
        if(qn>=real(lenja0)) go to 30
c
c       figure out drop tolerance in case of violation of maxfil
c
        jatrgt=int(qn+0.5e0_rknd)
        do i=1,nbin
            ibin(i)=0
        enddo
        do i=1,n
            do j=ja(i),ja(i+1)-1
                tt=eltest(a(i),a(j),a(j+amtx),a(ja(j)))/rtol
                if(tt>=1.0e0_rknd) then
                    it=min(nbin,1+int(log(tt)/alf))
                    ibin(it)=ibin(it)+1
                endif
            enddo
        enddo
        ibin(nbin)=ibin(nbin)+n+1
        do i=nbin-1,1,-1
            ibin(i)=ibin(i+1)+ibin(i)
        enddo
        if(ibin(1)<jatrgt) go to 50
        do i=1,nbin-1
            if(ibin(i+1)<=jatrgt.and.ibin(i)>jatrgt) go to 20
        enddo
        i=nbin
   20   rtol=rtol*fact**i
        write(6,*) 'sfac0',n,i,rtol
        go to 50
c
c       figure out drop tolerance in case of violation of minfil
c
   30   jatrgt=minfil*n+n+1
        if(jatrgt<lenja0) go to 50
        kount=n+1
        stol=rtol*1.0e-2_rknd
        do i=1,nbin
            ibin(i)=0
        enddo
        do i=1,n
            do j=ja(i),ja(i+1)-1
                tx=eltest(a(i),a(j),a(j+amtx),a(ja(j)))
                if(tx/rtol>=1.0e0_rknd) kount=kount+1
                tt=tx/stol
                if(tt>=1.0e0_rknd) then
                    it=min(nbin,1+int(log(tt)/alf))
                    ibin(it)=ibin(it)+1
                endif
            enddo
        enddo
        if(kount>=jatrgt) go to 50
        ibin(nbin)=ibin(nbin)+n+1
        do i=nbin-1,1,-1
            ibin(i)=ibin(i+1)+ibin(i)
        enddo
        if(ibin(1)<jatrgt) go to 50
        do i=1,nbin-1
            if(ibin(i+1)<=jatrgt.and.ibin(i)>jatrgt) go to 40
        enddo
        i=nbin
   40   rtol=stol*fact**i
        write(6,*) 'sfac1',n,i,rtol
c
c       now do it for real
c
   50   jai=ja(1)
        do i=1,n
            next=ja(i)
            do j=jai,ja(i+1)-1
                tt=eltest(a(i),a(j),a(j+amtx),a(ja(j)))/rtol
                if(tt>=1.0e0_rknd) then
                    ja(next)=ja(j)
                    a(next)=a(j)
                    a(next+amtx)=a(j+amtx)
                    next=next+1
                endif
            enddo
            jai=ja(i+1)
            ja(i+1)=next
        enddo
        lenja=ja(n+1)-1
        if(ispd/=1) then
            nnz=lenja-(n+1)
            do i=1,nnz
                a(lenja+i)=a(lenja0+i)
            enddo
            lena=lenja+nnz
        else
            lena=lenja
        endif
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine crsncr(n,nc,p,q,jf,vtype,qc,ncfact,
     +      ispd,ja,a,ju,u,nblock,bindx,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: p,q,jf,vtype,qc,
     +          bindx,ja,ju
            real(kind=rknd), dimension(*) :: a,u
            real(kind=rknd), dimension(2) :: coeff,coef0
            integer(kind=iknd), dimension(n) :: vtype0
            real(kind=rknd), dimension(n) :: dv,dw,z0,z,az
cy
        iflag=0
        nswap=-1
        itmax=2
        thresh=0.1e0_rknd
        eps=1.0e-10_rknd
        bias=1.0e1
        itnum=0
        nctrgt=max(int(real(n)/real(ncfact)),1_iknd)
c
        call ceig(n,ispd,ja,a,ju,u,dv,dw,0_iknd,vtype0,coef0)
        do i=1,n
            z0(i)=abs(dv(i))+abs(dw(i))
            if(vtype(i)/=1) then
                vtype(i)=0
                if(jf(i)==jf(i+1)) vtype(i)=-1
            endif
        enddo
c
c
   10   itnum=itnum+1
        if(itnum>itmax.or.nswap==0) go to 30
        call ceig(n,ispd,ja,a,ju,u,dv,dw,nc,vtype,coeff)
c
        do i=1,nblock
            az(i)=0.0e0_rknd
        enddo
        cw=10.0e0_rknd**(coef0(1)-coeff(1))
        cv=10.0e0_rknd**(coef0(2)-coeff(2))
        do i=1,n
            z(i)=abs(dw(i))*cw+abs(dv(i))*cv
            az(bindx(i))=max(az(bindx(i)),z0(i))
        enddo
        do i=1,n
            if(vtype(i)==1) then
                ss=z0(i)
                do jj=jf(i),jf(i+1)-1
                    ss=max(ss,z0(jf(jj)))
                enddo
                z(i)=bias*max(ss/az(bindx(i)),eps)
            else if(vtype(i)==-1) then
                z(i)=0.0e0_rknd
            else
                z(i)=-z(i)/az(bindx(i))
            endif
c
c       p/q share space with dv/dw
c
            p(i)=i
            q(i)=i
        enddo
c
        nn=n/2
        do k=nn,1,-1
            call updhp(k,n,p,q,z,0_iknd)
        enddo
        do k=n,2,-1
            kk=p(1)
            p(1)=p(k)
            p(k)=kk
            q(p(1))=1
            q(p(k))=k
            call updhp(1_iknd,k-1_iknd,p,q,z,0_iknd)
        enddo
c
c       1 -- (nf) are fine grid ordered by decreasing size
c       (nf+1) -- n are coarse grid ordered by increasing size
c
c       too many coarse points
c
        nf=n-nc
        if(nc>nctrgt) then
            do ii=nf+1,n-nctrgt
                i=p(ii)
                vtype(i)=0
                nc=nc-1
            enddo
            go to 10
c
c       too few coarse points
c
        else if(nc<nctrgt) then
            do ii=nf+1,n
               i=p(ii)
               if(abs(z(i))>thresh) go to 15
               vtype(i)=0
               nc=nc-1
            enddo
   15       do ii=1,nctrgt-nc
               i=p(ii)
               if(abs(z(i))<=thresh) go to 10
               vtype(i)=1
               nc=nc+1
            enddo
        else
c
c       simple swaps
c
            nn=min(nf,nc)
            nswap=0
            do i=1,nn
                kc=p(nf+i)
                kf=p(i)
                if(z(kc)>=abs(z(kf))) go to 10
                if(vtype(kf)==-1) go to 10
                nswap=nswap+1
                vtype(kf)=1
                vtype(kc)=0
            enddo
        endif
c
   30   ncc=nc
        nc=0
        do i=1,n
            if(vtype(i)>0) then
                nc=nc+1
                qc(nc)=i
            endif
        enddo
        if(nc/=ncc) stop 9551
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine crsncm(n,nc,p,q,jf,maxjf,vtype,qc,ncfact,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: p,q,jf,vtype,qc
            integer(kind=iknd), dimension(n) :: jr
cy
        iflag=0
        jtmax=1
        jtnum=0
        nctrgt=max(int(real(n)/real(ncfact)),1_iknd)
        ncmax=(nctrgt*3)/2
c
   10   do i=1,n
            p(i)=i
            q(i)=i
            vtype(i)=0
            if(jf(i)==jf(i+1)) vtype(i)=-1
        enddo
c
c       sort out regions
c
        itmax=2
        nr=0
        iptr=1
        next=1
   20   k=p(next)
        if(next>=iptr) then
            nr=nr+1
            jr(nr)=next
            iptr=iptr+1
        endif
        next=next+1
        do j=jf(k),jf(k+1)-1
            m=q(jf(j))
            if(m>=iptr.and.m<=n) then
                p(m)=p(iptr)
                p(iptr)=jf(j)
                q(p(m))=m
                q(p(iptr))=iptr
                iptr=iptr+1
            endif
        enddo
        if(next<=n) go to 20
        if(nr>=n) then
            iflag=1
            return
        endif
        jr(nr+1)=n+1
c
c       order more or less lengthwise using rcm
c
        do itnum=1,itmax
            do j=1,nr
                i1=jr(j)
                i2=jr(j+1)-1
                iseed=p(i2)
                p(i2)=p(i1)
                p(i1)=iseed
                q(p(i1))=i1
                q(p(i2))=i2
            enddo
            iptr=1
            next=1
   30       k=p(next)
            if(next>=iptr) iptr=iptr+1
            next=next+1
            do j=jf(k),jf(k+1)-1
                m=q(jf(j))
                if(m>=iptr.and.m<=n) then
                    p(m)=p(iptr)
                    p(iptr)=jf(j)
                    q(p(m))=m
                    q(p(iptr))=iptr
                    iptr=iptr+1
                endif
            enddo
            if(next<=n) go to 30
        enddo
c
c       mark coarse graph vertices
c
        nc=0
        do k=1,n
            i=p(k)
            if(vtype(i)==0) then
                nc=nc+1
                qc(nc)=i
                vtype(i)=1
                do j=jf(i),jf(i+1)-1
                    vtype(jf(j))=-1
                enddo
            endif
        enddo
c
        jtnum=jtnum+1
        if(jtnum>jtmax.or.nc<=ncmax) return
c
c       compute length of new jf array
c
        do i=1,n
            q(i)=0
        enddo
        num=0
        iptr=jf(n+1)
        jf(iptr)=iptr+n+1
        next=jf(iptr)
        do i=1,n
            call cp2(i,jf,len,p,q)
            if(next+len+jf(i+1)-jf(i)>maxjf) return
            num=num+len
            do j=jf(i),jf(i+1)-1
                jf(next)=jf(j)
                next=next+1
            enddo
            do j=1,len
                jf(next)=p(j)
                next=next+1
            enddo
            iptr=iptr+1
            jf(iptr)=next
        enddo
        ishift=jf(n+1)-1
        do i=1,jf(iptr)-1
            jf(i)=jf(i+ishift)
        enddo
        do i=1,n+1
            jf(i)=jf(i)-ishift
        enddo
        go to 10
c
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine cp2(i,jf,len,list,q)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: jf,list,q
cy
c       compute lenght two paths from vertex i (assume q init to 0)
c
        do k=jf(i),jf(i+1)-1
            q(jf(k))=-k
        enddo
        q(i)=-i
        len=0
        do jj=jf(i),jf(i+1)-1
            j=jf(jj)
            do kk=jf(j),jf(j+1)-1
                k=jf(kk)
                if(q(k)==0) then
c
c       first length 2 path to k
c
                    len=len+1
                    list(len)=-k
                    q(k)=len
c
c       subsequent length 2 paths to k
c
                else if(q(k)>0) then
                    list(q(k))=k
                endif
            enddo
        enddo
c
c
        do k=jf(i),jf(i+1)-1
            q(jf(k))=0
        enddo
        q(i)=0
        do k=1,len
            q(abs(list(k)))=0
        enddo
c
c       shorten list
c
        len0=len
        len=0
        do k=1,len0
            if(list(k)>0) then
                len=len+1
                list(len)=list(k)
            endif
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine cvf(nf,nc,ispd,jf,maxvf,vf,maxwt,qc,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: vf,jf,qc
            integer(kind=iknd), dimension(nf) :: vtype,q
cy
        iflag=0
        if(nf+2>maxvf) then
            iflag=20
            return
        endif
        do i=1,nf
            vtype(i)=-1
            vf(i+1)=0
        enddo
        do i=1,nc
            q(qc(i))=i
            vtype(qc(i))=1
        enddo
c
c       compute pointers
c
        vf(1)=nf+2
        do i=1,nf
            if(vtype(i)>=0) then
                len=1
            else
                len=0
                do jj=jf(i),jf(i+1)-1
                    j=jf(jj)
                    if(vtype(j)>=0) len=len+1
                enddo
            endif
            vf(i+1)=vf(i)+len
        enddo
c
        if(vf(nf+1)-1>maxvf) then
            iflag=20
            return
        endif
c
c       fill out
c
        do i=1,nf
            k=vf(i)
            if(vtype(i)>=0) then
                vf(k)=q(i)
            else
                do jj=jf(i),jf(i+1)-1
                    j=jf(jj)
                    if(vtype(j)>=0) then
                        vf(k)=q(j)
                        k=k+1
                    endif
                enddo
            endif
        enddo
c
        maxvf=vf(nf+1)-1
        if(ispd==1) then
            lenwt=vf(nf+1)-vf(1)
        else
            lenwt=2*(vf(nf+1)-vf(1))
        endif
c
        if(lenwt>maxwt) then
            iflag=20
            return
        endif
        maxwt=lenwt
c
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine cwt(nf,nc,ispd,vf,wt,vtype,ja,a,qc,dv,dw)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,vf,vtype,qc
            integer(kind=iknd) :: vmtx,wmtx,amtx
            real(kind=rknd), dimension(*) :: a,wt,dv,dw
cy
c
        anorm=a(nf+1)
        if(ispd==1) then
            wmtx=1-vf(1)
            vmtx=wmtx
            amtx=0
        else
            wmtx=1-vf(1)
            vmtx=wmtx+vf(nf+1)-vf(1)
            amtx=ja(nf+1)-ja(1)
        endif
        do i=1,nf
            vtype(i)=-1
        enddo
        do i=1,nc
            vtype(qc(i))=1
        enddo
c
        do i=1,nf
            if(vtype(i)>=0) then
                wt(vf(i)+vmtx)=1.0e0_rknd
                wt(vf(i)+wmtx)=1.0e0_rknd
            else if(vtype(i)<0) then
                if(abs(a(i))<=anorm) then
                    ainv=(a(i)/anorm)/anorm
                else
                    ainv=1.0e0_rknd/a(i)
                endif
                do jj=vf(i),vf(i+1)-1
                    j=qc(vf(jj))
                    call jacmap(i,j,ij,ji,ja,amtx)
                    if(ij>0) then
                        wt(jj+wmtx)=-a(ij)*ainv
                        wt(jj+vmtx)=-a(ji)*ainv
                    else
                        wt(jj+wmtx)=0.0e0_rknd
                        wt(jj+vmtx)=0.0e0_rknd
                    endif
                enddo
                call crw(i,wmtx,vf,wt,qc,dw)
                if(ispd/=1) call crw(i,vmtx,vf,wt,qc,dv)
            endif
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine crw(i,wmtx,vf,wt,qc,dw)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: vf,qc
            integer(kind=iknd) :: wmtx
            real(kind=rknd), dimension(*) :: wt,dw
cy
c
        if(vf(i+1)==vf(i)) return
        eps=1.0e-4_rknd
        eps1=1.0e-1_rknd
c
c       check simple scaling option
c
        aa=0.0e0_rknd
        do jj=vf(i),vf(i+1)-1
            aa=aa+abs(wt(jj+wmtx))
        enddo
c
c      make row exact for constant
c
        if(aa<=eps) then
            c=1.0e0_rknd/real(vf(i+1)-vf(i))
            do jj=vf(i),vf(i+1)-1
                wt(jj+wmtx)=c
            enddo
            aa=1.0e0_rknd
        else
            do jj=vf(i),vf(i+1)-1
                wt(jj+wmtx)=wt(jj+wmtx)/aa
            enddo
        endif
c
c       check smooth vector
c
        a11=real(vf(i+1)-vf(i))
        a12=0.0e0_rknd
        a22=0.0e0_rknd
        b2=0.0e0_rknd
        do jj=vf(i),vf(i+1)-1
            j=qc(vf(jj))
            a22=a22+dw(j)**2
            b2=b2+wt(jj+wmtx)*dw(j)
            if(wt(jj+wmtx)>0.0e0_rknd) then
                a12=a12+dw(j)
            else
                a12=a12-dw(j)
            endif
        enddo
        b2=dw(i)-b2
        det=a11*a22-a12**2
        if(det<=(a11+a22)*eps1) return
c
c       chose to interpolate two vectors correctly
c
        x1=-a12*b2/det
        x2=a11*b2/det
        do jj=vf(i),vf(i+1)-1
             j=qc(vf(jj))
             if(wt(jj+wmtx)>0.0e0_rknd) then
                 wt(jj+wmtx)=wt(jj+wmtx)+x1+x2*dw(j)
             else
                 wt(jj+wmtx)=wt(jj+wmtx)-x1+x2*dw(j)
             endif
        enddo
c
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ceig(n,ispd,ja,a,ju,u,dv,dw,nc,vtype,coeff)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ju,vtype
            real(kind=rknd), dimension(*) :: a,u,dv,dw
            real(kind=rknd), dimension(n) :: z,az
            real(kind=rknd), dimension(2) :: coeff
cy
c
        coeff(1)=0.0e0_rknd
        coeff(2)=0.0e0_rknd
        itmax=3
        eps=1.0e-3_rknd
        if(nc>=n) return
        if(nc<=0) then
            do i=1,n
                vtype(i)=0
            enddo
        endif
        do i=1,n
            if(vtype(i)>0) then
                dv(i)=0.0e0_rknd
                dw(i)=0.0e0_rknd
            else
                dv(i)=1.0e0_rknd
                dw(i)=1.0e0_rknd
            endif
        enddo
        do itnum=1,itmax
            call mtxmlt(n,ja,a,dw,az,ispd)
            call snsilu(n,ju,u,z,az,ispd)
            wnorm=0.0e0_rknd
            do i=1,n
                if(vtype(i)>0) then
                    az(i)=0.0e0_rknd
                else
                    az(i)=dw(i)-z(i)
                    wnorm=max(az(i),wnorm)
                endif
            enddo
            if(wnorm<eps) go to 50
            coeff(1)=coeff(1)-log10(wnorm)
            do i=1,n
                dw(i)=az(i)/wnorm
            enddo
        enddo
   50   if(ispd==1) then
            do i=1,n
                dv(i)=dw(i)
            enddo
            coeff(2)=coeff(1)
            return
        endif
        jspd=-(1+ispd)
        do itnum=1,itmax
            call mtxmlt(n,ja,a,dv,az,jspd)
            call snsilu(n,ju,u,z,az,jspd)
            vnorm=0.0e0_rknd
            do i=1,n
                if(vtype(i)>0) then
                    az(i)=0.0e0_rknd
                else
                    az(i)=dv(i)-z(i)
                    vnorm=max(az(i),vnorm)
                endif
            enddo
            if(vnorm<eps) return
            coeff(2)=coeff(2)-log10(vnorm)
            do i=1,n
                dv(i)=az(i)/vnorm
            enddo
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine updhp(i,len,p,q,qual,isw)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: p,q
            real(kind=rknd), dimension(*) :: qual
cy
c       this routine makes a heap with root at vertex i, assuming its
c       sons are already roots of heaps
c
        if(len<=0) return
        k=i
        if(isw==0.or.k==1) go to 10
        kfath=k/2
        if(qual(p(k))>qual(p(kfath))) go to 60
c
c       push
c
   10   kson=2*k
        if(kson>len) return
        if(kson<len) then
            if(qual(p(kson+1))>qual(p(kson))) kson=kson+1
        endif
        if(qual(p(k))>=qual(p(kson))) return
        itemp=p(k)
        p(k)=p(kson)
        p(kson)=itemp
        q(p(kson))=kson
        q(p(k))=k
        k=kson
        go to 10
c
c       pull
c
   50   kfath=k/2
        if(kfath==0) return
        if(qual(p(kfath))>qual(p(k))) return
   60   itemp=p(k)
        p(k)=p(kfath)
        p(kfath)=itemp
        q(p(kfath))=kfath
        q(p(k))=k
        k=kfath
        go to 50
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine mg(ispd,lvl,mxcg,eps1,ja,a,dr,br,ka,relerr,
     +      iflag,hist)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd), dimension(10,*) :: ka
            real(kind=rknd), dimension(*) :: a,dr,br,hist
cy
c       lenz=6*n+2*ns (ispd==1)
c       lenz=11*n+2*ns (ispd/=1)
c
c       compute initial norm of b
c
        iflag=0
        eps2=epsilon(1.0e0_rknd)*8.0e0_rknd
        eps=max(eps1,eps2)
        epsi=1.0e0_rknd/min(eps,eps2)
        n=ka(1,1)
        call getptr(lvl,lvl,nf,nptr,japtr,iaptr,
     +      juptr,iuptr,jvptr,ivptr,iqptr,ibptr,nc,ncptr,ka)
c
c
c
        do i=1,n
            dr(ja(iqptr+i-1))=br(i)
        enddo
        if(ispd==1) then
            call cscg(n,ispd,lvl,mxcg,eps,epsi,ja,a,br,dr,
     +          hist,ka,relerr,iflag)
c
        else
            call csbcg(n,ispd,lvl,mxcg,eps,epsi,ja,a,br,dr,
     +          hist,ka,relerr,iflag)
c
        endif
        if(iflag==0) then
            do i=1,n
                dr(i)=br(ja(iqptr+i-1))
            enddo
        else
            do i=1,n
                dr(i)=0.0e0_rknd
            enddo
        endif
c
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine perm(n,x,ja,isw)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            real(kind=rknd), dimension(*) :: x
            real(kind=rknd), dimension(n) :: z
cy
c       reorder x
c
        iqptr=ja(n+1)+n
        if(isw==1) then
            do i=1,n
                z(ja(iqptr+i-1))=x(i)
            enddo
        else
            do i=1,n
                z(i)=x(ja(iqptr+i-1))
            enddo
        endif
        do i=1,n
            x(i)=z(i)
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine cycle(ispd,lvl,ja,a,x,b,ka)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd), dimension(101) :: kount
            integer(kind=iknd), dimension(10,*) :: ka
            real(kind=rknd), dimension(*) :: x,a,b
            real(kind=rknd), allocatable, dimension(:) :: r,y,dy
cy
c       compute initial norm of b
c
        ivwcy=1
        n=ka(1,1)
        ns=ka(2,lvl+1)-1
        allocate(r(ns),y(ns),dy(n))
        level=lvl
        call getptr(level,lvl,nf,nptr,japtr,iaptr,
     +      juptr,iuptr,jvptr,ivptr,iqptr,ibptr,nc,ncptr,ka)
        do i=1,nf
            r(i+nptr-1)=b(i)
            y(i+nptr-1)=0.0e0_rknd
        enddo
        kount(lvl)=-1
c
c
c       the smoothing iterations
c
   10   if(level==1) then
            kount(level)=ivwcy+1
        else
            kount(level)=kount(level)+1
        endif
        call snsilu(nf,ja(juptr),a(iuptr),dy,r(nptr),ispd)
        call resid(nf,ispd,ja(japtr),a(iaptr),y(nptr),
     +      r(nptr),dy)
c
c
        if(level==lvl.and.kount(lvl)>=1) then
            do i=1,nf
                x(i)=y(i+nptr-1)
            enddo
            deallocate(r,y,dy)
            return
        endif
c
        if(kount(level)>ivwcy) then
c
c       increase level,  go to finer grid
c
            level=level+1
            call getptr(level,lvl,nf,nptr,japtr,iaptr,
     +          juptr,iuptr,jvptr,ivptr,iqptr,ibptr,nc,ncptr,ka)
            call cr2fn(nf,nc,ispd,dy,y(ncptr),
     +          ja(jvptr),a(ivptr))
            call resid(nf,ispd,ja(japtr),a(iaptr),y(nptr),
     +          r(nptr),dy)
        else
c
c       decrease level, go to coarse grid
c
            call fn2cr(nf,nc,ispd,r(nptr),r(ncptr),
     +          ja(jvptr),a(ivptr))
            do i=1,nc
                y(i+ncptr-1)=0.0e0_rknd
            enddo
            level=level-1
            kount(level)=0
            call getptr(level,lvl,nf,nptr,japtr,iaptr,
     +          juptr,iuptr,jvptr,ivptr,iqptr,ibptr,nc,ncptr,ka)
        endif
        go to 10
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine cscg(n,ispd,lvl,mxcg,eps,epsi,ja,a,dr,br,
     +      hist,ka,relerr,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd), dimension(10,*) :: ka
            real(kind=rknd), dimension(*) :: a,dr,br,hist
            real(kind=rknd), dimension(n) :: pr,apr,zr,azr
cy
c       initialize
c
        iflag=0
        epsmin=0.5e0_rknd
        relerr=0.0e0_rknd
c
c       compute initial norm of b
c
        do i=1,n
            dr(i)=0.0e0_rknd
        enddo
        brnorm=rl2nrm(n,br)
        call hist1(hist,0_iknd,brnorm)
        if(brnorm<=0.0e0_rknd) return
        rrnorm=brnorm
c
c       compute initial pr and apr
c
        call cycle(ispd,lvl,ja,a,pr,br,ka)
        call mtxmlt(n,ja,a,pr,apr,ispd)
        bp=rl2ip(n,pr,br)
        if(bp==0.0e0_rknd) return
c
c       the main loop
c
        do 100 itnum=1,mxcg
c
c       compute sigma, the next 'psuedo residual' and precondition
c
            pap=rl2ip(n,pr,apr)
            do i=1,n
                azr(i)=pap*br(i)-bp*apr(i)
            enddo
            zscale=rl2nrm(n,azr)
            if(zscale>0.0e0_rknd) then
                do i=1,n
                    azr(i)=azr(i)/zscale
                enddo
            endif
            call cycle(ispd,lvl,ja,a,zr,azr,ka)
c
c       compute alphas
c
            bz=rl2ip(n,zr,azr)*zscale/pap
            zap=-bz/bp
            do i=1,n
                zr(i)=zr(i)-zap*pr(i)
            enddo
            call mtxmlt(n,ja,a,zr,azr,ispd)
            zaz=rl2ip(n,zr,azr)
c
c       decide on pivoting strategy
c
            if(abs(pap)*rrnorm<zscale) then
                qscale=tstpiv(n,bp,bz,pap,zaz,br,apr,azr)
                if(qscale<abs(zscale*zaz)) go to  50
            endif
c
c       the case of a 1 x 1 pivot
c
            alpha=bp/pap
            bp=bz
            do i=1,n
                dr(i)=dr(i)+alpha*pr(i)
                br(i)=br(i)-alpha*apr(i)
                pr(i)=zr(i)
                apr(i)=azr(i)
            enddo
c
c       convergence test
c
            rrnorm=rl2nrm(n,br)
            call hist1(hist,itnum,rrnorm)
            relerr=rrnorm/brnorm
cc          write(6,*) itnum,relerr
            if(relerr<=eps.or.bp==0.0e0_rknd) return
            if(relerr>epsi) go to 200
            go to 100
c
c       the case of a 2 x 2 pivot
c
   50       alphap=bp/pap
            alphaz=bz/zaz
            do i=1,n
                dr(i)=dr(i)+(alphap*pr(i)+alphaz*zr(i))
                br(i)=br(i)-(alphap*apr(i)+alphaz*azr(i))
            enddo
c
c       convergence test
c
            rrnorm=rl2nrm(n,br)
            call hist1(hist,itnum,-rrnorm)
            relerr=rrnorm/brnorm
cc          write(6,*) -itnum,relerr
            if(relerr<=eps) return
            if(relerr>epsi) go to 200
c
c       compute next direction
c
            call cycle(ispd,lvl,ja,a,apr,br,ka)
            bp=rl2ip(n,apr,br)
            betaz=bp/bz
            do i=1,n
                pr(i)=apr(i)+betaz*zr(i)
            enddo
            call mtxmlt(n,ja,a,pr,apr,ispd)
  100   continue
        if(relerr>epsmin) iflag=12
c
        return
  200   iflag=-12
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine csbcg(n,ispd,lvl,mxcg,eps,epsi,ja,a,dr,br,
     +      hist,ka,relerr,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd), dimension(10,*) :: ka
            real(kind=rknd), dimension(*) :: a,dr,br,hist
            real(kind=rknd), dimension(n) :: pr,apr,zr,azr
            real(kind=rknd), dimension(n) :: pl,apl,zl,azl,bl
cy
c       initialize
c
        iflag=0
        epsmin=0.5e0_rknd
        relerr=0.0e0_rknd
c
c       compute initial norm of b
c
        do i=1,n
            dr(i)=0.0e0_rknd
            bl(i)=br(i)
        enddo
        jspd=-(1+ispd)
        brnorm=rl2nrm(n,br)
        blnorm=rl2nrm(n,bl)
        call hist1(hist,0_iknd,brnorm)
        if(min(brnorm,blnorm)<=0.0e0_rknd) return
        rrnorm=brnorm
c
c       compute initial pr and apr
c
        call cycle(ispd,lvl,ja,a,pr,br,ka)
        call cycle(jspd,lvl,ja,a,pl,bl,ka)
        call mtxmlt(n,ja,a,pr,apr,ispd)
        call mtxmlt(n,ja,a,pl,apl,jspd)
        bp=rl2ip(n,pl,br)
        if(bp==0.0e0_rknd) return
c
c       the main loop
c
        do 100 itnum=1,mxcg
c
c       compute sigma, the next 'psuedo residual' and precondition
c
            pap=rl2ip(n,pl,apr)
            do i=1,n
                azr(i)=pap*br(i)-bp*apr(i)
                azl(i)=pap*bl(i)-bp*apl(i)
            enddo
            zscale=rl2nrm(n,azr)
            if(zscale>0.0e0_rknd) then
                do i=1,n
                    azr(i)=azr(i)/zscale
                    azl(i)=azl(i)/zscale
                enddo
            endif
            call cycle(ispd,lvl,ja,a,zr,azr,ka)
            call cycle(jspd,lvl,ja,a,zl,azl,ka)
c
c       compute alphas
c
            bz=rl2ip(n,zl,azr)*zscale/pap
            zap=-bz/bp
            do i=1,n
                zr(i)=zr(i)-zap*pr(i)
                zl(i)=zl(i)-zap*pl(i)
            enddo
            call mtxmlt(n,ja,a,zr,azr,ispd)
            call mtxmlt(n,ja,a,zl,azl,jspd)
            zaz=rl2ip(n,zl,azr)
c
c       decide on pivoting strategy
c
            if(abs(pap)*rrnorm<zscale) then
                qscale=tstpiv(n,bp,bz,pap,zaz,br,apr,azr)
                if(qscale<abs(zscale*zaz)) go to 50
            endif
c
c       the case of a 1 x 1 pivot
c
            alpha=bp/pap
            bp=bz
            do i=1,n
                dr(i)=dr(i)+alpha*pr(i)
                br(i)=br(i)-alpha*apr(i)
                bl(i)=bl(i)-alpha*apl(i)
                pr(i)=zr(i)
                pl(i)=zl(i)
                apr(i)=azr(i)
                apl(i)=azl(i)
            enddo
c
c       convergence test
c
            rrnorm=rl2nrm(n,br)
cc          rlnorm=rl2nrm(n,bl)
            call hist1(hist,itnum,rrnorm)
            relerr=rrnorm/brnorm
cc          write(6,*) itnum,relerr
            if(relerr<=eps) return
            if(relerr>epsi) go to 200
            go to 100
c
c       the case of a 2 x 2 pivot
c
   50       alphap=bp/pap
            alphaz=bz/zaz
            do i=1,n
                dr(i)=dr(i)+(alphap*pr(i)+alphaz*zr(i))
                br(i)=br(i)-(alphap*apr(i)+alphaz*azr(i))
                bl(i)=bl(i)-(alphap*apl(i)+alphaz*azl(i))
            enddo
c
c       convergence test
c
            rrnorm=rl2nrm(n,br)
cc          rlnorm=rl2nrm(n,bl)
            call hist1(hist,itnum,-rrnorm)
            relerr=rrnorm/brnorm
cc          write(6,*) -itnum,relerr
            if(relerr<=eps) return
            if(relerr>epsi) go to 200
c
c       compute next direction
c
            call cycle(ispd,lvl,ja,a,apr,br,ka)
            call cycle(jspd,lvl,ja,a,apl,bl,ka)
            bp=rl2ip(n,apl,br)
            betaz=bp/bz
            do i=1,n
                pr(i)=apr(i)+betaz*zr(i)
                pl(i)=apl(i)+betaz*zl(i)
            enddo
            call mtxmlt(n,ja,a,pr,apr,ispd)
            call mtxmlt(n,ja,a,pl,apl,jspd)
  100   continue
        if(relerr>epsmin) iflag=12
c
        return
  200   iflag=-12
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        function tstpiv(n,bp,bz,pap,zaz,br,apr,azr)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: br,apr,azr
            real(kind=rknd) :: tstpiv
cy
c       compute norm to decide between 1x1 and 2x2 pivoting
c
        alphap=bp*zaz
        alphaz=bz*pap
        alpha=zaz*pap
        qscale=0.0e0_rknd
        qmax=0.0e0_rknd
        do i=1,n
            dq=alpha*br(i)-(alphap*apr(i)+alphaz*azr(i))
            if(abs(dq)<qmax) then
                qscale=qscale+(dq/qmax)**2
            else if(dq/=0.0e0_rknd) then
                qscale=1.0e0_rknd+qscale*(qmax/dq)**2
                qmax=abs(dq)
            endif
        enddo
        tstpiv=sqrt(qscale)*qmax
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine resid(n,ispd,ja,a,x,b,p)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            real(kind=rknd), dimension(*) :: a,b,x,p
            real(kind=rknd), dimension(n) :: ap
cy
c       residual update
c
        call mtxmlt(n,ja,a,p,ap,ispd)
        do i=1,n
            x(i)=x(i)+p(i)
            b(i)=b(i)-ap(i)
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine getptr(level,lvl,nf,nptr,japtr,iaptr,
     +      juptr,iuptr,jvptr,ivptr,iqptr,ibptr,nc,ncptr,ka)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(10,*) :: ka
cy
c       get pointers for level
c
        k=lvl+1-level
        nf   =ka(1,k)
        nptr =ka(2,k)
        japtr=ka(3,k)
        iaptr=ka(4,k)
        juptr=ka(5,k)
        iuptr=ka(6,k)
        jvptr=ka(7,k)
        ivptr=ka(8,k)
        iqptr=ka(9,k)
        ibptr=ka(10,k)
        if(level>1) then
            nc=ka(1,k+1)
            ncptr=ka(2,k+1)
        else
            nc=0
            ncptr=0
        endif
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine fn2cr(nf,nc,ispd,rf,rc,vf,wt)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: vf
            integer(kind=iknd) :: vmtx
            real(kind=rknd), dimension(*) :: rf,rc,wt
cy
c       fine to coarse transfer
c
        vmtx=1-vf(1)
        if(ispd==0) vmtx=vmtx+vf(nf+1)-vf(1)
        do i=1,nc
            rc(i)=0.0e0_rknd
        enddo
c
        do i=1,nf
            do j=vf(i),vf(i+1)-1
                rc(vf(j))=rc(vf(j))+wt(j+vmtx)*rf(i)
            enddo
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine cr2fn(nf,nc,ispd,xf,xc,vf,wt)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: vf
            integer(kind=iknd) :: wmtx
            real(kind=rknd), dimension(*) :: xf,xc,wt
cy
c       coarse to fine transfer
c
        wmtx=1-vf(1)
        if(ispd==-1) wmtx=wmtx+vf(nf+1)-vf(1)
        do i=1,nf
            xf(i)=0.0e0_rknd
            do j=vf(i),vf(i+1)-1
                xf(i)=xf(i)+wt(j+wmtx)*xc(vf(j))
            enddo
        enddo
c
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine mtxmlt(n,ja,a,x,b,ispd)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd) :: umtx,lmtx
            real(kind=rknd), dimension(*) :: a,x,b
cy
c       ispd = 1   symmetric
c            = 0   non-symmetric
c            =-1   non-symmetric for a-transpose
c
c       compute b=a*x
c
        lmtx=0
        umtx=0
        if(ispd==0) lmtx=ja(n+1)-ja(1)
        if(ispd==-1) umtx=ja(n+1)-ja(1)
c
        do i=1,n
            b(i)=a(i)*x(i)
        enddo
c
        do i=1,n
            do jj=ja(i),ja(i+1)-1
                j=ja(jj)
                b(i)=b(i)+a(jj+umtx)*x(j)
                b(j)=b(j)+a(jj+lmtx)*x(i)
            enddo
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine sfilu(n,ja,a,maxju,ju,maxu,u,ispd,dtol,maxfil,iflag)
c
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ju
            integer(kind=iknd), dimension(400) :: ibin
            integer(kind=iknd), dimension(n) :: list,mark,indx
            integer(kind=iknd)  :: amtx,umtx
            real(kind=rknd), dimension(*) :: a,u
            real(kind=rknd), dimension(n) :: tl,tu
cy
c       sparse numeric factorization
c
        nbin=400
        if(min(maxju,maxu)<n+1) then
            iflag=20
            return
        endif
        if(ispd/=1) then
            lenju=min(maxju,(maxu+n+1)/2)
        else
            lenju=min(maxju,maxu)
        endif
        qn=real(max(0,maxfil))*real(n)+real(n+1)
        if(qn<real(lenju)) lenju=int(qn+0.5e0_rknd)
        if(ispd/=1) then
            amtx=ja(n+1)-ja(1)
            umtx=lenju-(n+1)
        else
            amtx=0
            umtx=0
        endif
        fact=1.0e4**(1.0e0_rknd/real(nbin))
        alf=log(fact)
        eps=epsilon(1.0e0_rknd)
        rtol=max(eps,dtol)
        u(n+1)=a(n+1)
c
    5   ju(1)=n+2
        unorm=u(n+1)
        do i=1,n
            mark(i)=0
            list(i)=0
            indx(i)=0
        enddo
c
        num=0
        do i=1,nbin
            ibin(i)=0
        enddo
c
        do i=1,n
c
c       initialize row i and col i in tu and tl
c
            mark(i)=i
            len=0
            tu(i)=a(i)
            tl(i)=a(i)
            do jj=ja(i),ja(i+1)-1
                j=ja(jj)
                tu(j)=a(jj)
                tl(j)=a(jj+amtx)
                mark(j)=mark(i)
                mark(i)=j
                len=len+1
            enddo
c
c       do outer product updates
c
            lk=list(i)
   10       if(lk>0) then
                k=lk
                lk=list(k)
                j1=indx(k)
                j2=ju(k+1)-1
                if(abs(u(k))<=unorm) then
                    uinv=(u(k)/unorm)/unorm
                else
                    uinv=1.0e0_rknd/u(k)
                endif
                su=u(j1)*uinv
                sl=u(j1+umtx)*uinv
c
                do jj=j1,j2
                    j=ju(jj)
                    if(mark(j)/=0) then
                        tu(j)=tu(j)-sl*u(jj)
                        tl(j)=tl(j)-su*u(jj+umtx)
                    else
                        tu(j)=-sl*u(jj)
                        tl(j)=-su*u(jj+umtx)
                        mark(j)=mark(i)
                        mark(i)=j
                        len=len+1
                    endif
                enddo
                if(j1<j2) then
                    j=ju(j1+1)
                    list(k)=list(j)
                    list(j)=k
                    indx(k)=j1+1
                endif
                go to 10
            endif
c
c       check diagonal entry
c
            u(i)=tu(i)
c
c       make ju for this row
c
            next=ju(i)
            do j=1,len
                k=mark(i)
                tt=eltest(tu(i),tu(k),tl(k),a(k))/rtol
                if(tt>=1.0e0_rknd) then
                    it=min(nbin,1_iknd+int(log(tt)/alf))
                    ibin(it)=ibin(it)+1
                    if(next<lenju) then
                        ju(next)=k
                        next=next+1
                    else
                        num=num+1
                    endif
                endif
                mark(i)=mark(k)
                mark(k)=0
            enddo
            mark(i)=0
            ju(i+1)=next
            len=next-ju(i)
            if(len>1) call ihp(ju(ju(i)),len)
c
c       move tl, tu to u
c
            do jj=ju(i),ju(i+1)-1
                j=ju(jj)
                u(jj)=tu(j)
                u(jj+umtx)=tl(j)
            enddo
c
            if(ju(i)<ju(i+1)) then
                j=ju(ju(i))
                list(i)=list(j)
                list(j)=i
                indx(i)=ju(i)
            endif
        enddo
        if(num>0) then
            do i=nbin-1,1,-1
                ibin(i)=ibin(i+1)+ibin(i)
            enddo
            theta=1.01e0_rknd+0.4e0_rknd*real(num)
     +          /real(num+lenju-(n+1))
            nnz=int(real(lenju-(n+1))/theta)
            do i=1,nbin-1
                if(ibin(i+1)<=nnz.and.ibin(i)>nnz) go to 20
            enddo
            i=nbin
   20       rtol=rtol*fact**i
            write(6,*) 'sfilu',n,i,num,
     +          real(num)/real(num+lenju-(n+1)),rtol
            go to 5
        endif
        iflag=0
c
c       shift u for non symmetric case
c
        maxju=ju(n+1)-1
        if(ispd/=1) then
            nnz=maxju-(n+1)
            do i=1,nnz
                u(maxju+i)=u(lenju+i)
            enddo
            maxu=maxju+nnz
        else
            maxu=maxju
        endif
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        function eltest(a11,a12,a21,a22)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd) :: eltest
cy
        bd=a11*a22
        if(bd==0.0e0_rknd) then
            eltest=1.0e0_rknd
        else
            eltest=max(abs(a21),abs(a12))/sqrt(abs(bd))
        endif
c
cc      eltest=1.0e0_rknd
cc      aa=max(abs(a11),abs(a22),abs(a12),abs(a21))
cc      if(aa==0.0e0_rknd) return
cc      bd=(abs(a11)/aa)*(abs(a22)/aa)
cc      if(bd==0.0e0_rknd) return
cc      eltest=max(abs(a21),abs(a12))/(aa*sqrt(bd))
c
cc      bn=(abs(a21)/aa)*(abs(a12)/aa)
cc      eltest=sqrt(bn/bd)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine snfilu(n,ja,a,ju,u,ispd)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ju
            integer(kind=iknd), dimension(n) :: list,mark,indx
            integer(kind=iknd) :: amtx,umtx
            real(kind=rknd), dimension(*) :: a,u
            real(kind=rknd), dimension(n) :: tl,tu
cy
c       sparse numeric factorization
c
        if(ispd/=1) then
            amtx=ja(n+1)-ja(1)
            umtx=ju(n+1)-ju(1)
        else
            amtx=0
            umtx=0
        endif
        u(n+1)=a(n+1)
        unorm=u(n+1)
        do i=1,n
            mark(i)=0
            list(i)=0
            indx(i)=0
        enddo
c
        do i=1,n
c
c       initialize row i and col i in tu and tl
c
            mark(i)=1
            do jj=ju(i),ju(i+1)-1
                j=ju(jj)
                tu(j)=0.0e0_rknd
                tl(j)=0.0e0_rknd
                mark(j)=1
            enddo
            tu(i)=a(i)
            tl(i)=a(i)
            do jj=ja(i),ja(i+1)-1
                j=ja(jj)
                tu(j)=a(jj)
                tl(j)=a(jj+amtx)
            enddo
c
c       do outer product updates
c
            lk=list(i)
   10       if(lk>0) then
                k=lk
                lk=list(k)
                j1=indx(k)
                j2=ju(k+1)-1
                if(abs(u(k))<=unorm) then
                    uinv=(u(k)/unorm)/unorm
                else
                    uinv=1.0e0_rknd/u(k)
                endif
                su=u(j1)*uinv
                sl=u(j1+umtx)*uinv
c
                do jj=j1,j2
                    j=ju(jj)
                    if(mark(j)/=0) then
                        tu(j)=tu(j)-sl*u(jj)
                        tl(j)=tl(j)-su*u(jj+umtx)
                    endif
                enddo
                if(j1<j2) then
                    j=ju(j1+1)
                    list(k)=list(j)
                    list(j)=k
                    indx(k)=j1+1
                endif
                go to 10
            endif
c
c       move tl, tu to u
c
            u(i)=tu(i)
            mark(i)=0
            do jj=ju(i),ju(i+1)-1
                j=ju(jj)
                u(jj)=tu(j)
                u(jj+umtx)=tl(j)
                mark(j)=0
            enddo
c
            if(ju(i)<ju(i+1)) then
                j=ju(ju(i))
                list(i)=list(j)
                list(j)=i
                indx(i)=ju(i)
            endif
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine snsilu(n,ju,u,x,b,ispd)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ju
            integer(kind=iknd) :: lmtx,umtx
            real(kind=rknd), dimension(*) :: u,x,b
cy
c       ispd = 1   symmetric
c            = 0   non-symmetric
c            =-1   non-symmetric for a-transpose
c
c       solve a*x=b
c
        lmtx=0
        umtx=0
        if(ispd==0) lmtx=ju(n+1)-ju(1)
        if(ispd==-1) umtx=ju(n+1)-ju(1)
c
        do i=1,n
            x(i)=b(i)
        enddo
c
c       lower triangular system
c
        do i=1,n
            x(i)=x(i)/u(i)
            do jj=ju(i),ju(i+1)-1
                j=ju(jj)
                x(j)=x(j)-u(jj+lmtx)*x(i)
            enddo
        enddo
c
c       upper triangular system
c
        do i=n,1,-1
            s=0.0e0_rknd
            do jj=ju(i),ju(i+1)-1
                j=ju(jj)
                s=s+u(jj+umtx)*x(j)
            enddo
            x(i)=x(i)-s/u(i)
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine md(n,lenja,ja,p,lenu)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,p
            integer(kind=iknd), dimension(n) :: mark,equiv,befor,after
            integer(kind=iknd), dimension(n) :: list
            integer(kind=iknd), dimension(2*lenja-n) :: jc
cy
c
c       minimum degree algorithm
c
c       list        = linked list of equivalent vertices (v,e)
c                   = (temp) ptr to equiv vertex with clique imin (c)
c       equiv       = number of equivalent vertices (v)
c                   = ptr to equivant vertex (e)
c                   = size of clique (c)
c       befor/after = doubly linked list of verts. by degree (v)
c                     (temp) nvert/ncliq for verts in imin
c                     (temp) marker for outmatched verts in imin
c                   = (temp) switch/intersection size with imin (c)
c       mark        = temp linked list
c       p           = order (tail is head ptrs into befor/after)
c
        lenu=n+1
        mndeg=n+1
        iempty=0
        next=1
        do i=1,n
            p(i)=0
            equiv(i)=1
            list(i)=i
            befor(i)=0
            after(i)=0
            mark(i)=0
        enddo
        call ja2jc(n,ja,jc)
        do i=1,n
            ideg=jc(i+1)-jc(i)
            if(ideg<=0) then
                p(next)=i
                next=next+1
            else
                id=n+1-ideg
                if(p(id)/=0) befor(p(id))=i
                after(i)=p(id)
                p(id)=i
                befor(i)=-id
                mndeg=min(mndeg,ideg)
            endif
        enddo
        if(next>n) go to 100
c
c       order vertex of min degree
c
   10   id=n+1-mndeg
        if(p(id)==0) then
            mndeg=mndeg+1
            go to 10
        endif
        imin=p(id)
        if(after(imin)>0) befor(after(imin))=-id
        p(id)=after(imin)
        befor(imin)=0
        after(imin)=0
c
c       build the current clique (imin)
c
        call mkcliq(imin,jc,mark,equiv,ilen,imndeg,iempty)
c
        numequ=equiv(imin)
        i=imin
        do ii=1,numequ
            p(next)=i
            next=next+1
            equiv(i)=0
            lenu=lenu+imndeg+numequ-ii
            i=list(i)
        enddo
        if(next>n) go to 100
c
c       if the fillin will create a dense matrix....
c
        if(next+imndeg>n) then
            i=imin
            numequ=0
            do ii=1,ilen
                i=mark(i)
                inum=equiv(i)
                m=i
                do mm=1,inum
                    p(next)=m
                    next=next+1
                    equiv(m)=0
                    numequ=numequ+1
                    lenu=lenu+imndeg-numequ
                    m=list(m)
                enddo
            enddo
            go to 100
        endif
c
c       eliminate redundant vertices from adjacency lists of clique
c       members...this allows simple elimination of equivalent vertices
c
        i=imin
        numequ=0
        jx=imin
        jlen=0
        do ii=1,ilen
            i=mark(i)
            if(after(i)>0) befor(after(i))=befor(i)
            if(befor(i)<0) then
                id=-befor(i)
                if(id>=next)  p(id)=after(i)
            else
                after(befor(i))=after(i)
            endif
            befor(i)=0
            after(i)=0
c
c       update adjacency list
c
            call jcupdt(imin,i,jc,mark,equiv,befor,after,
     +          nvert,ncliq,ideg)
c
c       test for equivalence
c
            if(nvert==0.and.ncliq==1) then
                inum=equiv(i)
                m=i
                do mm=1,inum
                    p(next)=m
                    next=next+1
                    equiv(m)=0
                    numequ=numequ+1
                    lenu=lenu+imndeg-numequ
                    m=list(m)
                enddo
            endif
c
c       look for equivalent vertices
c
            if(nvert==0.and.ncliq==2) then
                jcj=-jc(jc(i))
                if(mark(jcj)==0) then
                    mark(jcj)=jx
                    jx=jcj
                    jlen=jlen+1
                    list(jcj)=i
                else
                    ieq=list(jcj)
                    inum=equiv(i)
                    equiv(ieq)=equiv(ieq)+inum
                    m=list(i)
                    do mm=1,inum
                        mnext=list(m)
                        list(m)=list(ieq)
                        list(ieq)=m
                        equiv(m)=-ieq
                        m=mnext
                    enddo
                endif
            endif
c
c       save partial degree (imin is not counted yet)
c
            if(equiv(i)>0) after(i)=ideg
        enddo
        if(next>n) go to 100
c
c       update degrees
c
        equiv(imin)=imndeg-numequ
        i=imin
        do ii=1,ilen
            i=mark(i)
            if(equiv(i)<=0) cycle
c
c       overcounting with three cliques requires this
c
            id=n+1-min(after(i)+equiv(imin)-1,n-next)
            if(p(id)/=0) befor(p(id))=i
            after(i)=p(id)
            p(id)=i
            befor(i)=-id
        enddo
c
c       clean up  mark, move clique to jc
c
        call svcliq(imin,jc,mark,equiv,ilen,iempty)
c
c     update cliques
c
        do jj=1,jlen
            jnext=mark(jx)
            call clqupd(jx,jc,mark,equiv,iempty)
            jx=jnext
        enddo
c
        mndeg=max(1,equiv(imin))
        if(next<=n)  go to 10
c
c       compute inverse permutation
c
  100   do i=1,n
            mark(p(i))=i
        enddo
        do i=1,n
            p(i)=mark(i)
        enddo
c
c       reversing order is specific to bank/smith bordering algorithm
c
cc      nn=n/2
cc      do i=1,nn
cc          ii=p(i)
cc          p(i)=p(n+1-i)
cc          p(n+1-i)=ii
cc      enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine mkcliq(imin,jc,mark,equiv,ilen,imndeg,iempty)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: jc,mark,equiv
cy
        mark(imin)=imin
        imndeg=0
        ilen=0
        do j=jc(imin),jc(imin+1)-1
            jcj=abs(jc(j))
            if(jcj==0) return
            if(jc(j)>0) then
c
c       merge a normal vertex
c
                if(mark(jcj)==0) then
                    mark(jcj)=mark(imin)
                    mark(imin)=jcj
                    imndeg=imndeg+equiv(jcj)
                    ilen=ilen+1
                endif
c
c       merge a clique
c
            else
   10           equiv(jcj)=0
                mark(jcj)=iempty
                iempty=jcj
                do m=jc(jcj),jc(jcj+1)-1
                    jcj=abs(jc(m))
                    if(jc(m)<0) go to 10
                    if(jc(m)==0) exit
                    if(mark(jcj)/=0) cycle
                    mark(jcj)=mark(imin)
                    mark(imin)=jcj
                    imndeg=imndeg+equiv(jcj)
                    ilen=ilen+1
                enddo
            endif
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine jcupdt(imin,i,jc,mark,equiv,befor,after,nvert,
     +      ncliq,ideg)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: jc,mark,equiv,
     +          befor,after
cy
c       update jc for vertex i
c
        iptr=jc(i)
        nvert=0
        ncliq=1
        ideg=0
        do j=jc(i),jc(i+1)-1
            jcj=abs(jc(j))
            if(jcj==0) exit
            if(jc(j)>0) then
c
c       check a normal vertex
c
                if(mark(jcj)==0) then
                    jc(iptr)=jcj
                    iptr=iptr+1
                    nvert=nvert+1
                    ideg=ideg+equiv(jcj)
                endif
            else
c
c       this loop overestimates degrees for vertices
c       connected to three or more cliques
c       on the first encounter, compute the intersection
c
                if(equiv(jcj)<=0) cycle
                if(befor(jcj)/=-imin) then
                    befor(jcj)=-imin
                    after(jcj)=0
                    jck=jcj
   10               do k=jc(jck),jc(jck+1)-1
                        jck=abs(jc(k))
                        if(jc(k)<0) go to 10
                        if(jc(k)==0) exit
                        if(mark(jck)<=0)
     +                      after(jcj)=after(jcj)+equiv(jck)
                    enddo
                endif
                if(after(jcj)>0) then
                    jc(iptr)=-jcj
                    ncliq=ncliq+1
                    iptr=iptr+1
                    ideg=ideg+after(jcj)
                endif
            endif
        enddo
        jc(iptr)=-imin
        if(iptr+1<jc(i+1)) jc(iptr+1)=0
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine clqupd(imin,jc,mark,equiv,iempty)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: jc,mark,equiv
cy
c       delete equivalent vertices from clique list
c
        jcj=imin
        jcnext=jc(jcj)
        jclast=jc(jcj+1)-1
   10   jcur=jc(jcj)
        jend=jc(jcj+1)-1
   20   jcj=abs(jc(jcur))
        if(jcj==0) go to 40
        if(jc(jcur)<0) then
            equiv(jcj)=0
            mark(jcj)=iempty
            iempty=jcj
            go to 10
        endif
        if(equiv(jcj)>0) then
            if(jcnext>jclast) then
                locsv=jclast
   30           if(mark(iempty)==0) then
                    next=iempty
                    iempty=mark(next)
                else
                    next=mark(iempty)
                    mark(iempty)=mark(next)
                endif
                jcnext=jc(next)+1
                jclast=jc(next+1)-1
                if(jcnext>jclast) go to 30
                jc(jcnext-1)=jc(locsv)
                jc(locsv)=-next
            endif
c
            jc(jcnext)=jcj
            jcnext=jcnext+1
        endif
        jcur=jcur+1
        if(jcur<=jend) go to 20
   40   if(jcnext<=jclast) jc(jcnext)=0
        mark(imin)=0
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine svcliq(imin,jc,mark,equiv,ilen,iempty)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: jc,mark,equiv
cy
c       save clique imin in jc
c
        jcnext=jc(imin)
        jclast=jc(imin+1)-1
        i=imin
        do ii=1,ilen
            is=i
            i=mark(i)
            mark(is)=0
            if(equiv(i)<=0) cycle
c
c       pop the stack if necessary
c
            if(jcnext>jclast) then
                locsv=jclast
   10           next=iempty
                iempty=mark(next)
                jcnext=jc(next)+1
                jclast=jc(next+1)-1
                if(jcnext>jclast) go to 10
                jc(jcnext-1)=jc(locsv)
                jc(locsv)=-next
            endif
c
            jc(jcnext)=i
            jcnext=jcnext+1
c
        enddo
        mark(i)=0
        if(jcnext<=jclast) jc(jcnext)=0
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine jamap(i,j,ij,ji,ja,amtx)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd) :: amtx
cy
c       compute location of a(i,j) and a(j,i)
c
        if(i<j) then
            do ij=ja(i),ja(i+1)-1
                if(ja(ij)==j) then
                    ji=ij+amtx
                    return
                endif
            enddo
c
        else
            do ji=ja(j),ja(j+1)-1
                if(ja(ji)==i) then
                    ij=ji+amtx
                    return
                endif
            enddo
        endif
        ij=0
        ji=0
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine jamap0(i,j,n,ispd,ij,ji,ja)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd) :: amtx
cy
c       find location in ja array of entry (i,j) in original ordering
c
        iqptr=ja(n+1)-1+n
        iq=ja(iqptr+i)
        if(i==j) then
            ij=iq
            ji=iq
        else
            jq=ja(iqptr+j)
            if(ispd/=1) then
                amtx=ja(n+1)-ja(1)
            else
                amtx=0
            endif
            call jamap(iq,jq,ij,ji,ja,amtx)
        endif
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        function rl2nrm(n,b)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: b
            real(kind=rknd) :: rl2nrm
cy
c       compute norm of b and update history
c
        bnorm=0.0e0_rknd
        bmax=0.0e0_rknd
        do i=1,n
            if(abs(b(i))<bmax) then
                bnorm=bnorm+(b(i)/bmax)**2
            else if(b(i)/=0.0e0_rknd) then
                bnorm=1.0e0_rknd+bnorm*(bmax/b(i))**2
                bmax=abs(b(i))
            endif
        enddo
        rl2nrm=sqrt(bnorm)*bmax
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        function rl2ip(n,x,y)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: x,y
            real(kind=rknd) :: rl2ip
cy
c       compute dot product
c
        rl2ip=0.0e0_rknd
        spmax=0.0e0_rknd
        snmax=0.0e0_rknd
        sp=0.0e0_rknd
        sn=0.0e0_rknd
        do i=1,n
            t=x(i)*y(i)
            if(t>=0.0e0_rknd) then
                if(t<spmax) then
                    sp=sp+t/spmax
                else if(t/=0.0e0_rknd) then
                    sp=1.0e0_rknd+sp*(spmax/t)
                    spmax=t
                endif
            else
                if(-t<snmax) then
                    sn=sn+t/snmax
                else
                    sn=-(1.0e0_rknd+sn*(snmax/t))
                    snmax=-t
                endif
            endif
        enddo
        rl2ip=sp*spmax+sn*snmax
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine hist1(hist,itnum,bnorm)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(22) :: hist
cy
c       update history array
c
        mxhist=20
        if(itnum<=0) then
            hist(mxhist+2)=bnorm
        else if(itnum>mxhist) then
            do i=1,mxhist-1
               hist(i)=hist(i+1)
            enddo
            hist(mxhist)=bnorm
        else
            hist(itnum)=bnorm
        endif
        if(itnum>=0) hist(mxhist+1)=real(itnum)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine chkja(n,ja,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd), dimension(n) :: mark
cy
c       check consistency of ja data structure
c
        iflag=0
        do i=1,n
            mark(i)=0
        enddo
c
c       check pointers
c
        if(ja(1)<=0) then
            iflag=-1
            return
        endif
        do i=1,n
            if(ja(i+1)<ja(i)) then
                iflag=-1
                return
            endif
        enddo
c
c       check column indices
c
        do i=1,n
            do jj=ja(i),ja(i+1)-1
                j=ja(jj)
c
c       j not in upper triangle
c
                if(j<=i) then
                    iflag=-2
                    return
                else if(j>n) then
                    iflag=-3
                    return
c
c       duplicate entry
c
                else if(mark(j)==1) then
                    iflag=-4
                    return
                else
                    mark(j)=1
                endif
            enddo
c
c
            do j=ja(i),ja(i+1)-1
                mark(ja(j))=0
            enddo
        enddo
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine chkib(n,nblock,ib,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ib
cy
c       check consistency of ja data structure
c
        iflag=0
c
c       check ib pointers
c
        if(ib(1)/=1) then
            iflag=-5
            return
        endif
        do i=1,nblock
            if(ib(i+1)<=ib(i)) then
                iflag=-5
                return
            endif
        enddo
        if(ib(nblock+1)/=n+1) then
            iflag=-5
            return
        endif
        return
        end
