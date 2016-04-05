c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine gphplt(ip,rp,sp,hist,ka,time)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(100) :: ip
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd), dimension(10,*) :: ka
            real(kind=rknd), dimension(100) :: rp
            real(kind=rknd), dimension(10) :: red,green,blue
            real(kind=rknd), dimension(*) :: hist,time
            character(len=80), dimension(100) :: sp
cy
c
        ip(26)=0
c
        do i=1,25
            jp(i)=0
        enddo
        mxcolr=max(2_iknd,ip(51))
        igrsw=ip(54)
        if(igrsw<-3.or.igrsw>3) igrsw=0
c
        jp(4)=1
        jp(5)=6
        jp(17)=mxcolr
        jp(18)=min(mxcolr,jp(5)+2)
        jp(10)=igrsw
        jp(11)=ip(2)
c
        jp(13)=ip(64)
        jp(14)=ip(65)
        jp(15)=ip(66)
c
        call clrmap(red,green,blue,jp)
c
        call pltutl(jp(18),red,green,blue)
c
c       convergence history
c
        if(igrsw==0) then
            call pframe(4_iknd)
            call title0(sp(3),0_iknd)
            call hbplt(hist,1_iknd,1_iknd,jp)
            call pframe(-4_iknd)
            call pframe(2_iknd)
            call kaplt(ka,jp)
            call pframe(-2_iknd)
            call pframe(3_iknd)
            call pieplt(time,jp)
            call pframe(-3_iknd)
c
c       level plot
c
        else if(igrsw==1) then
            call pframe(4_iknd)
            call title0(sp(3),0_iknd)
            call kaplt(ka,jp)
            call pframe(-4_iknd)
            call pframe(2_iknd)
            call hbplt(hist,1_iknd,1_iknd,jp)
            call pframe(-2_iknd)
            call pframe(3_iknd)
            call pieplt(time,jp)
            call pframe(-3_iknd)
c
c       time plot
c
        else if(igrsw==-1) then
            call pframe(4_iknd)
            call title0(sp(3),0_iknd)
            call pieplt(time,jp)
            call pframe(-4_iknd)
            call pframe(2_iknd)
            call hbplt(hist,1_iknd,1_iknd,jp)
            call pframe(-2_iknd)
            call pframe(3_iknd)
            call kaplt(ka,jp)
            call pframe(-3_iknd)
c
c       ip array
c
        else if(igrsw==2) then
            call pframe(1_iknd)
            call title0(sp(3),1_iknd)
            call prtip(ip,jp)
            call pframe(-1_iknd)
c
c       rp array
c
        else if(igrsw==-2) then
            call pframe(1_iknd)
            call title0(sp(3),1_iknd)
            call prtrp(rp,jp)
            call pframe(-1_iknd)
c
c       ka array
c
        else if(igrsw==3) then
            call pframe(1_iknd)
            call title0(sp(3),1_iknd)
            lvl=ip(2)
            call prtka(lvl,ka,jp)
            call pframe(-1_iknd)
c
c       sp array
c
        else if(igrsw==-3) then
            call pframe(1_iknd)
            call title0(sp(3),1_iknd)
            call prtsp(sp,jp)
            call pframe(-1_iknd)
        endif
c
        call pltutl(-1_iknd,red,green,blue)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine kaplt(ka,jp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd), dimension(4) :: icolor
            integer(kind=iknd), dimension(10,*) :: ka
            integer(kind=iknd) :: ccolor
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            real(kind=rknd), dimension(3,100) :: e
            real(kind=rknd), dimension(100) :: rn
            real(kind=rknd), dimension(4) :: ratio
            character(len=80) :: ichr,jchr
            character(len=80), save, dimension(4) :: label
cy
            data label/'ja','ju','vf','n'/
c
c     graph error
c
        lvl=jp(11)
        if(lvl<=0) return
        call linit(t,q)
        size=t(14)
        xshift=t(15)-size/2.0e0_rknd
        yshift=t(16)-size/2.0e0_rknd
        t(1)=xshift
        t(2)=yshift
        t(3)=size
c
cc      icolor(1)=ccolor(6_iknd,0_iknd,jp)
cc      icolor(2)=ccolor(2_iknd,0_iknd,jp)
        icolor(1)=ccolor(5_iknd,0_iknd,jp)
        icolor(2)=ccolor(3_iknd,0_iknd,jp)
        icolor(3)=ccolor(1_iknd,0_iknd,jp)
        icolor(4)=1
        if(lvl<=0) return
c
c       set up input arrays
c
        nsum=0
        jasum=0
        jusum=0
        jfsum=0
        do ii=1,lvl
            i=lvl+1-ii
            n=ka(1,i)
            lenja=ka(10,i)-ka(3,i)
            if(ka(5,i)==ka(3,i)) then
                lenju=lenja
            else
                lenju=ka(7,i)-ka(5,i)
            endif
            if(i<lvl) then
                lenvf=ka(3,i+1)-ka(7,i)
            else
                lenvf=n+1
            endif
            rn(ii)=log10(real(n))
            e(1,ii)=real(lenja-n-1)/real(n)
            e(2,ii)=real(lenju-n-1)/real(n)
            e(3,ii)=real(lenvf-n-1)/real(n)
            nsum=nsum+n
            jasum=jasum+2*(lenja-n-1)+n
            jusum=jusum+2*(lenju-n-1)+n
            jfsum=jfsum+lenvf-n-1
        enddo
        ja0=2*(lenja-n-1)+n
        ratio(1)=real(jasum)/real(ja0)
        ratio(2)=real(jusum)/real(ja0)
        ratio(3)=real(jfsum)/real(ja0)
        ratio(4)=real(nsum)/real(n)
c
        h=0.025e0_rknd
        h2=h/2.0e0_rknd
        xl=3.0e0_rknd*h
        xr=1.0e0_rknd-xl
        yl=xl
        yr=xr
        jmin=0
        jmax=jmin+int(rn(lvl))+1
        numx=jmax+1
        emx=e(1,1)
        emn=emx
        do i=1,lvl
            emx=max(e(1,i),e(2,i),emx)
            emn=min(e(1,i),e(2,i),emn)
        enddo
cc      imin=int(emn)
cc      if(emn<real(imin)) imin=imin-1
        imin=0
        imax=int(emx)
        if(emx>real(imax)) imax=imax+1
        if(imax-imin<4) then
            imin=imin-(4+imin-imax)/2
            imax=imin+4
        endif
        if(imax-imin<=6) then
            numy=imax-imin+1
            iy=1
        else if(imax-imin<=40) then
            imax=imin+((imax-imin-1)/4)*4+4
            numy=(imax-imin)/4+1
            iy=4
        else
            iy=((imax-imin-1)/100+1)*10
            imax=imin+((imax-imin-1)/iy)*iy+iy
            numy=(imax-imin)/iy+1
        endif
c
c       banner
c
        yyl=yr+1.2e0_rknd*h
        yyr=yyl+h
        ym=yyl+h2
        hx=(xr-xl)/3.5e0_rknd
        do j=1,4
            call fstr(ichr,nchr,label(j),0_iknd)
            ichr(nchr+1:nchr+1)=' '
            ii=3
            if(ratio(j)>=10.0e0_rknd) ii=4
            if(ratio(j)>=100.0e0_rknd) ii=5
            call sreal(jchr,mchr,ratio(j),ii,0_iknd)
            ichr(nchr+2:nchr+mchr+1)=jchr(1:mchr)
            nchr=nchr+mchr+1
            xxl=xl+real(j-1)*hx
            xxr=xxl+real(nchr)*h
            xm=xxr+h2
cc          xxm=(xxl+xxr)/2.0e0_rknd
cc          dxm=xxr-xxl
cc          dym=2.0e0_rknd*h
cc          call symbl(xxm,ym,dxm,dym,1_iknd,icolor(j),t)
            if(j/=4) call symbl(xm,ym,h,h,1_iknd,icolor(j),t)
            call htext(xxl,yyl,xxr,yyr,nchr,ichr,0_iknd,q,t,2_iknd)
        enddo
c
c       axis
c
        call xyaxis(xl,xr,yl,yr,h,t,q,numx,jmin,1_iknd,numy,imin,iy)
c
c        graph
c
        dx=(xr-xl)/real(jmax-jmin)
        hx=dx/12.0e0_rknd
        dy=(yr-yl)/real(imax-imin)
        do i=1,lvl
            xs=xl+dx*rn(i)-hx
            do j=1,3
                xm=xs+real(j-1)*hx
                hy=dy*(e(j,i)-real(imin))
                ym=yl+hy/2.0e0_rknd
                call symbl(xm,ym,hx,hy,1_iknd,icolor(j),t)
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
        subroutine prtip(ip,jp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(100) :: ip,ic
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd) :: ccolor
            integer(kind=iknd), dimension(3) :: icolor
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            character(len=15), dimension(300) :: name0
            character(len=15), dimension(100) :: name
            character(len=80) :: ichr
cy
c       print ip array
c
        call linit(t,q)
        mxcolr=jp(17)
        if(mxcolr>=8) then
            icolor(1)=2
            icolor(2)=ccolor(2_iknd,0_iknd,jp)
            icolor(3)=ccolor(6_iknd,0_iknd,jp)
         else
            icolor(1)=2
            icolor(2)=2
            icolor(3)=2
        endif
        call getnam(name0,nlen)
        do i=1,100
            name(i)=' '
            ic(i)=icolor(1)
            call sint(ichr,length,i)
            name(i)(4-length:3)=ichr(1:length)
        enddo
        do i=1,nlen
            if(name0(i)(15:15)=='i') then
                call cint(name0(i),3_iknd,indx,jerr)
                name(indx)(4:10)=name0(i)(4:10)
                if(name0(i)(12:13)==' ') then
                    ic(indx)=icolor(2)
                else
                    ic(indx)=icolor(3)
                endif
            endif
        enddo
c
        size=t(14)
        dy=size/25.0e0_rknd
        dx=(size+0.5e0_rknd)/4.0e0_rknd
        h=min(dy*0.9e0_rknd,dx/20.0e0_rknd)
c
        do i=1,25
            do j=1,4
                k=(j-1)*25+i
                xl=real(j-1,rknd)*dx+.05e0_rknd+dx/10.0e0_rknd
                xr=xl+dx/2.0e0_rknd
                yl=0.95e0_rknd-(real(i,rknd)*dy)
                yr=yl+h
                call htext(xl,yl,xr,yr,10_iknd,name(k),
     +              -1_iknd,q,t,ic(k))
                xl=xl+dx/2.0e0_rknd
                xr=xl+3.0e0_rknd*dx/10.0e0_rknd
                ichr=' '
                call sint(ichr(6:6),nchr,ip(k))
                m=min(nchr,6)
                nchr=max(6,nchr)
                call htext(xl,yl,xr,yr,nchr,ichr(m:m),
     +              1_iknd,q,t,ic(k))
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
        subroutine prtrp(rp,jp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd) :: ccolor
            integer(kind=iknd), dimension(3) :: icolor
            integer(kind=iknd), dimension(100) :: ic
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            real(kind=rknd), dimension(100) :: rp
            character(len=15), dimension(300) :: name0
            character(len=15), dimension(100) :: name
            character(len=80) :: ichr
cy
c       print rp array
c
        call linit(t,q)
        mxcolr=jp(17)
        if(mxcolr>=8) then
            icolor(1)=2
            icolor(2)=ccolor(2_iknd,0_iknd,jp)
            icolor(3)=ccolor(6_iknd,0_iknd,jp)
         else
            icolor(1)=2
            icolor(2)=2
            icolor(3)=2
        endif
        call getnam(name0,nlen)
        do i=1,100
            name(i)=' '
            ic(i)=icolor(1)
            call sint(ichr,length,i)
            name(i)(4-length:3)=ichr(1:length)
        enddo
        do i=1,nlen
            if(name0(i)(15:15)/='r') cycle
            call cint(name0(i),3_iknd,indx,jerr)
            name(indx)(4:10)=name0(i)(4:10)
            if(name0(i)(12:13)==' ') then
                ic(indx)=icolor(2)
            else
                ic(indx)=icolor(3)
            endif
        enddo
c
        size=t(14)
        dy=size/25.0e0_rknd
        dx=(size+0.5e0_rknd)/4.0e0_rknd
        h=min(dy*0.9e0_rknd,dx/20.0e0_rknd)
c
        do i=1,25
            do j=1,4
                k=(j-1)*25+i
                xl=real(j-1,rknd)*dx+.05e0_rknd+dx/22.0e0_rknd
                xr=xl+10.0e0_rknd*dx/22.0e0_rknd
                yl=0.95e0_rknd-(real(i,rknd)*dy)
                yr=yl+h
                call htext(xl,yl,xr,yr,10_iknd,name(k),
     +              -1_iknd,q,t,ic(k))
                xl=xr+dx/22.0e0_rknd
                xr=xl+9.0e0_rknd*dx/22.0e0_rknd
                ichr=' '
                call sreal(ichr(9:9),nchr,rp(k),3_iknd,0_iknd)
                m=min(nchr,9)
                nchr=max(9,nchr)
                call htext(xl,yl,xr,yr,nchr,ichr(m:m),
     +              1_iknd,q,t,ic(k))
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
        subroutine prtsp(sp,jp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd) :: ccolor
            integer(kind=iknd), dimension(3) :: icolor
            integer(kind=iknd), dimension(100) :: ic
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            character(len=15), dimension(300) :: name0
            character(len=15), dimension(100) :: name
            character(len=80) :: ichr
            character(len=80), dimension(100) :: sp
cy
c       print sp array
c
        call linit(t,q)
        mxcolr=jp(17)
        if(mxcolr>=8) then
            icolor(1)=2
            icolor(2)=ccolor(2_iknd,0_iknd,jp)
            icolor(3)=ccolor(6_iknd,0_iknd,jp)
         else
            icolor(1)=2
            icolor(2)=2
            icolor(3)=2
        endif
        call getnam(name0,nlen)
        do i=1,100
            name(i)=' '
            ic(i)=icolor(1)
            call sint(ichr,length,i)
            name(i)(4-length:3)=ichr(1:length)
        enddo
        do i=1,nlen
            isw=1
            if(name0(i)(15:15)=='r') isw=0
            if(name0(i)(15:15)=='i') isw=0
            if(isw==1) then
                call cint(name0(i),3_iknd,indx,jerr)
                name(indx)(4:10)=name0(i)(4:10)
                if(name0(i)(12:13)==' ') then
                    ic(indx)=icolor(2)
                else
                    ic(indx)=icolor(3)
                endif
            endif
        enddo
c
        size=t(14)
        dy=size/25.0e0_rknd
        dx=(size+0.5e0_rknd)/4.0e0_rknd
        h=min(dy*0.9e0_rknd,dx/20.0e0_rknd)
c
        do i=1,25
            do j=1,2
                k=(j-1)*25+i
                xl=real(j-1,rknd)*dx*2.0e0_rknd
     +              +.05e0_rknd+dx/10.0e0_rknd
                xr=xl+dx/2.0e0_rknd
                yl=0.95e0_rknd-(real(i,rknd)*dy)
                yr=yl+h
                call htext(xl,yl,xr,yr,10_iknd,name(k),
     +              -1_iknd,q,t,ic(k))
                xl=xl+dx/2.0e0_rknd
                xr=xl+1.5e0_rknd*dx
                call fstr(ichr,nchr,sp(k),0_iknd)
                if(nchr>0) call htext(xl,yl,xr,yr,nchr,ichr,
     +              -1_iknd,q,t,ic(k))
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
        subroutine prtka(lvl,ka,jp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(10,*) :: ka
            integer(kind=iknd), dimension(8) :: list,total
            integer(kind=iknd), dimension(25) :: jp
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            character(len=80) :: ichr
            character(len=80), save, dimension(8) :: name
cy
c
            data name/'level','  n  ','lenja','lena ','lenju','lenu ',
     +          'lenvf','lenwt'/
c
c       print ip array
c
        call linit(t,q)
        dy=0.9e0/25.0e0_rknd
        dx=1.4e0/8.0e0_rknd
        h=min(dy*0.6e0_rknd,dx/8.0e0_rknd)
c
        xl=.05e0_rknd
        xr=xl+dx
        yl=0.95e0_rknd-dy
        yr=yl+h
        do k=1,8
           total(k)=0
           call fstr(ichr,nchr,name(k),0_iknd)
           call htext(xl,yl,xr,yr,nchr,ichr,0_iknd,q,t,2_iknd)
           xl=xr
           xr=xr+dx
        enddo
        do i=1,lvl
            xl=.05e0_rknd
            xr=xl+dx
            yl=0.95e0_rknd-(real(i+1)*dy)
            yr=yl+h
            list(1)=lvl+1-i
            list(2)=ka(1,i)
            list(3)=ka(10,i)-ka(3,i)
            list(4)=ka(6,i)-ka(4,i)
            if(ka(5,i)==ka(3,i)) then
                list(5)=list(3)
            else
                list(5)=ka(7,i)-ka(5,i)
            endif
            list(6)=ka(8,i)-ka(6,i)
            if(ka(6,i)==ka(4,i)) then
                list(4)=list(6)
            endif
            kk=6
            if(i<lvl) then
                list(7)=ka(3,i+1)-ka(7,i)
                list(8)=ka(4,i+1)-ka(8,i)
                kk=8
            else
                list(7)=0
                list(8)=0
            endif
            do k=2,8
                total(k)=total(k)+list(k)
            enddo
            do k=1,kk
                call sint(ichr,nchr,list(k))
                call htext(xl,yl,xr,yr,nchr,ichr,0_iknd,q,t,2_iknd)
                xl=xr
                xr=xr+dx
            enddo
        enddo
c
c       totals
c
        xl=.05e0_rknd+2.0e0_rknd*dx
        xr=xl+dx
        yl=0.95e0_rknd-(real(lvl+2)*dy)
        yr=yl+h
        icolor=jp(18)
        do k=3,8
            call sint(ichr,nchr,total(k))
            call htext(xl,yl,xr,yr,nchr,ichr,0_iknd,q,t,icolor)
            xl=xr
            xr=xr+dx
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
        subroutine hbplt(hist,numhst,lab,jp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(22,3) :: icolor
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd), dimension(3) :: num,num1
            integer(kind=iknd), dimension(2,3) :: jc
            integer(kind=iknd) :: ccolor
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            real(kind=rknd), dimension(22,*) :: hist
            real(kind=rknd), dimension(22,3) :: e
            real(kind=rknd), dimension(25) :: x,y,z
            real(kind=rknd), dimension(3) :: ave,sfact
            character(len=80) :: ichr
            character(len=80), save, dimension(2) :: rate
cy
            data rate/'multigraph iteration','singular vector'/
c
c     graph error
c
        call linit(t,q)
        size=t(14)
        xshift=t(15)-size/2.0e0_rknd
        yshift=t(16)-size/2.0e0_rknd
        zshift=t(5)
        mxhist=20
        t(1)=xshift
        t(2)=yshift
        t(3)=size
        jc(1,1)=ccolor(4_iknd,0_iknd,jp)
        jc(2,1)=ccolor(6_iknd,0_iknd,jp)
        jc(1,2)=ccolor(2_iknd,0_iknd,jp)
        jc(2,2)=ccolor(5_iknd,0_iknd,jp)
        jc(1,3)=ccolor(3_iknd,0_iknd,jp)
        jc(2,3)=ccolor(1_iknd,0_iknd,jp)
        sfact(1)=0.8e0_rknd/sqrt(2.0e0_rknd)
        sfact(2)=0.8e0_rknd
        sfact(3)=0.6e0_rknd
        do j=1,3
            num(j)=0
            num1(j)=0
        enddo
c
        do j=1,numhst
            num1(j)=int(hist(mxhist+1,j))
            num(j)=min(num1(j),mxhist)
            if(num(j)>0) then
                e1=abs(hist(mxhist+2,j))
                if(e1>0.0e0_rknd) e1=1.0e0_rknd/e1
                do i=1,num(j)
                    qq=abs(hist(i,j))*e1
                    e(i,j)=0.0e0_rknd
                    if(qq>0.0e0_rknd) e(i,j)=log10(qq)
                    ee=e(i,j)
                    if(hist(i,j)>=0.0e0_rknd) then
                        icolor(i,j)=jc(1,j)
                    else
                        icolor(i,j)=jc(2,j)
                    endif
                enddo
                ave(j)=10.0e0_rknd**(e(num(j),j)/real(num1(j)))
            endif
        enddo
        n1max=num1(1)
        n1min=num1(1)
        do j=1,numhst
            if(num1(j)>0) then
                n1max=max(num1(j),n1max)
                n1min=min(num1(j),n1min)
            endif
        enddo
        if(n1max==0) return
c
        h=0.025e0_rknd
        h2=h/2.0e0_rknd
        xl=3.0e0_rknd*h
        xr=1.0e0_rknd-xl
        yl=xl
        yr=xr
        if(n1max-n1min+4<=mxhist) then
            jmin=max(n1max-mxhist,0)
c*          jmax=jmin+max(((n1max-jmin-1)/4)*4+4,8)
            jmax=jmin+mxhist
        else
            jmin=max(n1min-4,0)
            jmax=jmin+((n1max-jmin-1)/4)*4+4
        endif
        if(jmax-jmin==8) then
            numx=5
            is=2
        else if(jmax-jmin<=40) then
            numx=(jmax-jmin)/4+1
            is=4
        else
            jmax=jmin+((n1max-jmin-1)/10)*10+10
            numx=(jmax-jmin)/10+1
            is=10
        endif
        emx=ee
        emn=ee
        do j=1,numhst
            if(num(j)>0) then
                do i=1,num(j)
                    emx=max(e(i,j),emx)
                    emn=min(e(i,j),emn)
                enddo
            endif
        enddo
        imin=int(emn)
        if(emn<real(imin)) imin=imin-1
        imax=int(emx)
        if(emx>real(imax)) imax=imax+1
        if(imax-imin<4) then
            imin=imin-(4+imin-imax)/2
            imax=imin+4
        endif
        numy=imax-imin+1
c
c       banners
c
        yyl=yr+1.8e0_rknd*h
        yyr=yyl+h
        ym=yyl+h2
        call fstr(ichr,nchr,rate(lab),0_iknd)
        xxl=2.0e0_rknd*h
        xxl=h
        xxr=xxl+20.0e0_rknd*h
        call htext(xxl,yyl,xxr,yyr,nchr,ichr,-1_iknd,q,t,2_iknd)
        do j=1,numhst
            if(num(j)>0) then
                call sreal(ichr,nchr,ave(j),2_iknd,0_iknd)
                xxl=real(10+8*j)*h
                xxr=xxl+7.0e0_rknd*h
                xm=xxl-h
                call htext(xxl,yyl,xxr,yyr,nchr,ichr,-1_iknd,q,t,2_iknd)
                call symbl(xm,ym,h,h,j,jc(1,j),t)
            endif
        enddo
c
c       axis
c
        call xyaxis(xl,xr,yl,yr,h,t,q,numx,jmin,
     +      is,numy,imin,1_iknd)
c
c        graph
c
        dy=(yr-yl)/real(numy-1)
        dx=(xr-xl)/real(jmax-jmin)
        do j=1,numhst
            hh=h*sfact(j)
            ishift=max(num1(j)-mxhist,0)-jmin
            i0=max(1,-ishift)
            do i=i0,num(j)
                xs=xl+dx*real(i+ishift)
                ys=yl+dy*(e(i,j)-real(imin))
                x(i)=xs*size+xshift
                y(i)=ys*size+yshift
                z(i)=zshift
                call symbl(xs,ys,hh,hh,j,icolor(i,j),t)
            enddo
            nn=num(j)-i0+1
            if(nn>1) call pline(x(i0),y(i0),z(i0),nn,2_iknd)
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
        subroutine pieplt(time,jp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd), dimension(2) :: icolor
            integer(kind=iknd) :: ccolor
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            real(kind=rknd), dimension(*) :: time
            real(kind=rknd), dimension(21) :: th,dt
            real(kind=rknd), dimension(90) :: x,y,z
            character(len=80) :: ichr
            character(len=80), save, dimension(2) :: label
cy
            data   label/'mginit','mg'/
c
c     graph times
c
        call linit(t,q)
        size=t(14)
        xshift=t(15)-size/2.0e0_rknd
        yshift=t(16)-size/2.0e0_rknd
        t(1)=xshift
        t(2)=yshift
        t(3)=size
        zshift=t(5)
        scale=t(3)
c
        num=2
        icolor(1)=ccolor(6_iknd,0_iknd,jp)
        icolor(2)=ccolor(2_iknd,0_iknd,jp)
c
c       set up input arrays
c
        tot=0.0e0_rknd
        do i=1,num
            tot=tot+time(i)
        enddo
        if(tot<=0.0e0_rknd) return
        pi=3.141592653589793e0_rknd
        th(1)=pi/2.0e0_rknd
        do i=1,num
            fr=time(i)/tot
            dt(i)=fr*2.0e0_rknd*pi
            th(i+1)=th(i)+dt(i)
        enddo
c
c       make pie chart
c
        xcen=0.5e0_rknd
        ycen=0.35e0_rknd
        rad=0.3e0_rknd
        dd=pi/32.0e0_rknd
        do i=1,num
            m=int(dt(i)/dd)
            x(1)=xcen*scale+xshift
            y(1)=ycen*scale+yshift
            z(1)=zshift
            dtheta=dt(i)/real(m+1)
            theta=th(i)
            do j=1,m+2
                ang=theta+dtheta*real(j-1)
                xx=xcen+rad*cos(ang)
                yy=ycen+rad*sin(ang)
                x(j+1)=xx*scale+xshift
                y(j+1)=yy*scale+yshift
                z(j+1)=zshift
            enddo
            x(m+4)=x(1)
            y(m+4)=y(1)
            z(m+4)=z(1)
            call pfill(x,y,z,m+3_iknd,icolor(i))
            call pline(x,y,z,m+4_iknd,2_iknd)
        enddo
c
        xl=t(15)-size/2.0e0_rknd
        yt=t(16)+size/2.0e0_rknd
        h=1.0e0_rknd/20.0e0_rknd
        h2=h/2.0e0_rknd
        do i=1,num
            xxl=xl
            yyl=yt-h*real(i)*2.0e0_rknd
            xm=xxl+h
            ym=yyl+h2
            call symbl(xm,ym,h,h,1_iknd,icolor(i),t)
c
            call fstr(ichr,nchr,label(i),0_iknd)
            xxl=xl+3.0e0_rknd*h
            xxr=xxl+real(nchr)*h
            yyr=yyl+h
            call htext(xxl,yyl,xxr,yyr,nchr,ichr,-1_iknd,q,t,2_iknd)
            ii=3
            if(time(i)>=10.0e0_rknd) ii=4
            if(time(i)>=100.0e0_rknd) ii=5
            call sreal(ichr,nchr,time(i),ii,0_iknd)
            xxl=xl+8.0e0_rknd*h
            xxr=xl+14.0e0_rknd*h
            call htext(xxl,yyl,xxr,yyr,nchr,ichr,1_iknd,q,t,2_iknd)
            fr=time(i)/tot*100.0e0_rknd
            call sreal(ichr,nchr,fr,3_iknd,0_iknd)
            xxl=xl+15.0e0_rknd*h
            xxr=xl+19.0e0_rknd*h
            call htext(xxl,yyl,xxr,yyr,nchr,ichr,1_iknd,q,t,2_iknd)
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
        subroutine xyaxis(xl,xr,yl,yr,h,t,q,numx,iminx,
     +      incx,numy,iminy,incy)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            real(kind=rknd), dimension(2) :: x,y,z
            character(len=80) :: ichr
cy
        dx=(xr-xl)/(numx-1)
        dy=(yr-yl)/(numy-1)
        xshift=t(1)
        yshift=t(2)
        zshift=t(5)
        scale=t(3)
        h2=h/2.0e0_rknd
c
c    x - axis
c
        z(1)=zshift
        z(2)=zshift
        x(1)=xl*scale+xshift
        y(1)=yl*scale+yshift
        x(2)=xr*scale+xshift
        y(2)=y(1)
        call pline(x,y,z,2_iknd,2_iknd)
        do i=1,numx
            k=iminx+(i-1)*incx
            call sint(ichr,nchr,k)
            xx=xl+real(i-1,rknd)*dx
            x(1)=xx*scale+xshift
            y(1)=(yl+h2)*scale+yshift
            x(2)=x(1)
            y(2)=yl*scale+yshift
            call pline(x,y,z,2_iknd,2_iknd)
            xxl=xx-real(nchr,rknd)*h2
            xxr=xx+real(nchr,rknd)*h2
            yyl=yl-2.25e0_rknd*h
            yyr=yyl+h
            call htext(xxl,yyl,xxr,yyr,nchr,ichr,0_iknd,q,t,2_iknd)
        enddo
c
c    y-axis
c
        x(1)=xl*scale+xshift
        y(1)=yl*scale+yshift
        x(2)=x(1)
        y(2)=yr*scale+yshift
        call pline(x,y,z,2_iknd,2_iknd)
        do i=1,numy
            k=iminy+(i-1)*incy
            call sint(ichr,nchr,k)
            yy=yl+real(i-1,rknd)*dy
            x(1)=(xl+h2)*scale+xshift
            y(1)=yy*scale+yshift
            x(2)=xl*scale+xshift
            y(2)=y(1)
            call pline(x,y,z,2_iknd,2_iknd)
            xxl=xl-real(2*nchr+1,rknd)*h2
            xxr=xl-h2
            yyl=yy-h2
            yyr=yy+h2
            call htext(xxl,yyl,xxr,yyr,nchr,ichr,0_iknd,q,t,2_iknd)
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
        subroutine symbl(xm,ym,hx,hy,itype,icolor,t)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), save, dimension(5) :: iptr
            integer(kind=iknd), save, dimension(14) :: px,py
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(5) :: x,y,z
cy
            data iptr/1,5,9,12,15/
            data px/-1,1,1,-1,1,0,-1,0,-1,1,0,0,1,-1/
            data py/-1,-1,1,1,0,1,0,-1,-1,-1,1,-1,1,1/
c
c       itype = 1 box   itype = 2 diamond  itype = 3,4 triangle
c
        xshift=t(1)
        yshift=t(2)
        scale=t(3)
        zshift=t(5)
        istart=iptr(itype)
        num=iptr(itype+1)-istart
        do i=1,num
            px0=real(px(i+istart-1),rknd)/2.0e0_rknd
            py0=real(py(i+istart-1),rknd)/2.0e0_rknd
            x(i)=(xm+hx*px0)*scale+xshift
            y(i)=(ym+hy*py0)*scale+yshift
            z(i)=zshift
        enddo
        x(num+1)=x(1)
        y(num+1)=y(1)
        z(num+1)=z(1)
        call pfill(x,y,z,num,icolor)
        call pline(x,y,z,num+1_iknd,2_iknd)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine htext(xl,yl,xr,yr,nchr,cchr,ijust,q,t,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(80) :: ichr
            integer(kind=iknd), save, dimension(640) :: symbcd
            integer(kind=iknd), save, dimension(94) :: istart
            integer(kind=iknd), save, dimension(128) :: map
            real(kind=rknd), save, dimension(94) :: width
            real(kind=rknd), dimension(2) :: x,y,z
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            character(len=*) :: cchr
            character(len=1) :: cc
cy
c       writes text given in cchr array in the rectangle defined by its
c       lower left corner of world coordinates xl,yl and its upper right
c       corner of world coordinates xr,yr.
c
c       ijust=-1  for justification on the left
c       ijust= 0  for centered text
c       ijust=+1  for justification on the right
c
c        the symbol numbers are
c        1-26   upper case roman simplex
c       27-52   lower case roman simplex
c       53-62   simplex numbers
c       63-78   symbols + - ( ) , . / = * $ < > { } @ ^
c       79-94   symbols [ ] # : ; ! ? % & ~ " ' _ \ | `
c
c
c       symbol parameters taken from n.m.wolcott, fortran iv enhanced
c       character graphics, nbs
c       ichr(j) contains the symbol number of the jth symbol
c       everything outside this range is considered a space
c
        data (symbcd(i),i=1,60)/
     +  443556555,443557579,432612882,        0,433070987,433071584,
     1  323987166,328083226,325854871,317404054,317400725,325723922,
     2  327657165,323364299,298156032,462268125,321889760,309339231,
     3  300852123,296493907,298329038,304489675,317040204,325527312,
     4          0,433070987,433071456,319792797,325953304,327788240,
     5  323429900,312845195,        0,433070987,433071840,432743830,
     6  432383691,        0,433070987,433071840,432743830,        0,
     7  462268125,321889760,309339231,300852123,296493907,298329038,
     8  304489675,317040204,325527312,327792083,327778304,433070987,
     9  462432011,432744214,        0,433070987,        0,449848720/
        data (symbcd(i),i=61,120)/
     +  312911116,306553867,298197837,294134546,        0,433070987,
     1  462431122,443262731,        0,433070987,432383627,        0,
     2  433070987,433071499,466625931,466626443,        0,433070987,
     3  433071883,462432011,        0,443556959,300852123,296493907,
     4  298329038,304489675,317040204,325527312,329885528,328050397,
     5  321889760,309329920,433070987,433071584,323987166,328083225,
     6  325822102,317367189,        0,443556959,300852123,296493907,
     7  298329038,304489675,317040204,325527312,329885528,328050397,
     8  321889760,309343631,327450624,433070987,433071584,323987166,
     9  328083226,325854871,317399958,447424267,        0,460236383/
        data (symbcd(i),i=121,180)/
     +  315630752,300917597,296592281,300688471,317367892,323593937,
     1  325527116,314942603,300294990,        0,441459851,426780256,
     2          0,433070993,300360780,310748555,321267406,327722784,
     3          0,426779851,460334283,        0,428876875,449848395,
     4  449849035,470820555,        0,430974667,460333899,        0,
     5  426779862,308655840,309002240,460333899,430974688,430286539,
     6          0,455910987,455812568,313304217,302785430,296330065,
     7  298263564,306554187,317072974,        0,433070987,432743448,
     8  307012953,317466198,323593873,321332684,312845451,302392206,
     9          0,455812568,313304217,302785430,296330065,298263564/
        data (symbcd(i),i=181,240)/
     +  306554187,317072974,        0,456140363,455812568,313304217,
     1  302785430,296330065,298263564,306554187,317072974,        0,
     2  430548563,321562135,317465945,307012632,298525523,296264590,
     3  302392459,312845772,321323008,445654176,303014876,300266265,
     4  309100544,455910985,318973381,312616068,302167638,317465945,
     5  307012632,298525523,296264590,302392459,312845772,321323008,
     6  433070987,432710744,309110169,319563349,321224704,430973855,
     7  300950433,296760217,298156032,435168287,305144865,300954649,
     8  302261189,295838404,        0,433070987,453813135,441034315,
     9          0,433070987,        0,432841611,432710744,309110169/
        data (symbcd(i),i=241,300)/
     +  319563349,321238613,327952281,338471128,344631563,        0,
     1  432841611,432710744,309110169,319563349,321224704,441230360,
     2  298525523,296264590,302392459,312845772,321332881,323593814,
     3  317465945,307003392,432841604,432743448,307012953,317466198,
     4  323593873,321332684,312845451,302392206,        0,455910980,
     5  455812568,313304217,302785430,296330065,298263564,306554187,
     6  317072974,        0,432841611,432645078,304882905,315392000,
     7  453715416,311207001,298591062,298460179,313075153,319268366,
     8  317072651,304456588,296157184,435168207,302392459,310752025,
     9  309100544,432841615,300295243,310748556,321369689,321224704/
        data (symbcd(i),i=301,360)/
     +  428647563,453813387,        0,430744651,447521867,447522379,
     1  464299595,        0,430745099,453813067,        0,428647563,
     2  453813387,302228357,293741252,        0,453813067,430745113,
     3  430286347,        0,443556895,298722135,296362895,302392523,
     4  312845836,323462868,325822108,319792480,309329920,437134493,
     5  313533771,        0,432907164,300885023,307242400,319792734,
     6  323888794,321660373,296068811,        0,435168928,311174616,
     7  321627798,325691089,323429900,312845451,300295053,296189952,
     8  451945298,327759328,317030400,456139744,298558424,307012953,
     9  319563414,325691089,323429900,312845451,300295053,296189952/
        data (symbcd(i),i=361,420)/
     +  458139231,315630880,305112028,298558354,300360780,310748491,
     1  319170190,325625554,323659287,313271576,304849877,298385408,
     2  460334155,430974688,        0,441459679,298754971,300721240,
     3  313239062,323626706,325559949,321267083,306553804,298230607,
     4  296297364,302720215,317466201,323856029,321889696,307232768,
     5  458008150,317334803,308913172,298525529,296559517,303015136,
     6  311436767,321824409,323626575,317072651,306553804,298254336,
     7  451847627,432678932,        0,432678932,        0,447882466,
     8  305112027,298525586,300328009,308487492,        0,431104994,
     9  305112283,311108882,308716617,300098372,        0,436609995/
        data (symbcd(i),i=421,480)/
     +  298197965,302392330,300163975,        0,434545548,300262412,
     1  300318720,466756356,        0,432777239,432580625,        0,
     2  441263246,430679505,451650385,        0,441590919,449979783,
     3  460236383,315630752,300917597,296592281,300688471,317367892,
     4  323593937,325527116,314942603,300294990,        0,466527124,
     5  331710464,432973716,298156032,443688035,303113184,300885020,
     6  304981145,306947093,439460897,303015005,307111130,309077142,
     7  298460306,308815054,306586699,302294023,304264211,306750607,
     8  304522252,300229576,302195781,308412416,435299427,307307744,
     9  309273756,304981017,302752917,439461025,307209309,302916570/
        data (symbcd(i),i=481,540)/
     +  300688406,311043090,300426190,302392395,306488455,304264339,
     1  302556175,304522380,308618440,306390085,300023808,462169818,
     2  321758619,311239897,306914451,308847952,319301265,325694875,
     3  311207126,308913425,313014043,325691089,329787344,338241685,
     4  340502618,336471966,328181344,315630815,305079260,298656599,
     5  296362897,300393549,308684171,321234700,331786190,464365331,
     6  327722832,        0,426321109,325661394,309012178,        0,
     7  433202052,435299268,433202532,432153924,        0,443688132,
     8  445785348,431105316,430056708,        0,447751044,460334340,
     9  432711445,430417615,        0,434938776,300655640,300725197/
        data (symbcd(i),i=541,600)/
     +  298197963,302392269,        0,434938776,300655640,300725195,
     1  298197965,302392330,300163975,        0,435168158,300491806,
     2  300954590,300692429,298197963,302392269,        0,432939995,
     3  298656603,296625054,300917856,311436767,319759964,321725976,
     4  317433045,308884768,315598302,319694362,317465942,442934412,
     5  308651276,308707328,468722507,441459998,311305434,304915417,
     6  296592221,298820640,307242271,317662878,330278880,459875921,
     7  319268365,323331851,331753422,333981522,325648384,468461463,
     8  334178327,336340953,332179288,327886481,319235468,310748235,
     9  298197838,296264595,311141785,317564381,315598112,307209309/
        data (symbcd(i),i=601,640)/
     +  304981144,311076430,325461899,333817868,335983691,300295054,
     1  298361811,304788571,307013262,327559051,        0,430482259,
     2  298525719,306947350,319399570,327755667,334148435,298492950,
     3  306914581,319366801,327722898,334145495,        0,435168153,
     4  437265305,451945881,454043033,        0,443557017,445654169,
     5          0,432351242,        0,429008772,        0,439493700,
     6          0,430973849,428876697,        0/
c
        data istart/
     +     1,   5,  16,  26,  34,  39,  43,  54,  58,  60,  66,  70,
     1    73,  78,  82,  93, 100, 112, 120, 131, 134, 140, 143, 148,
     2   151, 154, 158, 167, 176, 184, 193, 202, 206, 217, 222, 226,
     3   232, 236, 238, 247, 252, 261, 270, 279, 283, 292, 296, 301,
     4   304, 309, 312, 317, 321, 330, 333, 341, 349, 352, 361, 373,
     5   376, 391, 403, 406, 408, 414, 420, 425, 428, 430, 433, 437,
     6   450, 452, 454, 473, 492, 519, 523, 528, 533, 538, 544, 551,
     7   558, 573, 588, 612, 624, 629, 632, 634, 636, 638/
c
        data (width(i),i=1,40)/
     +  18.0e0_rknd,21.0e0_rknd,21.0e0_rknd,21.0e0_rknd,
     1  19.0e0_rknd,18.0e0_rknd,21.0e0_rknd,22.0e0_rknd,
     2   8.0e0_rknd,16.0e0_rknd,21.0e0_rknd,17.0e0_rknd,
     3  24.0e0_rknd,22.0e0_rknd,22.0e0_rknd,21.0e0_rknd,
     4  22.0e0_rknd,21.0e0_rknd,20.0e0_rknd,16.0e0_rknd,
     5  22.0e0_rknd,18.0e0_rknd,24.0e0_rknd,20.0e0_rknd,
     6  18.0e0_rknd,20.0e0_rknd,19.0e0_rknd,19.0e0_rknd,
     7  18.0e0_rknd,19.0e0_rknd,18.0e0_rknd,12.0e0_rknd,
     8  19.0e0_rknd,19.0e0_rknd, 8.0e0_rknd,10.0e0_rknd,
     9  17.0e0_rknd, 8.0e0_rknd,30.0e0_rknd,19.0e0_rknd/
        data (width(i),i=41,80)/
     +  19.0e0_rknd,19.0e0_rknd,19.0e0_rknd,13.0e0_rknd,
     1  17.0e0_rknd,12.0e0_rknd,19.0e0_rknd,16.0e0_rknd,
     2  22.0e0_rknd,17.0e0_rknd,16.0e0_rknd,17.0e0_rknd,
     3  20.0e0_rknd,20.0e0_rknd,20.0e0_rknd,20.0e0_rknd,
     4  20.0e0_rknd,20.0e0_rknd,20.0e0_rknd,20.0e0_rknd,
     5  20.0e0_rknd,20.0e0_rknd,26.0e0_rknd,26.0e0_rknd,
     6  14.0e0_rknd,14.0e0_rknd,10.0e0_rknd,10.0e0_rknd,
     7  22.0e0_rknd,26.0e0_rknd,16.0e0_rknd,20.0e0_rknd,
     8  24.0e0_rknd,24.0e0_rknd,14.0e0_rknd,14.0e0_rknd,
     9  27.0e0_rknd,22.0e0_rknd,14.0e0_rknd,14.0e0_rknd/
        data (width(i),i=81,94)/
     +  21.0e0_rknd,10.0e0_rknd,10.0e0_rknd,10.0e0_rknd,
     1  18.0e0_rknd,24.0e0_rknd,25.0e0_rknd,24.0e0_rknd,
     2  16.0e0_rknd, 8.0e0_rknd,26.0e0_rknd,22.0e0_rknd,
     3  14.0e0_rknd, 8.0e0_rknd/
c
            data map/
     +       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     1       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     2       0,84,89,81,72,86,87,90,65,66,71,63,67,64,68,69,
     3      53,54,55,56,57,58,59,60,61,62,82,83,73,70,74,85,
     4      77, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,
     5      16,17,18,19,20,21,22,23,24,25,26,79,92,80,78,91,
     6      94,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,
     7      42,43,44,45,46,47,48,49,50,51,52,75,93,76,88, 0/
c
c       ixtrct gets nbits from iword starting at the nstart
c       bit from the right
c
        ixtrct(nstart,nbits,iword)=mod(iword/(2**(nstart-nbits)),
     +      2**nbits)+((1-sign(1_iknd,iword))/2)*
     1      (2**nbits-min(1,mod(-iword,2**(nstart-nbits))))
c
        if(nchr<=0) return
        if(xl>=xr) return
        if(yl>=yr) return
c
        do i=1,nchr
            cc=cchr(i:i)
            ii=ichar(cc)
            ichr(i)=map(ii+1)
        enddo
        dx=xr-xl
        dy=yr-yl
c
c       find width of strings to be plotted
c
        wid=0.0e0_rknd
        do i=1,nchr
            ic=ichr(i)
            if(ic<1.or.ic>94) then
                wid=wid+20.0e0_rknd
            else
                wid=wid+width(ic)
            endif
        enddo
        wid=wid/21.0e0_rknd
c
        height=min(dx/wid,dy)
        if(height<dy) then
            x0=xl
            y0=yl+(dy-height)/2.0e0_rknd
        else
c
c       justification
c
            y0=yl
            if(ijust==-1) then
                x0=xl
            else if(ijust==0) then
                x0=xl+(dx-wid*height)/2.0e0_rknd
            else if(ijust==1) then
                x0=xr-wid*height
            endif
        endif
c
        scale=t(3)
        xshift=t(1)
        yshift=t(2)
        zshift=t(5)
c
        rscale=height/21.0e0_rknd
        xi=x0
        yi=y0
c
        do i=1,nchr
            ic=ichr(i)
            if(ic<=0.or.ic>94) then
c
c        plot a space
c
                xi=xi+20.0e0_rknd*rscale
            else
c
c       plot a single symbol
c
                is=istart(ic)
                ib=30
   70           ipen=ixtrct(ib,3_iknd,symbcd(is))
                if(ipen==0)then
                    xi=xi+rscale*width(ic)
                    cycle
                endif
                ix=ixtrct(ib-3_iknd,6_iknd,symbcd(is))
                iy=ixtrct(ib-9_iknd,6_iknd,symbcd(is))
                xx=xi+(ix-10)*rscale
                yy=yi+(iy-11)*rscale
                xm=xx*q(1,1)+yy*q(2,1)
                ym=xx*q(1,2)+yy*q(2,2)
                zm=xx*q(1,3)+yy*q(2,3)
                xx=xm*scale+xshift
                yy=ym*scale+yshift
                zz=zm*scale+zshift
                if(ipen==2) then
                    x(2)=xx
                    y(2)=yy
                    z(2)=zz
                    call lwindw(x,y,z,2_iknd,t,icolor)
                endif
                x(1)=xx
                y(1)=yy
                z(1)=zz
                ib=45-ib
                if(ib==30)is=is+1
                go to 70
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
        subroutine pwindw(x,y,z,llen,t,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: x,y,z
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(22) :: xn,yn,zn,x0,y0,z0
            real(kind=rknd), dimension(4) :: cx,cy,cc
cy
c       map a polygon onto the current window
c
        rmag=t(12)
        if(rmag<=1.0e0_rknd) then
            call pfill(x,y,z,llen,icolor)
            return
        endif
c
        nmax=22
        eps=t(7)/rmag
        shift=(1.0e0_rknd-t(14))/2.0e0_rknd
c
        cx(1)=1.0e0_rknd
        cx(2)=-cx(1)
        cx(3)=0.0e0_rknd
        cx(4)=cx(3)
c
        cy(1)=cx(3)
        cy(2)=cx(4)
        cy(3)=cx(1)
        cy(4)=cx(2)
c
        cc(1)=-t(8)
        cc(2)=t(9)
        cc(3)=-t(10)
        cc(4)=t(11)
c
        do i=1,llen
            xn(i)=x(i)
            yn(i)=y(i)
            zn(i)=z(i)
        enddo
        num=llen
c
        do k=1,4
            len=num
            num=0
            do i=1,len
                x0(i)=xn(i)
                y0(i)=yn(i)
                z0(i)=zn(i)
            enddo
            do i=1,len
                si=x0(i)*cx(k)+y0(i)*cy(k)+cc(k)
                if(si>=eps) then
                    num=num+1
                    xn(num)=x0(i)
                    yn(num)=y0(i)
                    zn(num)=z0(i)
                else
                    ibef=i-1
                    if(i==1) ibef=len
                    iaft=i+1
                    if(i==len) iaft=1
                    j=ibef
                    do jj=1,2
                        s=x0(j)*cx(k)+y0(j)*cy(k)+cc(k)
                        if(s>eps) then
                            num=num+1
                            f=s/(s-si)
                            xn(num)=x0(i)*f+x0(j)*(1.0e0_rknd-f)
                            yn(num)=y0(i)*f+y0(j)*(1.0e0_rknd-f)
                            zn(num)=z0(i)*f+z0(j)*(1.0e0_rknd-f)
                        endif
                        j=iaft
                    enddo
                endif
            enddo
            if(num<=2) return
            if(num>=nmax-2) stop 7577
        enddo
        do i=1,num
            xn(i)=(xn(i)+cc(1))*rmag+shift
            yn(i)=(yn(i)+cc(3))*rmag+shift
cc          zn(i)=zn(i)*rmag
        enddo
        call pfill(xn,yn,zn,num,icolor)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine lwindw(x,y,z,n,t,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: x,y,z
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(2) :: xx,yy,zz
cy
c       draw the part of the picture within the current window
c
        rmag=t(12)
        if(rmag<=1.0e0_rknd) then
            call pline(x,y,z,n,icolor)
            return
        endif
c
        xl=t(8)
        xr=t(9)
        yb=t(10)
        yt=t(11)
        shift=(1.0e0_rknd-t(14))/2.0e0_rknd
c
c       the main loop
c
        do i=2,n
            xx(1)=x(i-1)
            yy(1)=y(i-1)
            zz(1)=z(i-1)
            xx(2)=x(i)
            yy(2)=y(i)
            zz(2)=z(i)
c
c       fit line into window in x direction
c
            jl=1
            if(xx(2)<xx(1)) jl=2
            jr=3-jl
            if(xx(jr)<=xl.or.xx(jl)>=xr) cycle
c
            if(xx(jl)<xl) then
               f=(xx(jr)-xl)/(xx(jr)-xx(jl))
               xx(jl)=xl
               yy(jl)=yy(jl)*f+yy(jr)*(1.0e0_rknd-f)
               zz(jl)=zz(jl)*f+zz(jr)*(1.0e0_rknd-f)
            endif
c
            if(xx(jr)>xr) then
                f=(xr-xx(jl))/(xx(jr)-xx(jl))
                xx(jr)=xr
                yy(jr)=yy(jr)*f+yy(jl)*(1.0e0_rknd-f)
                zz(jr)=zz(jr)*f+zz(jl)*(1.0e0_rknd-f)
            endif
c
c       fit line into window in y direction
c
            jb=1
            if(yy(2)<yy(1)) jb=2
            jt=3-jb
            if(yy(jt)<=yb.or.yy(jb)>=yt) cycle
c
            if(yy(jb)<yb) then
                f=(yy(jt)-yb)/(yy(jt)-yy(jb))
                yy(jb)=yb
                xx(jb)=xx(jb)*f+xx(jt)*(1.0e0_rknd-f)
                zz(jb)=zz(jb)*f+zz(jt)*(1.0e0_rknd-f)
            endif
c
            if(yy(jt)>yt) then
                f=(yt-yy(jb))/(yy(jt)-yy(jb))
                yy(jt)=yt
                xx(jt)=xx(jt)*f+xx(jb)*(1.0e0_rknd-f)
                zz(jt)=zz(jt)*f+zz(jb)*(1.0e0_rknd-f)
            endif
c
c       rescale and then draw
c
            do j=1,2
                xx(j)=(xx(j)-xl)*rmag+shift
                yy(j)=(yy(j)-yb)*rmag+shift
cc              zz(j)=zz(j)*rmag
            enddo
            call pline(xx,yy,zz,2_iknd,icolor)
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
        subroutine title0(title,isw)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            character(len=80) :: title,ichr
cy
c       draw the title for the picture
c
        call linit(t,q)
        size=t(14)
        xl=t(15)-size/2.0e0_rknd
        xr=t(15)+size/2.0e0_rknd
        if(isw==1) xr=xr+0.5e0_rknd
        yb=t(16)+size/2.0e0_rknd
        yt=t(16)+t(3)/2.0e0_rknd
        yl=yb+(yt-yb)*0.25e0_rknd
        yr=yb+(yt-yb)*0.75e0_rknd
        call fstr(ichr,nchr,title,0_iknd)
        call htext(xl,yl,xr,yr,nchr,ichr,0_iknd,q,t,2_iknd)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine mtxplt(ip,rp,sp,ja,a,ka)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            integer(kind=iknd), dimension(10,*) :: ka
            integer(kind=iknd), dimension(100) :: ip
            integer(kind=iknd), dimension(25) :: jp,jpl
            real(kind=rknd), dimension(*) :: a
            real(kind=rknd), dimension(100) :: rp
            real(kind=rknd), dimension(25) :: t,tl
            real(kind=rknd), dimension(3,3) :: q,ql
            real(kind=rknd), dimension(256) :: red,green,blue
            character(len=80), dimension(100) :: sp
            real(kind=rknd), allocatable, dimension(:) :: z
            integer(kind=iknd), allocatable, dimension(:) :: jz,color
cy
c       user specified ip variables
c
        ispd=ip(8)
        imtxsw=abs(ip(55))
        if(imtxsw<=0.or.imtxsw>6) imtxsw=1
        if(ip(55)<0) then
            ip(55)=-imtxsw
        else
            ip(55)=imtxsw
        endif
        level=ip(67)
        lvl=ip(2)
        if(level<=0.or.level>lvl) level=lvl
        call getptr(level,lvl,nf,nptr,japtr,iaptr,juptr,
     +      iuptr,jvptr,ivptr,iqptr,ibptr,nc,ncptr,ka)
c
c       error flags
c
        ip(26)=0
c
        lenja=ja(japtr+nf)-1
        lenju=ja(juptr+nf)-1
        if(ispd==1) then
            lena=lenja
            lenu=lenju
        else
            lena=2*lenja-(nf+1)
            lenu=2*lenju-(nf+1)
        endif
        maxjz=max(lenja,lenju)
        maxz=max(lena,lenu)
        allocate(jz(maxjz),color(maxz),z(maxz))
        if(imtxsw>=5) then
            call sferr(nf,ispd,jz,z,ja(japtr),a(iaptr),
     +          ja(juptr),a(iuptr),maxjz,maxz,iflag)
            if(iflag/=0) then
                ip(26)=20
                go to 10
            endif
        else
            if(imtxsw<=2) then
                nnz=ja(juptr+nf)-ja(juptr)
                ll=nf+1+nnz
                if(ispd/=1) ll=ll+nnz
                do i=1,ll
                    z(i)=a(iuptr+i-1)
                enddo
            else
                nnz=ja(japtr+nf)-ja(japtr)
                ll=nf+1+nnz
                if(ispd/=1) ll=ll+nnz
                do i=1,ll
                    z(i)=a(iaptr+i-1)
                enddo
            endif
        endif
c
c
        call linit(t,q)
        call linit(tl,ql)
        call minit(ip,rp,nf,ja(japtr),ja(juptr),jz,z,
     +      color,jp,jpl,t,tl,q,ql)
        jtype=imtxsw-(imtxsw/2)*2
c
        call clrmap(red,green,blue,jp)
c
        call pltutl(jp(18),red,green,blue)
c
c       main plot
c
        call pframe(4_iknd)
        call title0(sp(4),0_iknd)
        call pframe(-4_iknd)
        call pframe(5_iknd)
        if(imtxsw>=5) then
            call mplot1(jp,t,q,jz,z,color)
        else if(imtxsw>=3) then
            call mplot1(jp,t,q,ja(japtr),z,color)
        else
            call mplot1(jp,t,q,ja(juptr),z,color)
        endif
        call pframe(-5_iknd)
c
c       legend plot
c
        call pframe(2_iknd)
        if(jtype==1) then
            call legnd5(jp,t)
        else
            call legnd4(jp,tl,z)
        endif
        call pframe(-2_iknd)
c
c       small plot
c
        call pframe(3_iknd)
        if(imtxsw>=5) then
            call mplot1(jpl,tl,ql,jz,z,color)
        else if(imtxsw>=3) then
            call mplot1(jpl,tl,ql,ja(japtr),z,color)
        else
            call mplot1(jpl,tl,ql,ja(juptr),z,color)
        endif
        call legnd0(t)
        call pframe(-3_iknd)
c
        call pltutl(-1_iknd,red,green,blue)
        ip(26)=0
   10   deallocate(z,jz,color)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine sferr(n,ispd,je,e,ja,a,ju,u,maxje,maxe,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ju,je
            integer(kind=iknd), dimension(n) :: mark,list,indx
            integer(kind=iknd) :: amtx,umtx,emtx
            real(kind=rknd), dimension(*) :: a,u,e
            real(kind=rknd), dimension(n) :: tl,tu
cy
c
c       compute sparsity structure for error matrix
c
        if(ispd/=1) then
            lenje=min(maxje,(maxe+n+1)/2)
            amtx=ja(n+1)-ja(1)
            umtx=ju(n+1)-ju(1)
            emtx=lenje-(n+1)
        else
            lenje=min(maxje,maxe)
            amtx=0
            umtx=0
            emtx=0
        endif
c
c
c
        je(1)=n+2
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
                    uinv=1.0e0/u(k)
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
            tu(i)=tu(i)-u(i)
            tl(i)=tu(i)-u(i)
            do jj=ju(i),ju(i+1)-1
                j=ju(jj)
                tu(j)=0.0e0_rknd
                tl(j)=0.0e0_rknd
c               tu(j)=tu(j)-u(jj)
c               tl(j)=tl(j)-u(jj+umtx)
            enddo
c
c       make je for this row
c
            next=je(i)
            do j=1,len
                k=mark(i)
                tt=max(abs(tl(k)),abs(tu(k)))
                if(tt>0.0e0_rknd) then
                    if(next<lenje) then
                        je(next)=k
                        next=next+1
                    else
                        iflag=i
                        return
                    endif
                endif
                mark(i)=mark(k)
                mark(k)=0
            enddo
            mark(i)=0
            je(i+1)=next
            len=next-je(i)
            if(len>1) call ihp(je(je(i)),len)
c
c       move tl, tu to e
c
            e(i)=-tu(i)
            do jj=je(i),je(i+1)-1
                j=je(jj)
                e(jj)=-tu(j)
                e(jj+emtx)=-tl(j)
            enddo
c
            if(ju(i)<ju(i+1)) then
                j=ju(ju(i))
                list(i)=list(j)
                list(j)=i
                indx(i)=ju(i)
            endif
        enddo
        iflag=0
c
c       shift u for non symmetric case
c
        maxje=je(n+1)-1
        if(ispd/=1) then
            nnz=maxje-(n+1)
            do i=1,nnz
                e(maxje+i)=e(lenje+i)
            enddo
            maxe=maxje+nnz
        else
            maxe=maxje
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
        subroutine mplot1(jp,t,q,ju,u,color)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ju,color
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd) :: ccolor
            real(kind=rknd), dimension(*) :: u
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
cy
        ispd=jp(7)
        n=jp(2)
        isw=1
c
        lshift=0
        if(ispd/=1) lshift=ju(n+1)-ju(1)
        if(q(1,3)+q(2,3)>=0.0e0_rknd) then
c
c       border for lower triangle
c
            call mtxbrd(t,q,0_iknd)
c
c       lower triangle
c
            if(q(1,3)<0.0e0_rknd) then
                n1=n
                n2=1
                ns=-1
            else
                n1=1
                n2=n
                ns=1
            endif
            if(q(2,3)<0.0e0_rknd) then
                i1=0
                i2=1
                ks=1
            else
                i1=1
                i2=0
                ks=-1
            endif
            do i=n1,n2,ns
                do k=ju(i+i1)-i1,ju(i+i2)-i2,ks
                    icolor=ccolor(color(k+lshift),0_iknd,jp)
                    ix=i
                    call centry(ix,ju(k),u(k+lshift),t,q,jp,icolor,isw)
                enddo
            enddo
c
c       diagonal
c
            do i=n1,n2,ns
                icolor=ccolor(color(i),0_iknd,jp)
                ix=i
                iy=i
                call centry(ix,iy,u(i),t,q,jp,icolor,isw)
            enddo
c
c       upper triangle
c
            if(q(2,3)>0.0e0_rknd) then
                n1=n
                n2=1
                ns=-1
            else
                n1=1
                n2=n
                ns=1
            endif
            if(q(1,3)>0.0e0_rknd) then
                i1=0
                i2=1
                ks=1
            else
                i1=1
                i2=0
                ks=-1
            endif
            do i=n1,n2,ns
                do k=ju(i+i1)-i1,ju(i+i2)-i2,ks
                    icolor=ccolor(color(k),0_iknd,jp)
                    iy=i
                    call centry(ju(k),iy,u(k),t,q,jp,icolor,isw)
                enddo
            enddo
c
c       border for upper triangle
c
            call mtxbrd(t,q,1_iknd)
        else
c
c       border for upper triangle
c
            call mtxbrd(t,q,1_iknd)
c
c       upper triangle
c
            if(q(2,3)>0.0e0_rknd) then
                n1=n
                n2=1
                ns=-1
            else
                n1=1
                n2=n
                ns=1
            endif
            if(q(1,3)>0.0e0_rknd) then
                i1=0
                i2=1
                ks=1
            else
                i1=1
                i2=0
                ks=-1
            endif
            do i=n1,n2,ns
                do k=ju(i+i1)-i1,ju(i+i2)-i2,ks
                    icolor=ccolor(color(k),0_iknd,jp)
                    iy=i
                    call centry(ju(k),iy,u(k),t,q,jp,icolor,isw)
                enddo
            enddo
c
c       diagonal
c
            do i=n1,n2,ns
                icolor=ccolor(color(i),0_iknd,jp)
                ix=i
                iy=i
                call centry(ix,iy,u(i),t,q,jp,icolor,isw)
            enddo
c
c       lower triangle
c
            if(q(1,3)<0.0e0_rknd) then
                n1=n
                n2=1
                ns=-1
            else
                n1=1
                n2=n
                ns=1
            endif
            if(q(2,3)<0.0e0_rknd) then
                i1=0
                i2=1
                ks=1
            else
                i1=1
                i2=0
                ks=-1
            endif
            do i=n1,n2,ns
                do k=ju(i+i1)-i1,ju(i+i2)-i2,ks
                    icolor=ccolor(color(k+lshift),0_iknd,jp)
                    ix=i
                    call centry(ix,ju(k),u(k+lshift),t,q,jp,icolor,isw)
                enddo
            enddo
c
c       border for lower triangle
c
            call mtxbrd(t,q,0_iknd)
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
        subroutine centry(ix,iy,val,t,q,jp,icolor,isw)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd), save, dimension(6) :: order
            integer(kind=iknd), save, dimension(3,3) :: index
            integer(kind=iknd), save, dimension(4,6) :: face
            integer(kind=iknd), save  :: n,istrt
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
            real(kind=rknd), save :: h,hl,hr,dz
            real(kind=rknd), save, dimension(8) :: px,py,pz
            real(kind=rknd), dimension(5) :: xn,yn,zn
            character(len=80) :: ichr
cy
        data px/0.0e0_rknd,1.0e0_rknd,1.0e0_rknd,0.0e0_rknd,
     +      0.0e0_rknd,1.0e0_rknd,1.0e0_rknd,0.0e0_rknd/
        data py/0.0e0_rknd,0.0e0_rknd,1.0e0_rknd,1.0e0_rknd,
     +      0.0e0_rknd,0.0e0_rknd,1.0e0_rknd,1.0e0_rknd/
        data pz/0.0e0_rknd,0.0e0_rknd,0.0e0_rknd,0.0e0_rknd,
     +      1.0e0_rknd,1.0e0_rknd,1.0e0_rknd,1.0e0_rknd/
        data face/4,1,5,8,2,3,7,6,1,2,6,5,3,4,8,7,4,3,2,1,5,6,7,8/
        data index/1,2,3,2,3,1,3,1,2/
c
        if(isw==1) then
            isw=0
            n=jp(2)
            h=1.0e0_rknd/real(n)
            hl=h/10.0e0_rknd
            hr=h-hl
c
c       compute order
c
            kmin=1
            if(abs(q(kmin,3))>abs(q(2,3))) kmin=2
            if(abs(q(kmin,3))>abs(q(3,3))) kmin=3
            kmid=index(2,kmin)
            kmax=index(3,kmin)
            if(abs(q(kmid,3))>abs(q(kmax,3))) kmid=kmax
            kmax=6-kmin-kmid
c
            if(q(kmax,3)>0.0e0_rknd) then
                order(1)=2*kmax-1
                order(6)=2*kmax
            else
                order(6)=2*kmax-1
                order(1)=2*kmax
            endif
            if(q(kmid,3)>0.0e0_rknd) then
                order(2)=2*kmid-1
                order(5)=2*kmid
            else
                order(5)=2*kmid-1
                order(2)=2*kmid
            endif
            if(q(kmin,3)>0.0e0_rknd) then
                order(3)=2*kmin-1
                order(4)=2*kmin
            else
                order(4)=2*kmin-1
                order(3)=2*kmin
            endif
c
            tol=1.e-3_rknd
            istrt=6
            if(abs(q(kmin,3))>tol) then
                istrt=4
            else if(abs(q(kmid,3))>tol) then
                istrt=5
            endif
cc          istrt=1
            zmin=t(24)
            zmax=t(25)
            if(zmax>zmin) then
                dz=1.0e0_rknd/(zmax-zmin)
            else
                dz=0.0e0_rknd
            endif
        endif
c
        lines=jp(20)
        numbrs=jp(21)
        i3d=jp(22)
        xshift=t(1)
        yshift=t(2)
        zshift=t(5)
        scale=t(3)
        zl=t(23)
        zmin=t(24)
        zmax=t(25)
c
        x=real(ix-1)*h
        y=real(n-iy)*h
        if(i3d/=0) then
            zz=(val-zmin)*dz
            if(zz>zl) then
                z=zz
            else
                z=zl
                zl=zz
            endif
        else
            z=zl
        endif
        do i=istrt,6
            ii=order(i)
            do j=1,4
                xx=x+h*px(face(j,ii))
                yy=y+h*py(face(j,ii))
                zz=zl+(z-zl)*pz(face(j,ii))
                xn(j)=(xx*q(1,1)+yy*q(2,1))*scale+xshift
                yn(j)=(xx*q(1,2)+yy*q(2,2)+zz*q(3,2))*scale+yshift
                zn(j)=(xx*q(1,3)+yy*q(2,3)+zz*q(3,3))*scale+zshift
            enddo
            xn(5)=xn(1)
            yn(5)=yn(1)
            zn(5)=zn(1)
            call pwindw(xn,yn,zn,4_iknd,t,icolor)
            if(lines==-2) call lwindw(xn,yn,zn,5_iknd,t,2_iknd)
        enddo
c
c
c
        if(numbrs>=0) return
        if(numbrs==-1) then
            call sreal(ichr,nn,val,3_iknd,1_iknd)
        else if(numbrs==-2) then
            ichr(1:1)='('
            call sint(ichr(2:2),iylen,iy)
            ichr(iylen+2:iylen+2)=','
            call sint(ichr(iylen+3:iylen+3),ixlen,ix)
            nn=3+ixlen+iylen
            ichr(nn:nn)=')'
        endif
        xl=x+hl
        xr=x+hr
        yb=y+hl
        yt=y+hr
        call htext(xl,yb,xr,yt,nn,ichr,0_iknd,q,t,2_iknd)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine mtxbrd(t,q,isw)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3) :: x,y,z
            real(kind=rknd), dimension(3,3) :: q
cy
        xshift=t(1)
        yshift=t(2)
        zshift=t(5)
        scale=t(3)
        zl=t(23)
c
c       border for lower triangle
c
        x(1)=q(2,1)*scale+xshift
        y(1)=(q(2,2)+zl*q(3,2))*scale+yshift
        z(1)=(q(2,3)+zl*q(3,3))*scale+zshift
        if(isw==0) then
            x(2)=xshift
            y(2)=zl*q(3,2)*scale+yshift
            z(2)=zl*q(3,3)*scale+zshift
        else
c
c       border for upper triangle
c
            x(2)=(q(1,1)+q(2,1))*scale+xshift
            y(2)=(q(1,2)+q(2,2)+zl*q(3,2))*scale+yshift
            z(2)=(q(1,3)+q(2,3)+zl*q(3,3))*scale+zshift
        endif
        x(3)=q(1,1)*scale+xshift
        y(3)=(q(1,2)+zl*q(3,2))*scale+yshift
        z(3)=(q(1,3)+zl*q(3,3))*scale+zshift
c
        call lwindw(x,y,z,3_iknd,t,2_iknd)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine mtxclr(ja,ju,je,jp,t,u,color)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ju,color,je
            integer(kind=iknd), dimension(25) :: jp
            real(kind=rknd), dimension(*) :: u
            real(kind=rknd), dimension(25) :: t
cy
c       compute types
c
        ispd=jp(7)
        n=jp(2)
        imtxsw=abs(jp(3))
        ncolor=jp(5)
        iscale=jp(19)
        if(imtxsw>=5) then
            len=je(n+1)-1
            icolor=6
        else if(imtxsw>=3) then
            len=ja(n+1)-1
            icolor=4
        else
            len=ju(n+1)-1
            icolor=2
        endif
        if(ispd==1) then
            lenu=len
            lshift=0
        else
            lenu=2*len-(n+1)
            lshift=len-(n+1)
        endif
c
        ity=imtxsw-(imtxsw/2)*2
        if(ity==0) go to 10
c
c       color by type
c
c       type = 2 fillin       (blue)
c       type = 4 original     (green)
c       type = 5 diagonal     (yellow)
c       type = 6 neglected    (red)
c
        do i=1,lenu
            color(i)=icolor
        enddo
c
        if(imtxsw>=5) then
            do i=1,n
                color(i)=5
                do j=ju(i),ju(i+1)-1
                    call jamap(i,ju(j),ij,ji,je,lshift)
                    color(ij)=2
                    color(ji)=2
                enddo
                do j=ja(i),ja(i+1)-1
                    call jamap(i,ja(j),ij,ji,je,lshift)
                    color(ij)=4
                    color(ji)=4
                enddo
            enddo
        else if(imtxsw>=3) then
            do i=1,n
                color(i)=5
            enddo
        else
            do i=1,n
                color(i)=5
                do j=ja(i),ja(i+1)-1
                    call jamap(i,ja(j),ij,ji,ju,lshift)
                    color(ij)=4
                    color(ji)=4
                enddo
c
            enddo
        endif
        return
c
c
   10   umin=t(19)
        umax=t(20)
        zmin=fscale(umin,iscale,0_iknd)
        zmax=fscale(umax,iscale,0_iknd)
        eps=t(7)
        if(zmax>zmin) then
            zscale=(1.0e0_rknd-eps)*real(ncolor)/(zmax-zmin)
        else
            zscale=0.0e0_rknd
        endif
c
        do i=1,lenu
            zz=(fscale(u(i),iscale,0_iknd)-zmin)*zscale
            color(i)=max(0,int(zz))+1
            color(i)=min(color(i),ncolor)
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
        subroutine minit(ip,rp,n,ja,ju,je,u,color,jp,jpl,t,tl,q,ql)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(100) :: ip
            integer(kind=iknd), dimension(25) :: jp,jpl
            integer(kind=iknd), dimension(*) :: ja,ju,color,je
            real(kind=rknd), dimension(3,3) :: q,ql
            real(kind=rknd), dimension(100) :: rp
            real(kind=rknd), dimension(25) :: t,tl
            real(kind=rknd), dimension(*) :: u
cy
c       initialize for mtxplt
c
        do i=1,25
            jp(i)=0
        enddo
        call zoombx(rp,t)
        rmag=t(12)
c
c       check control parameters in ip and rp
c
        ispd=ip(8)
        iscale=ip(58)
        if(iscale<0.or.iscale>2) iscale=0
        lines=ip(59)
        if(lines/=-2) lines=0
        numbrs=ip(60)
        if(numbrs>0.or.numbrs<-2) numbrs=0
        mxcolr=max(2,ip(51))
        imtxsw=abs(ip(55))
        jtype=imtxsw-(imtxsw/2)*2
        if(jtype==1) then
            ncolor=6
        else
            ncon=ip(56)
            ncolor=max(1,ncon)
        endif
        lenja=ja(n+1)-1
        lenju=ju(n+1)-1
        lenje=je(n+1)-1
        if(imtxsw>=5) then
            len=lenje
        else if(imtxsw>=3) then
            len=lenja
        else
            len=lenju
        endif
        if(ispd/=1) len=2*len-(n+1)
        u(n+1)=u(1)
        if(ip(55)>=0) then
            do i=1,len
                u(i)=abs(u(i))
            enddo
        endif
        umin=u(1)
        umax=u(1)
        do i=1,len
            umin=min(umin,u(i))
            umax=max(umax,u(i))
        enddo
        if(iscale==1.and.umin<=0.0e0_rknd) iscale=2
c
c       set up rotated coordinate system
c
        nx=ip(64)
        ny=ip(65)
        nz=ip(66)
        i3d=1
        if(numbrs/=0) i3d=0
cc      if(nx==0.and.ny==0) i3d=0
c
        call mkrot(nx,ny,nz,q)
c
        xmin=min(0.0e0_rknd,q(1,1))+min(0.0e0_rknd,q(2,1))
        xmax=max(0.0e0_rknd,q(1,1))+max(0.0e0_rknd,q(2,1))
        ymin=min(0.0e0_rknd,q(1,2))+min(0.0e0_rknd,q(2,2))
        ymax=max(0.0e0_rknd,q(1,2))+max(0.0e0_rknd,q(2,2))
        zmin=min(0.0e0_rknd,q(1,3))+min(0.0e0_rknd,q(2,3))
        zmax=max(0.0e0_rknd,q(1,3))+max(0.0e0_rknd,q(2,3))
        if(i3d==1) then
            ymax=ymax+q(3,2)
            zmin=zmin+min(0.0e0_rknd,q(3,3))
            zmax=zmax+max(0.0e0_rknd,q(3,3))
        endif
        size=t(14)
        xs=t(15)
        ys=t(16)
        zs=t(17)
        scale=size/max(xmax-xmin,ymax-ymin)
        xshift=xs-scale*(xmax+xmin)/2.0e0_rknd
        yshift=ys-scale*(ymax+ymin)/2.0e0_rknd
        zshift=zs-scale*(zmax+zmin)/2.0e0_rknd
c
c       set up jp
c
        jp(1)=len
        jp(2)=n
        jp(3)=imtxsw
        jp(4)=1
        jp(5)=ncolor
        if(jtype==1) then
            jp(6)=0
        else
            jp(6)=1
        endif
        jp(7)=ispd
c
        jp(13)=ip(64)
        jp(14)=ip(65)
        jp(15)=ip(66)
        jp(16)=0
c
        jp(17)=mxcolr
        jp(18)=min(ncolor+2,mxcolr)
        jp(19)=iscale
        jp(20)=lines
        jp(21)=numbrs
        jp(22)=i3d
c
        t(1)=xshift
        t(2)=yshift
        t(3)=scale
        t(5)=zshift
c
        if(rp(8)<rp(9)) then
            t(19)=rp(8)
            t(20)=rp(9)
        else
            t(19)=umin
            t(20)=umax
        endif
        if(i3d==1.and.umin<min(0.0e0_rknd,umax)) then
            t(23)=-umin/(umax-umin)
        else
            t(23)=0.0e0_rknd
        endif
        t(24)=umin
        t(25)=umax
c
c       parameters for legend plot
c
        do i=1,25
            tl(i)=t(i)
            jpl(i)=jp(i)
        enddo
        tl(12)=1.0e0_rknd
c
        jpl(20)=0
        jpl(21)=0
        if(rmag<=1.0e0_rknd) jpl(22)=0
c
c       set q0, scale,xshift, yshift correctly for picture
c
        if(rmag/=1.0e0_rknd) then
            do i=1,3
                do j=1,3
                    ql(i,j)=q(i,j)
                enddo
            enddo
        else
            tl(1)=xs-size/2.0e0_rknd
            tl(2)=ys-size/2.0e0_rknd
            tl(5)=zs
            tl(3)=size
            jpl(13)=0
            jpl(14)=0
            jpl(15)=1
        endif
c
c       set colors
c
        call mtxclr(ja,ju,je,jp,t,u,color)
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
        subroutine legnd5(jp,t)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd), save, dimension(4) ::mcic
            integer(kind=iknd) :: ccolor
            real(kind=rknd), dimension(5) :: x,y,z
            real(kind=rknd), dimension(25) :: t,tt
            real(kind=rknd), dimension(3,3) :: q
            character(len=80) :: ichr
            character(len=80), save :: title
            character(len=80), save, dimension(4) :: label
            character(len=80), save, dimension(2) :: mnmx
c
            data title/'element types'/
            data label/'diagonal','original','fillin','neglected'/
            data mnmx/'min','max'/
            data mcic/5,4,2,6/
cy
        call linit(tt,q)
        size=tt(14)
        zshift=tt(5)
        xs=tt(15)
        ys=tt(16)
c
        xl=xs-size/2.0e0_rknd
        xr=xs+size/2.0e0_rknd
        yb=ys-size/2.0e0_rknd
        yt=ys+size/2.0e0_rknd
        dx=(xr-xl)/14.5e0_rknd
        dy=(yt-yb)/7.0e0_rknd
        h=min(0.9e0_rknd*dy,dx)
c
        call fstr(ichr,nchr,title,0_iknd)
        xxl=xl+2.25e0_rknd*dx
        xxr=xxl+15.0e0_rknd*dx
        yyl=yt-dy
        yyr=yyl+h
        call htext(xxl,yyl,xxr,yyr,15_iknd,ichr,-1_iknd,q,tt,2_iknd)
c
        do i=1,4
            yy=yt-real(i+1)*dy
c
c       square icon
c
            x(1)=xl+0.25e0_rknd*dx
            x(2)=x(1)+h
            x(3)=x(2)
            x(4)=x(1)
            x(5)=x(1)
            y(1)=yy
            y(2)=yy
            y(3)=yy+h
            y(4)=y(3)
            y(5)=y(1)
            do j=1,5
                z(j)=zshift
            enddo
            ii=ccolor(mcic(i),0_iknd,jp)
            call pfill(x,y,z,4_iknd,ii)
            call pline(x,y,z,5_iknd,2_iknd)
c
c     label
c
            call fstr(ichr,nchr,label(i),0_iknd)
            xxl=xl+2.25e0_rknd*dx
            xxr=xxl+15.0e0_rknd*dx
            yyl=yy
            yyr=yyl+h
            call htext(xxl,yyl,xxr,yyr,15_iknd,ichr,-1_iknd,q,tt,2_iknd)
        enddo
c
c       min-max values
c
        do i=1,2
            yy=yt-real(i+5)*dy
c
c     label
c
            call fstr(ichr,nchr,mnmx(i),0_iknd)
            xxl=xl+0.25e0_rknd*dx
            xxr=xxl+3.0e0_rknd*dx
            yyl=yy
            yyr=yyl+h
            call htext(xxl,yyl,xxr,yyr,3_iknd,ichr,-1_iknd,q,tt,2_iknd)
c
c      value
c
            ichr=' '
            call sreal(ichr(4:4),nchr,t(23+i),3_iknd,1_iknd)
            if(nchr<7) then
                ii=nchr-3
                nchr=7
            else
                ii=4
            endif
            xxl=xl+4.25e0_rknd*dx
            xxr=xr
            call htext(xxl,yyl,xxr,yyr,nchr,ichr(ii:ii),
     +          1_iknd,q,tt,2_iknd)
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
        function fscale(f,iscale,invrse)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
cy
c       set scaling function
c
        if(iscale==0) then
c
c       linear scale
c
            fscale=f
            return
        else if(iscale==1) then
c
c       log scale
c
            if(invrse==0) then
                fscale=log(f)
                return
            else
                fscale=exp(f)
                return
            endif
        else
c
c       arcsinh scale
c
            if(invrse==0) then
                af=abs(f)
                if(af<1.0e0_rknd) then
                    q=sqrt(1.0e0_rknd+f*f)+af
                    fx=log(q)
                    fscale=fx+(af-sinh(fx))/cosh(fx)
                else
                    q=1.0e0_rknd/f
                    q=sqrt(1.0e0_rknd+q*q)+1.0e0_rknd
                    fscale=log(q)+log(af)
                endif
                if(f<0.0e0_rknd) fscale=-fscale
                return
            else
                fscale=sinh(f)
                return
            endif
        endif
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine clrmap(red,green,blue,jp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd), save, dimension(7) :: ir,ig,ib
            real(kind=rknd), dimension(*) :: red,green,blue
cy
            data ir/1,1,1,0,0,0,1/
            data ig/0,0,1,1,1,0,0/
            data ib/1,0,0,0,1,1,1/
c
c       set up a color map
c
        ncolor=jp(5)
        icplt=jp(4)
        nshade=jp(16)
        mxcolr=jp(17)
        maplen=jp(18)
        gamma=0.7e0_rknd
        theta=1.0e0_rknd
c
c       background color (white)
c
        red(1)=1.0e0_rknd
        green(1)=1.0e0_rknd
        blue(1)=1.0e0_rknd
c
c       line-drawing color (black)
c
        red(2)=0.0e0_rknd
        green(2)=0.0e0_rknd
        blue(2)=0.0e0_rknd
c
        if(maplen<=2) return
c
        if(ncolor>=mxcolr-2) then
            jcolor=mxcolr-2
        else
            jcolor=ncolor
        endif
c
c       the primary set of colors
c
        red(3)=real(ir(7),rknd)
        green(3)=real(ig(7),rknd)
        blue(3)=real(ib(7),rknd)
        if(jcolor==1) go to 20
        if(icplt/=0) then
            h=5.0e0_rknd/real(jcolor-1,rknd)
        else
            h=6.0e0_rknd/real(jcolor,rknd)
        endif
        do ii=2,jcolor
            i=ii+2
            x=6.0e0_rknd-h*real(ii-1,rknd)
            k=1+int(x)
            dl=real(k,rknd)-x
            dr=1.0e0_rknd-dl
            red(i)=dl*real(ir(k),rknd)+dr*real(ir(k+1),rknd)
            red(i)=max(0.0e0_rknd,red(i))**gamma
            red(i)=min(1.0e0_rknd,red(i))
            green(i)=dl*real(ig(k),rknd)+dr*real(ig(k+1),rknd)
            green(i)=max(0.0e0_rknd,green(i))**gamma
            green(i)=min(1.0e0_rknd,green(i))
            blue(i)=dl*real(ib(k),rknd)+dr*real(ib(k+1),rknd)
            blue(i)=max(0.0e0_rknd,blue(i))**gamma
            blue(i)=min(1.0e0_rknd,blue(i))
        enddo
c
c       shading
c
   20   if(nshade==0) return
        if(icplt/=0) then
            bmax=0.5e0_rknd/real(nshade,rknd)
            wmax=0.5e0_rknd/real(nshade,rknd)
        else
            bmax=0.45e0_rknd/real(nshade,rknd)
            wmax=0.75e0_rknd/real(nshade,rknd)
        endif
        do j=1,nshade
            jplus=j*ncolor+2
            jminus=jplus+nshade*ncolor
            fb=(1.0e0_rknd-real(j,rknd)*bmax)**theta
            fw=(1.0e0_rknd-real(j,rknd)*wmax)**theta
            w=1.0e0_rknd-fw
            do i=1,ncolor
                k=i+jplus
                red(k)=red(i+2)*fw+w
                red(k)=max(red(k),0.0e0_rknd)
                red(k)=min(red(k),1.0e0_rknd)
                green(k)=green(i+2)*fw+w
                green(k)=max(green(k),0.0e0_rknd)
                green(k)=min(green(k),1.0e0_rknd)
                blue(k)=blue(i+2)*fw+w
                blue(k)=max(blue(k),0.0e0_rknd)
                blue(k)=min(blue(k),1.0e0_rknd)
                k=i+jminus
                red(k)=red(i+2)*fb
                red(k)=max(red(k),0.0e0_rknd)
                red(k)=min(red(k),1.0e0_rknd)
                green(k)=green(i+2)*fb
                green(k)=max(green(k),0.0e0_rknd)
                green(k)=min(green(k),1.0e0_rknd)
                blue(k)=blue(i+2)*fb
                blue(k)=max(blue(k),0.0e0_rknd)
                blue(k)=min(blue(k),1.0e0_rknd)
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
        function ccolor(icolor,ishade,jp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd) :: ccolor
cy
c       compute the color index
c
        ncolor=jp(5)
        nshade=jp(16)
        mxcolr=jp(17)
        if(icolor<=0.or.icolor>ncolor
     +      .or.abs(ishade)>nshade) then
            ccolor=1
        else if(ishade==0) then
            ccolor=icolor+2-((icolor-1)/(mxcolr-1))*(mxcolr-1)
            if(ccolor>mxcolr) ccolor=1
        else
            if(ishade>0) then
                ccolor=icolor+2+ncolor*ishade
            else
                ccolor=icolor+2+ncolor*(nshade-ishade)
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
        subroutine legnd0(t)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(6) :: x,y,z
            real(kind=rknd), dimension(25) :: t,tt
            real(kind=rknd), dimension(3,3) :: q
cy
c       helps locating current window (draw boundary in small window)
c
        call linit(tt,q)
        zshift=tt(5)
        scale=tt(3)
        size=tt(14)
        dd=(scale+size)/4.0e0_rknd
        x0=tt(15)-dd
        x1=tt(15)+dd
        y0=tt(16)-dd
        y1=tt(16)+dd
c
c       mark magnified area
c
        do i=1,6
            z(i)=zshift
        enddo
        if(t(12)>1.0e0_rknd) then
            xl=max(x0,t(8))
            xr=min(x1,t(9))
            yb=max(y0,t(10))
            yt=min(y1,t(11))
c
c           mark the box in the window
c
            x(1)=(xl+xr)/2.0e0_rknd
            x(2)=x(1)
            y(1)=y0
            y(2)=yb
            call pline(x,y,z,2_iknd,2_iknd)
            y(1)=yt
            y(2)=y1
            call pline(x,y,z,2_iknd,2_iknd)
            x(1)=x0
            x(2)=xl
            y(1)=(yb+yt)/2.0e0_rknd
            y(2)=y(1)
            call pline(x,y,z,2_iknd,2_iknd)
            x(1)=xr
            x(2)=x1
            call pline(x,y,z,2_iknd,2_iknd)
            x(1)=xl
            y(1)=yb
            x(2)=xr
            y(2)=y(1)
            x(3)=x(2)
            y(3)=yt
            x(4)=x(1)
            y(4)=y(3)
            x(5)=x(1)
            y(5)=y(1)
            call pline(x,y,z,5_iknd,2_iknd)
        endif
c
        x(1)=x0
        y(1)=y0
        x(2)=x1
        y(2)=y(1)
        x(3)=x(2)
        y(3)=y1
        x(4)=x(1)
        y(4)=y(3)
        x(5)=x(1)
        y(5)=y(1)
        call pline(x,y,z,5_iknd,2_iknd)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine legnd4(jp,t,e)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(20) :: nchr
            integer(kind=iknd), dimension(25) :: jp
            integer(kind=iknd), dimension(22) :: kd
            integer(kind=iknd) :: ccolor
            real(kind=rknd), dimension(12) :: f
            real(kind=rknd), dimension(44) :: x,y,z
            real(kind=rknd), dimension(25) :: t,tt
            real(kind=rknd), dimension(3,3) :: qq
            real(kind=rknd), dimension(*) :: e
            character(len=80), dimension(15) :: ichr
cy
        call linit(tt,qq)
        size=tt(14)
        zshift=tt(5)
        xs=tt(15)
        ys=tt(16)
c
        ntf=jp(1)
        ierrsw=jp(6)
        icolor=jp(5)
        if(icolor<=0) return
        ncolor=min(icolor,11)
        iscale=jp(19)
c
c       set function values
c
        zmin=fscale(t(19),iscale,0_iknd)
        zmax=fscale(t(20),iscale,0_iknd)
        df=(zmax-zmin)/real(ncolor)
        do i=1,ncolor+1
            zz=zmin+df*real(i-1)
            f(i)=fscale(zz,iscale,1_iknd)
        enddo
c
c       make boxes for each color
c
        xf=xs
        xi=xf-size*0.45e0_rknd
        xc=xf+0.04e0_rknd*size
        xx=xc+0.4e0_rknd*size
        yi=ys-size*0.45e0_rknd
        yf=ys+size*0.45e0_rknd
        yinc=0.04e0_rknd*size
        tic=0.02e0_rknd*size
        if(icolor==ncolor) yf=yi+(yf-yi)*ncolor/11.0e0_rknd
c
        do i=1,5
            z(i)=zshift
        enddo
        x(1)=xi
        x(2)=xf
        x(3)=xf
        x(4)=xi
        x(5)=xi
        dy=(yf-yi)/real(icolor)
        do i=1,icolor
            y(1)=yi+dy*real(i)
            y(2)=y(1)
            y(3)=yi+dy*real(i-1)
            y(4)=y(3)
            ii=ccolor(i,0_iknd,jp)
            call pfill(x,y,z,4_iknd,ii)
        enddo
c
c       draw the border and tick marks
c
        y(1)=yi
        y(2)=yi
        y(3)=yf
        y(4)=yf
        y(5)=yi
        call pline(x,y,z,5_iknd,2_iknd)
c
c
        x(1)=xf
        scale=(yf-yi)/real(ncolor)
        do i=0,ncolor
            yp=yi+scale*i
            x(2)=xf+tic
            y(1)=yp
            y(2)=yp
            call pline(x,y,z,2_iknd,2_iknd)
        enddo
c
c       compute error distribution
c
        if(ierrsw==1.and.df/=0.0e0_rknd) then
            num=2*ncolor
            do i=1,num
                kd(i)=0
            enddo
            dd=2.0e0_rknd/df
            do i=1,ntf
                ff=(fscale(e(i),iscale,0_iknd)-zmin)*dd
                iq=max(1,int(ff)+1)
                iq=min(num,iq)
                kd(iq)=kd(iq)+1
            enddo
            kdm=0
            do i=1,num
                kdm=max(kdm,kd(i))
            enddo
            ddy=(yf-yi)/real(num)
            xxi=xi+0.05e0_rknd*(xf-xi)
            ddx=0.9e0_rknd*(xf-xi)
            do i=1,num
                j=2*i-1
                x(j)=xxi+ddx*(real(kd(i))/real(kdm))
                x(j+1)=x(j)
                y(j+1)=yi+ddy*real(i)
                y(j)=yi+ddy*real(i-1)
                z(j)=zshift
                z(j+1)=zshift
            enddo
            num=2*num
            call pline(x,y,z,num,2_iknd)
        endif
c
c       label the tick marks
c
        mxchr=0
        do i=1,ncolor+1
            ichr(i)=' '
            zc=f(i)
            if(zc<0.0e0_rknd) then
                call sreal(ichr(i)(2:2),nn,zc,3_iknd,1_iknd)
                nchr(i)=nn+1
            else
                call sreal(ichr(i)(3:3),nn,zc,3_iknd,1_iknd)
                nchr(i)=nn+2
            endif
            mxchr=max(mxchr,nchr(i))
        enddo
        do i=1,ncolor+1
            yc=yi+scale*real(i-1)-yinc/2.0e0_rknd
            yf=yc+yinc
            call htext(xc,yc,xx,yf,mxchr,ichr(i),-1_iknd,qq,tt,2_iknd)
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
        subroutine linit(t,q)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(25) :: t
            real(kind=rknd), dimension(3,3) :: q
cy
c
c       initial for legends and graphs
c
        size=0.9e0_rknd
        do i=1,25
            t(i)=0.0e0_rknd
        enddo
        t(3)=1.0e0_rknd
        t(5)=0.5e0_rknd
        t(7)=1.0e1_rknd*epsilon(1.0e0_rknd)
        t(12)=1.0e0_rknd
        t(14)=size
        t(15)=0.5e0_rknd
        t(16)=0.5e0_rknd
        t(17)=0.5e0_rknd
        do i=1,3
            do j=1,3
                q(i,j)=0.0e0_rknd
            enddo
            q(i,i)=1.0e0_rknd
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
        subroutine zoombx(rp,t)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(100) :: rp
            real(kind=rknd), dimension(25) :: t
cy
c       compute the zoom-in window
c
        size=t(14)
        xs=t(15)
        ys=t(16)
        rmag=max(1.0e0_rknd,rp(8))
        cenx=max(0.0e0_rknd,rp(9))
        cenx=min(1.0e0_rknd,cenx)
        ceny=max(0.0e0_rknd,rp(10))
        ceny=min(1.0e0_rknd,ceny)
        h=1.0e0_rknd/(2.0e0_rknd*rmag)
        hx=xs-size/2.0e0_rknd
        hy=ys-size/2.0e0_rknd
        t(8)=size*(cenx-h)+hx
        t(9)=size*(cenx+h)+hx
        t(10)=size*(ceny-h)+hy
        t(11)=size*(ceny+h)+hy
        t(12)=rmag
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
        subroutine mkrot(nx,ny,nz,q)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(3,3) :: q
            real(kind=rknd), dimension(3) :: d
cy
c       compute rotation matrix
c
        d(1)=real(nx,rknd)
        d(2)=real(ny,rknd)
        d(3)=real(nz,rknd)
        do i=1,3
            do j=1,3
                q(j,i)=0.0e0_rknd
            enddo
            q(i,i)=1.0e0_rknd
        enddo
        dl=sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
        if(dl>0.0e0_rknd) then
            do i=1,3
                q(i,3)=d(i)/dl
            enddo
        endif
        dl=sqrt(q(1,3)*q(1,3)+q(2,3)*q(2,3))
        if(dl>0.0e0_rknd) then
            q(1,1)=-q(2,3)/dl
            q(2,1)=q(1,3)/dl
            q(1,2)=-q(2,1)*q(3,3)
            q(2,2)=q(1,1)*q(3,3)
            q(3,2)=dl
        else
            if(q(3,3)<0.0e0_rknd) q(1,1)=-q(1,1)
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
        subroutine setgr(filnam,n,ispd,lenja,ja,lena,a,lenb,b,
     +      nblock,ib,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd),  dimension(*) :: ja,ib
            integer(kind=iknd),  save :: iunit=8
            real(kind=rknd), dimension(*) :: a,b
            integer(kind=iknd), allocatable, dimension(:,:) :: irc
            real(kind=rknd), allocatable, dimension(:) :: a0
            character(len=80) :: filnam
cy
c
        allocate(a0(lena),irc(2,lena))
c
        nnz=0
        open(unit=iunit,form='formatted',status='old',
     +      file=filnam,access='sequential',err=100)
        read(iunit,*) n,ispd,nblock
        if(n>lenb) then
            iflag=20
            go to 90
        endif
        do i=1,nblock+1
            read(iunit,*) k,ibk
            ib(k)=ibk
        enddo
        do i=1,n
            read(iunit,*) k, bk
            b(k)=bk
        enddo
        nnz=0
   10   read(iunit,*,end=20) i,j,aij
        if(nnz>lena) then
            iflag=21
            go to 90
        endif
        nnz=nnz+1
        irc(1,nnz)=i
        irc(2,nnz)=j
        a0(nnz)=aij
        go to 10
   20   close(unit=iunit)
c
c       make ja
c
        if(lenja<n+1) then
            iflag=22
            go to 90
        endif
        do i=1,n+1
            ja(i)=0
        enddo
c
c       count entries for each row
c       redundancy is built into this algorithm to be safe
c       (eg. if a nonsymmetric zero/nonzero structure is specified for ispd=0
c       or both u and l are provided for the case ispd=1)
c
        do i=1,nnz
            irow=min(irc(1,i),irc(2,i))
            icol=max(irc(1,i),irc(2,i))
            if(irow<icol) ja(irow+1)=ja(irow+1)+1
        enddo
c
c       compute pointers
c
        ja(1)=n+2
        do i=1,n
            ja(i+1)=ja(i)+ja(i+1)
        enddo
        if(ja(n+1)-1>lenja) then
            iflag=23
            go to 90
        endif
        do i=ja(1),ja(n+1)-1
            ja(i)=0
        enddo
c
c       compute col indices
c
        do i=1,nnz
            irow=min(irc(1,i),irc(2,i))
            icol=max(irc(1,i),irc(2,i))
            if(irow<icol) then
                k=ja(irow)
                ja(irow)=ja(irow)+1
                ja(k)=icol
            endif
        enddo
        do i=n+1,2,-1
            ja(i)=ja(i-1)
        enddo
        ja(1)=n+2
c
c       remove duplicate entries
c
        jai=ja(1)
        do i=1,n
           len=ja(i+1)-jai
           call ihp(ja(jai),len)
           is=ja(i)
           last=0
           do j=jai,ja(i+1)-1
               if(ja(j)/=last) then
                   ja(is)=ja(j)
                   is=is+1
                   last=ja(j)
               endif
           enddo
           jai=ja(i+1)
           ja(i+1)=is
        enddo
c
c       now compute a
c
        if(ispd==1) then
            lmtx=0
        else
            lmtx=ja(n+1)-ja(1)
        endif
c
        do i=1,ja(n+1)-1+lmtx
            a(i)=0.0e0
        enddo
        do i=1,nnz
            irow=irc(1,i)
            icol=irc(2,i)
            if(irow==icol) then
                a(irow)=a0(i)
            else
                call jamap(irow,icol,jrc,jcr,ja,lmtx)
                a(jrc)=a0(i)
            endif
        enddo
c
        iflag=0
   90   deallocate(a0,irc)
        return
  100   iflag=30
        go to 90
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine setmtx(ip,ja0,a0,b0,ja,a,b,ib,ka)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ja0,ib
            integer(kind=iknd), dimension(100) :: ip
            integer(kind=iknd), dimension(2,*) :: ka(10,*)
            real(kind=rknd), dimension(*) :: a,b,a0,b0
cy
        n=ip(1)
        nblock=ip(3)
        maxja0=ip(71)
        do i=1,maxja0
            ja(i)=ja0(i)
        enddo
        maxa0=ip(72)
        do i=1,maxa0
            a(i)=a0(i)
        enddo
        do i=1,n
            b(i)=b0(i)
        enddo
        lvl=1
        ip(2)=lvl
c
        ispd=ip(8)
        lenja=ja(n+1)-1
        if(ispd==1) then
            lena=lenja
        else
            lena=2*lenja-(n+1)
        endif
        nptr=1
        japtr=1
        iaptr=1
        ibptr=japtr+lenja
        iqptr=ibptr+n
        juptr=japtr
        iuptr=iaptr
c
        ka(1,lvl)=n
        ka(2,lvl)=nptr
        ka(3,lvl)=japtr
        ka(4,lvl)=iaptr
        ka(5,lvl)=juptr
        ka(6,lvl)=iuptr
        ka(7,lvl)=iqptr+n
        ka(8,lvl)=iuptr+lena
        ka(9,lvl)=iqptr
        ka(10,lvl)=ibptr
c
        do i=1,nblock
           do j=ib(i),ib(i+1)-1
              ja(lenja+j)=i
           enddo
        enddo
        do i=1,n
            ja(iqptr+i-1)=i
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
        subroutine linsys(ngrid,mtxtyp,name,n,ispd,lenja,ja,
     +      lena,a,lenb,b,nblock,ib,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ib
            real(kind=rknd), dimension(*) :: a,b
            character(len=80) :: name
cy
        if(ngrid<=1) then
            write(unit=name,fmt='(a10)') 'ngrid<=1'
            iflag=4
        else if(mtxtyp==0) then
            write(unit=name,fmt='(a2,i4)') 'gg',ngrid
            diag=4.0e0_rknd
            off=-1.0e0_rknd
            call star5(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +          nblock,ib,diag,off,iflag)
        else if(mtxtyp==1) then
            write(unit=name,fmt='(a2,i4)') 'mm',ngrid
            diag=4.0e0_rknd
            off=1.0e0_rknd
            call star5(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +          nblock,ib,diag,off,iflag)
        else if(mtxtyp==2) then
            write(unit=name,fmt='(a2,i4)') 'dd',ngrid
            diag=6.0e0_rknd
            off=-1.0e0_rknd
            call star7(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +          nblock,ib,diag,off,iflag)
        else if(mtxtyp==3) then
            write(unit=name,fmt='(a2,i4)') 'ss',ngrid
            call stoke(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +          nblock,ib,iflag)
        else if(mtxtyp==4) then
            write(unit=name,fmt='(a2,i4)') 'g9',ngrid
            diag=8.0e0_rknd
            off0=-1.0e0_rknd
            off1=-1.0e0_rknd
            call star9(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +          nblock,ib,diag,off0,off1,iflag)
        else if(mtxtyp==5) then
            write(unit=name,fmt='(a2,i4)') 'm9',ngrid
            diag=8.0e0_rknd
            off0=1.0e0_rknd
            off1=1.0e0_rknd
            call star9(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +          nblock,ib,diag,off0,off1,iflag)
        else
            write(unit=name,fmt='(a6)') 'notype'
            iflag=4
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
        subroutine star7(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +      nblock,ib,diag,off,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ib
            real(kind=rknd), dimension(*) :: a,b
cy
        n=ngrid**3
        ispd=1
        nblock=1
        ib(1)=1
        ib(2)=n+1
c
        n1=ngrid
        n2=ngrid**2
c
        rhs=real(n2)
c
        if(lenb<n) then
            iflag=1
            return
        endif
        if(lenja<4*n) then
            iflag=2
            return
        endif
        if(lena<4*n) then
            iflag=3
            return
        endif
c
        iflag=0
c
        ja(1)=n+2
        do i=1,n1
           do j=1,n1
                do k=1,n1
                    m=(i-1)*n2+(j-1)*n1+k
                    b(m)=rhs
                    a(m)=diag
                    next=ja(m)
                    if(k<n1) then
                         ja(next)=m+1
                         a(next)=off
                         next=next+1
                    endif
                    if(j<n1) then
                         ja(next)=m+n1
                         a(next)=off
                         next=next+1
                    endif
                    if(i<n1) then
                         ja(next)=m+n2
                         a(next)=off
                         next=next+1
                    endif
                    ja(m+1)=next
                enddo
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
        subroutine star5(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +      nblock,ib,diag,off,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ib
            real(kind=rknd), dimension(*) :: a,b
cy
        n=ngrid**2
        ispd=1
        nblock=1
        ib(1)=1
        ib(2)=n+1
c
        n1=ngrid
c
        rhs=real(n)
c
        if(lenb<n) then
            iflag=1
            return
        endif
        if(lenja<3*n) then
            iflag=2
            return
        endif
        if(lena<3*n) then
            iflag=3
            return
        endif
c
        iflag=0
c
        ja(1)=n+2
        do i=1,n1
           do j=1,n1
                m=(i-1)*n1+j
                b(m)=rhs
                a(m)=diag
                next=ja(m)
                if(j<n1) then
                     ja(next)=m+1
                     a(next)=off
                     next=next+1
                endif
                if(i<n1) then
                     ja(next)=m+n1
                     a(next)=off
                     next=next+1
                endif
                ja(m+1)=next
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
        subroutine star9(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +      nblock,ib,diag,off0,off1,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ib
            real(kind=rknd), dimension(*) :: a,b
cy
        n=ngrid**2
        ispd=1
        nblock=1
        ib(1)=1
        ib(2)=n+1
c
        n1=ngrid
c
        rhs=real(n)
c
        if(lenb<n) then
            iflag=1
            return
        endif
        if(lenja<5*n) then
            iflag=2
            return
        endif
        if(lena<5*n) then
            iflag=3
            return
        endif
c
        iflag=0
c
        ja(1)=n+2
        do i=1,n1
           do j=1,n1
                m=(i-1)*n1+j
                b(m)=rhs
                a(m)=diag
                next=ja(m)
                if(j<n1) then
                     ja(next)=m+1
                     a(next)=off0
                     next=next+1
                endif
                if(i<n1) then
                     if(j>1) then
                         ja(next)=m+n1-1
                         a(next)=off1
                         next=next+1
                     endif
                     ja(next)=m+n1
                     a(next)=off0
                     next=next+1
                     if(j<n1) then
                         ja(next)=m+n1+1
                         a(next)=off1
                         next=next+1
                     endif
                 endif
                ja(m+1)=next
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
        subroutine stoke(ngrid,n,ispd,lenja,ja,lena,a,lenb,b,
     +      nblock,ib,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,ib
            integer(kind=iknd), dimension(9) :: iv
            real(kind=rknd), dimension(*) :: a,b
            real(kind=rknd), dimension(9) :: be
            real(kind=rknd), dimension(9,9) :: ae
cy
c
        n1=ngrid
        n2=n1**2
        n=3*n2
        ispd=1
        nblock=3
        ib(1)=1
        ib(2)=n2+1
        ib(3)=2*n2+1
        ib(4)=n+1
c
        if(lenb<n) then
            iflag=1
            return
        endif
        if(lenja<9*n) then
            iflag=2
            return
        endif
        if(lena<9*n) then
            iflag=3
            return
        endif
        iflag=0
c
        call setgr0(ngrid,ja)
        ll=ja(n+1)-1
        do i=1,n
            b(i)=0.0e0_rknd
        enddo
        do i=1,ll
            a(i)=0.0e0_rknd
        enddo
        ishift=0
        do ix=1,n1-1
        do iy=1,n1-1
           call ap(ix,iy,ngrid,iv,ae,be)
           do i=1,9
                ii=iv(i)
                b(ii)=b(ii)+be(i)
                a(ii)=a(ii)+ae(i,i)
                do j=i+1,9
                    jj=iv(j)
                    if(ae(i,j)/=0) then
                       call jamap(ii,jj,ij,ji,ja,ishift)
                       a(ij)=a(ij)+ae(i,j)
                    endif
                enddo
           enddo
        enddo
        enddo
        do ix=2,ngrid
        do iy=2,ngrid
           call am(ix,iy,ngrid,iv,ae,be)
           do i=1,9
                ii=iv(i)
                b(ii)=b(ii)+be(i)
                a(ii)=a(ii)+ae(i,i)
                do j=i+1,9
                    jj=iv(j)
                    if(ae(i,j)/=0) then
                       call jamap(ii,jj,ij,ji,ja,ishift)
                       a(ij)=a(ij)+ae(i,j)
                    endif
                enddo
           enddo
        enddo
        enddo
c
c       compress
c
        jai=ja(1)
        do i=1,n
           k=ja(i)
           do j=jai,ja(i+1)-1
                if(a(j)/=0.0e0_rknd) then
                    ja(k)=ja(j)
                    a(k)=a(j)
                    k=k+1
               endif
           enddo
           jai=ja(i+1)
           ja(i+1)=k
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
        subroutine setgr0(n,ja)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
cy
c
c       first block
c
        nn=n**2
        ja(1)=3*nn+2
        n3=2*nn
        do j=1,n
            do i=1,n
                k=i+(j-1)*n
                next=ja(k)
                if(i/=n) then
                    ja(next)=k+1
                    next=next+1
                endif
                if(i/=1.and.j/=n) then
                    ja(next)=k-1+n
                    next=next+1
                endif
                if(j/=n) then
                    ja(next)=k+n
                    next=next+1
                endif
                if(j/=1) then
                    ja(next)=k-n+n3
                    next=next+1
                endif
                if(j/=1.and.i/=n) then
                    ja(next)=k-n+1+n3
                    next=next+1
                endif
                if(i/=1) then
                    ja(next)=k-1+n3
                    next=next+1
                endif
                ja(next)=k+n3
                next=next+1
                if(i/=n) then
                    ja(next)=k+1+n3
                    next=next+1
                endif
                if(i/=1.and.j/=n) then
                    ja(next)=k-1+n+n3
                    next=next+1
                endif
                if(j/=n) then
                    ja(next)=k+n+n3
                    next=next+1
                endif
                ja(k+1)=next
            enddo
        enddo
c
c       second block
c
        n3=nn
        do j=1,n
            do i=1,n
                k=i+(j-1)*n+nn
                next=ja(k)
                if(i/=n) then
                    ja(next)=k+1
                    next=next+1
                endif
                if(i/=1.and.j/=n) then
                    ja(next)=k-1+n
                    next=next+1
                endif
                if(j/=n) then
                    ja(next)=k+n
                    next=next+1
                endif
                if(j/=1) then
                    ja(next)=k-n+n3
                    next=next+1
                endif
                if(j/=1.and.i/=n) then
                    ja(next)=k-n+1+n3
                    next=next+1
                endif
                if(i/=1) then
                    ja(next)=k-1+n3
                    next=next+1
                endif
                ja(next)=k+n3
                next=next+1
                if(i/=n) then
                    ja(next)=k+1+n3
                    next=next+1
                endif
                if(i/=1.and.j/=n) then
                    ja(next)=k-1+n+n3
                    next=next+1
                endif
                if(j/=n) then
                    ja(next)=k+n+n3
                    next=next+1
                endif
                ja(k+1)=next
            enddo
        enddo
c
c       third block
c
        do j=1,n
            do i=1,n
                k=i+(j-1)*n+nn*2
                next=ja(k)
                if(i/=n) then
                    ja(next)=k+1
                    next=next+1
                endif
                if(i/=1.and.j/=n) then
                    ja(next)=k-1+n
                    next=next+1
                endif
                if(j/=n) then
                    ja(next)=k+n
                    next=next+1
                endif
                ja(k+1)=next
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
        subroutine ap(ix,iy,n,iv,a,b)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(9) :: iv
            real(kind=rknd), dimension(9) :: b
            real(kind=rknd), dimension(9,9) :: a
cy
c
c
        h=1.0e0_rknd/real(n-1)
        cl=1.0e0_rknd/2.0e0_rknd
        cs=-h**2/80.0e0_rknd
        cb=-h/6.0e0_rknd
c
c       compute iv
c
        nn=n**2
        k=ix+(iy-1)*n
        iv(1)=k
        iv(2)=k+1
        iv(3)=k+n
        iv(4)=iv(1)+nn
        iv(5)=iv(2)+nn
        iv(6)=iv(3)+nn
        iv(7)=iv(4)+nn
        iv(8)=iv(5)+nn
        iv(9)=iv(6)+nn
c
c       initialize
c
        do i=1,9
            if(i<=6) then
                b(i)=h**2
            else
                b(i)=0.0e0_rknd
            endif
            do j=1,9
                a(i,j)=0.0e0_rknd
            enddo
        enddo
        a(1,1)=2.0e0_rknd*cl
        a(1,2)=-cl
        a(1,3)=-cl
        a(2,1)=-cl
        a(2,2)=cl
        a(3,1)=-cl
        a(3,3)=cl
c
        a(4,4)=2.0e0_rknd*cl
        a(4,5)=-cl
        a(4,6)=-cl
        a(5,4)=-cl
        a(5,5)=cl
        a(6,4)=-cl
        a(6,6)=cl
c
        a(7,7)=2.0e0_rknd*cs
        a(7,8)=-cs
        a(7,9)=-cs
        a(8,7)=-cs
        a(8,8)=cs
        a(9,7)=-cs
        a(9,9)=cs
c
        a(1,7)=-cb
        a(2,7)=-cb
        a(3,7)=-cb
        a(1,8)=cb
        a(2,8)=cb
        a(3,8)=cb
c
        a(7,1)=-cb
        a(7,2)=-cb
        a(7,3)=-cb
        a(8,1)=cb
        a(8,2)=cb
        a(8,3)=cb
c
        a(4,7)=-cb
        a(5,7)=-cb
        a(6,7)=-cb
        a(4,9)=cb
        a(5,9)=cb
        a(6,9)=cb
c
        a(7,4)=-cb
        a(7,5)=-cb
        a(7,6)=-cb
        a(9,4)=cb
        a(9,5)=cb
        a(9,6)=cb
c
c       dirichlet bc
c
        if(ix==1.or.iy==1) then
            b(1)=0.0e0_rknd
            b(4)=0.0e0_rknd
            do i=1,9
                if(i/=1) then
                    a(1,i)=0.0e0_rknd
                    a(i,1)=0.0e0_rknd
                 endif
                if(i/=4) then
                    a(4,i)=0.0e0_rknd
                    a(i,4)=0.0e0_rknd
                 endif
            enddo
        endif
        if(iy==1.or.ix==n-1) then
            b(2)=0.0e0_rknd
            b(5)=0.0e0_rknd
            do i=1,9
                if(i/=2) then
                    a(2,i)=0.0e0_rknd
                    a(i,2)=0.0e0_rknd
                 endif
                if(i/=5) then
                    a(5,i)=0.0e0_rknd
                    a(i,5)=0.0e0_rknd
                 endif
            enddo
        endif
        if(ix==1.or.iy==n-1) then
            b(3)=0.0e0_rknd
            b(6)=0.0e0_rknd
            do i=1,9
                if(i/=3) then
                    a(3,i)=0.0e0_rknd
                    a(i,3)=0.0e0_rknd
                 endif
                if(i/=6) then
                    a(6,i)=0.0e0_rknd
                    a(i,6)=0.0e0_rknd
                 endif
            enddo
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
        subroutine am(ix,iy,n,iv,a,b)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(9) :: iv
            real(kind=rknd), dimension(9) :: b
            real(kind=rknd), dimension(9,9) :: a
cy
c
c
        h=1.0e0_rknd/real(n-1)
        cl=1.0e0_rknd/2.0e0_rknd
        cs=-h**2/80.0e0_rknd
        cb=h/6.0e0_rknd
c
c       compute iv
c
        nn=n**2
        k=ix+(iy-1)*n
        iv(1)=k
        iv(2)=k-1
        iv(3)=k-n
        iv(4)=iv(1)+nn
        iv(5)=iv(2)+nn
        iv(6)=iv(3)+nn
        iv(7)=iv(4)+nn
        iv(8)=iv(5)+nn
        iv(9)=iv(6)+nn
c
c       initialize
c
        do i=1,9
            if(i<=6) then
                b(i)=h**2
            else
                b(i)=0.0e0_rknd
            endif
            do j=1,9
                a(i,j)=0.0e0_rknd
            enddo
        enddo
        a(1,1)=2.0e0_rknd*cl
        a(1,2)=-cl
        a(1,3)=-cl
        a(2,1)=-cl
        a(2,2)=cl
        a(3,1)=-cl
        a(3,3)=cl
c
        a(4,4)=2.0e0_rknd*cl
        a(4,5)=-cl
        a(4,6)=-cl
        a(5,4)=-cl
        a(5,5)=cl
        a(6,4)=-cl
        a(6,6)=cl
c
        a(7,7)=2.0e0_rknd*cs
        a(7,8)=-cs
        a(7,9)=-cs
        a(8,7)=-cs
        a(8,8)=cs
        a(9,7)=-cs
        a(9,9)=cs
c
        a(1,7)=-cb
        a(2,7)=-cb
        a(3,7)=-cb
        a(1,8)=cb
        a(2,8)=cb
        a(3,8)=cb
c
        a(7,1)=-cb
        a(7,2)=-cb
        a(7,3)=-cb
        a(8,1)=cb
        a(8,2)=cb
        a(8,3)=cb
c
        a(4,7)=-cb
        a(5,7)=-cb
        a(6,7)=-cb
        a(4,9)=cb
        a(5,9)=cb
        a(6,9)=cb
c
        a(7,4)=-cb
        a(7,5)=-cb
        a(7,6)=-cb
        a(9,4)=cb
        a(9,5)=cb
        a(9,6)=cb
c
c       dirichlet bc
c
        if(ix==n.or.iy==n) then
            b(1)=0.0e0_rknd
            b(4)=0.0e0_rknd
            do i=1,9
                if(i/=1) then
                    a(1,i)=0.0e0_rknd
                    a(i,1)=0.0e0_rknd
                 endif
                if(i/=4) then
                    a(4,i)=0.0e0_rknd
                    a(i,4)=0.0e0_rknd
                 endif
            enddo
        endif
        if(iy==n.or.ix==2) then
            b(2)=0.0e0_rknd
            b(5)=0.0e0_rknd
            do i=1,9
                if(i/=2) then
                    a(2,i)=0.0e0_rknd
                    a(i,2)=0.0e0_rknd
                 endif
                if(i/=5) then
                    a(5,i)=0.0e0_rknd
                    a(i,5)=0.0e0_rknd
                 endif
            enddo
        endif
        if(ix==n.or.iy==2) then
            b(3)=0.0e0_rknd
            b(6)=0.0e0_rknd
            do i=1,9
                if(i/=3) then
                    a(3,i)=0.0e0_rknd
                    a(i,3)=0.0e0_rknd
                 endif
                if(i/=6) then
                    a(6,i)=0.0e0_rknd
                    a(i,6)=0.0e0_rknd
                 endif
            enddo
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
        subroutine scale0(n,ispd,ja,a,s,b)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            real(kind=rknd), dimension(*) :: a,b,s
cy
c       nonsymmetric scaling with positive diagonal
c
        call s2ns(n,ispd,ja,a)
        ishift=0
        if(ispd/=1) ishift=ja(n+1)-ja(1)
c
c
        itmax=5
        do itnum=1,itmax
c
c       row
c
            do i=1,n
                s(i)=abs(a(i))
            enddo
            do i=1,n
                do jj=ja(i),ja(i+1)-1
                    j=ja(jj)
                    s(i)=max(s(i),abs(a(jj)))
                    s(j)=max(s(j),abs(a(jj+ishift)))
                enddo
            enddo
            do i=1,n
                if(a(i)<0.0e0_rknd) s(i)=-s(i)
                s(i)=1.0e0_rknd/s(i)
                b(i)=b(i)*s(i)
            enddo
            do i=1,n
                a(i)=a(i)*s(i)
                do jj=ja(i),ja(i+1)-1
                    j=ja(jj)
                    aa=a(jj)
                    a(jj+ishift)=a(jj+ishift)*s(j)
                    a(jj)=aa*s(i)
                enddo
            enddo
c
c       column
c
            do i=1,n
                s(i)=abs(a(i))
            enddo
            do i=1,n
                do jj=ja(i),ja(i+1)-1
                    j=ja(jj)
                    s(j)=max(s(j),abs(a(jj)))
                    s(i)=max(s(i),abs(a(jj+ishift)))
                enddo
            enddo
            do i=1,n
                if(a(i)<0.0e0_rknd) s(i)=-s(i)
                s(i)=1.0e0_rknd/s(i)
            enddo
            do i=1,n
                a(i)=a(i)*s(i)
                do jj=ja(i),ja(i+1)-1
                    j=ja(jj)
                    aa=a(jj)
                    a(jj+ishift)=a(jj+ishift)*s(i)
                    a(jj)=aa*s(j)
                enddo
            enddo
c
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
        subroutine scale1(n,ispd,ja,a,s,b)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            real(kind=rknd), dimension(*) :: a,b,s
cy
c       symmetric scaling
c
        ishift=0
        if(ispd/=1) ishift=ja(n+1)-ja(1)
c
        itmax=5
        do itnum=1,itmax
c
            do i=1,n
                s(i)=sqrt(abs(a(i)))
            enddo
            do i=1,n
                do jj=ja(i),ja(i+1)-1
                    j=ja(jj)
                    s(j)=max(abs(a(jj)),s(j))
                    s(i)=max(abs(a(jj+ishift)),s(i))
                enddo
            enddo
            do i=n,1,-1
                do jj=ja(i),ja(i+1)-1
                    j=ja(jj)
                    au=a(jj)*s(j)
                    al=a(jj+ishift)*s(j)
                    a(jj)=au
                    a(jj+ishift)=al
                    s(i)=max(s(i),abs(au),abs(al))
                enddo
                if (s(i)>0.0e0_rknd) then
                   s(i)=1.0e0_rknd/s(i)
                else
                   s(i)=1.0e0_rknd
                end if
                b(i)=s(i)*b(i)
                a(i)=s(i)*a(i)*s(i)
                do jj=ja(i),ja(i+1)-1
                    au=a(jj)*s(i)
                    al=a(jj+ishift)*s(i)
                    a(jj)=au
                    a(jj+ishift)=al
                enddo
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
        subroutine s2ns(n,ispd,ja,a)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja
            real(kind=rknd), dimension(*) :: a
cy
c       convert symmetric to non symmetric data structure
c
        if (ispd/=1) return
        ishift=ja(n+1)-ja(1)
        ispd=0
c
        do k=ja(1),ja(n+1)-1
           a(k+ishift)=a(k)
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
        subroutine setls(n,ispd,maxja,ja,jc,maxa,a,c,b,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: jc,ja
            real(kind=rknd), dimension(*) :: a,b,c
            real(kind=rknd), dimension(n) :: z
cy
c       assumes ja/jc and a/c are same size (as ja0/ja a0/a in atest)
c
        iflag=0
        if(ispd==1) then
            jspd=1
        else
            jspd=-(1+ispd)
        endif
        lenjc=ja(n+1)-1+ja(n+1)-ja(1)
        if(lenjc>maxja.or.lenjc>maxa) then
            iflag=1
            return
        endif
c
c       form a*t b and stor in b
c
        call mtxmlt(n,ja,a,b,z,jspd)
        do i=1,n
            b(i)=z(i)
        enddo
c
c       form jc/c from ja/a
c
        call ja2jc(n,ja,jc)
        call a2c(n,ispd,ja,jc,a,c)
c
c       form ja2/a2 from jc/c
c
        call jc2ja2(n,jc,ja,maxja,iflag)
        if(iflag/=0) return
        call c2a2(n,ja,jc,a,c)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine jc2ja2(n,jc,ja2,maxja,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: jc,ja2
            integer(kind=iknd), dimension(n) :: mark,list
cy
c       make ja2 data structure from ja data structure
c
        iflag=0
        do i=1,n
            mark(i)=0
        enddo
        ja2(1)=n+2
        do i=1,n
c
c       make list of paths of length 1 and 2
c
            len=0
            do jj=jc(i),jc(i+1)-1
                j=jc(jj)
                if(mark(j)==0) then
                    len=len+1
                    list(len)=j
                    mark(j)=len
                endif
                do kk=jc(j),jc(j+1)-1
                    k=jc(kk)
                    if(mark(k)==0) then
                        len=len+1
                        list(len)=k
                        mark(k)=len
                    endif
                enddo
            enddo
c
c       put strict upper triangle indices into ja2
c
            next=ja2(i)
            do j=1,len
                if(list(j)>i) then
                    if(next>maxja) then
                        iflag=1
                        return
                    endif
                    ja2(next)=list(j)
                    next=next+1
                endif
                mark(list(j))=0
            enddo
            ja2(i+1)=next
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
        subroutine jcmap(i,j,ij,ji,jc)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: jc
cy
c       compute location of a(i,j) and a(j,i) in jc
c       this assumes rowwise storage
c
        do ij=jc(i),jc(i+1)-1
            if(jc(ij)==j) go to 10
        enddo
        stop 2233
   10   do ji=jc(j),jc(j+1)-1
            if(jc(ji)==i) return
        enddo
        stop 2244
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine a2c(n,ispd,ja,jc,a,c)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja,jc
            real(kind=rknd), dimension(*) :: a,c
cy
c       make rowise c/jc  matrix from a/ja
c
        if(ispd/=1) then
            lmtx=ja(n+1)-ja(1)
        else
            lmtx=0
        endif
        do i=1,jc(n+1)-1
            c(i)=0.0e0_rknd
        enddo
        do i=1,n
            c(i)=a(i)
            do jj=ja(i),ja(i+1)-1
                j=ja(jj)
                call jcmap(i,j,ij,ji,jc)
                c(ij)=a(jj)
                c(ji)=a(jj+lmtx)
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
        subroutine c2a2(n,ja2,jc,a2,c)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ja2,jc
            real(kind=rknd), dimension(*) :: a2,c
cy
c       make a*a/ja2  matrix from c/jc
c
        do i=1,ja2(n+1)-1
            a2(i)=0.0e0_rknd
        enddo
        do i=1,n
            a2(i)=a2(i)+c(i)*c(i)
            do jj=jc(i),jc(i+1)-1
                j=jc(jj)
                call jamap(i,j,ij,ji,ja2,0)
                if(j>i) then
                    a2(ij)=a2(ij)+c(i)*c(jj)
                else
                    a2(ij)=a2(ij)+c(j)*c(jj)
                endif
                a2(j)=a2(j)+c(jj)*c(jj)
                do kk=jj+1,jc(i+1)-1
                    k=jc(kk)
                    call jamap(k,j,kj,jk,ja2,0)
                    a2(kj)=a2(kj)+c(kk)*c(jj)
               enddo
            enddo
        enddo
        return
        end
c*********************** machine dependent routine *********************
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine pltutl(ncolor,red,green,blue)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: red,green,blue
            character(len=80) :: fname,fname0,sp
            common /atest1/ip(100),rp(100),sp(100)
            common /atest5/idevce
cy
c       ncolor > 0 -- initialize graphics using ncolor colors
c       ncolor <= 0 -- exit graphics
c
c       socket graphics
c
        if(idevce>=0.and.idevce<=3) then
            isock=idevce
            call fstr(fname,length,sp(21),0_iknd)
            call vutl(ncolor,red,green,blue,isock,fname)
            if(ncolor<0) call cpause()
c
c       bh file
c
        else if(idevce==4) then
            if(ncolor>0) then
                call mkname(fname0,sp(20))
                call stfile(fname,fname0)
            endif
            call vutl(ncolor,red,green,blue,-1_iknd,fname)
c
c       postscript file
c
        else if(idevce==5) then
            if(ncolor>0) then
                call mkname(fname0,sp(18))
                call stfile(fname,fname0)
            endif
            call psutl(ncolor,red,green,blue,fname)
c
c       xpm file
c
        else if(idevce==6) then
            if(ncolor>0) then
                call mkname(fname0,sp(19))
                call stfile(fname,fname0)
            endif
            call xpmutl(ncolor,red,green,blue,fname)
c
c       classic x graphics
c
        else if(idevce>=7.and.idevce<=10) then
            isock=idevce-7
            call xutl(ncolor,red,green,blue,isock)
            if(ncolor<0) call cpause()
c
c       hf file
c
        else if(idevce==11) then
            call hfutl(ncolor,red,green,blue,sp(17))
        endif
        return
        end
c*********************** machine dependent routine *********************
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine pframe(iframe)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            common /atest5/idevce
cy
c       frame/list equivalence table
c        ___ ___ ___       ___ ___ ___
c       |       |   |     |           |
c       |       | 2 |     |           |
c       |   4   |___|     |     1     |
c       |       |   |     |           |
c       |       | 3 |     |           |
c       |___ ___|___|     |___ ___ ___|
c
c        list    frame        type
c
c          1       1          non-rotating, non-lighted
c
c          2       2          non-rotating, non-lighted
c
c          3       3          non-rotating, non-lighted
c
c          4       4          non-rotating, non-lighted
c          5       4              rotating, non-lighted
c          6       4              rotating, non-lighted
c          7       4              rotating,     lighted
c          8       4              rotating,     lighted
c          9       4          non-rotating,     lighted
c
c
        if(idevce>=0.and.idevce<=4) then
            call vframe(iframe)
        else if(idevce==5) then
            call sframe(iframe)
        else if(idevce>=6.and.idevce<=10) then
            call xframe(iframe)
        else if(idevce==11) then
            call hframe(iframe)
        endif
        return
        end
c*********************** machine dependent routine *********************
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine pline(x,y,z,n,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: x,y,z
            common /atest5/idevce
cy
c       subroutine pline moves the pen (or whatever)
c       to the point (x(1),y(1)), and then draws the
c       n-1 line segments (x(i-1),y(i-1)) to (x(i),y(i)),
c       i=2,3,....n.
c
        if(idevce>=0.and.idevce<=4) then
            call vline(x,y,z,n,icolor)
        else if(idevce==5) then
            call pspath(x,y,z,n,icolor,0_iknd)
        else if(idevce>=6.and.idevce<=10) then
            call xline(x,y,z,n,icolor)
        else if(idevce==11) then
            call hline(x,y,z,n,icolor)
        endif
        return
        end
c*********************** machine dependent routine *********************
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine pfill(x,y,z,n,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: x,y,z
            common /atest5/idevce
cy
c       subroutine pfill fills the n-sided polygon with
c       vertices (x(i),y(i)) with the indicated color
c
        if(idevce>=0.and.idevce<=4) then
            call vfill(x,y,z,n,icolor)
        else if(idevce==5) then
            call pspath(x,y,z,n,icolor,1_iknd)
        else if(idevce>=6.and.idevce<=10) then
            call xfill(x,y,z,n,icolor)
        else if(idevce==11) then
            call hfill(x,y,z,n,icolor)
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
        subroutine psutl(ncolor,red,green,blue,fname)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), save :: length
            real(kind=rknd), dimension(*) :: red,green,blue
            character(len=16), save :: hex='0123456789abcdef'
            character(len=80) :: msg,fname
            character(len=80), save :: sname
            common /ps0/id
            common /ps1/scale,fscale,xshift,yshift
cy
c       postscript graphics implementation for pltutl
c       this version is based on suggestions of klas samuelsson for
c       reducing the size of the postscript files
c
c       postscript limited to 256 colors
c
c       print picture
c
        if(ncolor<=0) then
            msg='showpage'
            call ascstr(id,msg,8_iknd,iflag)
            call ascutl(id,sname,'c',iflag)
            return
        endif
c
c       ipl = 1 (0) is portrait (landscape) mode
c       center in 8.5 x 11 inch paper
c       picture is 8 (10.5) inches wide in portrait (landscape)
c       note there are 72 points per inch
c
        ipl=1
c
c       scale factor is 5.e3 (about 4 digits of resolution)
c
        scale=5.0e3_rknd
        fscale=1.0e0_rknd
        xshift=0.0e0_rknd
        yshift=0.0e0_rknd
c
        call fstr(sname,length,fname,0_iknd)
        call ascutl(id,sname,'w',iflag)
c
c       set main definitions
c
        msg='%!'
        call ascstr(id,msg,2_iknd,iflag)
c
        if(ipl==1) then
c***        msg='%%BoundingBox: 18 204 402 588'
            msg='%%BoundingBox: 18 204 594 588'
            call ascstr(id,msg,29_iknd,iflag)
            msg='[384 0 0 384 18 204] concat'
            call ascstr(id,msg,27_iknd,iflag)
        else
            msg='%%BoundingBox: 54 18 558 774'
            call ascstr(id,msg,28_iknd,iflag)
            msg='[0 504 -504 0 558 18] concat'
            call ascstr(id,msg,28_iknd,iflag)
        endif
c
        si=1.0e0_rknd/scale
        write(unit=msg,fmt='(2(f8.6,1x),a5)') si,si,'scale'
        call ascstr(id,msg,23_iknd,iflag)
c
        msg='1 setlinewidth'
        call ascstr(id,msg,14_iknd,iflag)
        msg='2 setlinejoin'
        call ascstr(id,msg,13_iknd,iflag)
        msg='/s {setrgbcolor newpath moveto} def'
        call ascstr(id,msg,35_iknd,iflag)
        msg='/r {count 2 idiv {rlineto} repeat} def'
        call ascstr(id,msg,38_iknd,iflag)
        msg='/f {s r closepath fill} def'
        call ascstr(id,msg,27_iknd,iflag)
        msg='/g {s r stroke} def'
        call ascstr(id,msg,19_iknd,iflag)
c
c       define colors
c
        do i=1,ncolor
            i1=(i-1)/16
            i0=i-1-i1*16
c
            write(unit=msg,fmt='(a2,a1,a1,a2,3(f4.2,1x),a6)')
     +          '/b',hex(i1+1:i1+1),hex(i0+1:i0+1),' {',
     1           red(i),green(i),blue(i),'g} def'
            call ascstr(id,msg,27_iknd,iflag)
c
            write(unit=msg,fmt='(a2,a1,a1,a2,3(f4.2,1x),a6)')
     +          '/c',hex(i1+1:i1+1),hex(i0+1:i0+1),' {',
     1           red(i),green(i),blue(i),'f} def'
            call ascstr(id,msg,27_iknd,iflag)
c
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
        subroutine sframe(iframe)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=80) :: msg
            common /ps0/id
            common /ps1/scale,fscale,xshift,yshift
cy
        write(unit=msg,fmt='(a3,i3)') '%%l',iframe
        call ascstr(id,msg,6_iknd,iflag)
c
        if(iframe==2) then
            fscale=scale/2.0e0_rknd
            xshift=scale
            yshift=scale/2.0e0_rknd
        else if(iframe==3) then
            fscale=scale/2.0e0_rknd
            xshift=scale
            yshift=0.0e0_rknd
        else
            fscale=scale
            xshift=0.0e0_rknd
            yshift=0.0e0_rknd
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
        subroutine pspath(x,y,z,n,icolor,itype)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: x,y,z
            character(len=100) :: list
            character(len=16), save :: hex='0123456789abcdef'
            common /ps0/id
            common /ps1/scale,fscale,xshift,yshift
cy
c       print a path in compact integer form
c
c       look for first nontrivial entry
c
c***    if(scale/=fscale) return
        length=0
        npts=0
        do i=n-1,1,-1
            ix=nint((x(i+1)-x(i))*fscale)
            iy=nint((y(i+1)-y(i))*fscale)
            if(ix==0.and.iy==0) cycle
            npts=npts+1
            call sint(list(length+1:length+1),lenx,ix)
            length=length+lenx+1
            list(length:length)=' '
            call sint(list(length+1:length+1),leny,iy)
            length=length+leny+1
            list(length:length)=' '
c
            if(length<=60) cycle
            call ascstr(id,list,length-1_iknd,iflag)
            length=0
        enddo
c
c       first point
c
        if(npts==0) return
        ix=nint(x(1)*fscale+xshift)
        iy=nint(y(1)*fscale+yshift)
        call sint(list(length+1:length+1),lenx,ix)
        length=length+lenx+1
        list(length:length)=' '
        call sint(list(length+1:length+1),leny,iy)
        length=length+leny+1
        list(length:length)=' '
c
c       set color, and line/fill
c
        if(itype==1) then
            list(length+1:length+1)='c'
        else
            list(length+1:length+1)='b'
        endif
        i1=(icolor-1)/16
        i0=icolor-1-i1*16
        list(length+2:length+2)=hex(i1+1:i1+1)
        list(length+3:length+3)=hex(i0+1:i0+1)
        length=length+3
        call ascstr(id,list,length,iflag)
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
        subroutine xutl(ncolor,red,green,blue,id)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=isngl) image
            integer(kind=iknd), save :: nx,ny
            real(kind=rknd), dimension(*) :: red,green,blue
            common /xpm0/iscale,jscale,ishift,image(540000)
            common /xpm1/scale,fscale,xshift,yshift
            common /atest3/mode,jnlsw,jnlr,jnlw,ibatch
cy
c       xwindows graphics implementation for pltutl
c
        if(mode/=0) return
        if(ncolor<=0) then
            call xgdisp(nx,ny,ishift,image)
            return
        endif
c
c       initialize bitmap
c
        do i=1,ncolor
           image(3*i-2)=int(red(i)*65535.0e0_rknd)
           image(3*i-1)=int(green(i)*65535.0e0_rknd)
           image(3*i)=int(blue(i)*65535.0e0_rknd)
        enddo
        call xginit(ncolor,image,id,ix,iy)
        ny=min(600,iy)
        nx=ny*3/2
        scale=real(ny,rknd)
        iscale=nx
        jscale=ny
        ishift=4096
        do k=1,nx*ny
            image(k)=0
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
        subroutine xpmutl(ncolor,red,green,blue,fname)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), save :: length,nc,nx,ny,lenc,id=0
            integer(kind=isngl) image
            real(kind=rknd), dimension(*) :: red,green,blue
            character(len=1), save, dimension(92) :: cdef
            character(len=2), save, dimension(4096) :: cmap
            character(len=2) :: cs
            character(len=80) :: fname
            character(len=80), save :: sname
            character(len=4000) :: msg
            common /xpm0/iscale,jscale,ishift,image(540000)
            common /xpm1/scale,fscale,xshift,yshift
cy
            data (cdef(i),i=1,92)/
     +          ' ','.','X','o','O','+','@','#','$','%',
     1          '&','*','=','-',';',':','>',',','<','1',
     2          '2','3','4','5','6','7','8','9','0','q',
     3          'w','e','r','t','y','u','i','p','a','s',
     4          'd','f','g','h','j','k','l','z','x','c',
     5          'v','b','n','m','M','N','B','V','C','Z',
     6          'A','S','D','F','G','H','J','K','L','P',
     7          'I','U','Y','T','R','E','W','Q','!','~',
     8          '^','/','(',')','_','`','|',']','[','{',
     9          '}','|'/
c
c       xpm graphics implementation for pltutl
c
c       xpm limited to 4096 colors
c
        if(ncolor<=0) go to 10
c
        ny=600
c***    ny=260
        nx=ny*3/2
c***    nx=ny
        scale=real(ny,rknd)
        iscale=nx
        jscale=ny
        ishift=4096
        nc=1
        lenc=91
        if(ncolor>lenc) nc=2
c
c       initialize bitmap
c
        do k=1,nx*ny
            image(k)=0
        enddo
c
        call fstr(sname,length,fname,0_iknd)
        call ascutl(id,sname,'w',iflag)
c
c       set main definitions
c
        msg='/* XPM */'
        call ascstr(id,msg,9_iknd,iflag)
        msg(1:14)='static char * '
        if(sname(length-3:length)=='.xpm') then
            msg(15:10+length)=sname(1:length-4)
            ll=10+length
        else
            msg(15:14+length)=sname(1:length)
            ll=14+length
        endif
        msg(ll+1:ll+10)='_xpm[] = {'
        call ascstr(id,msg,ll+10_iknd,iflag)
c
        write(unit=msg,fmt='(a1,i4,1x,i4,1x,i4,1x,i1,a2)')
     +      '"',nx,ny,ncolor,nc,'",'
        call ascstr(id,msg,19_iknd,iflag)
c
c       define colors
c
        do i=1,ncolor
            msg='"       c #ffffffffffff",'
            i2=(i-1)/lenc
            i1=i-1-lenc*i2
            cs(1:1)=cdef(i1+1)
            cs(2:2)=cdef(i2+1)
            msg(2:3)=cs
            cmap(i)=cs
            call hexclr(red(i),green(i),blue(i),msg(12:23))
            call ascstr(id,msg,25_iknd,iflag)
        enddo
        return
c
c       print bitmap
c
   10   do j=ny,1,-1
            msg(1:1)='"'
            if(nc==1) then
                do i=1,nx
                    idx=i+(j-1)*iscale
                    ic=image(idx)-(image(idx)/ishift)*ishift+1
                    msg(i+1:i+2)=cmap(ic)
                enddo
            else
                do i=1,nx
                    idx=i+(j-1)*iscale
                    ic=image(idx)-(image(idx)/ishift)*ishift+1
                    msg(2*i:2*i+1)=cmap(ic)
                enddo
            endif
            if(j/=1) then
                msg(nc*nx+2:nc*nx+3)='",'
                call ascstr(id,msg,nc*nx+3_iknd,iflag)
            else
                msg(nc*nx+2:nc*nx+2)='"'
                call ascstr(id,msg,nc*nx+2_iknd,iflag)
            endif
        enddo
        msg(1:2)='};'
        call ascstr(id,msg,2_iknd,iflag)
        call ascutl(id,sname,'c',iflag)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine hexclr(r,g,b,color)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(3) :: ic
            character(len=12) :: color
            character(len=16), save :: hex='0123456789abcdef'
cy
c       translate (r,g,b) colors to hexidecimal
c
        ic(1)=int(r*65535.0e0_rknd)
        ic(2)=int(g*65535.0e0_rknd)
        ic(3)=int(b*65535.0e0_rknd)
        do i=1,3
            jj=max(0,ic(i))
            jj=min(65535,jj)
            do j=1,4
                kk=jj/16
                ii=jj-kk*16
                color(4*i+1-j:4*i+1-j)=hex(ii+1:ii+1)
                jj=kk
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
        subroutine hfutl(ncolor,red,green,blue,fname)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: red,green,blue
            character(len=80) :: fname,fname0,msg
            character(len=80), save  :: sname
            common /hf0/id
cy
        if(ncolor<=0) then
            msg='putl'
            call ascstr(id,msg,4_iknd,iflag)
            call sint(msg,length,ncolor)
            call ascstr(id,msg,length,iflag)
            call ascutl(id,sname,'c',iflag)
            return
        endif
c
        call mkname(fname0,fname)
        call stfile(sname,fname0)
        call ascutl(id,sname,'w',iflag)
        msg='putl'
        call ascstr(id,msg,4_iknd,iflag)
        call sint(msg,length,ncolor)
        call ascstr(id,msg,length,iflag)
        do i=1,ncolor
            write(unit=msg,fmt='(f17.15,1x,f17.15,1x,f17.15)')
     +          red(i),green(i),blue(i)
            call ascstr(id,msg,53_iknd,iflag)
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
        subroutine hframe(iframe)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=80) :: msg
            common /hf0/id
cy
        msg='list'
        call ascstr(id,msg,4_iknd,iflag)
        call sint(msg,length,iframe)
        call ascstr(id,msg,length,iflag)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine hline(x,y,z,n,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(n) :: x,y,z
            character(len=80) :: msg
            common /hf0/id
cy
        msg='line'
        call ascstr(id,msg,4_iknd,iflag)
        call sint(msg,length,icolor)
        call ascstr(id,msg,length,iflag)
        call sint(msg,length,n)
        call ascstr(id,msg,length,iflag)
c
        do i=1,n
            write(unit=msg,fmt='(f17.15,1x,f17.15,1x,f17.15)')
     +          x(i),y(i),z(i)
            call ascstr(id,msg,53_iknd,iflag)
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
        subroutine hfill(x,y,z,n,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(n) :: x,y,z
            character(len=80) :: msg
            common /hf0/id
cy
        msg='fill'
        call ascstr(id,msg,4_iknd,iflag)
        call sint(msg,length,icolor)
        call ascstr(id,msg,length,iflag)
        call sint(msg,length,n)
        call ascstr(id,msg,length,iflag)
        do i=1,n
            write(unit=msg,fmt='(f17.15,1x,f17.15,1x,f17.15)')
     +          x(i),y(i),z(i)
            call ascstr(id,msg,53_iknd,iflag)
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
        subroutine xframe(iframe)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            common /xpm1/scale,fscale,xshift,yshift
cy
c
        if(iframe==2) then
            fscale=scale/2.0e0_rknd
            xshift=scale
            yshift=scale/2.0e0_rknd
        else if(iframe==3) then
            fscale=scale/2.0e0_rknd
            xshift=scale
            yshift=0.0e0_rknd
        else
            fscale=scale
            xshift=0.0e0_rknd
            yshift=0.0e0_rknd
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
        subroutine xline(x,y,z,n,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: x,y,z
            common /xpm1/scale,fscale,xshift,yshift
cy
c       pline for xpm graphics
c
c***    if(scale/=fscale) return
        zshift=fscale*0.01e0_rknd
        ix=int(x(1)*fscale+xshift+0.5e0_rknd)
        iy=int(y(1)*fscale+yshift+0.5e0_rknd)
        iz=int(z(1)*fscale+zshift+0.5e0_rknd)
        do i=2,n
            jx=ix
            jy=iy
            jz=iz
            ix=int(x(i)*fscale+xshift+0.5e0_rknd)
            iy=int(y(i)*fscale+yshift+0.5e0_rknd)
            iz=int(z(i)*fscale+zshift+0.5e0_rknd)
            call iline(ix,iy,iz,jx,jy,jz,icolor)
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
        subroutine iline(nix,niy,niz,njx,njy,njz,ic)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=isngl) image
            common /xpm0/iscale,jscale,ishift,image(540000)
cy
c       update bitmap for a line segment
c
        ix=max(nix,1)
        ix=min(ix,iscale)
        jx=max(njx,1)
        jx=min(jx,iscale)
        iy=max(niy,1)
        iy=min(iy,jscale)
        jy=max(njy,1)
        jy=min(jy,jscale)
        iz=niz
        jz=njz
c
        if(ix/=jx) then
            kmin=min(ix,jx)
            kmax=max(ix,jx)
            do k=kmin,kmax
                x=real((k-ix)*jx+(jx-k)*ix,rknd)/real(jx-ix,rknd)
                y=real((k-ix)*jy+(jx-k)*iy,rknd)/real(jx-ix,rknd)
                z=real((k-ix)*jz+(jx-k)*iz,rknd)/real(jx-ix,rknd)
                kx=max(int(x+0.5e0_rknd),1_iknd)
                kx=min(kx,iscale)
                ky=max(int(y+0.5e0_rknd),1_iknd)
                ky=min(ky,jscale)
                kz=int(z+0.5e0_rknd)
                idx=kx+(ky-1)*iscale
                if(kz>=image(idx)/ishift) image(idx)=kz*ishift+ic-1
            enddo
        endif
        if(iy/=jy) then
            kmin=min(iy,jy)
            kmax=max(iy,jy)
            do k=kmin,kmax
                x=real((k-iy)*jx+(jy-k)*ix,rknd)/real(jy-iy,rknd)
                y=real((k-iy)*jy+(jy-k)*iy,rknd)/real(jy-iy,rknd)
                z=real((k-iy)*jz+(jy-k)*iz,rknd)/real(jy-iy,rknd)
                kx=max(int(x+0.5e0_rknd),1_iknd)
                kx=min(kx,iscale)
                ky=max(int(y+0.5e0_rknd),1_iknd)
                ky=min(ky,jscale)
                kz=int(z+0.5e0_rknd)
                idx=kx+(ky-1)*iscale
                if(kz>=image(idx)/ishift) image(idx)=kz*ishift+ic-1
            enddo
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
        subroutine xfill(x,y,z,n,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: x,y,z
            real(kind=rknd), dimension(200) :: rm,rz
            common /xpm1/scale,fscale,xshift,yshift
cy
c       pfill for xpm graphics
c
c***    if(scale/=fscale) return
        ixmin=int(x(1)*fscale+xshift+0.5e0_rknd)
        ixmax=ixmin
        iymin=int(y(1)*fscale+yshift+0.5e0_rknd)
        iymax=iymin
        do i=2,n
            ix=int(x(i)*fscale+xshift+0.5e0_rknd)
            iy=int(y(i)*fscale+yshift+0.5e0_rknd)
            ixmin=min(ixmin,ix)
            ixmax=max(ixmax,ix)
            iymin=min(iymin,iy)
            iymax=max(iymax,iy)
        enddo
        if(ixmax-ixmin<iymax-iymin) then
c
c       scan by row index
c
            do k=ixmin,ixmax
c
c       find intersections
c
                xx=(real(k,rknd)-xshift)/fscale
                np=0
                nm=0
                num=0
                j=n
                do i=1,n
                    if(x(i)>xx.and.x(j)<=xx) then
                        np=np+1
                    else if(x(i)<=xx.and.x(j)>xx) then
                        nm=nm+1
                    else
                        go to 5
                    endif
                    num=num+1
                    rm(num)=((xx-x(j))*y(i)+(x(i)-xx)*y(j))/(x(i)-x(j))
                    rz(num)=((xx-x(j))*z(i)+(x(i)-xx)*z(j))/(x(i)-x(j))
                    do m=num-1,1,-1
                       if(rm(m)<rm(m+1)) go to 5
                       rr=rm(m)
                       rm(m)=rm(m+1)
                       rm(m+1)=rr
                       rr=rz(m)
                       rz(m)=rz(m+1)
                       rz(m+1)=rr
                    enddo
    5               j=i
                enddo
                if(nm/=np) stop 6123
c
c       update bitmap along line k
c
                do j=1,num,2
                    iy=int(rm(j  )*fscale+yshift+0.5e0_rknd)
                    iz=int(rz(j  )*fscale       +0.5e0_rknd)
                    jy=int(rm(j+1)*fscale+yshift+0.5e0_rknd)
                    jz=int(rz(j+1)*fscale       +0.5e0_rknd)
                    call iline(k,iy,iz,k,jy,jz,icolor)
                enddo
            enddo
        else
c
c       scan by column index
c
            do k=iymin,iymax
c
c       find intersections
c
                yy=(real(k,rknd)-yshift)/fscale
                np=0
                nm=0
                num=0
                j=n
                do i=1,n
                    if(y(i)>yy.and.y(j)<=yy) then
                        np=np+1
                    else if(y(i)<=yy.and.y(j)>yy) then
                        nm=nm+1
                    else
                        go to 10
                    endif
                    num=num+1
                    rm(num)=((yy-y(j))*x(i)+(y(i)-yy)*x(j))/(y(i)-y(j))
                    rz(num)=((yy-y(j))*z(i)+(y(i)-yy)*z(j))/(y(i)-y(j))
                    do m=num-1,1,-1
                       if(rm(m)<rm(m+1)) go to 10
                       rr=rm(m)
                       rm(m)=rm(m+1)
                       rm(m+1)=rr
                       rr=rz(m)
                       rz(m)=rz(m+1)
                       rz(m+1)=rr
                    enddo
   10               j=i
                enddo
                if(nm/=np) stop 6124
c
c       update bitmap along line k
c
                do j=1,num,2
                    ix=int(rm(j  )*fscale+xshift+0.5e0_rknd)
                    iz=int(rz(j  )*fscale       +0.5e0_rknd)
                    jx=int(rm(j+1)*fscale+xshift+0.5e0_rknd)
                    jz=int(rz(j+1)*fscale       +0.5e0_rknd)
                    call iline(ix,k,iz,jx,k,jz,icolor)
                enddo
            enddo
        endif
c
c       trace boundary
c
        ix=int(x(n)*fscale+xshift+0.5e0_rknd)
        iy=int(y(n)*fscale+yshift+0.5e0_rknd)
        iz=int(z(n)*fscale       +0.5e0_rknd)
        do i=1,n
            jx=ix
            jy=iy
            jz=iz
            ix=int(x(i)*fscale+xshift+0.5e0_rknd)
            iy=int(y(i)*fscale+yshift+0.5e0_rknd)
            iz=int(z(i)*fscale       +0.5e0_rknd)
            call iline(ix,iy,iz,jx,jy,jz,icolor)
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
        subroutine ascutl(id,fname,mode,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), save, dimension(10) :: istack
            integer(kind=iknd), save :: length,ifirst=1,next
            character(len=1) :: mode
            character(len=80), save :: sname
            character(len=80) :: fname
            common /asc/maxid,irw(10),iunit(10)
cy
c       iflag= 0  ok
c              1  error on open
c              2  bad mode (not c,r, or w)
c              3  exceed maxid id's
c              4  invalid id
c              5  file not open
c              6  read error
c              7  write error
c             -1  end of file
c
        if(ifirst==1) then
            maxid=10
            do i=1,maxid
                iunit(i)=20+i
                irw(i)=0
                istack(i)=i+1
            enddo
            istack(maxid)=-1
            next=1
            ifirst=0
        endif
        iflag=0
c
c       close
c
        if(mode=='c') then
c
c       ckeck for valid id
c
            if(id<=0.or.id>maxid) then
                iflag=4
                return
            endif
            if(irw(id)==0) then
                iflag=5
                return
            endif
            irw(id)=0
            istack(id)=next
            next=id
            close(unit=iunit(id))
            return
        endif
c
c       get next available id
c
        if(next>0) then
            id=next
            next=istack(id)
        else
c
c       too many files open
c
            iflag=3
            return
        endif
c
c       process filename
c
        call fstr(sname,length,fname,0_iknd)
c
c       open for write
c
        if(mode=='w') then
            open(unit=iunit(id),form='formatted',status='unknown',
     +          file=sname,access='sequential',err=10)
            irw(id)=1
        else if(mode=='r') then
c
c       open for read
c
            open(unit=iunit(id),form='formatted',status='old',
     +          file=sname,access='sequential',err=10)
            irw(id)=-1
        else
            iflag=2
            go to 20
        endif
        return
c
c       if open failed, restore id to available stack
c
   10   iflag=1
   20   irw(id)=0
        istack(id)=next
        next=id
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ascstr(id,sval,length,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=*) :: sval
            common /asc/maxid,irw(10),iunit(10)
cy
c       write a character string
c
c       the long formats are to accomodate xpm files
c       normally should be (80a1)
c
        iflag =0
        if(id<=0.or.id>maxid) then
            iflag=4
            return
        endif
        if(irw(id)==0) then
            iflag=5
            return
        endif
        if(irw(id)<0) then
            read(iunit(id),fmt='(2000a1)',end=10,err=20)
     +          (sval(i:i),i=1,length)
        else
            write(iunit(id),fmt='(2000a1)',err=30)
     +          (sval(i:i),i=1,length)
        endif
        flush(iunit(id))
        return
   10   iflag=-1
        return
   20   iflag=6
        return
   30   iflag=7
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ascint(id,ival,length,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ival
            common /asc/maxid,irw(10),iunit(10)
cy
c       write an integer array
c
        iflag =0
        if(id<=0.or.id>maxid) then
            iflag=4
            return
        endif
        if(irw(id)==0) then
            iflag=5
            return
        endif
        if(irw(id)<0) then
            read(iunit(id),fmt='(6(2x,i11))',end=10,err=20)
     +          (ival(i),i=1,length)
        else
            write(iunit(id),fmt='(6(2x,i11))',err=30)
     +          (ival(i),i=1,length)
        endif
        flush(iunit(id))
        return
   10   iflag=-1
        return
   20   iflag=6
        return
   30   iflag=7
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine ascflt(id,rval,length,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: rval
            common /asc/maxid,irw(10),iunit(10)
cy
c       write a real array
c
        iflag =0
        if(id<=0.or.id>maxid) then
            iflag=4
            return
        endif
        if(irw(id)==0) then
            iflag=5
            return
        endif
        if(irw(id)<0) then
            read(iunit(id),fmt='(3(2x,e23.15))',end=10,err=20)
     +          (rval(i),i=1,length)
        else
            write(iunit(id),fmt='(3(2x,e23.15))',err=30)
     +          (rval(i),i=1,length)
        endif
        flush(iunit(id))
        return
   10   iflag=-1
        return
   20   iflag=6
        return
   30   iflag=7
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine shrtnm(outnam,innam)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=80) :: innam,outnam,temp
cy
c       set file name
c
        call fstr(temp,length,innam,0_iknd)
        do k=length,1,-1
            if(temp(k:k)=='/') go to 10
        enddo
        outnam=temp
        return
   10   outnam=' '
        outnam(1:length-k)=temp(k+1:length)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine reset(num,name,nptr,labels,values,ip,rp,sp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ip,nptr
            integer(kind=isngl), dimension(101) :: sptr
            real(kind=rknd), dimension(*) :: rp
            character(len=15), dimension(*) :: name
            character(len=80), dimension(*) :: sp,labels,values
            character(len=80), dimension(100) :: sval
            character(len=80) :: list,ss,msg
            character(len=1) :: cmd
            character(len=1), dimension(100) :: typ,mark
            character(len=6) :: cmdtyp
            character(len=9), dimension(100) :: tval
            common /atest3/mode,jnlsw,jnlr,jnlw,ibatch
            common /atest4/jcmd,cmdtyp,list
cy
c       reset user paremeters
c
        cmd=list(1:1)
        call lookup(name,num,ip,rp,sp,list,ierr,length)
c
c       print  parameters
c
        if(mode==-1) call disply(name,num,ip,rp,sp)
c
        if(ierr/=0) then
            ss='command error'
            call filutl(ss,0_iknd)
        endif
        if(length>1.and.ierr==0) return
c
c       x-windows display
c
        if(jnlsw==0) then
            do i=1,num
                mark(i)='f'
                sptr(i)=nptr(i)
                call cint(name(i),3_iknd,indx,jerr)
                tval(i)(1:9)=name(i)(5:13)
                if(tval(i)(9:9)==' ') then
                    tval(i)(9:9)=tval(i)(8:8)
                    tval(i)(8:8)=' '
                endif
                typ(i)=name(i)(15:15)
                sval(i)=' '
                if(name(i)(15:15)=='i') then
                    call sint(sval(i),length,ip(indx))
                else if(name(i)(15:15)=='r') then
                    call sreal(sval(i),length,rp(indx),5_iknd,0_iknd)
                else
                    sval(i)=sp(indx)
                endif
            enddo
            sptr(num+1)=nptr(num+1)
c
            if(num==1.and.typ(1)=='f') then
                call xfile(list,sval,tval,jcmd)
                if(sp(indx)/=sval(1)) mark(1)='t'
            else
                call xreset(list,num,typ,sval,mark,tval,
     +              sptr,labels,values,jcmd)
            endif
c
            do i=1,num
                if(mark(i)=='t') then
                    call cint(name(i),3_iknd,indx,jerr)
                    if(name(i)(15:15)=='i') then
                        call cint(sval(i),80_iknd,ival,jerr)
                        if(jerr==0) ip(indx)=ival
                    else if(name(i)(15:15)=='r') then
                        call creal(sval(i),80_iknd,rval,jerr)
                        if(jerr==0) rp(indx)=rval
                    else
                        jerr=0
                        sp(indx)=sval(i)
                    endif
                    if(jerr==0) then
                        ss=' '
                        if(name(i)(15:15)=='l') then
                            call fstr(ss,length,sval(i),1_iknd)
                        else
                            call fstr(ss,length,sval(i),0_iknd)
                        endif
                        write(unit=msg,fmt='(a1,a6,a1,a72)')
     +                      cmd,name(i)(5:10),'=',ss
                        call star0(msg)
                        call filutl(msg,1_iknd)
                    endif
                endif
            enddo
c
        else if(jnlsw/=-2) then
            call getcmd(list)
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
        subroutine menu(ip,rp,sp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), save, dimension(24) :: st
            integer(kind=iknd), dimension(*) :: ip
            integer(kind=iknd), save, dimension(525) :: iptr
            integer(kind=iknd), save, dimension(301) :: nptr
            integer(kind=iknd), dimension(101) :: sptr
            integer(kind=iknd), save :: iustat,lowera=97,lowerz=122
            integer(kind=iknd), save :: mpibtn=0,ifirst=-1,ncmd
            real(kind=rknd), dimension(*) :: rp
            character(len=1), save, dimension(24) :: sty
            character(len=6) :: cmdtyp
            character(len=15), save, dimension(300) :: name
            character(len=15), dimension(100) :: sname
            character(len=15), save, dimension(24) :: ctable
            character(len=80), dimension(100) :: sp
            character(len=80) :: list,filnam
            character(len=80), save :: ulist
            character(len=80), dimension(500) :: file
            character(len=80), save, dimension(500) :: labels,values
            character(len=80), dimension(200) :: slabel,svalue
            common /atest3/mode,jnlsw,jnlr,jnlw,ibatch
            common /atest4/jcmd,cmdtyp,list
            common /atest6/nproc,myid,mpisw,mpiint,mpiflt
cy
c
        if(ifirst==-1) then
            call mpiutl(1_iknd)
            mode=0
            do i=1,100
                ip(i)=0
                rp(i)=0.0e0_rknd
                sp(i)=' '
            enddo
            call gtfile(file,len)
            call mkcmd(file,len,name,nlen,nptr,labels,values,
     +          ncmd,ctable,st,sty,iptr,ip,rp,sp)
c*************************
cc      call prtfl(file,len)
cc      call prtfl0(file,len)
c**********************
            ip(42)=mode
            ip(48)=mpisw
            ip(49)=nproc
            ip(50)=myid+1
            ifirst=1
            return
        endif
        sp(12)=' '
        if(ifirst==1) then
            jcmd=ncmd
            sp(12)(1:6)='quit  '
            mode=ip(42)
            jnlsw=mode
            if(mode>1.or.mode<-1) then
                sp(11)='menu: bad value for mode'
                return
            endif
            if(myid/=0) then
                if(mode==1) then
                    call mkjnl(sp,kflag)
                    if(kflag/=0) go to 40
                    call stfile(filnam,sp(10))
                    call ascutl(jnlr,filnam,'r',kflag)
                    if(kflag/=0) go to 40
                    jnlsw=2
                    mode=-2
                else
                    mode=-2
                    jnlsw=mode
                endif
            else if(mode==0) then
                call xwinit(ncmd,ctable,sp(13))
                call grinit(ip(43))
                do i=1,ncmd
                    if(ctable(i)(10:15)=='mpicmd')  then
                        mpibtn=i
                        call xmpi(mpisw,mpibtn)
                    endif
                enddo
            else if(mode==1) then
                call mkjnl(sp,kflag)
                if(kflag/=0) go to 40
                call stfile(filnam,sp(10))
                call ascutl(jnlr,filnam,'r',kflag)
                if(kflag/=0) go to 40
            endif
            call stfile(filnam,sp(8))
            call ascutl(jnlw,filnam,'w',kflag)
            if(kflag/=0) go to 40
            call stfile(filnam,sp(9))
            call ascutl(ibatch,filnam,'w',kflag)
            if(kflag/=0) go to 40
c
            ulist=' '
            list=' '
            iustat=0
            ifirst=0
        endif
c
        ierr=0
    5   if(ierr>0) sp(11)='command error'
        if(sp(11)/=' ') call filutl(sp(11),0_iknd)
        if(iustat==1.and.ulist==list) iustat=0
        if(iustat==0) then
            if(jnlsw==0) then
                call xgtcmd(list)
            else if(jnlsw/=-2) then
                call getcmd(list)
            endif
        endif
c
c       mpi communication
c
        if(mode==-2.and.jnlsw==-2) then
            call star0(list)
            call parcmd(ncmd,ctable,list,length,nequal,
     +          jcmd,cmdtyp,ierr)
        else
            call parcmd(ncmd,ctable,list,length,nequal,
     +          jcmd,cmdtyp,ierr)
            call star0(list)
        endif
c
        iustat=0
        sp(11)=' '
        if(ierr/=0) go to 5
        if(length==0) then
            if(mode==-1) call discmd(ncmd,ctable)
            ierr=0
            go to 5
        endif
c
c       quit and mpicmd are always executed by all processors
c
        if(mode==-2.and.jnlsw>=1.and.mpisw==-1) then
            if(cmdtyp/='mpicmd'.and.cmdtyp/='quit  ') then
                ii=ichar(list(1:1))
                if(ii>=lowera.and.ii<=lowerz) list(1:1)=char(ii-32)
                if(length<=1) go to 5
            endif
        endif
        if(list(1:1)==ctable(jcmd)(8:8)) go to 30
c
c       reset parameters with display
c
        iustat=1
        ulist=list
        if(nequal==0.and.st(jcmd)>0.and.length>1) then
            call shrtfm(ip,rp,sp,length,sty,st,ierr)
        else
            num=iptr(jcmd+1)-iptr(jcmd)
            call mktabl(jcmd,name,iptr,sname,
     +          nptr,labels,values,sptr,slabel,svalue)
            call reset(num,sname,sptr,slabel,svalue,ip,rp,sp)
            ierr=0
        endif
        if(ctable(jcmd)(1:6)=='mpicmd') ip(48)=mpisw
        sp(11)=' '
        go to 5
c
   30   sp(12)(1:6)=ctable(jcmd)(1:6)
        if(length==1) go to 40
c
c       short form of command
c
        if(nequal==0.and.st(jcmd)>0) then
            call shrtfm(ip,rp,sp,length,sty,st,ierr)
c
c       long form of command
c
        else
            num=iptr(jcmd+1)-iptr(jcmd)
            call mktabl(jcmd,name,iptr,sname,
     +          nptr,labels,values,sptr,slabel,svalue)
            call lookup(sname,num,ip,rp,sp,list,ierr,length)
        endif
c
        if(ierr/=0) go to 5
c
c       quit command
c
   40   if(sp(12)(1:6)=='quit  '.or.cmdtyp=='quit  ') then
            call mpiutl(-1_iknd)
            jcmd=-1
            if(mode==0) call xwinit(jcmd,ctable,sp(13))
            if(jnlsw>=1) call ascutl(jnlr,filnam,'c',kflag)
            call ascutl(jnlw,filnam,'c',kflag)
            call ascutl(ibatch,filnam,'c',kflag)
c
c       journal command
c
        else if(cmdtyp=='journl') then
            ierr=0
            if(jnlsw<=0) then
                call mkjnl(sp,kflag)
                if(kflag/=0) go to 5
                call stfile(filnam,sp(10))
                call ascutl(jnlr,filnam,'r',kflag)
                if(kflag/=0) then
                    sp(11)='journl: cannot open file'
                else
                    sp(11)='journl: ok'
                    jnlsw=1
                endif
                go to 5
            else
                go to 5
            endif
c
c       user command
c
        else if(cmdtyp=='usrcmd') then
            iustat=1
            ulist=list
            sp(11)='usrcmd: ok'
c
c       mpi command
c
        else if(cmdtyp=='mpicmd') then
            if(length==1) then
                mpisw=-mpisw
                ip(48)=mpisw
            else
                if(ip(48)/=1) ip(48)=-1
                mpisw=ip(48)
            endif
            if(mpisw==1) then
                sp(11)='mpi is on'
            else
                sp(11)='mpi is off'
            endif
            ierr=0
            if(mode==0) call xmpi(mpisw,mpibtn)
            go to 5
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
        subroutine parcmd(ncmd,ctable,list,length,nequal,
     +      jcmd,cmdtyp,ierr)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(24) :: lequal,lcomma
            character(len=1) :: lcmd,ucmd
            character(len=6) :: cmdtyp
            character(len=15), dimension(*) :: ctable
            character(len=80) :: list
cy
        call fxcase(list,length,ncomma,nequal,ndbleq,
     +      lcomma,lequal,icomnt)
c
c       obvious errors
c
        call filutl(list,1_iknd)
        if(length==0) then
            ierr=0
            return
        endif
        ierr=1
        jcmd=0
        cmdtyp=' '
        if(icomnt==1) then
            ierr=-1
            return
        endif
        if(nequal>0) then
            if(ncomma/=nequal-1) return
        else
            if(ncomma>0) return
        endif
        if((ndbleq/2)*2/=ndbleq) return
c
c       find command code
c
        do icmd=1,ncmd
            lcmd=ctable(icmd)(8:8)
            ii=ichar(lcmd)-32
            ucmd=char(ii)
            if(lcmd==list(1:1).or.ucmd==list(1:1)) go to 20
        enddo
        return
   20   if(lcmd==list(1:1)) cmdtyp=ctable(icmd)(10:15)
        jcmd=icmd
        ierr=0
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine shrtfm(ip,rp,sp,length,sty,st,ierr)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: st,ip
            real(kind=rknd), dimension(*) :: rp
            character(len=1), dimension(*) :: sty
            character(len=6) :: cmdtyp
            character(len=80), dimension(100) :: sp
            character(len=80) :: list
            common /atest4/jcmd,cmdtyp,list
cy
c       short form of command
c
        ierr=0
        ll=length-1
        if(sty(jcmd)=='i') then
            call cint(list(2:2),ll,ival,ierr)
            if(ierr==0) ip(st(jcmd))=ival
        else if(sty(jcmd)=='r') then
            call creal(list(2:2),ll,rval,ierr)
            if(ierr==0) rp(st(jcmd))=rval
        else if(sty(jcmd)=='l') then
            sp(st(jcmd))=' '
            sp(st(jcmd))(1:ll-2)=list(3:length-1)
        else
            sp(st(jcmd))=' '
            sp(st(jcmd))(1:ll)=list(2:length)
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
        subroutine prtfl(file,len)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(24) :: lcomma,lequal
            character(len=80), dimension(*) :: file
            character(len=80) :: lstr,line
            character(len=10), save :: mark='+123456789'
            character(len=1), save, dimension(4) :: cc
cy
            data cc/'n','c','r','s'/
c
c       get rid of comments, blank lines and spaces
c
        ishift=0
        do i=1,len
            lstr=file(i)
            call fxcase(lstr,length,ncomma,nequal,ndbleq,
     +          lcomma,lequal,icomnt)
            if(icomnt==1.or.length==0) then
                ishift=ishift+1
            else
                file(i-ishift)=' '
                file(i-ishift)(1:length)=lstr(1:length)
            endif
        enddo
        len=len-ishift
c
c
c
        ii=1
        k=1
        is=8
        do m=1,4
            do i=1,len
                if(file(i)(1:1)/=cc(m)) cycle
                call fstr(lstr,length,file(i),0_iknd)
                line=' '
                line(6:6)=mark(ii:ii)
                line(is+1:is+1)=char(39)
                ll=is+1+length
                line(is+2:ll)=lstr(1:length)
                line(ll+1:ll+1)=char(39)
                if(ii/=10.and.k<len) then
                    line(ll+2:ll+2)=','
                else
                    line(ll+2:ll+2)='/'
                endif
                if(ii==1) then
                    k9=min(k+9,len)
                    write(unit=10,fmt='(12x,a17,i3,a1,i3,a2)')
     +              'data (file0(i),i=',k,',',k9,')/'
                endif
                write(unit=10,fmt='(a80)') line
                k=k+1
                ii=ii+1
                if(ii>10) ii=1
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
        subroutine prtfl0(file,len)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(24) :: lcomma,lequal
            character(len=80), dimension(*) :: file
            character(len=80) :: lstr
            character(len=1), save, dimension(4) :: cc
cy
            data cc/'n','c','r','s'/
c
c       get rid of comments, blank lines and spaces
c
        ishift=0
        do i=1,len
            lstr=file(i)
            call fxcase(lstr,length,ncomma,nequal,ndbleq,
     +          lcomma,lequal,icomnt)
            if(icomnt==1.or.length==0) then
                ishift=ishift+1
            else
                file(i-ishift)=' '
                file(i-ishift)(1:length)=lstr(1:length)
            endif
        enddo
        len=len-ishift
c
c
        do m=1,4
            do i=1,len
                if(file(i)(1:1)/=cc(m)) cycle
                call fstr(lstr,length,file(i),0_iknd)
                write(unit=11,fmt='(a80)') lstr
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
        subroutine getnam(name,nlen)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(24) :: lcomma,lequal
            integer(kind=iknd), dimension(2) :: ig
            real(kind=rknd), dimension(2) :: rg
            character(len=15), dimension(*) :: name
            character(len=15), save, dimension(20) :: name0
            character(len=80) :: lstr
            character(len=80), dimension(500) :: file
            character(len=80), dimension(5) :: sg
cy
            data (name0(i),i= 1, 5)/
     +      '  1 index  i  s','  2 vname  n  s','  3 alias  a  s',
     1      '  4 vtype  t  s','  5 deflt  d  l'/
c
c
        call gtfile(file,len)
        nlen=0
        do i=1,len
            lstr=file(i)
            call fxcase(lstr,length,ncomma,nequal,ndbleq,
     +          lcomma,lequal,icomnt)
            if(icomnt==1.or.length==0) cycle
            if(lstr(1:1)/='n') cycle
c
            nlen=nlen+1
            name(nlen)=' '
            do j=1,5
                sg(j)=' '
            enddo
            call lookup(name0,5_iknd,ig,rg,sg,lstr,ierr,length)
            name(nlen)(1:3)=sg(1)(1:3)
            name(nlen)(5:10)=sg(2)(1:6)
            name(nlen)(12:13)=sg(3)(1:2)
            name(nlen)(15:15)=sg(4)(1:1)
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
        subroutine mkcmd(file,len,name,nlen,nptr,labels,values,
     +      ncmd,ctable,st,sty,iptr,ip,rp,sp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: iptr,st,nptr
            integer(kind=iknd), dimension(24) :: lcomma,lequal,inum,snum
            integer(kind=iknd), dimension(300) :: num
            integer(kind=iknd), dimension(2) :: ig
            integer(kind=iknd), dimension(500) :: jv
            integer(kind=iknd), dimension(200) :: iv,ic
            integer(kind=iknd), dimension(100) :: ip
            integer(kind=iknd), save :: mxnam=300,mxcmd=24,
     +          mxvar=500,mxlst=500
            real(kind=rknd), dimension(2) :: rg
            real(kind=rknd), dimension(100) :: rp
            character(len=1) :: typ,jtyp,uppera,upperz,cc
            character(len=1), dimension(*) :: sty
            character(len=15), dimension(*) :: name,ctable
            character(len=15), save, dimension(20) :: name0
            character(len=15) :: ntemp
            character(len=80) :: lstr
            character(len=80), dimension(*) :: labels,values,file
            character(len=80), dimension(5) :: sg
            character(len=80), dimension(100) :: sp
            character(len=80), dimension(500) :: l0,v0
cy
            data (name0(i),i= 1, 14)/
     +      '  1 index  i  s','  2 vname  n  s','  3 alias  a  s',
     1      '  4 vtype  t  s','  5 deflt  d  l','  1 cname  c  s',
     2      '  2 cmdkey k  s','  3 ctype  t  s','  1 cname  c  s',
     3      '  2 vname  n  s','  3 short  s  s','  1 vname  n  s',
     4      '  2 value  v  s','  3 label  l  l'/
c
c       get rid of comments, blank lines and spaces
c
        ishift=0
        uppera=char(65)
        upperz=char(90)
        do i=1,len
            lstr=file(i)
            call fxcase(lstr,length,ncomma,nequal,ndbleq,
     +          lcomma,lequal,icomnt)
            if(icomnt==1.or.length==0) then
                ishift=ishift+1
            else
                cc=lstr(1:1)
                if(cc>=uppera.and.cc<=upperz) then
                    ii=ichar(cc)+32
                    lstr(1:1)=char(ii)
                endif
                file(i-ishift)=' '
                file(i-ishift)(1:length)=lstr(1:length)
            endif
        enddo
        len=len-ishift
c
c       name and ctable
c
        ncmd=0
        ilen=0
        nlen=0
        do i=1,len
c
c       name
c
            if(file(i)(1:1)=='n') then
                nlen=nlen+1
                if(nlen>mxnam) stop 3001
                name(nlen)=' '
                do j=1,5
                    sg(j)=' '
                enddo
                call lookup(name0(1),5_iknd,ig,rg,sg,file(i),
     +              ierr,length)
                name(nlen)(1:3)=sg(1)(1:3)
                name(nlen)(5:10)=sg(2)(1:6)
                name(nlen)(12:13)=sg(3)(1:2)
                typ=sg(4)(1:1)
                name(nlen)(15:15)=typ
                if(typ=='i'.or.typ=='r'.or.typ=='s') ilen=ilen+1
                if(sg(5)/=' ') then
                    call cint(sg(1),3_iknd,indx,ierr)
                    call fstr(lstr,length,sg(5),0_iknd)
                    if(typ=='i') then
                        call cint(lstr,length,ip(indx),ierr)
                    else if(typ=='r') then
                        call creal(lstr,length,rp(indx),ierr)
                    else
                        sp(indx)=' '
                        sp(indx)(1:length)=lstr(1:length)
                    endif
                endif
c
c       command
c
            else if(file(i)(1:1)=='c') then
                ncmd=ncmd+1
                if(ncmd>mxcmd) stop 3002
                ctable(ncmd)=' '
                do j=1,3
                    sg(j)=' '
                enddo
                call lookup(name0(6),3_iknd,ig,rg,sg,file(i),
     +              ierr,length)
                ctable(ncmd)(1:6)=sg(1)(1:6)
                ctable(ncmd)(8:8)=sg(2)(1:1)
                ctable(ncmd)(10:15)=sg(3)(1:6)
            endif
        enddo
c
c       sort
c
        nn=ilen+1
        do 5 i=1,ilen
            typ=name(i)(15:15)
            if(typ=='i'.or.typ=='r'.or.typ=='s') cycle
            do j=nn,nlen
                jtyp=name(j)(15:15)
                if(jtyp=='i'.or.jtyp=='r'.or.jtyp=='s') then
                    ntemp=name(i)
                    name(i)=name(j)
                    name(j)=ntemp
                    nn=j+1
                    go to 5
                endif
            enddo
            stop 9413
    5   continue
c
c       iptr, nptr
c
        do i=1,nlen
            num(i)=0
        enddo
        do i=1,ncmd
            inum(i)=0
            snum(i)=0
        enddo
c
        ilen=0
        jlen=0
        do i=1,len
c
c       reset variable
c
            if(file(i)(1:1)=='r') then
                ilen=ilen+1
                if(ilen>mxvar) stop 3003
                do j=1,3
                    sg(j)=' '
                enddo
                call lookup(name0(9),3_iknd,ig,rg,sg,file(i),
     +              ierr,length)
                do j=1,ncmd
                    if(sg(1)(1:6)==ctable(j)(1:6)) go to 10
                enddo
                stop 1001
   10           ic(ilen)=j
                do j=1,nlen
                    if(sg(2)(1:6)==name(j)(5:10)) go to 20
                enddo
cc              write(6,*) file(i)
                stop 1002
   20           iv(ilen)=j
                if(sg(3)(1:1)=='1') iv(ilen)=-j
                typ=name(j)(15:15)
                if(typ=='i'.or.typ=='r'.or.typ=='s') then
                    inum(ic(ilen))=inum(ic(ilen))+1
                else
                    snum(ic(ilen))=snum(ic(ilen))+1
                endif
c
c       switch
c
            else if(file(i)(1:1)=='s') then
                jlen=jlen+1
                if(jlen>mxlst) stop 3004
                do j=1,3
                    sg(j)=' '
                enddo
                call lookup(name0(12),3_iknd,ig,rg,sg,file(i),
     +              ierr,length)
                do j=1,nlen
                    if(sg(1)(1:6)==name(j)(5:10)) go to 30
                enddo
                stop 1003
   30           jv(jlen)=j
                v0(jlen)=sg(2)
                l0(jlen)=sg(3)
                num(j)=num(j)+1
            endif
        enddo
c
c       compute start of iptr
c
        iptr(1)=ncmd+2
        do i=1,ncmd
            iptr(i+1)=iptr(i)+inum(i)+snum(i)
            snum(i)=iptr(i)+inum(i)
            inum(i)=iptr(i)
            st(i)=0
        enddo
c
c       compute the rest of iptr
c
        do i=1,ilen
            icmd=ic(i)
            ivar=abs(iv(i))
            typ=name(ivar)(15:15)
            if(typ=='i'.or.typ=='r'.or.typ=='s') then
                k=inum(icmd)
                inum(icmd)=k+1
            else
                k=snum(icmd)
                snum(icmd)=k+1
            endif
            iptr(k)=ivar
            if(iv(i)<0) then
                call cint(name(ivar),3_iknd,indx,jerr)
                st(icmd)=indx
                sty(icmd)=typ
            endif
        enddo
c
c       compute nptr
c
        nptr(1)=1
        do i=1,nlen
            nptr(i+1)=nptr(i)+num(i)
            num(i)=nptr(i)
        enddo
c
c       compute labels and values
c
        do i=1,jlen
            ivar=jv(i)
            k=num(ivar)
            num(ivar)=k+1
            labels(k)=l0(i)
            values(k)=v0(i)
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
        subroutine usrset(file,len,ip,rp,sp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(5) :: iptr,st
            integer(kind=iknd), dimension(301) :: nptr
            integer(kind=iknd), dimension(100) :: ip
            real(kind=rknd), dimension(100) :: rp
            character(len=1), dimension(5) :: sty
            character(len=15), dimension(300) :: name
            character(len=15), dimension(5) :: ctable
            character(len=80), dimension(500) :: labels,values
            character(len=80), dimension(*) :: file
            character(len=80), dimension(100) :: sp
cy
c       mkcmd interface for usrcmd
c
        if(len>500) return
        call mkcmd(file,len,name,nlen,nptr,labels,values,
     +      ncmd,ctable,st,sty,iptr,ip,rp,sp)
        call reset(nlen,name,nptr,labels,values,ip,rp,sp)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine gtfile(file,len)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=80), save, dimension(500) :: file0
            character(len=80), dimension(*) :: file
cy
            data (file0(i),i=  1, 10)/
     +  'ni=1,n=n,t=i',
     1  'ni=2,n=lvl,t=i',
     2  'ni=3,n=nblock,t=i,d="1"',
     3  'ni=8,n=ispd,t=i,a=i,d="0"',
     4  'ni=9,n=method,t=i,a=m,d="0"',
     5  'ni=10,n=mxcg,t=i,a=c,d="20"',
     6  'ni=12,n=maxlvl,t=i,a=ml,d="20"',
     7  'ni=13,n=maxfil,t=i,a=mf,d="20"',
     8  '#ni=20,n=lenw,t=i',
     9  'ni=21,n=maxn,t=i'/
            data (file0(i),i= 11, 20)/
     +  'ni=22,n=maxja,t=i',
     1  'ni=23,n=maxa,t=i',
     2  'ni=26,n=iflag,t=i',
     3  'ni=42,n=mode,t=i,d="0"',
     4  'ni=43,n=ngraph,t=i,d="0"',
     5  'ni=45,n=gdevce,t=i,a=d,d="1"',
     6  'ni=47,n=mdevce,t=i,a=d,d="2"',
     7  'ni=51,n=mxcolr,t=i,a=mc,d="100"',
     8  'ni=54,n=igrsw,t=i,a=i,d="0"',
     9  'ni=55,n=imtxsw,t=i,a=i,d="2"'/
            data (file0(i),i= 21, 30)/
     +  'ni=56,n=ncon,t=i,a=c,d="11"',
     1  'ni=58,n=iscale,t=i,a=s,d="0"',
     2  'ni=59,n=lines,t=i,a=l,d="0"',
     3  'ni=60,n=numbrs,t=i,a=n,d="0"',
     4  'ni=64,n=mx,t=i,a=mx,d="1"',
     5  'ni=65,n=my,t=i,a=my,d="1"',
     6  'ni=66,n=mz,t=i,a=mz,d="1"',
     7  'ni=67,n=level,t=i,a=l,d="0"',
     8  'ni=71,n=lenja0,t=i',
     9  'ni=72,n=lena0,t=i'/
            data (file0(i),i= 31, 40)/
     +  'ni=73,n=fil0,t=i',
     1  'ni=74,n=ngrid,t=i,a=n,d="10"',
     2  'ni=75,n=mtxtyp,t=i,a=t,d="0"',
     3  'ni=1,n=eps,t=r,a=e,d="1.0e-6"',
     4  'ni=2,n=relerr,t=r',
     5  'ni=3,n=dtol,t=r,a=d,d="1.0e-2"',
     6  'ni=6,n=smin,t=r,a=sn,d="0.0e0"',
     7  'ni=7,n=smax,t=r,a=sx,d="0.0e0"',
     8  'ni=8,n=rmag,t=r,a=m,d="1.0e0"',
     9  'ni=9,n=cenx,t=r,a=cx,d="0.5e0"'/
            data (file0(i),i= 41, 50)/
     +  'ni=10,n=ceny,t=r,a=cy,d="0.5e0"',
     1  'ni=26,n=mginit,t=r',
     2  'ni=27,n=mg,t=r',
     3  'ni=3,n=gtitle,t=l,a=t,d="gphplt"',
     4  'ni=4,n=mtitle,t=l,a=t,d="mtxplt"',
     5  'ni=5,n=shcmd,t=l,a=c',
     6  'ni=6,n=rwfile,t=f,a=f,d="mtx001"',
     7  'ni=7,n=jrfile,t=f,a=f,d="mgraph.jnl"',
     8  'ni=8,n=jwfile,t=f,d="journl.jnl"',
     9  'ni=9,n=bfile,t=f,d="output.out"'/
            data (file0(i),i= 51, 60)/
     +  'ni=10,n=jtfile,t=f,d="jnltmp_mpixxx.jnl"',
     1  'ni=11,n=iomsg,t=l',
     2  'ni=12,n=cmd,t=s',
     3  'ni=13,n=logo,t=l,d="multigraph"',
     4  'ni=14,n=bgclr,t=l,d="gray85"',
     5  'ni=15,n=btnbg,t=l,d="gray30"',
     6  'ni=18,n=psfile,t=f,d="figxxx.ps"',
     7  'ni=19,n=xpfile,t=f,d="figxxx.xpm"',
     8  'ni=20,n=bhfile,t=f,d="figxxx.bh"',
     9  'ni=48,n=mpisw,t=i,d="-1"'/
            data (file0(i),i= 61, 70)/
     +  'ni=49,n=nproc,t=i,d="1"',
     1  'ni=50,n=irgn,t=i,d="1"',
     2  'ni=21,n=sghost,t=f,d="localhost"',
     3  'ni=11,n=ncfact,t=i,a=nc,d="4"',
     4  'cc=factor,k=f,t=popup',
     5  'cc=solve,k=s,t=popup',
     6  'cc=gphplt,k=g,t=popup',
     7  'cc=mtxplt,k=m,t=popup',
     8  'cc=linsys,k=l,t=popup',
     9  'cc=read,k=r,t=file'/
            data (file0(i),i= 71, 80)/
     +  '#cc=lstsq,k=a,t=popup',
     1  'cc=journl,k=j,t=journl',
     2  '#cc=shell,k=k,t=popup',
     3  'cc=quit,k=q,t=quit',
     4  'rc=factor,n=maxlvl',
     5  'rc=factor,n=maxfil',
     6  'rc=factor,n=method',
     7  'rc=factor,n=dtol',
     8  'rc=solve,n=mxcg',
     9  'rc=solve,n=eps'/
            data (file0(i),i= 81, 90)/
     +  'rc=gphplt,n=igrsw,s=1',
     1  'rc=gphplt,n=gdevce',
     2  'rc=gphplt,n=mxcolr',
     3  'rc=gphplt,n=gtitle',
     4  'rc=mtxplt,n=imtxsw,s=1',
     5  'rc=mtxplt,n=iscale',
     6  'rc=mtxplt,n=lines',
     7  'rc=mtxplt,n=numbrs',
     8  'rc=mtxplt,n=mdevce',
     9  'rc=mtxplt,n=mx'/
            data (file0(i),i= 91,100)/
     +  'rc=mtxplt,n=my',
     1  'rc=mtxplt,n=mz',
     2  'rc=mtxplt,n=ncon',
     3  'rc=mtxplt,n=level',
     4  'rc=mtxplt,n=mxcolr',
     5  'rc=mtxplt,n=smin',
     6  'rc=mtxplt,n=smax',
     7  'rc=mtxplt,n=rmag',
     8  'rc=mtxplt,n=cenx',
     9  'rc=mtxplt,n=ceny'/
            data (file0(i),i=101,110)/
     +  'rc=mtxplt,n=mtitle',
     1  'rc=read,n=rwfile,s=1',
     2  'rc=journl,n=jrfile,s=1',
     3  '#rc=shell,n=shcmd,s=1',
     4  'rc=linsys,n=ngrid',
     5  'rc=linsys,n=mtxtyp',
     6  'rc=factor,n=ncfact',
     7  'sn=method,v=0,l="drop tol"',
     8  'sn=method,v=1,l="ilu"',
     9  'sn=method,v=2,l="sgs"'/
            data (file0(i),i=111,120)/
     +  'sn=gdevce,v=0,l="socket 0"',
     1  'sn=gdevce,v=1,l="socket 1"',
     2  'sn=gdevce,v=2,l="socket 2"',
     3  'sn=gdevce,v=3,l="socket 3"',
     4  'sn=gdevce,v=4,l="bh file"',
     5  'sn=gdevce,v=5,l="ps file"',
     6  'sn=gdevce,v=6,l="xpm file"',
     7  'sn=gdevce,v=7,l="popup 0"',
     8  'sn=gdevce,v=8,l="popup 1"',
     9  'sn=gdevce,v=9,l="popup 2"'/
            data (file0(i),i=121,130)/
     +  'sn=gdevce,v=10,l="popup 3"',
     1  'sn=mdevce,v=0,l="socket 0"',
     2  'sn=mdevce,v=1,l="socket 1"',
     3  'sn=mdevce,v=2,l="socket 2"',
     4  'sn=mdevce,v=3,l="socket 3"',
     5  'sn=mdevce,v=4,l="bh file"',
     6  'sn=mdevce,v=5,l="ps file"',
     7  'sn=mdevce,v=6,l="xpm file"',
     8  'sn=mdevce,v=7,l="popup 0"',
     9  'sn=mdevce,v=8,l="popup 1"'/
            data (file0(i),i=131,140)/
     +  'sn=mdevce,v=9,l="popup 2"',
     1  'sn=mdevce,v=10,l="popup 3"',
     2  'sn=igrsw,v=0,l="convergence history"',
     3  'sn=igrsw,v=1,l="levels"',
     4  'sn=igrsw,v=-1,l="time pie chart"',
     5  'sn=igrsw,v=2,l="ip array"',
     6  'sn=igrsw,v=-2,l="rp array"',
     7  'sn=igrsw,v=-3,l="sp array"',
     8  'sn=igrsw,v=3,l="ka array"',
     9  'sn=imtxsw,v=1,l="|ilu| by type"'/
            data (file0(i),i=141,150)/
     +  'sn=imtxsw,v=-1,l="ilu by type"',
     1  'sn=imtxsw,v=2,l="|ilu| by value"',
     2  'sn=imtxsw,v=-2,l="ilu by value"',
     3  'sn=imtxsw,v=3,l="|mtx a| by type"',
     4  'sn=imtxsw,v=-3,l="mtx a by type"',
     5  'sn=imtxsw,v=4,l="|mtx a| by value"',
     6  'sn=imtxsw,v=-4,l="mtx a by value"',
     7  'sn=imtxsw,v=5,l="|error| by type"',
     8  'sn=imtxsw,v=-5,l="error by type"',
     9  'sn=imtxsw,v=6,l="|error| by value"'/
            data (file0(i),i=151,160)/
     +  'sn=imtxsw,v=-6,l="error by value"',
     1  'sn=iscale,v=0,l="linear"',
     2  'sn=iscale,v=1,l="log"',
     3  'sn=iscale,v=2,l="arcsinh"',
     4  'sn=lines,v=0,l="no lines"',
     5  'sn=lines,v=-2,l="matrix elements"',
     6  'sn=numbrs,v=0,l="none"',
     7  'sn=numbrs,v=-1,l="mtx value"',
     8  'sn=numbrs,v=-2,l="mtx index"',
     9  'sn=mtxtyp,v=0,l="star 5"'/
            data (file0(i),i=161,165)/
     +  'sn=mtxtyp,v=1,l="|star 5|"',
     1  'sn=mtxtyp,v=2,l="star 7"',
     2  'sn=mtxtyp,v=3,l="stokes"',
     3  'sn=mtxtyp,v=4,l="star 9"',
     4  'sn=mtxtyp,v=5,l="|star 9|"'/
c
        data len0/165/
c
        len=len0
        do i=1,len
            file(i)=file0(i)
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
        subroutine mktabl(icmd,name,iptr,sname,
     +      nptr,labels,values,sptr,slabel,svalue)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: iptr,nptr,sptr
            character(len=15), dimension(*) :: name,sname
            character(len=80), dimension(*) :: labels,values,slabel
            character(len=80), dimension(*) :: svalue
cy
c       compute sname, sptr, slabel, svalue
c
        sptr(1)=1
        do i=iptr(icmd),iptr(icmd+1)-1
            k=i+1-iptr(icmd)
            nl=iptr(i)
            sname(k)=name(nl)
            ii=sptr(k)
            do j=nptr(nl),nptr(nl+1)-1
                slabel(ii)=labels(j)
                svalue(ii)=values(j)
                ii=ii+1
            enddo
            sptr(k+1)=ii
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
        subroutine getcmd(list)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=80) :: list
            common /atest3/mode,jnlsw,jnlr,jnlw,ibatch
cy
c       get the next command from
c       the tty or the command file
c
c       jnlsw  > 0    get command from journal file
c              = 0    get command from x-windows interface
c              < 0    get command for terminal window
c
c
c       print a prompt symbol
c
        if(jnlsw<0) then
            call crtutl(list,'r','command:')
        else if(jnlsw>0) then
            call ascstr(jnlr,list,80_iknd,kflag)
            if(kflag/=0) then
                call ascutl(jnlr,list,'c',kflag)
                if(mode==1.or.jnlsw==2) then
                    list='q'
                else
                    list=' '
                    jnlsw=mode
                endif
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
        subroutine lookup(name,num,ip,rp,sp,list,ierr,length)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(24) :: ival,lequal,lcomma,iptr
            integer(kind=iknd), dimension(*) :: ip
            real(kind=rknd), dimension(*) :: rp
            real(kind=rknd), dimension(24) :: rval
            character(len=80), dimension(*) :: sp
            character(len=80), dimension(24) :: sval
            character(len=15), dimension(*) :: name
            character(len=6) :: lname
            character(len=2) :: sname
            character(len=*) :: list
            character(len=1), save :: dbleq='"'
cy
c       determine number of entries
c
        ierr=0
        call fxcase(list,length,ncomma,nequal,ndbleq,
     +      lcomma,lequal,icomnt)
        if(num<=0.or.length==1) return
        ierr=1
        if(icomnt==1.or.(ndbleq/2)*2/=ndbleq) return
        if(nequal==0.or.ncomma/=nequal-1) return
        if(ncomma>0) then
            do i=1,ncomma
                if(lcomma(i)<lequal(i).or.
     +              lcomma(i)>lequal(i+1)) return
            enddo
        endif
c
c       the main loop
c
        do ii=1,nequal
            lname=' '
            sname=' '
            istart=2
            imid=lequal(ii)
            iend=length
            if(ii>1) istart=lcomma(ii-1)+1
            if(ii<nequal) iend=lcomma(ii)-1
            if(iend<=imid) return
            if(istart>=imid) return
            if(istart+6<imid) return
c
c       search name array for the variable
c
            do i=istart,imid-1
                j=i+1-istart
                if(imid<=istart+2) sname(j:j)=list(i:i)
                lname(j:j)=list(i:i)
            enddo
            do i=1,num
                if(name(i)(12:13)/='  '.and.
     +               name(i)(12:13)==sname) go to 9
                if(name(i)(5:10)==lname) go to 9
            enddo
            return
c
c       compute the value
c
    9       iptr(ii)=i
            ll=iend-imid
            if(name(i)(15:15)=='i') then
                call cint(list(imid+1:imid+1),ll,ival(ii),jerr)
                if(jerr/=0) return
            else if(name(i)(15:15)=='r') then
                call creal(list(imid+1:imid+1),ll,rval(ii),jerr)
                if(jerr/=0) return
            else if(name(i)(15:15)=='s') then
                sval(ii)=' '
                do j=imid+1,iend
                    k=j-imid
                    sval(ii)(k:k)=list(j:j)
                enddo
            else if(name(i)(15:15)=='f') then
                sval(ii)=' '
                do j=imid+1,iend
                    k=j-imid
                    sval(ii)(k:k)=list(j:j)
                enddo
            else
                if(list(iend:iend)/=dbleq) return
                if(list(imid+1:imid+1)/=dbleq) return
                sval(ii)=' '
                do j=imid+2,iend-1
                    k=j-imid-1
                    sval(ii)(k:k)=list(j:j)
                enddo
            endif
        enddo
c
c       update ip and rp arrays
c
        do ii=1,nequal
            i=iptr(ii)
            call cint(name(i),3_iknd,indx,jerr)
            if(name(i)(15:15)=='i') then
                ip(indx)=ival(ii)
            else if(name(i)(15:15)=='r') then
                rp(indx)=rval(ii)
            else
                sp(indx)=sval(ii)
            endif
        enddo
        ierr=0
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine disply(name,num,ip,rp,sp)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: ip
            real(kind=rknd), dimension(*) :: rp
            character(len=15), dimension(*) :: name
            character(len=80), dimension(*) :: sp
            character(len=80), dimension(100) :: sval
            character(len=80) :: stemp
            character(len=100) :: msg
cy
c       print reset paremeters
c
        nn=0
        do i=1,num
            call cint(name(i),3_iknd,indx,jerr)
            if(name(i)(15:15)=='i') then
                sval(i)=' '
                call sint(sval(i),ll,ip(indx))
                nn=nn+1
            else if(name(i)(15:15)=='r') then
                sval(i)=' '
                nn=nn+1
                call sreal(sval(i),ll,rp(indx),3_iknd,0_iknd)
            else
                go to 10
            endif
        enddo
c
   10   do i=1,nn,4
            write(unit=msg,fmt='(4(a9,1x,a10))')
     +          (name(j)(5:13),sval(j)(1:10),
     1           j=i,min(i+3,nn))
            call crtutl(msg,'w','prompt: ')
        enddo
c
        do i=nn+1,num
            call cint(name(i),3_iknd,indx,jerr)
            if(name(i)(15:15)=='s') then
                call fstr(stemp,length,sp(indx),0_iknd)
            else if(name(i)(15:15)=='f') then
                call fstr(stemp,length,sp(indx),0_iknd)
            else
                call fstr(stemp,length,sp(indx),1_iknd)
            endif
            write(unit=msg,fmt='(a9,1x,80a1)')
     +          name(i)(5:13),(stemp(k:k),k=1,length)
            call crtutl(msg,'w','prompt: ')
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
        subroutine discmd(ncmd,ctable)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=15), dimension(*) :: ctable
            character(len=100) :: msg
cy
c       print command list
c
        if(ncmd<=6) then
            nstep=ncmd
        else
            nstep=min((ncmd+1)/2,6)
        endif
        do k=1,ncmd,nstep
            write(unit=msg,fmt='(6(a8,4x))')
     +          (ctable(j)(1:8),j=k,min(k+nstep-1,ncmd))
            call crtutl(msg,'w','prompt: ')
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
        subroutine fxcase(list,length,ncomma,nequal,ndbleq,
     +      lcomma,lequal,icomnt)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(*) :: lequal,lcomma
            character(len=*) :: list
            character(len=1), save :: blank=' ',comma=',',
     +          equal='=',dbleq='"',commnt='#'
            character(len=1) :: cc,uppera,upperz
cy
        llist=80
        length=0
        ncomma=0
        nequal=0
        ndbleq=0
        if(list(1:1)==commnt) then
            icomnt=1
            do i=llist,1,-1
                if(list(i:i)/=blank) then
                    length=i
                    return
                endif
            enddo
            return
        else
            icomnt=0
        endif
        uppera=char(65)
        upperz=char(90)
c
c       delete blanks, find equal, commas, and double quotes
c       convert upper to lower case except for command code
c
        do i=1,llist
            cc=list(i:i)
            list(i:i)=blank
            if(ndbleq-(ndbleq/2)*2==0) then
                if(cc/=blank) then
                    length=length+1
                    if(cc==comma) then
                        ncomma=ncomma+1
                        lcomma(ncomma)=length
                    else if(cc==equal) then
                        nequal=nequal+1
                        lequal(nequal)=length
                    else if(cc==dbleq) then
                        ndbleq=ndbleq+1
                    else if(cc>=uppera.and.cc<=upperz) then
                        ii=ichar(cc)+32
                        if(length>1.and.nequal==ncomma)
     +                      cc=char(ii)
                    endif
                    list(length:length)=cc
                endif
            else
                length=length+1
                list(length:length)=cc
                if(cc==dbleq) ndbleq=ndbleq+1
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
        subroutine sreal(list,length,val,ndig,just)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd) :: elen,mlen,just
            character(len=*) :: list
            character(len=1), save :: zero='0',minus='-',
     +          e='e',dot='.'
            character(len=100) :: ex,mant
cy
c       compute character string for floating point number
c
        if(val==0.0e0_rknd) then
            length=3
            list(1:1)=zero
            list(2:2)=dot
            list(3:3)=zero
        else
            zc=abs(val)
            zz=log10(zc)
            iex=int(zz)
            ratio=10.0e0_rknd**(zz-real(iex,rknd))
c**         ratio=zc*(10.0e0_rknd**(-iex))
            if(iex==-1) then
                h=0.5e0_rknd*10.0e0_rknd**(2-ndig)
            else
                h=0.5e0_rknd*10.0e0_rknd**(1-ndig)
            endif
            if(ratio+h<1.0e0_rknd) then
                ratio=ratio*10.0e0_rknd
                iex=iex-1
            else if(ratio+h>=10.0e0_rknd) then
                ratio=ratio/10.0e0_rknd
                iex=iex+1
            endif
c
c       exponent field
c
            call sint(ex,elen,iex)
c
c       mantissa field
c
            if(iex==-1) then
                n=int(ratio*10.0e0_rknd**(ndig-2)+0.5e0_rknd)
            else
                n=int(ratio*10.0e0_rknd**(ndig-1)+0.5e0_rknd)
            endif
c
            if(just/=1) then
   90           k=n/10
                j=n-10*k
                if(j==0) then
                    n=k
                    go to 90
                endif
            endif
            call sint(mant,mlen,n)
            if(val>0) then
                is=0
            else
                is=1
                list(1:1)=minus
            endif
            if(iex==-1) then
                list(is+1:is+1)=zero
                list(is+2:is+2)=dot
                do i=1,mlen
                    list(is+i+2:is+i+2)=mant(i:i)
                enddo
                mlen=mlen+1
                iex=0
            else if(iex==1) then
                list(is+1:is+1)=mant(1:1)
                list(is+2:is+2)=zero
                list(is+3:is+3)=dot
                list(is+4:is+4)=zero
                if(mlen<=2) then
                    if(mlen==2) list(is+2:is+2)=mant(2:2)
                    mlen=3
                else
                    list(is+2:is+2)=mant(2:2)
                    do i=3,mlen
                        list(is+i+1:is+i+1)=mant(i:i)
                    enddo
                endif
                iex=0
            else
                list(is+1:is+1)=mant(1:1)
                list(is+2:is+2)=dot
                if(mlen==1) then
                    list(is+3:is+3)=zero
                    mlen=mlen+1
                else
                    do i=2,mlen
                        list(is+i+1:is+i+1)=mant(i:i)
                    enddo
                endif
            endif
            if(iex/=0) then
                length=elen+mlen+2+is
                list(is+mlen+2:is+mlen+2)=e
                do i=1,elen
                    list(is+mlen+2+i:is+mlen+2+i)=ex(i:i)
                enddo
            else
                length=mlen+1+is
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
        subroutine sint(list,length,ival)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(100) :: temp
            character(len=*) :: list
            character(len=10), save :: num='0123456789'
            character(len=1), save :: minus='-'
cy
c       compute character string for integer
c
        if(ival==0) then
            length=1
            list(1:1)=num(1:1)
        else
            length=0
            n=abs(ival)
   10       j=n/10
            i=n-j*10
            length=length+1
            temp(length)=i+1
            n=j
            if(n>0) go to 10
            if(ival<0) then
                list(1:1)=minus
                do i=1,length
                    j=temp(length+1-i)
                    list(i+1:i+1)=num(j:j)
                enddo
                length=length+1
            else
                do i=1,length
                    j=temp(length+1-i)
                    list(i:i)=num(j:j)
                enddo
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
        subroutine fstr(ss,length,sval,iquote)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), save :: mxchar=80
            character(len=1), save :: blank=' ',dbleq='"'
            character(len=80) :: ss,sval
cy
        istart=mxchar+1
        istop=0
        ss=' '
        do i=1,mxchar
            if(sval(i:i)==blank) cycle
            istart=min(istart,i)
            istop=max(istop,i)
        enddo
        if(iquote==1) then
            ss(1:1)=dbleq
            if(istart>istop) then
                length=3
            else
                length=istop-istart+3
                if(length>mxchar) then
                    istop=istop-(length-mxchar)
                    length=mxchar
                endif
                ss(2:length-1)=sval(istart:istop)
            endif
            ss(length:length)=dbleq
        else
            if(istart>istop) then
                length=1
            else
                length=istop-istart+1
                ss(1:length)=sval(istart:istop)
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
        subroutine mkname(outnam,innam)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), save :: num=0
            character(len=80) :: innam,outnam,temp
            common /atest6/nproc,myid,mpisw,mpiint,mpiflt
cy
c       look for key string and insert number
c
        num=num+1
cccc    if(mpisw==1) call exnum(num)
        call fstr(outnam,length,innam,0_iknd)
        do i=6,length
            if(outnam(i-5:i)=='figxxx') then
                outnam(i-2:i)='000'
                call sint(temp,len,num)
                outnam(i+1-len:i)=temp(1:len)
                return
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
        subroutine creal(list,length,val,ierr)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd) :: zero
            character(len=*) :: list
            character(len=1), save :: dot='.',lce='e',blank=' ',
     +          plus='+',minus='-',lcd='d'
            character(len=1) :: uce,cc,ucd
            character(len=1), dimension(80) :: temp
cy
c       compute a real number from a-format input
c
        ii=ichar(lce)-32
        uce=char(ii)
        ii=ichar(lcd)-32
        ucd=char(ii)
        val=0.0e0_rknd
        ierr=1
        newlen=0
        idot=length+1
        iee=length+1
        do i=1,length
            cc=list(i:i)
            list(i:i)=blank
            if(cc==blank) cycle
            newlen=newlen+1
            temp(newlen)=cc
            list(newlen:newlen)=cc
            if(temp(newlen)==lce) iee=newlen
            if(temp(newlen)==uce) iee=newlen
            if(temp(newlen)==lcd) iee=newlen
            if(temp(newlen)==ucd) iee=newlen
            if(temp(newlen)==dot) idot=newlen
        enddo
        if(newlen==0) return
c
c       exponent
c
        if(iee<=newlen) then
            if(iee==1.or.iee==newlen) return
            ll=newlen-iee
            call cint(temp(iee+1),ll,ix,jerr)
            if(jerr/=0) return
            newlen=iee-1
        else
            ix=0
        endif
c
c       mantissa
c
        if(idot<=newlen) then
            if(newlen==1) return
            ix=ix+idot-newlen
            newlen=newlen-1
            if(idot<=newlen) then
                do i=idot,newlen
                    temp(i)=temp(i+1)
                enddo
            endif
        endif
c
c       sign
c
        if(temp(1)==minus.or.temp(1)==plus) then
            if(newlen==1) return
            ii=2
        else
            ii=1
        endif
c
        zero=ichar('0')
        value=0.0e0_rknd
        do i=ii,newlen
            kx=ichar(temp(i))-zero
            if(kx<0.or.kx>9) return
            value=10.0e0_rknd*value+real(kx,rknd)
        enddo
        if(temp(1)==minus) then
            val=-value*(10.0e0_rknd**ix)
        else
            val=value*(10.0e0_rknd**ix)
        endif
        ierr=0
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine cint(list,length,ival,ierr)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd) :: zero
            character(len=*) :: list
            character(len=1), save :: blank=' ',plus='+',minus='-'
            character(len=1), dimension(80) :: temp
            character(len=1) :: cc
cy
c       compute an integer from a-format input
c
        ierr=1
        ival=0
        newlen=0
        do i=1,length
            cc=list(i:i)
            list(i:i)=blank
            if(cc==blank) cycle
            newlen=newlen+1
            temp(newlen)=cc
            list(newlen:newlen)=cc
        enddo
        if(newlen==0) return
c
c       sign
c
        if(temp(1)==minus.or.temp(1)==plus) then
            if(newlen==1) return
            ii=2
        else
            ii=1
        endif
c
c
        zero=ichar('0')
        do i=ii,newlen
            ix=ichar(temp(i))-zero
            if(ix<0.or.ix>9) return
            ival=10*ival+ix
        enddo
        if(temp(1)==minus) ival=-ival
        ierr=0
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine cpause()
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=80) :: cc
            common /atest3/mode,jnlsw,jnlr,jnlw,ibatch
cy
c       wait for user to view picture
c
        if(mode==0.and.jnlsw==1) then
            call xpause()
        else if(mode==-1.and.jnlsw==1) then
            call crtutl(cc,'r','pause:  ')
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
        subroutine crtutl(list,mode,prompt)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), save :: icrtr=5,icrtw=6
            character(len=1) :: mode
            character(len=8) :: prompt
            character(len=80) :: list
cy
c       print a prompt symbol
c
        if(mode=='r') then
            write(icrtw,fmt='(/ a8 $)') prompt
            flush(icrtw)
            read(icrtr,fmt='(a80)') list
            flush(icrtr)
        else if(mode=='w') then
            write(icrtw,fmt='(a80)') list
            flush(icrtw)
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
        subroutine filutl(list,isw)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=1), save :: lowerj='j'
            character(len=80) :: list
            character(len=80), save :: blank=' '
            character(len=100) :: msg
            common /atest3/mode,jnlsw,jnlr,jnlw,ibatch
cy
        if(isw==1) then
            if(list(1:1)==lowerj) then
                write(unit=msg,fmt='(a1,a80)') '#',list
            else
                write(unit=msg,fmt='(a80)') list
            endif
            len=1
            do i=2,80
                if(msg(i:i)/=' ') len=i
            enddo
            call ascstr(jnlw,msg,len,iflag)
            write(unit=msg,fmt='(a8,a80)') 'command:',list
c
            call ascstr(ibatch,blank,1_iknd,iflag)
            call ascstr(ibatch,msg,80_iknd,iflag)
c
            if(mode==0) then
                call xtext(blank)
                call xtext(msg)
            endif
c
            if(mode==-1.and.jnlsw==1) then
                call crtutl(blank,'w','prompt: ')
                call crtutl(msg,'w','prompt: ')
            endif
        else
            call ascstr(ibatch,list,80_iknd,iflag)
            if(mode==0) call xtext(list)
            if(mode==-1) call crtutl(list,'w','prompt: ')
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
        subroutine mkjnl(sp,iflag)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd), dimension(11) :: jnlr
            integer(kind=iknd), dimension(24) :: lequal,lcomma
            integer(kind=iknd), save :: maxd=10
            character(len=1) :: lowerj,upperj
            character(len=80), dimension(100) :: sp
            character(len=80), dimension(11) :: name
            character(len=80) :: list,filnam
cy
c       make journal file
c
        lowerj=char(106)
        upperj=char(74)
        do i=1,maxd
            name(i)=' '
            jnlr(i)=-1
        enddo
        if(sp(10)==' ') then
            sp(11)='journl: bad filename'
            go to 50
        endif
        call stfile(filnam,sp(10))
        call ascutl(jnlr(maxd+1),filnam,'w',kflag)
        if(kflag/=0) then
            sp(11)='journl: cannot open file'
            go to 50
        endif
        iflag=0
        sp(11)='journl: ok'
        level=1
        name(1)=sp(7)
c
c       open file
c
   10   if(name(level)==' ') then
            sp(11)='journl: bad filename'
            go to 50
        endif
        if(level>=maxd) then
            sp(11)='journl: too many levels'
            go to 50
        endif
        do i=1,level-1
            if(name(level)==name(i)) then
                sp(11)='journl: bad filename'
                go to 50
            endif
        enddo
        call stfile(filnam,name(level))
        call ascutl(jnlr(level),filnam,'r',kflag)
        if(kflag/=0) then
            sp(11)='journl: cannot open file'
            go to 50
        endif
c
c       read next command
c
   20   call ascstr(jnlr(level),list,80_iknd,kflag)
        if(kflag>0) then
            sp(11)='journl: read error'
            go to 50
        endif
        if(kflag==-1) then
c
c       close current file, reduce level
c
            call ascutl(jnlr(level),filnam,'c',jflag)
            if(jflag/=0) then
                sp(11)='journl: cannot close file'
                return
            endif
            jnlr(level)=-1
            level=level-1
            if(level>=1) go to 20
            call ascutl(jnlr(maxd+1),filnam,'c',jflag)
            return
        endif
c
c       process this command
c
        call fxcase(list,length,ncomma,nequal,ndbleq,
     +      lcomma,lequal,icomnt)
        if(length<=0) then
            go to 20
c
c       check for journal commands
c
        else if(list(1:1)==lowerj.or.list(1:1)==upperj) then
            if(ncomma>0.or.ndbleq>0.or.nequal>=2) then
                sp(11)='journl: command error'
                go to 50
            endif
            if(nequal==1) then
                ll=length-lequal(1)
                name(level+1)=' '
                name(level+1)(1:ll)=list(lequal(1)+1:length)
            endif
            if(list(1:1)==lowerj) then
                level=level+1
                go to 10
            else
                go to 20
            endif
        else
c
c       print this command
c
            call ascstr(jnlr(maxd+1),list,length,kflag)
            if(kflag>0) then
                sp(11)='journl: write error'
                go to 50
            endif
            go to 20
        endif
c
c       close all open files
c
   50   do i=1,maxd+1
            if(jnlr(i)==-1) cycle
            call ascutl(jnlr(i),filnam,'c',kflag)
        enddo
        iflag=-7
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine stfile(outnam,innam)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=80) :: innam,outnam,temp
            common /atest6/nproc,myid,mpisw,mpiint,mpiflt
cy
c       look for key strng and replace with proc number
c
        call fstr(outnam,length,innam,0_iknd)
        do i=6,length
            if(outnam(i-5:i)=='mpixxx') then
                outnam(i-2:i)='000'
                call sint(temp,len,myid+1)
                outnam(i+1-len:i)=temp(1:len)
                return
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
        subroutine mpiutl(isw)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            common /atest6/nproc,myid,mpisw,mpiint,mpiflt
cy
c       initialize
c
        if(isw==1) then
            nproc=1
            myid=0
            mpisw=-1
            mpiint=0
            mpiflt=0
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
        subroutine star0(cmd)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            character(len=80) :: cmd
cy
        return
        end
