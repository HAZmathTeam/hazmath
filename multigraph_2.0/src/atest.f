c*************************** file: atest.f *****************************
c-----------------------------------------------------------------------
c
c             piecewise linear triangle multi grid package
c
c                   edition 8.5 - - - december, 2000
c
c-----------------------------------------------------------------------
        program atest
c
c       storage allocation
c
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
c
c       main array sizes
c
            integer(kind=iknd), parameter :: maxn=1000000
            integer(kind=iknd), parameter :: maxja=32*maxn
            integer(kind=iknd), parameter :: maxa=2*maxja
c
c
c       main  array declarations
c
            integer(kind=iknd), dimension(maxja) :: ja,ja0
            integer(kind=iknd), dimension(10,100) :: ka
            integer(kind=iknd), dimension(100) :: ib
c
            real(kind=rknd), dimension(maxa) :: a,a0
            real(kind=rknd), dimension(maxn) :: br,b0,dr
            real(kind=rknd), dimension(10) :: time
            real(kind=rknd), dimension(22) :: hist
            real(kind=rknd), dimension(2) :: temp
            real(kind=rknd) :: etime
c
            character(len=80) :: sp,nam0
c
            common /atest1/ip(100),rp(100),sp(100)
            common /atest5/idevce
            common /atest6/nproc,myid,mpisw,mpiint,mpiflt
            data mxhist/20/
c
c       mode  =  1    run in batch mode
c             =  0    use x-windows interface
c             = -1    use terminal window interface
c
        mode=0
c
c       initialize the ip, rp, sp arrays
c
        call menu(ip,rp,sp)
c
c       storage parameters
c
        ip(21)=maxn
        ip(22)=maxja
        ip(23)=maxa
c
c       parameters for atest
c
        ip(41)=1
        ip(42)=mode
        ip(43)=0
        ip(48)=mpisw
        ip(49)=nproc
        ip(50)=myid+1
        sp(21)='localhost'
c
c       get command
c
   50   call menu(ip,rp,sp)
        write(unit=sp(11),fmt='(a6,a4)') sp(12)(1:6),': ok'
c
c       factor
c
        if(sp(12)(1:6).eq.'factor') then
            n=ip(1)
            nblock=ip(3)
            ispd=ip(8)
            method=ip(9)
            ncfact=ip(11)
            maxlvl=ip(12)
            maxfil=ip(13)
            dtol=rp(3)
c
c       move original matrix from ja0/a0 to ja/a
c
            call setmtx(ip,ja0,a0,b0,ja,a,br,ib,ka)
c
c       setup phase
c
            call cpu_time(tx)
            call mginit(n,ispd,nblock,ib,maxja,ja,maxa,a,ncfact,
     +          maxlvl,maxfil,ka,lvl,dtol,method,iflag)
            call cpu_time(ty)
            time(1)=ty-tx
c
c       uncomment for debugging mgilu
c
c*          call mgilu(ja,a,lvl,ka)
c
            rp(26)=time(1)
            ip(2)=lvl
            ip(26)=iflag
            call shrtnm(nam0,sp(6))
            if(iflag.eq.0) then
                write(unit=sp(11),
     +              fmt='(a15,a5,e10.3,a9,i3,a6,i2,a9,i6,a9,i2)')
     1              nam0(1:15),'dtol:',dtol,'  maxlvl:',maxlvl,
     2              '  lvl:',lvl,'  maxfil:',maxfil,
     3              '  ncfact:',ncfact
            else
                write(unit=sp(11),fmt='(a12,i3)') 'mginit error',iflag
            endif
            write(unit=6,fmt='(a76)') sp(11)
c
c       solve
c
        else if(sp(12)(1:6).eq.'solve ') then
            n=ip(1)
            lvl=ip(2)
            ispd=ip(8)
            mxcg=ip(10)
            eps=rp(1)
c
c       move right hand side from b0 to br
c
            do i=1,n
                br(i)=b0(i)
            enddo
c
c       solve phase
c
            call cpu_time(tx)
            call mg(ispd,lvl,mxcg,eps,ja,a,dr,br,ka,relerr,
     +          iflag,hist)
            call cpu_time(ty)
            time(2)=ty-tx
c
            rp(27)=time(2)
            ip(26)=iflag
            rp(2)=relerr
            ii=int(hist(mxhist+1))
            digits=-log10(relerr)
c
            call shrtnm(nam0,sp(6))
            if(iflag.eq.0) then
                write(unit=sp(11),
     +              fmt='(a15,a7,f5.2,a8,i3,a8,e10.3,a9,e10.3)')
     1              nam0(1:15),'digits:',digits,'   iter:',ii,
     2              '   init:',time(1),'   solve:',time(2)
            else
                write(unit=sp(11),fmt='(a8,i3)') 'mg error',iflag
            endif
            write(unit=6,fmt='(a75)') sp(11)
c
c       graph output data
c
        else if(sp(12)(1:6).eq.'gphplt') then
            idevce=ip(45)
            call gphplt(ip,rp,sp,hist,ka,time)
c
c       plot matrix
c
        else if(sp(12)(1:6).eq.'mtxplt') then
            idevce=ip(47)
            call mtxplt(ip,rp,sp,ja,a,ka)
c
c       generate matrix
c
        else if(sp(12)(1:6).eq.'linsys') then
            ngrid=ip(74)
            mtxtyp=ip(75)
            call linsys(ngrid,mtxtyp,sp(6),n,ispd,maxja,ja0,
     +          maxa,a0,maxn,b0,nblock,ib,iflag)
            call shrtnm(nam0,sp(6))
c
            ip(26)=iflag
            if(iflag.eq.0) then
                ip(1)=n
                ip(3)=nblock
                ip(8)=ispd
                ip(71)=ja0(n+1)-1
                ip(72)=ip(71)
                ip(73)=(ja0(n+1)-ja0(1))/n+1
                if(ispd.ne.1) ip(72)=2*ip(72)-(n+1)
                call setmtx(ip,ja0,a0,b0,ja,a,br,ib,ka)
                write(unit=sp(11),fmt='(a15,a3,i6,a10,i2)')
     +              nam0(1:15),'n =',n,'    ispd =',ispd
            else
                write(unit=sp(11),fmt='(a12,i3)') 'linsys error',iflag
            endif
            write(unit=6,fmt='(a36)') sp(11)
c
c       read file
c
        else if(sp(12)(1:6).eq.'read  ') then
            call setgr(sp(6),n,ispd,maxja,ja0,maxa,a0,maxn,b0,
     +          nblock,ib,iflag)
c
            call shrtnm(nam0,sp(6))
            ip(26)=iflag
            if(iflag.eq.0) then
                ip(1)=n
                ip(3)=nblock
                ip(8)=ispd
                ip(71)=ja0(n+1)-1
                ip(72)=ip(71)
                ip(73)=(ja0(n+1)-ja0(1))/n+1
                if(ispd.ne.1) ip(72)=2*ip(72)-(n+1)
                call setmtx(ip,ja0,a0,b0,ja,a,br,ib,ka)
                write(unit=sp(11),fmt='(a15,a3,i6,a10,i2)')
     +              nam0(1:15),'n =',n,'    ispd =',ispd
            else
                write(unit=sp(11),fmt='(a11,i3)') 'setgr error',iflag
            endif
            write(unit=6,fmt='(a36)') sp(11)
c
c       form least squares system
c
c*      else if(sp(12)(1:6).eq.'lstsq ') then
c*          call setls(n,ispd,maxja,ja0,ja,maxa,a0,a,br,iflag)
c*          call shrtnm(nam0,sp(6))
c*          ip(26)=iflag
c*          if(iflag.eq.0) then
c*              ispd=1
c*              ip(8)=ispd
c*              ip(71)=ja0(n+1)-1
c*              ip(72)=ip(71)
c*              ip(73)=(ja0(n+1)-ja0(1))/n+1
c*              if(ispd.ne.1) ip(72)=2*ip(72)-(n+1)
c*              call setmtx(ip,ja0,a0,b0,ja,a,br,ib,ka)
c*              write(unit=sp(11),fmt='(a15,a3,i6,a10,i2)')
c*   +              nam0(1:15),'n =',n,'    ispd =',ispd
c*          else
c*              write(unit=sp(11),fmt='(a11,i3)') 'setls error',iflag
c*          endif
c*          write(unit=6,fmt='(a36)') sp(11)
c
c       shell
c
c*      else if(sp(12)(1:6).eq.'shell ') then
c*          call cshex(sp(5))
c
c       quit
c
        else if(sp(12)(1:6).eq.'quit  ') then
            stop
        endif
        go to 50
c
        end





