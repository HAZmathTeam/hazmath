C=====================================================================
      subroutine init_guess(ia,ja,a,b,sol,n,nx,ny,h)
C=====================================================================
C
      integer ia(1),ja(1),k,n,p1,p2, iseed(4),nx,ny
      real*8 a(1),sol(1),b(1),
     >     solmax,solmin,xxx,yyy,x0,y0,pi,h
      real time(2),t1,t2,t3
C--------------------------------------------------------------------
C... Initial guess:
C
      pi = 3.141592653589793d0
      t1 = etime(time)
      call nullv(sol,n)
C
      do k = 1 , 4
         iseed(k) = 1
      end do
C... Get a uniformly distributed random numbers as initial guess
cccccccccccc      call dlaruv( iseed, n, u )
cccccccccccc      call dlarnv( 2, iseed, n, sol )
      write(*,*) 'NX,NY', nx,ny,h
      do i = 1 , nx
         x0 = (i-1)*h
         do j = 1 , ny
            y0 = (j-1)*h
            node = nomxy(i,j,nx,ny)
cccccccccccccccc            write(*,*) i,j,node
            yyy = j*pi*0.05
            xxx = i*pi*0.05
            sol(node)=dsin(xxx)*dsin(yyy)+sol(node)
         end do
      end do
C
      do k = 1 , n
         p1 = ia(k)
         p2 = ia(k+1)-1
         if(p2-p1 .lt. 0) stop 10
         if(p2-p1 .eq. 0) then
            sol(k) = b(k)
         end if
      end do
 222  continue
      return
      end

CLALALLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      subroutine blockit(levelf,iblock,jblock,
     >     mgend,ia,ja,lsub,nx,ny,n)
      implicit real*8 (a-h,o-z)
      dimension iblock(1),jblock(1),mask(1),ia(1),ja(1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      kh = 2
      istep = 2**(levelf-mgend)
      iistep = 2**(levelf-mgend)
      i0 = istep+1
      j0 = istep+1
      jn = ny-istep
      in = nx-istep
      nccx = 2**mgend+1
      nccy = nccx
      lvl = levelf - mgend + 1
      lsub = 0
      iblock(1) = 1
      ip=1
      do ii = 1, nccx
         i = (ii-1)*iistep+1
         nxbegin = max0(i - istep,1)
         nxend = min0(i + istep,nx)
         if(ii .ne. 1) nxbegin = nxbegin + 1
         if(ii .ne. nccx) nxend = nxend - 1
         do jj = 1, nccy
            j = (jj-1)*iistep+1
            nybegin = max0(j - istep,1)
            nyend = min0(j + istep,ny)
            if(jj .ne. 1) nybegin = nybegin + 1
            if(jj .ne. nccy) nyend = nyend - 1
c     write(*,*) ii,jj,nxbegin,nxend,nybegin,nyend
c     read(*,*)
            lsub = lsub + 1
            do i = nxbegin, nxend
               do j = nybegin, nyend
                  node = nomxy(i,j,nx,ny)
                  jblock(ip) = node
cccccccccwrite(*,*) ip,node
                  ip = ip + 1
               end do
            end do
            iblock(lsub+1) = ip
         end do
      end do
      return
      end
C________________________________________________________________
      subroutine subd(n,mask,iblock,jblock,lsub,
     >     ia,ja,a,x,b,m,r,ijkl,mnp)
C________________________________________________________________
C
C   y=BPXmaybe* x
C________________________________________________________________
      implicit real*8 (a-h,o-z)
      dimension mask(1),iblock(1),jblock(1)
      dimension ia(1),ja(1),a(1),x(1),b(1)
      dimension m(1),r(1)
      call inullv(mask,n)
      ibeg = 1
      iend = lsub
      if (ijkl .ge. mnp) then 
         ibeg=1
         iend=ibeg
      end if 
      do k = ibeg,iend
         if(lsub .le. 81) then
            write(*,'(a$)')  '*'
         elseif(mod(k,81) .eq. 1) then
            write(*,'(a$)')  '*'
         end if
         ip =iblock(k)
         ip1 = iblock(k+1)
         nsub = ip1-ip
         kxsub = 1
         ifree = 1
         kfree = kxsub + nsub
C         write(*,*) ' kfree is: ', kfree
         do i = 1 , nsub
            ibig = jblock(ip+i-1)
            uu = b(ibig)
            do j=ia(ibig),ia(ibig+1)-1
               uu = uu - a(j)*x(ja(j))
            end do
            ismall = kxsub + i-1
            r(ismall)=uu
         end do
C
         call upper2(n,jblock(ip),mask,
     >        ia,ja,a,m(ifree),r(kfree),r(kxsub),
     >        nsub,nwk)
         write(*,*) k, ' has', nsub, ' unknowns and ', 
     >        nwk, ' nonzeroes in the envelope'
         do i = 1 , nsub
            i1 =  i - 1
            ibig=jblock(ip+i1)
            ismall = kxsub + i1
            x(ibig) = x(ibig)+r(ismall)
         end do
      end do
      return
      end
      subroutine upper2(nbig,jblock,mask,
     >     xadj,adjncy,a,maxa,au,b,n,nwk)
C=====================================================================
      integer*4 adjncy(1),xadj(1),maxa(1),jblock(1),mask(1)
      real*8 a(1),b(1),au(1)
C---------------------------------------------------------------------
C...  This subroutine transforms the adjacency graph structure into 
C...  envelope one. The goal is to invert them into AL and AU by gauss 
C...  elimination. 
C---------------------------------------------------------------------
c      write(*,*) 'Jblock: ',(jblock(kk),kk=1,n)
c      call lprr(xadj,adjncy,a,nbig)
      do k = 1 , n
         mask(jblock(k)) = k
      end do
cc      write(*,*) ' MASK: ', (mask(kk),kk=1,nbig)
cc      read(*,*)
      maxa(1) = 1
      do  k = 2, n+1
         nodei = k-1
         icol  = 0
         nodeb = jblock(nodei)
         do jk = xadj(nodeb), xadj(nodeb+1) - 1

            nodej=mask(adjncy(jk))

            if(nodej .ne. 0) then
               icol = max0(icol,k-nodej)
            end if
         end do
         maxa(k) = maxa(k-1)+icol
ccccccc         write(*,*) maxa(k)
      end do
C
      nnm = n+1
      nwk = maxa(nnm) - 1
C
      call nullv(au,nwk)
C
cc      write(*,*) ' WHATEVER has', n, ' unknowns and ', 
cc     >        nwk, ' nonzeroes in the envelope'
      kkz = 0
      do k = 1, n
         nodei = k
         nodeb = jblock(nodei)
         do jk = xadj(nodeb),xadj(nodeb+1)-1
            nodej = mask(adjncy(jk))
            idist = iabs(nodej - k)
            if(nodej-k) 50, 30, 40
 30         continue
            au(maxa(k)) = a(jk)
            go to 50
 40         continue
            au(maxa(nodej)+idist) = a(jk)
 50         continue
         end do
      end do
C
      nll = nwk
      call nullv(au(nwk+1),nwk)
C
      do k = 1 , n
         nodei = k
         nodeb = jblock(nodei)
         nxx = maxa(nodei) + nll
         do jk = xadj(nodeb),xadj(nodeb+1)-1
            nodej = mask(adjncy(jk))
            idist = iabs(nodej - k)
            if(nodej-k) 110, 120, 140
 110        continue
            if(nodej .ne. 0) then
               au(nxx+idist) = a(jk)
            end if
            go to 140
 120        continue
            au(nxx) = 1.0d00
 140        continue
         end do
 150     continue
      end do
      do k = 1 , n
         mask(jblock(k)) = 0
      end do
      call decons(au,au(nwk+1),maxa,n)
      call redbns(au,au(nwk+1),b,maxa,n)
      return
      end
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLC
