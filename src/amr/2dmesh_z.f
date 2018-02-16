      subroutine squtri(nop)
      dimension nop(3,2,2)
C
C... This forms the correspondance between the square numbering (4)
C...  nodes and the local triangle numbering (2 triangles in the
C...  square).  THUS: nop(i,j,1) is the GLOBAL number in the square of
C...  of the i-th vertex in the j-th triangle . last index is if we have
C...  diagonal SW-NE or NW-SE
C
!     first triangle:
      it=1 !SW-NE
      nop(1,1,it) = 1
      nop(2,1,it) = 3
      nop(3,1,it) = 4
!     second triangle:
      nop(1,2,it) = 4
      nop(2,2,it) = 2
      nop(3,2,it) = 1
!
      it=2                      !SW-NE
!     first triangle:
      nop(1,1,it) = 3
      nop(2,1,it) = 1
      nop(3,1,it) = 2
!     second triangle:
      nop(1,2,it) = 2
      nop(2,2,it) = 4
      nop(3,2,it) = 3
      return
      end
      subroutine xyloc(xy,jsqu,nx,ny,
     >     i,j,hx,hy,xmin,ymin)
      implicit real*8(a-h,o-z)
      dimension xy(4,2),jsqu(*)
c
      im1=i-1
      jm1=j-1
c
      xy(1,1) = (im1)*hx + xmin
      xy(1,2) = (jm1)*hy + ymin
c
      xy(3,1) = (i)*hx   + xmin
      xy(3,2) = (j-1)*hy + ymin
c
      xy(4,1) = (i)*hx   + xmin
      xy(4,2) = (j)*hy   + ymin
c
      xy(2,1) = (im1)*hx + xmin
      xy(2,2) = (j)*hy   + ymin
c
      ip1=i+1
      jp1=j+1
C (i,j)
      jsqu(1) = nomxy(i,j,nx,ny)
!      jsqu(1) = (j-1)*nx + i
C (i,j+1)
      jsqu(2) = nomxy(i,jp1,nx,ny)
!      jsqu(2) =  (j)*nx   + i 
C (i+1,j)
      jsqu(3) = nomxy(ip1,j,nx,ny)
!      jsqu(3) = (j-1)*nx + i+1 
C (i+1,j+1)
      jsqu(4) = nomxy(ip1,jp1,nx,ny)
!      jsqu(4) = (j)*nx + i + 1 
      return
      end

      subroutine bndry2(ib,nx,ny,
     >     mini,minj,maxi,maxj,inumb,
     >     minneu,maxneu)
      implicit real*8 (a-h,o-z)
      dimension ib(*),inumb(*)
      do j = 1 , ny
         do i = 1 , nx
            iiold = nomxy(i,j,nx,ny)
            ii = inumb(iiold)
            if(ii .ne. 0) then 
               ib(ii) = 0
            end if
C     ii = (k-1)*nx*ny + (j-1)*nx + i
            if(  i.eq.1 .or. i.eq.nx .or.
     >           j.eq.1 .or. j.eq.ny) then
c     write(*,*) i,j,k
               if(ii .eq. 0) stop 4
               if(i.eq. 1) ib(ii) = minneu ! neumann conditions are > 16; 
               if(j.eq. 1) ib(ii) = minneu+1 ! neumann conditions are > 16; 
               if(i.eq. nx) ib(ii) = minneu+2 ! neumann conditions are > 16; 
               if(j.eq. ny) ib(ii) = minneu+3 ! neumann conditions are > 16; 
            else 
               if((i.eq.mini .or. i.eq.maxi) .and.
     >              ((j.ge.minj .and. j.le.maxj))) then
c     write(*,*) i,j,k
                  if(ii .eq. 0) stop 5
                  ib(ii) = -1 
               else if ((j.eq.minj .or. j.eq.maxj) .and.
     >                 ((i.ge.mini .and. i.le.maxi))) then
c     write(*,*) i,j,k
                  if(ii .eq. 0) stop 6
                  ib(ii) = -1 
               end if
            end if
         end do
      end do
c     read(*,*)
      return
      end
      subroutine bess2(ib,nx,ny,inumb,ibess,jbess,
     >     ibcode,minneu,maxneu)
      implicit real*8 (a-h,o-z)
      dimension ib(*),inumb(*)
      logical doi,doj

!     sets boundary conditions: Natural everywhere, but essential on the
      ! come bundaries. Inefficiently done, but simple enough. ibess is
      ! the parameter that will tell us. Natural boundary conditions are
      ! imposed on boundarys of the form: (ibess,j) and (i,jbess) if we
      ! want essential conditions everywhere, it can be called 2 times
      ! with ibess=jbess=1 and ibess=jbess=(nx,ny).

!     We assume the boundary conditions are already set on all boundaries. 
!     ibcode() is an array specifying what code to be put on what
!     boundary. The boundaries are ordered as
!     ibcode(1:6) = [west,east,south,north]
!
!     the minneu and maxneu are the min and max codes used on Neumann
!     boundaries. Any code in the interval [minneu,maxneu] is ignored.
      
      if(ibcode .eq. 0) return
      if(ibcode .ge. minneu .and. ibcode .le. maxneu) return
      doi=(ibess .le. nx .and. ibess .ge. 1)
      doj=(jbess .le. ny .and. jbess .ge. 1)
      if(doi) then
         i = ibess
! only interior to the manifold "i=ibess"
         do j = 1 , ny
            iiold = nomxy(i,j,nx,ny)
            ii = inumb(iiold)
C     ii = (k-1)*nx*ny + (j-1)*nx + i
            ib(ii) = ibcode         ! <=16 essential condition;
         end do
      end if
      if(doj) then
         j=jbess
         do i = 1 , nx
            iiold = nomxy(i,j,nx,ny)
            ii = inumb(iiold)
            ib(ii) = ibcode         ! <= 16 essential condition;
         end do
      end if
      return
      end
c=====================================================================
      subroutine getm2(nx,ny,nvert,nel,
     >     xcoord,ycoord,
     >     je,iflags,ib,inumb,
     >     ibcode,minneu,maxneu)
      implicit real*8(a-h,o-z), integer(i-n)
      dimension xcoord(*),ycoord(*)
      dimension je(*),ib(*), inumb(*),iflags(*)
      dimension xy(4,2),jsqu(4),jcolo(3)
      dimension nop(3,2,2)
      dimension ibcode(*)
C
C... ndl = 3 ! number of degrees of freedom per element
C
      do k = 1 , nvert
         inumb(k) = 0
         ib(k) = 0
      end do
      ntold = nvert
      nvert = 0
      xmin = 0.0d00
      ymin = 0.0d00
!     
      hx = 1d0/dfloat(nx-1)
      hy = 1d0/dfloat(ny-1)
      ishift=0
      iiold = 0
      mini=nx+1
      maxi=0
      minj=ny+1
      maxj=0
      do j = 1 , ny
         do i = 1 , nx
            iiold = iiold + 1
            iichk = nomxy(i,j,nx,ny)
            if(iichk .ne. iiold) then
               write(*,*) ' ERROR: The function nomxy() '
               write(*,*) '  is not in accordance '
               write(*,*) ' with the lexicographical ordering: '
               write(*,*) ' first x; second y '
               stop 2
            end if
            if( (i .gt. mini .and.
     >           j .gt. minj) .and. 
     >           (i .lt. maxi .and.
     >           j .lt. maxj)) then
C...  this is a node outside of our domain, and so we increase the shift
               ishift = ishift + 1
C     write(*,*) 'ijk is in the hole', i,j,k
            else
               nvert = iiold - ishift
               inumb(iiold) = nvert
            end if
         end do
      end do
c      read(*,*) iii
      call squtri(nop) ! this is done only once. 
      jjkk = 1
      kl = 1
      isift=0 
      nxm1=nx-1
      nym1=ny-1
C
      do j = 1 , nym1 
         do i = 1 , nxm1
!     we choose how to orient the quadrilateral now:
            it=mod(iabs(i-j),2)+1
!!            write(*,*) 'it=', it
            if(  i .lt. mini .or.
     >           j .lt. minj .or.
     >           i .ge. maxi .or.
     >           j .ge. maxj) then
               call xyloc(xy,jsqu,nx,ny,
     >              i,j,hx,hy,xmin,ymin)
               do jk = 1 , 2
                  do mm = 1 , 3
                     inode = inumb(jsqu(nop(mm,jk,it)))
                     if(inode .eq. 0) then
                        write(*,*) 'ERROR:i,j and inode = 0', 
     >                       i,j,inode
                        stop 3
                     end if
                     xcoord(inode) = xy(nop(mm,jk,it),1)
                     ycoord(inode) = xy(nop(mm,jk,it),2)
                     je(kl) = inode
                     kl = kl + 1
                  end do
                  iflags(jjkk) = j ! this is the element flag changing
                                   ! in vertical direction
                  jjkk = jjkk + 1
               end do
            end if
         end do
      end do
      nel = jjkk  - 1
      call bndry2(ib,nx,ny,
     >     mini,minj,maxi,maxj,inumb,
     >     minneu,maxneu)
!     west(left)
      ibess=1
      jbess=0
      call bess2(ib,nx,ny,inumb,ibess,jbess,
     >     ibcode(1),minneu,maxneu)
!     east(right)
      ibess=0
      jbess=nx
      call bess2(ib,nx,ny,inumb,ibess,jbess,
     >     ibcode(2),minneu,maxneu)
!     south(front)
      ibess=0
      jbess=1
      call bess2(ib,nx,ny,inumb,ibess,jbess,
     >     ibcode(3),minneu,maxneu)
!     north(back)
      ibess=0
      jbess=ny
      call bess2(ib,nx,ny,inumb,ibess,jbess,
     >     ibcode(4),minneu,maxneu)
!      write(*,'(a,i12,a,i12)') ' Elements=', nel, '  Nodes=', nvert
!      write(*,'(a,2i7)') ' Essential conditions on (iess, jess): ',
!     >     ibess,jbess
      return
      end
      integer function nomxy(i,j,nx,ny)
      implicit real*8(a-h,o-z),integer (i-n)
      nomxy=0
      if(  i.gt.0 .and. i.le.nx .and.
     >     j.gt.0 .and. j.le.ny) then
         nomxy = (j-1)*nx + i
      end if
      return
      end
