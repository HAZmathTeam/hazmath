!     numbering is first wrt z second y third x, it follows the ordering
!     of the vertices of the unit cube if we consider their vertices as
!     binary numbers and order the corresponding binary numbers      
!... SHOULD BE REPLACED BY DIM INDEPENDENT. 
!=======================================================================
      subroutine cubtet(nop)
      dimension nop(4,6,4)
      dimension ip1(8),ip2(8),ip3(8),ip4(8)
      dimension ip(8,3)
!!,ip777(8,3)
!!      data ip1 /4,3,2,1,8,7,6,5/
!!      data ip2 /2,4,1,3,6,8,5,7/
!!      data ip3 /3,1,4,2,7,5,8,6/
C... should be replaced by dim independent. 
c$$$      data ip777 /
c$$$     >     3,1,4,2,7,5,8,6,
c$$$     >     2,4,1,3,6,8,5,7,
c$$$     >     4,3,2,1,8,7,6,5
c$$$     >     /
c$$$      data ip /
c$$$     >     5,1,7,3,6,2,8,4,
c$$$     >     3,7,1,5,4,8,2,6,
c$$$     >     7,5,3,1,8,6,4,2
c$$$  >     /
c$$$      data ip /
c$$$     >     5,6,7,8,1,2,3,4,
c$$$     >     3,4,1,2,7,8,5,6,
c$$$     >     2,1,4,3,6,5,8,7
c$$$     >     /
      data ip /
     >     5,6,7,8,1,2,3,4,
     >     3,4,7,8,1,2,5,6,
     >     2,4,6,8,1,3,5,7
     >     /
      
C
C... This forms the correspondance between the cube numbering (8) nodes
C...  and the local tetrahedra numbering (6 tetrahedra in the cube).
C... THUS:
C...  nop(i,j) is the GLOBAL number in the cube of of the i-th 
C...  vertex in the j-th tetrahedra. 
C
C
      it=1
      nop(1,1,it) = 1
      nop(2,1,it) = 2
      nop(3,1,it) = 4
      nop(4,1,it) = 8
c
c     second tetrahedra:::
c
      nop(1,2,it) = 1
      nop(2,2,it) = 2
      nop(3,2,it) = 6
      nop(4,2,it) = 8
c
c     third tetrahedra:::
c
      nop(1,3,it) = 1
      nop(2,3,it) = 3
      nop(3,3,it) = 4
      nop(4,3,it) = 8
c
c     fourth tetrahedra:::
c
      nop(1,4,it) = 1
      nop(2,4,it) = 3
      nop(3,4,it) = 7
      nop(4,4,it) = 8
c
c     fifth tetrahedra:::
c
      nop(1,5,it) = 1
      nop(2,5,it) = 5
      nop(3,5,it) = 6
      nop(4,5,it) = 8
c
c     sixth tetrahedra:::
      nop(1,6,it) = 1
      nop(2,6,it) = 5
      nop(3,6,it) = 7
      nop(4,6,it) = 8
C
!      write(*,'(18i3)') (ip(j,1),j=1,8)
!      write(*,'(18i3)') (ip(j,2),j=1,8)
!      write(*,'(18i3)') (ip(j,3),j=1,8)
!     read(*,*) iii
      do it=2,4                 ! cube type (2,3,4 depending on the
                                ! diagonal used; type 1 uses diagonal
                                ! 1,8; the rest are permutations of this
                                ! using ip
         do j = 1 , 6
            do l = 1 , 4
               np=nop(l,j,1)
               nop(l,j,it) = ip(np,it-1)
            end do
         end do
      end do
      return
      end
C================================================================
C      SUBROUTINE nullv(u,n)
CC================================================================
C      implicit real*8(a-h,o-z)
C      dimension u(n)
C      do i = 1 , n
C         u(i) = 0.0d00
C      end do
C      return
C      end
C================================================================

C... The next FOUR subroutines have already stored in 
C... the old version of the program, so kick 
C...  them out from  here before attaching this routine to the package.
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine xyzloc(xyz,jcub,nx,ny,nz,
     >     i,j,k,hx,hy,hz,xmin,ymin,zmin)
      implicit real*8(a-h,o-z)
      dimension xyz(8,3),jcub(*)
c
      xyz(1,1) = (i-1)*hx + xmin
      xyz(1,2) = (j-1)*hy + ymin
      xyz(1,3) = (k-1)*hz + zmin
c
      xyz(2,1) = (i-1)*hx + xmin
      xyz(2,2) = (j-1)*hy   + ymin
      xyz(2,3) = (k)*hz + zmin
c
      xyz(3,1) = (i-1)*hx   + xmin
      xyz(3,2) = (j)*hy + ymin
      xyz(3,3) = (k-1)*hz + zmin
c
      xyz(4,1) = (i-1)*hx   + xmin
      xyz(4,2) = (j)*hy   + ymin
      xyz(4,3) = (k)*hz + zmin
c
      xyz(5,1) = (i)*hx + xmin
      xyz(5,2) = (j-1)*hy + ymin
      xyz(5,3) = (k-1)*hz   + zmin
c
      xyz(6,1) = (i)*hx + xmin
      xyz(6,2) = (j-1)*hy   + ymin
      xyz(6,3) = (k)*hz   + zmin
C
      xyz(7,1) = (i)*hx   + xmin
      xyz(7,2) = (j)*hy + ymin
      xyz(7,3) = (k-1)*hz   + zmin
c
      xyz(8,1) = (i)*hx   + xmin
      xyz(8,2) = (j)*hy   + ymin
      xyz(8,3) = (k)*hz   + zmin
c
      ip1=i+1
      jp1=j+1
      kp1=k+1
      im1=i-1
      jm1=j-1
      km1=k-1
c
C (i,j,k)
      jcub(1) = nomxyz(i,j,k,nx,ny,nz)
!      jcub(1) = (k-1)*nx*ny + (j-1)*nx + i
C (i,j+1,k)
      jcub(2) = nomxyz(i,j,kp1,nx,ny,nz)
!      jcub(2) = (k-1)*nx*ny + (j)*nx   + i 
C (i+1,j,k)
      jcub(3) = nomxyz(i,jp1,k,nx,ny,nz)
!      jcub(3) = (k-1)*nx*ny + (j-1)*nx + i+1 
C (i+1,j+1,k)
      jcub(4) = nomxyz(i,jp1,kp1,nx,ny,nz)
!      jcub(4) = (k-1)*nx*ny + (j)*nx + i + 1 
C     
C (i,j,k+1)
      jcub(5) = nomxyz(ip1,j,k,nx,ny,nz)
!      jcub(5) = (k)*nx*ny + (j-1)*nx + i            
C (i,j+1,k+1)
      jcub(6) = nomxyz(ip1,j,kp1,nx,ny,nz)
!      jcub(6) = (k)*nx*ny + (j)*nx   + i            
C (i+1,j,k+1)
      jcub(7) = nomxyz(ip1,jp1,k,nx,ny,nz)
!      jcub(7) = (k)*nx*ny + (j-1)*nx + i+1 
C (i+1,j+1,k+1)
      jcub(8) = nomxyz(ip1,jp1,kp1,nx,ny,nz)
!     jcub(8) = (k)*nx*ny + (j)*nx + i + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$      write(*,*) i,j,k
c$$$      write(*,*) i,j,kp1
c$$$      write(*,*) i,jp1,k
c$$$      write(*,*) i,jp1,kp1
c$$$      write(*,*) ip1,j,k
c$$$      write(*,*) ip1,j,kp1
c$$$      write(*,*) ip1,jp1,k
c$$$      write(*,*) ip1,jp1,kp1
c$$$      read(*,*)
      return
      end
      subroutine bndry3(ib,nx,ny,nz,
     >     mini,minj,mink,maxi,maxj,maxk,inumb,
     >     minneu,maxneu)
      implicit real*8 (a-h,o-z)
      dimension ib(*),inumb(*)
      do k = 1 , nz
         do j = 1 , ny
            do i = 1 , nx
               iiold = nomxyz(i,j,k,nx,ny,nz)
               ii = inumb(iiold)
               if(ii .ne. 0) then 
                  ib(ii) = 0
               end if
C     ii = (k-1)*nx*ny + (j-1)*nx + i
               if(  i.eq.1 .or. i.eq.nx .or.
     >              j.eq.1 .or. j.eq.ny .or.
     >              k.eq.1 .or. k.eq.nz) then
c     write(*,*) i,j,k
                  if(ii .eq. 0) stop 4
                  if(i.eq. 1) ib(ii) = minneu ! neumann conditions are > 16; 
                  if(j.eq. 1) ib(ii) = minneu+1 ! neumann conditions are > 16; 
                  if(k.eq. 1) ib(ii) = minneu+2 ! neumann conditions are > 16; 
                  if(i.eq. nx) ib(ii) = minneu+3 ! neumann conditions are > 16 and < 33; 
                  if(j.eq. ny) ib(ii) = minneu+4 ! neumann conditions are > 16; 
                  if(k.eq. nz) ib(ii) = minneu+5 ! neumann conditions are > 16; 
               else 
                  if((i.eq.mini .or. i.eq.maxi) .and.
     >                 ((j.ge.minj .and. j.le.maxj) .and.
     >                 (k.ge.mink .and. k.le.maxk))) then
c                     write(*,*) i,j,k
                     if(ii .eq. 0) stop 5
                     ib(ii) = -1 
                  else if ((j.eq.minj .or. j.eq.maxj) .and.
     >                 ((i.ge.mini .and. i.le.maxi) .and.
     >                 (k.ge.mink .and. k.le.maxk))) then
c                     write(*,*) i,j,k
                     if(ii .eq. 0) stop 6
                     ib(ii) = -1 
                  else if((k.eq.mink .or. k.eq.maxk) .and.
     >                 ((j.ge.minj .and. j.le.maxj) .and.
     >                 (i.ge.mini .and. i.le.maxi))) then
c                     write(*,*) i,j,k
                     if(ii .eq. 0) stop 7
                     ib(ii) = -1 
                  end if
               end if
            end do
         end do
      end do
c     read(*,*)
      return
      end      
      subroutine bess3(ib,nx,ny,nz,inumb,ibess,jbess,kbess,
     >     ibcode,minneu,maxneu)
      implicit real*8 (a-h,o-z), integer (i-n)
      dimension ib(*),inumb(*)
      logical doi,doj,dok
!     sets boundary conditions: First, natural everywhere and then set
!     essential conditions on some specified boundaries. Inefficiently
!     done, but simple enough. ibess is the parameter that will tell
!     us. Natural boundary conditions are imposed on boundarys of the
!     form: (ibess,j) and (i,jbess)
      
!     Assume natural conditions are already set.  if we want essential
!     conditions everywhere, we need to call this twice with
!     ibess=jbess=kbess=1 and ibess=jbess=kbess=nx,ny,nz.
!     
!     essential conditions on the left: ibess=1, jbess=kbess=0.
!
!     ibcode() is an array specifying what code to be put on what
!     boundary. The boundaries are ordered as
!     ibcode(1:6) = [left,right,front,back, bottom, top]
!
!     the minneu and maxneu are the min and max codes used on Neumann
!     boundaries. any code in the interval [minneu,maxneu] is ignored.
      
      if(ibcode .eq. 0) return
      if(ibcode .ge. minneu .and. ibcode .le. maxneu) return
      doi=(ibess .le. nx .and. ibess .ge. 1)
      doj=(jbess .le. ny .and. jbess .ge. 1)
      dok=(kbess .le. nz .and. kbess .ge. 1)
      if(doi) then
         i=ibess
         do k = 1 , nz
            do j = 1 , ny
               iiold = nomxyz(i,j,k,nx,ny,nz)
               ii = inumb(iiold)
               if(ii .eq. 0) stop 4
               ib(ii) = ibcode
            end do
         end do
      end if
      if(doj) then
         j=jbess
         do k = 1 , nz
            do i = 1 , nx
               iiold = nomxyz(i,j,k,nx,ny,nz)
               ii = inumb(iiold)
               if(ii .eq. 0) stop 4
               ib(ii) = ibcode
            end do
         end do
      end if
      if(dok) then
         k=kbess
         do j = 1 , ny
            do i = 1 , nx
               iiold = nomxyz(i,j,k,nx,ny,nz)
               ii = inumb(iiold)
               if(ii .eq. 0) stop 4
               ib(ii) = ibcode
            end do
         end do
      end if
      return
      end
c=====================================================================
      subroutine getm3(nd,nvert,nel,
     >     xcoord,ycoord,zcoord,
     >     je,iflags,ib,inumb,
     >     ibcode,minneu,maxneu)
      implicit real*8(a-h,o-z)
      dimension xcoord(*),ycoord(*),zcoord(*)
      ! nd is the number of divisions in every direction
      dimension nd(*),je(*),ib(*), inumb(*),iflags(*)
      dimension xyz(8,3),jcub(8),nop(4,6,4)
      dimension ibcode(*)
C     
      nx=nd(1)
      ny=nd(2)
      nz=nd(3)
C...  Given nx,ny,nz, get the femesh: xcoord,ycoord,zcoord,ie,je.
C...  ndl = 4 ! number of degrees of freedom per element
C     
      do k = 1 , nvert
         inumb(k) = 0
         ib(k) = 0
      end do
      ntold = nvert
      nvert = 0
      xmin = 0.0d00
      ymin = 0.0d00
      zmin = 0.0d00
!
      hx = 1d0/dfloat(nx-1)
      hy = 1d0/dfloat(ny-1)
      hz = 1d0/dfloat(nz-1)
C... The loop below forms xcoord, ycoord and zcoord, and ie(),je()
C... mark all the points that are excluded and put the dof number in the 
C... array...
      ishift=0
      iiold = 0
      mini=nx+1
      maxi=0
      minj=ny+1
      maxj=0
      mink=nz+1
      maxk=0
      do k = 1 , nz
         do j = 1 , ny 
            do i = 1 , nx
               iiold = iiold + 1
               iichk = nomxyz(i,j,k,nx,ny,nz)
               if(iichk .ne. iiold) then
                  write(*,*) ' ERROR: The function nomxyz() '
                  write(*,*) '  is not in accordance '
                  write(*,*) ' with the lexicographical ordering: '
                  write(*,*) ' first x; second y; last z '
                  stop 2
               end if
               if( (i .gt. mini .and.
     >              j .gt. minj .and.
     >              k .gt. mink) .and. 
     >              (i .lt. maxi .and.
     >              j .lt. maxj .and.
     >              k .lt. maxk)) then
C... this is a node outside of our domain, and so we increase the shift
                  ishift = ishift + 1
C                  write(*,*) 'ijk is in the hole', i,j,k
               else
                  nvert = iiold - ishift
                  inumb(iiold) = nvert
c                  write(*,*) 'ijk is NOT in the hole', i,j,k
               end if
            end do
         end do
      end do
c      read(*,*) iii
      call cubtet(nop)
      jjkk = 1
      kl = 1
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
C
      do k = 1 , nzm1
         km=mod(k+1,2)+1
         do j = 1 , nym1 
            jm=2*mod(j+1,2)
            do i = 1 , nxm1
!!!!  jm=mod(j+1,2)+1
!!!!  im = mod(i+1,2)+1
!!!!  it = 2*(jm-1)+im
               it = jm + mod(i+1,2)+1
               if(km .ne. 1) then
                  it=5-it
               endif
!!               write(*,*) 'i,j,k,it',i,j,k,it
!               write(*,*) 'im,jm,km,it',im,jm,km,it
               if(  i .lt. mini .or.
     >              j .lt. minj .or.
     >              k .lt. mink .or. 
     >              i .ge. maxi .or.
     >              j .ge. maxj .or.
     >              k .ge. maxk) then
                  call xyzloc(xyz,jcub,nx,ny,nz,
     >                 i,j,k,hx,hy,hz,xmin,ymin,zmin)
!!!                  write(*,*) 'element= ', i,j,k,' type=',it
!!!                  write(*,'(a,8i5,a)') '(',(jcub(mm),mm=1,8),')'
                  do jk = 1 , 6
                     do mm = 1 , 4
                        inode = inumb(jcub(nop(mm,jk,it)))
!!!                        write(*,'(i7$)') inode
                        if(inode .eq. 0) then
                           write(*,*) 'ERROR:i,j,k and inode = 0', 
     >                          i,j,k,inode
                           stop 3
                        end if
                        xcoord(inode) = xyz(nop(mm,jk,it),1)
                        ycoord(inode) = xyz(nop(mm,jk,it),2)
                        zcoord(inode) = xyz(nop(mm,jk,it),3)
                        je(kl) = inode
                        kl = kl + 1
                     end do
!!                     write(*,*)
                     iflags(jjkk)=k
                     jjkk = jjkk + 1
                  end do
               end if
            end do
         end do
      end do
      nel = jjkk  - 1
      call bndry3(ib,nx,ny,nz,
     >     mini,minj,mink,maxi,maxj,maxk,inumb,
     >     minneu,maxneu)
      !Left, x=0 
      ibess=1
      jbess=0
      kbess=0
      call bess3(ib,nx,ny,nz,inumb,ibess,jbess,kbess,
     >     ibcode(1),minneu,maxneu)
      ! Right, x=xmax
      ibess=nx
      jbess=0
      kbess=0
      call bess3(ib,nx,ny,nz,inumb,ibess,jbess,kbess,
     >     ibcode(2),minneu,maxneu)
      ! Front, y=0
      ibess=0
      jbess=1
      kbess=0
      call bess3(ib,nx,ny,nz,inumb,ibess,jbess,kbess,
     >     ibcode(3),minneu,maxneu)
! Back, y=ymax
      ibess=0
      jbess=ny
      kbess=0
      call bess3(ib,nx,ny,nz,inumb,ibess,jbess,kbess,
     >     ibcode(4),minneu,maxneu)
! Bottom, nz=1
      ibess=0
      jbess=0
      kbess=1
      call bess3(ib,nx,ny,nz,inumb,ibess,jbess,kbess,
     >     ibcode(5),minneu,maxneu)
! Top, z=zmax
      ibess=0
      jbess=0
      kbess=nz
      call bess3(ib,nx,ny,nz,inumb,ibess,jbess,kbess,
     >     ibcode(6),minneu,maxneu)
      !!write(*,'(a,i12,a,i12)') ' Elements=', nel, '  Nodes=', nvert
      return
      end
C
      integer function nomxyz(i,j,k,nx,ny,nz)
      implicit real*8(a-h,o-z),integer (i-n)
      nomxyz=0
      if(  i.gt.0 .and. i.le.nx .and.
     >     j.gt.0 .and. j.le.ny .and.
     >     k.gt.0 .and. k.le.nz) then
         nomxyz = (k-1)*nx*ny + (j-1)*nx + i
      end if
      return
      end
