c=====================================================================
      program triangulation
c=====================================================================
cc      parameter(nlast=60 000 000, nrlast=30 000 000)
      parameter(nlast=8 000 000, nrlast=4 000 000)
      implicit real*8 (a-h,o-z), integer*4 (i-n)
      common /top/ m(nlast)              !To use the machine witt
      common /top1/ r(nrlast)            !To use the machine witt
cc      dimension m(nlast), r(nrlast)
      character*132 fname
C---------------------------------------------------------------------
C...  On the rectangle of size (XA,XB) x(YA,YB) this subroutine generates
C...  a triangular mesh with N1 grid points on the x-axis and N2 grid 
C...  points on the y-axis. Various boundary conditions are imposed on 
C...  the boundary.
C---------------------------------------------------------------------
      xa = 0.d0
      xb = 1.d0
      ya = 0.d0
      yb = 1.d0
!
      WRITE (*,*) ' File name in 6 characters ?'
      read(*,'(a)') fname(1:6)
!      WRITE (*,*) ' Xmin, Xmax, Ymin, Ymax in double precision ?'
!      read(*,*) xa,xb,ya,yb
      WRITE (*,*) ' Input n1,n2'
      read(*,*) n1,n2
      WRITE (*,*) ' Uniform mesh :1; jiggled one: 2 ? '
      read(*,*) irglr
      WRITE (*,*) ' Formatted:1; unformatted : 2 ? '
      read(*,*) iform
      WRITE (*,*) ' 1 === Pure Dirichlet conditions '
      WRITE (*,*) ' 2 === Neumann conditions on top '
      WRITE (*,*) ' 3 === Neumann conditions on top and right '
      WRITE (*,*) ' 4 === Neumann conditions on top and bottom '
      WRITE (*,*) ' 5 === Neumann conditions on top and left and right'
      WRITE (*,*) ' 6 === Neumann conditions everywhere '
      WRITE (*,*) ' 7 === Dirichlet BC only on [0,0.25] and [0.75,1]
     >     in both directions. '
      WRITE (*,'(a$)') ' B.C. :: 1,2,3,4,5 6,7?   '
      read(*,*) ibc
C
      nodes = n1*n2
      nel = n1*n2-n1-n2+1
      nel = 2*nel
      write(*,'(i14)') nodes,nel
C     
      if(iform .eq. 1) then
         open(16,file=fname,form='formatted',status='unknown')
      else
         open(16,file=fname,
     >        form='unformatted',status='unknown')
cc         open(16,file=fname,access='sequential',
cc     >        form='unformatted',status='unknown')
      end if
C
      itwp  = 2
      ndl   = 3
      nside = 3
C     
      ix     = 1
      iy     = ix + nodes       
      irad = iy + nodes  
      iphi = irad + maxr
      irlast = iphi + maxr
C
      inop = 1
      inona = inop + nel * ndl
      imlast = inona + nodes
C
      if (imlast .gt. nlast .or. irlast .gt. nrlast) then
         write(*,*) '***** Insuficient storage ***** real*8 memo.', 
     >        irlast,' -- used' , nrlast,' required'
         write(*,*) '***** Insuficient storage ***** integer*4 memo.', 
     >        imlast,' -- used' , nlast,' required'
         stop
      else
         write(*,*) 
     >        ' Storage for real*8 ',irlast, ' for integer*4 ',imlast
      end if
C
      do 60 i = 1,imlast
         m(i) = 0
 60   continue
C
      do 70 i = 1,irlast
         r(i) = 0
 70   continue
C
      call trian(iform,ibc,r(ix),r(iy),m(inop),nodes,nel,ndl,nside,
     >     xa,xb,ya,yb,n1,n2,iflag,m(inona),r(irad),r(iphi),maxr,irglr)
C
      endfile(16)
      close(16)
      stop
      end
C=====================================================================
      subroutine trian(iform,ibc,x,y,nop,nodes,nel,ndl,nside,
     >     xa,xb,ya,yb,n1,n2,iflag,nona,rad,phi,maxr,irglr)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension x(1),y(1),rad(1),phi(1),
     >          nop(ndl,nel),nona(1),nos3(2,3),nver(3)
      data nos3 /1,2,2,3,3,1/,nver/3,1,2/
C---------------------------------------------------------------------
C...  On the rectangle of size (XA,XB) x(YA,YB) this subroutine generates
C...  a triangular mesh with N1 grid points on the x-axis and N2 grid 
C...  points on the y-axis. Various boundary conditions are imposed on 
C...  the boundary.
C---------------------------------------------------------------------
      xslp = 0.5d0*(xb - xa) 
      yslp = 0.5d0*(yb - ya) 
      xint = 0.5d0*(xb + xa) 
      yint = 0.5d0*(yb + ya) 
C
      x0 = -1.d0
      y0 = -1.d0
      hx = 2.d0/dfloat(n1-1)
      hy = 2.d0/dfloat(n2-1)

      do 20 l = 1 , n1
         do 10 k = 1 , n2
            indx = nomxy(l,k,n1,n2)
            x(indx) = x0 + dfloat(l-1)*hx
            y(indx) = y0 + dfloat(k-1)*hy
 10      end do
 20   end do
C
      do 40 i = 1 , n1-1
         do 30 j = 1 , n2-1
            nmnel1 = 2*nomxy(i,j,n1-1,n2-1)-1
            nmnel2 = 2*nomxy(i,j,n1-1,n2-1)
            nop(1,nmnel1) = nomxy(i,j,n1,n2)
            nop(2,nmnel1) = nomxy(i+1,j,n1,n2)
            nop(3,nmnel1) = nomxy(i+1,j+1,n1,n2)
            nop(1,nmnel2) = nomxy(i,j,n1,n2)
            nop(2,nmnel2) = nomxy(i+1,j+1,n1,n2)
            nop(3,nmnel2) = nomxy(i,j+1,n1,n2)
 30      end do
 40   end do
C     
      if(irglr .ne. 1) then
         call jiji(x,y,nodes,n1,n2)     
         xmin = 1.d20
         ymin = 1.d20
         xmax = -1.d20
         ymax = -1.d20
         do 50 i = 1, nodes
            xmin = dmin1(xmin,x(i))
            xmax = dmax1(xmax,x(i))
            ymin = dmin1(ymin,y(i))
            ymax = dmax1(ymax,y(i))
 50      end do
C     
         if(dabs(xmin-xmax) .gt. 1.d-7) then 
            delx1 = 2.d0/(xmax-xmin)
            delx2 = -(xmin+xmax)*0.5*delx1
         else 
            delx1 = 0.0d00
            delx2 = 0.0d00
         end if
         if(dabs(ymin-ymax) .gt. 1.d-07) then
            dely1 = 2.d0/(ymax-ymin)
            dely2 = -(ymin+ymax)*0.5*dely1
         else
            dely1 = 0.0d00
            dely2 = 0.0d00
         end if
C
         do 60 i =1 , nodes
            x(i) =  delx1*x(i) + delx2
            y(i) =  dely1*y(i) + dely2
 60      end do
      end if
C
      do 70 k = 1 , nodes
         x(k) = xslp*x(k) + xint
         y(k) = yslp*y(k) + yint
 70   end do
C
C...  OUTPUT
C
      ned = 2*n1 + 2*n2 - 4
      if(iform .eq. 1) then
         write(16,*)n1,n2
         write(16,*)nel,nodes,ned
cc         write(16,*)((nop(i,j),i=1,ndl),j=1,nel)
         do j = 1, nel
            write(16,*)(nop(i,j),i=1,ndl)
         end do
cc         write(16,*)(x(ik),y(ik),ik=1,nodes)
         do ik = 1, nodes
            write(16,*) x(ik), y(ik)
         end do
      else
         write(16)n1,n2
         write(*,*) ' first record ok '
         write(16)nel,nodes,ned
         write(*,*) ' second record ok '
         write(16)((nop(i,j),i=1,ndl),j=1,nel)
         write(*,*) ' third record ok ',i,j
         write(16)(x(ik),y(ik),ik=1,nodes)
         write(*,*) ' forth record ok '
      end if
C
      go to (110,120,130,140,150,160,170),ibc
 110  continue
      nfirstx = n1
      nfirsty = n2
      nseconx = 1
      nsecony = 1
      go to 200
 120  continue   
      nfirstx = n1
      nfirsty = n2
      nseconx = n1
      nsecony = 1
      go to 200
 130  continue   
      nfirstx = n1
      nfirsty = n2
      nseconx = n1
      nsecony = n2
      go to 200
 140  continue   
      nfirstx = 1
      nfirsty = n2
      nseconx = n1
      nsecony = 1
      go to 200
 150  continue   
      nfirstx = n1
      nfirsty = 1
      nseconx = n1
      nsecony = n2
      go to 200
 160  continue
      nfirstx = 1
      nfirsty = 1
      nseconx = n1
      nsecony = n2
      go to 200
 170  continue
      nfirstx = 0
      nseconx = 0
      ky1 = 1
      kyn = n2
      do 180 kx = 1 , n1
         nx1 = nomxy(kx,ky1,n1,n2)
         nx2 = nomxy(kx,kyn,n1,n2)
         if(x(nx1) .le. 0.25000001) nfirstx = kx
         if(x(nx2) .lt. 0.7500001) nseconx = kx
 180  end do
      nfirsty = 0
      nsecony = 0
      kx1 = 1
      kxn = n1
      do 190 ky = 1 , n2
         ny1 = nomxy(kx1,ky,n1,n2)
         ny2 = nomxy(kxn,ky,n1,n2)
         if(y(ny1) .le. 0.25000001) nfirsty = ky
         if(y(ny2) .lt. 0.750001) nsecony = ky
 190  end do
C
 200  continue
C
      if (iform .eq. 1) then
         write(16,'(3i10)')
     >      (nomxy(1,j,n1,n2), nomxy(1,j+1,n1,n2),1,j=1,nfirsty-1),
     >      (nomxy(n1,j,n1,n2), nomxy(n1,j+1,n1,n2),2,j=nsecony,n2-1),
     >      (nomxy(i,1,n1,n2), nomxy(i+1,1,n1,n2),3,i=1,nfirstx-1),
     >      (nomxy(i,n2,n1,n2), nomxy(i+1,n2,n1,n2),4,i=nseconx,n1-1),
C...  THE REST IS NEUMANN
     >      (nomxy(1,j,n1,n2), nomxy(1,j+1,n1,n2),-1,j=nfirsty,n2-1),
     >      (nomxy(n1,j,n1,n2), nomxy(n1,j+1,n1,n2),-2,j=1,nsecony-1),
     >      (nomxy(i,1,n1,n2), nomxy(i+1,1,n1,n2),-3,i=nfirstx,n1-1),
     >      (nomxy(i,n2,n1,n2), nomxy(i+1,n2,n1,n2),-4,i=1,nseconx-1)
      else
         write(16)
     >      (nomxy(1,j,n1,n2), nomxy(1,j+1,n1,n2),1,j=1,nfirsty-1),
     >      (nomxy(n1,j,n1,n2), nomxy(n1,j+1,n1,n2),2,j=nsecony,n2-1),
     >      (nomxy(i,1,n1,n2), nomxy(i+1,1,n1,n2),3,i=1,nfirstx-1),
     >      (nomxy(i,n2,n1,n2), nomxy(i+1,n2,n1,n2),4,i=nseconx,n1-1),
C...  THE REST IS NEUMANN
     >      (nomxy(1,j,n1,n2), nomxy(1,j+1,n1,n2),-1,j=nfirsty,n2-1),
     >      (nomxy(n1,j,n1,n2), nomxy(n1,j+1,n1,n2),-2,j=1,nsecony-1),
     >      (nomxy(i,1,n1,n2), nomxy(i+1,1,n1,n2),-3,i=nfirstx,n1-1),
     >      (nomxy(i,n2,n1,n2), nomxy(i+1,n2,n1,n2),-4,i=1,nseconx-1)
      end if
      endfile(16)
      close(16)
      return
      end
C=====================================================================
      function nomxy(i,j,nx,ny)
C=====================================================================
C---------------------------------------------------------------------
C...  NOMXY gives the global number of the node (i,j)
C---------------------------------------------------------------------
      NOMXY = (i-1)*ny + j
      if(i .le. 0 .or. j .le. 0) then
         nomxy = 0
      end if
      return
      end
C=====================================================================
      subroutine jiji(x,y,nodes,n1,n2)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension x(1),y(1),xel(1625),yel(1625)
      integer modei(40),modej(40)
      data pi/3.141592653589793/
C---------------------------------------------------------------------
C...  
C---------------------------------------------------------------------
      hx = 2.d0/dfloat(n1-1)
      hy = 2.d0/dfloat(n2-1)
      xx0 = 1.d00
      yy0 = 1.d00
      modei(1) = 0
      modej(1) = 0
      nhh = 3
      if(n1 .gt. 30) nhh =n1/15
      mmmi = (n1-mod(n1,nhh))/nhh  
      mmmj = (n2-mod(n2,nhh))/nhh
C
      do 10 i = 1,mmmi-1
         modei(i+1) = modei(i) + nhh
 10   end do 
      do 20 i = 1,mmmj-1
         modej(i+1) = modej(i) + nhh
 20   end do 
      modei(mmmi+1) = n1
      modej(mmmj+1) = n2
C
      write(*,*)(modei(i),i=1,mmmi+1)
      write(*,*)(modej(i),i=1,mmmj+1)
c     read(*,*)
C
      smx = 1.1
      smy = 1.1
C
      if(mmmi+1 .ge. 9) smx = 2.25
      if(mmmj+1 .ge. 9) smy = 2.25
      do 40 i = 1 , mmmi+1
         do 30 j = 1 , mmmj+1
            nomm = nomxy(i,j,mmmi+1,mmmj+1)
            xel(nomm) = dsin((i+j)*pi*smx/dble(mmmi+3))+i
            yel(nomm) = dcos((i-j)*pi*smy/dble(mmmj+3))+j
 30      end do
 40   end do
C
      do 80 kk = 1 , mmmi
         do 70 ll = 1 , mmmj
            xel1 = xel(nomxy(kk,ll,mmmi+1,mmmj+1))
            yel1 = yel(nomxy(kk,ll,mmmi+1,mmmj+1))
            xel2 = xel(nomxy(kk+1,ll,mmmi+1,mmmj+1))
            yel2 = yel(nomxy(kk+1,ll,mmmi+1,mmmj+1))
            xel3 = xel(nomxy(kk,ll+1,mmmi+1,mmmj+1))
            yel3 = yel(nomxy(kk,ll+1,mmmi+1,mmmj+1))
            xel4 = xel(nomxy(kk+1,ll+1,mmmi+1,mmmj+1))
            yel4 = yel(nomxy(kk+1,ll+1,mmmi+1,mmmj+1))
            nomer1 = nomxy(modei(kk)+1,modej(ll)+1,n1,n2)
            x0 = x(nomer1)-hx
            y0 = y(nomer1)-hy
            xn = x0 + dble(modei(kk+1)-modei(kk))*hx
            yn = y0 + dble(modej(ll+1)-modej(ll))*hy
            delx1 = 2./(xn-x0)
            dely1 = 2./(yn-y0)
            delx2 = -(x0+xn)/(xn-x0)
            dely2 = -(y0+yn)/(yn-y0)
cc            write(*,*) '     news ****'
cc            write(*,*) xel1,xel2,xel3,xel4
cc            write(*,*) yel1,yel2,yel3,yel4
cc            write(*,*) x0,xn
cc            write(*,*) y0,yn
cc            read(*,*)
C
            do 60 ikk = modei(kk)+1,modei(kk+1)
               do 50 jkk = modej(ll)+1,modej(ll+1)
                  i = nomxy(ikk,jkk,n1,n2)
                  r = x(i)*delx1+delx2
                  s = y(i)*dely1+dely2
                  phi1 = 0.25*(1-r)*(1-s)-0.125*(1-r*r)*(1-s*s)
                  phi2 = 0.25*(1+r)*(1-s)+0.125*(1-r*r)*(1-s*s)
                  phi3 = 0.25*(1-r)*(1+s)+0.125*(1-r*r)*(1-s*s)
                  phi4 = 0.25*(1+r)*(1+s)-0.125*(1-r*r)*(1-s*s)
                  x(i) = xel1*phi1+xel2*phi2+xel3*phi3+xel4*phi4
                  y(i) = yel1*phi1+yel2*phi2+yel3*phi3+yel4*phi4
cc                  write(*,*) ikk,jkk
 50            end do
 60         end do
 70      end do
 80   end do
      return
      end
C=====================================================================




