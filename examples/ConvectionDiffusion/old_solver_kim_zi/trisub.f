c=====================================================================
      program triangulation
c=====================================================================
      parameter(nlast=27 200 000, nrlast=27 200 000)
      implicit real*8 (a-h,o-z), integer*4 (i-n)
      common /top/ m(nlast)              !To use the machine witt
      common /top1/ r(nrlast)            !To use the machine witt
cc      dimension m(nlast), r(nrlast)
      integer n1x(10),n2x(10),n1y(10),n2y(10)
C
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
      WRITE (*,*) ' File name in 5 characters ?'
      read(*,'(a)') fname(1:5)
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
      ie = 1
      je = ie + nel + 1
      inona = je + nel * ndl
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
      call trian(iform,ibc,r(ix),r(iy),
     >     m(je),nodes,nel,ndl,nside,
     >     xa,xb,ya,yb,n1,n2,iflag,m(inona),r(irad),r(iphi),maxr,irglr)
C
      m(ie) = 1
      do k = 1 , nel
         m(ie+k) = m(ie+k-1) + ndl
      end do

      nsub = 4
      n1x(1) = 1
      n2x(1) = 6
      n1y(1) = 1
      n2y(1) = 6
      n1x(2) = 1
      n2x(2) = 6
      n1y(2) = 4
      n2y(2) = 9
      n1x(3) = 4
      n2x(3) = 9
      n1y(3) = 1
      n2y(3) = 6
      n1x(4) = 4
      n2x(4) = 9
      n1y(4) = 4
      n2y(4) = 9
C
      nelsub = 0
      nodesub = 0
      nbdry = 0
      do k = 1 , nsub
         nelsub = nelsub + (n2x(k)-n1x(k))*(n2y(k)-n1y(k))*2
         nodesub = nodesub + (n2x(k)-n1x(k)+1)*(n2y(k)-n1y(k)+1)
         nbdry = nbdry + 2*((n2x(k)-n1x(k)+1)+(n2y(k)-n1y(k)+1))
      end do
      write(*,*) nelsub,nodesub,nbdry
      read(*,*)
C
      isub = inona
      jsub = isub + nsub + 1
      isub_bdry = jsub + nelsub
      jsub_bdry =  isub_bdry + nsub + 1
      isub_all =  jsub_bdry +  nbdry
      jsub_all = isub_all + nsub + 1
      iwork = jsub_all + nodesub
      imlast = iwork + max0(nel,nodes)
C
      call formsub(nel,nodes,n1,n2,
     >     nsub,n1x,n2x,n1y,n2y,m(isub),m(jsub),m(ie),m(je),
     >     m(isub_bdry),m(jsub_bdry),m(isub_all),
     >     m(jsub_all),m(iwork))
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
         write(16,*)((nop(i,j),i=1,ndl),j=1,nel)
         write(16,*)(x(ik),y(ik),ik=1,nodes)
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
      subroutine formsub(nel,n,n1,n2,
     >     nsub,n1x,n2x,n1y,n2y,isub,jsub,ie,je,
     >     isub_bdry,jsub_bdry,isub_all,jsub_all,iwork)
C=====================================================================
      integer n1x(nsub),n2x(nsub),n1y(nsub),n2y(nsub),isub(1),jsub(1)
      integer isub_bdry(1),jsub_bdry(1),isub_all(1),jsub_all(1)
      integer ie(1),je(1),iwork(1),nel,n
C---------------------------------------------------------------------
C...  Generate mesh information about the each subdomain.
C...  
C...  Parameters
C...    NEL, N, N1, N2, IE, JE - Mesh information on the global domain
C...    NSUB    -  # of subdomains
C...    N1X     -  global node numbering of left end point in the 
C...               x-direction of each subdomain 
C...    N2X     -  global node numbering of right end point in the
C...               x-direction of each subdomain 
C...    N1Y     -  global node numbering of bottom end point in the
C...               y-direction of each subdomain 
C...    N2Y     -  global node numbering of top end point in the
C...               y-direction of each subdomain 
C---------------------------------------------------------------------
C
C...  ERROR STOP
C
      if(n1*n2 - n .ne.  0) stop 
C
      call inullv(iwork,n)
C
      isub_all(1) = 1
      ipoint = 1
      do k = 1 , nsub
         k1x = n1x(k)
         k2x = n2x(k)
         k1y = n1y(k)
         k2y = n2y(k)
         do jx = k1x,k2x
            do jy = k1y,k2y
C...           Get the  node number 
               nodenum = nomxy(jx,jy,n1,n2)
               jsub_all(ipoint) = nodenum
               ipoint = ipoint + 1
            end do
         end do
         isub_all(k+1) = ipoint
      end do
C
      call lpri(isub_all,jsub_all,nsub)
C
      isub(1) = 1
      ipoint = 1
      do k = 1 , nsub
         do jk = isub_all(k),isub_all(k+1)-1
            iwork(jsub_all(jk)) = k
         end do
         do iel = 1 , nel
            do jel = ie(iel),ie(iel+1)-1
               nodein = je(jel)
c               read(*,*)
               if(iwork(nodein) .ne. k) go to 100
            end do
            jsub(ipoint) = iel
            ipoint = ipoint + 1
 100        continue
         end do
         isub(k+1) = ipoint
      end do
C
      write(*,*) '---------------------------------------------------',
     >     '--------------'
      call lpri(isub,jsub,nsub)
C
      isub_bdry(1) = 1
      ipoint = 1
      do k = 1 , nsub
         k1x = n1x(k)
         k2x = n2x(k)
         k1y = n1y(k)
         k2y = n2y(k)
         do jx = k1x,k2x
            if (jx .eq. k1x .or. jx .eq. k2x) then
               do jy = k1y,k2y
C...              Get the node number 
                  nodenum = nomxy(jx,jy,n1,n2)
                  jsub_bdry(ipoint) = nodenum
                  ipoint = ipoint + 1
               end do
            else 
C...           Get the node number 
               nodenum = nomxy(jx,k1y,n1,n2)
               jsub_bdry(ipoint) = nodenum
               ipoint = ipoint + 1
C...           Get the node number 
               nodenum = nomxy(jx,k2y,n1,n2)
               jsub_bdry(ipoint) = nodenum
               ipoint = ipoint + 1
            end if
         end do
         isub_bdry(k+1) = ipoint
      end do
C
      write(*,*) '---------------------------------------------------',
     >     '--------------'
      call lpri(isub_bdry,jsub_bdry,nsub)
C
      return
      end
C=====================================================================
      function nomxy(i,j,nx,ny)
C=====================================================================
C---------------------------------------------------------------------
C...  NOMXY gives the global number of the node (i,j)
C---------------------------------------------------------------------
      NOMXY = (i-1)*ny + j
C
      if(i .le. 0 .or. j .le. 0) then
         nomxy = 0
      end if
      return
      end
C====================================================================
      subroutine lpri(ip,jp,n)
C====================================================================
      integer*4 ip(1),jp(1),n,k,kk
C
      do k = 1, n
         write(*,*) ' subdomain ', k, ' entries ',ip(k+1)-ip(k)
         write(*,*)' elements: ', (jp(kk),kk=ip(k),ip(k+1)-1)
         if(mod(k,10) .eq. 1) read(*,*)
      end do
      return
      end
C====================================================================
      subroutine lprr(ip,jp,p,n)
C====================================================================
      real*8 p(1)
      integer*4 ip(1),jp(1),n,k,kk
C
      do k = 1, n
         write(*,*) ' subdomain: ', k
         write(*,*)(jp(kk),p(kk),kk=ip(k),ip(k+1)-1)
         if(mod(k,10) .eq. 1) read(*,*)
      end do
      return
      end
C====================================================================
      subroutine inullv(iu,n)
C====================================================================
      dimension iu(n)
C--------------------------------------------------------------------
C...  IU = 0
C--------------------------------------------------------------------
      do i = 1 , n
         iu(i) = 0
      end do
      return
      end
C====================================================================
