C
C compute the local stiffness matrix
C
      subroutine local(xx,yy,ndl,stiff,rhs)
      implicit real*8(a-h,o-z)
      dimension xx(3), yy(3), stiff(ndl,ndl), rhs(ndl)
      common /menu/ ich(200) 

      lmethod = ich(4)
      goto (10,20,30) lmethod
      print *,'Error: no such method in LOCAL '
      stop
C
 10   continue
C...  EAFE method
      call elt_stf_eafe(xx,yy,ndl,stiff,rhs)
      return
      
 20   continue
C
C     Standart Galerkin
C     
      call elt_stf(xx,yy,ndl,stiff,rhs)
      return
      
 30   continue
C
C... Mixture between standard Galerkin and EAFE.
C      
      rho1 = dsqrt(xx(1)*xx(1) + yy(1)*yy(1))
      rho2 = dsqrt(xx(2)*xx(2) + yy(2)*yy(2))
      rho3 = dsqrt(xx(3)*xx(3) + yy(3)*yy(3))
C
      h1 = 1.d00/16d00+0.00001
      h1 = -1.d00
      rho0 = 0.75

c      if(
c     >     (rho1 .le. rho0 - h1 .or. rho1 .ge. rho0 + h1) .and.
c     >     (rho2 .le. rho0 - h1 .or. rho2 .ge. rho0 + h1) .and.
c     >     (rho3 .le. rho0 - h1 .or. rho3 .ge. rho0 + h1)
c     >     ) then
      if(
     >     ((xx(1) .le. h1 .or. yy(1) .ge. 1.d00-h1)) .and. 
     >     ((xx(2) .le. h1 .or. yy(2) .ge. 1.d00-h1)) .and. 
     >     ((xx(3) .le. h1 .or. yy(3) .ge. 1.d00-h1)) 
C
     >     ) then
         call elt_stf_eafe(xx,yy,ndl,stiff,rhs)
ccccccccccc         write(*,*) ' EAFE '
      else
         call elt_stf_1(xx,yy,ndl,stiff,rhs)
ccccccccccc         write(*,*) ' Galerkin '
      end if
      return
C
      end
