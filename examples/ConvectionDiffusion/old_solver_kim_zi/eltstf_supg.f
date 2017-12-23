C======================================================================
      subroutine elt_stf_upwind(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2),wrk(3),aed(3,3)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C
      call nullv(ae,nd*nd)
      call nullv(aed,nd*nd)
      call nullv(be,nd)
C
      call advection_upwind_matrix(xe,ye,wrk,ae)
C     
      call diffusion_element_matrix(xe,ye,wrk,aed)
C     
      call source_element_vector(xe,ye,be)
C     
      call uupluv(ae,aed,nd*nd)
C
      return
      end
C======================================================================
      subroutine elt_stf_supg(xe,ye,nd,ae,be)
C======================================================================
      implicit real*8(a-h,o-z)
      dimension xe(nd),ye(nd),ae(nd,nd),be(nd)
      dimension phx(3),phy(3),b(2,2),c(2,2),wrk(3),aed(3,3)
      common /quad/ w(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C
      call nullv(ae,nd*nd)
      call nullv(aed,nd*nd)
      call nullv(be,nd)
C
      call advection_element_matrix(xe,ye,wrk,ae)
C     
      call diffusion_element_matrix(xe,ye,wrk,aed)
C     
      call source_element_vector(xe,ye,be)
C     
      call uupluv(ae,aed,nd*nd)
C
      return
      end
C======================================================================
      subroutine diffusion_element_matrix(x,y,phi,diffusion_matrix)
C======================================================================
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3),phi(3)
      dimension diffusion_matrix(3,3)
      dimension element_normal(3,2)
C---------------------------------------------------------------------- 
c  GALERKIN DISCRETIZATION OF DIFFUSION TERMS: Assume C^0 linear elements
c
c  Compute     \int_T  w( (c \phi_x)_x + (d \phi_y)_y) \ d\,x d\,y = 
c       = - \int_T  (c w_x \phi_x  + d w_y \phi_y) \ d\,x d\,y 
c         + \int_{\partial T} w ((c \phi_x , d \phi_y) \cdot {\bf n}) dl
c
c  with either w=0 or (c \phi_x , d \phi_y) \cdot {\bf n} = 0 on the
c  boundary.
c
c   The element matrix is symmetric but we will store the full 
c   3x3 matrix anyway.
C----------------------------------------------------------------------
c
c (1) calculate the element area
c
      element_area = area(x(1),y(1),x(2),y(2),x(3),y(3))
c
c (2) calculate the 3 scaled edge normals
c
      do n = 1,3
         np  = mod(n,3) + 1
         npp = mod(np,3) + 1
         element_normal(n,1) =   y(npp)-y(np)
         element_normal(n,2) = -(x(npp)-x(np))
      enddo

c (3) compute the integral averaged values of the diffusion coeffs
c     by numerical quadrature
c
c     c_ave =  1/A \int_T c(x,y) dA,  d_ave = 1/a \int_T d(x,y) dA 
c
      call coefficients(x(1),y(1),a1,b1,c1,d1,s1)
      call coefficients(x(2),y(2),a2,b2,c2,d2,s2)
      call coefficients(x(3),y(3),a3,b3,c3,d3,s3)

      c_ave = (c1+c2+c3)/3        
      d_ave = (d1+d2+d3)/3        
c
c (4) form the Galerkin intgral exploiting the simple formula for 
c     linear elements:
c                                             ->           ->           ->
c     grad(\phi) = -1/(2 element_area) (\phi_1 n_1 + \phi_2 n_2 + \phi_3 n_3)
c
c     diffusion_matrix(n,m): contribution to n-th row of matrix by m-th col
c
      earea = 0.25/element_area
      do n=1,3
         do m=1,3
            diffusion_matrix(n,m) = 
ccccc     >           -( element_normal(n,1)*element_normal(m,1)*c_ave
     >           (element_normal(n,1)*element_normal(m,1)*c_ave
     >           +element_normal(n,2)*element_normal(m,2)*d_ave)*earea
         enddo
      enddo

      return
      end
C======================================================================
      subroutine advection_upwind_matrix(x,y,phi,advection_matrix)
C======================================================================
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3),phi(3),ak(3),akp(3),akm(3)
      dimension advection_matrix(3,3)
      dimension element_normal(3,2)
      parameter (epsilon = 1.e-12)
c----------------------------------------------------------------------
c  Upwind approximation on triangles (produces M-matrix)
c----------------------------------------------------------------------
c
c (1) calculate the element area
c
      element_area = area(x(1),y(1),x(2),y(2),x(3),y(3))
      if(element_area .lt. epsilon) then
         print * , ' area is not positive '
         stop 99
      end if
c
c (2) calculate the 3 scaled edge normals. 
c
      do n = 1,3
         np  = mod(n,3) + 1
         npp = mod(np,3) + 1
         element_normal(n,1) = -(y(npp)-y(np))
         element_normal(n,2) =  (x(npp)-x(np))
      enddo
c
c (3) compute the integral averaged values of the diffusion coeffs
c     by numerical quadrature
c
c     a_ave =  1/A \int_T a(x,y) dA,  b_ave = 1/A \int_T b(x,y) dA 
c
      call coefficients(x(1),y(1),a1,b1,c1,d1,s1)
      call coefficients(x(2),y(2),a2,b2,c2,d2,s2)
      call coefficients(x(3),y(3),a3,b3,c3,d3,s3)

      a_ave = (a1+a2+a3)/3        
      b_ave = (b1+b2+b3)/3
c
      do n=1,3
         ak(n) = 0.5d0*
     >        (a_ave*element_normal(n,1)+b_ave*element_normal(n,2))
         akp(n) = dmax1(0.d0,ak(n))
         akm(n) = dmin1(0.d0,ak(n))
      enddo

      asum = akp(1)+akp(2)+akp(3)
      
      if(dabs(asum) .gt. 1.d-12) then
         do n=1,3
c           term_n = akp(n)*phi(n)+
c     >            akp(n)*akm(1)/asum*phi(1)
c     >            akp(n)*akm(2)/asum*phi(2)
c     >            akp(n)*akm(3)/asum*phi(3)
            
            advection_matrix(n,1) = akp(n)*akm(1)/asum     
            advection_matrix(n,2) = akp(n)*akm(2)/asum     
            advection_matrix(n,3) = akp(n)*akm(3)/asum   
C     
            advection_matrix(n,n) = advection_matrix(n,n) + akp(n)      
         enddo
      else
         do n=1,3            
            advection_matrix(n,1) = 0.0d00
            advection_matrix(n,2) = 0.0d00
            advection_matrix(n,3) = 0.0d00
C     
            advection_matrix(n,n) = advection_matrix(n,n) + akp(n)
         enddo
      end if
c
c      do n=1,3
c
c      summm = 0.0d00
c         summm =  summm +    advection_matrix(n,1) +
c     >     advection_matrix(n,2) +       advection_matrix(n,3) 
c       write(*,*) summm
c      enddo
c      read(*,*)

      return
      end
C======================================================================
      subroutine advection_element_matrix(x,y,phi,advection_matrix)
C======================================================================
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3),phi(3)
      dimension advection_matrix(3,3)
      dimension element_normal(3,2)
      parameter (epsilon = 1.e-12)
C----------------------------------------------------------------------
C...  SUPG
C----------------------------------------------------------------------
c
c *** GALERKIN DISCRETIZATION OF ADVECTION with LEAST-SQUARES STABILIZATION ***
c
c     strategy: compute the Galerkin contribution (step 4) then add the
c              least-square integral (step 5)
c----------------------------------------------------------------------
c (1) calculate the element area
c
      element_area = area(x(1),y(1),x(2),y(2),x(3),y(3))
      if(element_area .lt. epsilon) then
         print * , ' area is not positive '
         stop 99
      end if
c
c (2) calculate the 3 scaled edge normals. 
c
      do n = 1,3
         np  = mod(n,3) + 1
         npp = mod(np,3) + 1
         element_normal(n,1) =   -(y(npp)-y(np))
         element_normal(n,2) = (x(npp)-x(np))
      enddo
c
c (3) compute the integral averaged values of the diffusion coeffs
c     by numerical quadrature
c
c     a_ave =  1/A \int_T a(x,y) dA,  b_ave = 1/A \int_T b(x,y) dA 
c
      call coefficients(x(1),y(1),a1,b1,c1,d1,s1)
      call coefficients(x(2),y(2),a2,b2,c2,d2,s2)
      call coefficients(x(3),y(3),a3,b3,c3,d3,s3)

      a_ave = (a1+a2+a3)/3        
      b_ave = (b1+b2+b3)/3
c
c (4) form the Galerkin contribution exploiting the simple formula 
c     for linear elements:
c                                              ->           ->           ->
c     grad(\phi) = -1/(2 element_area) (\phi_1 n_1 + \phi_2 n_2 + \phi_3 n_3)
c
c     compute  \int_T     w (a \phi_x + b \phi_y) \ d\,x d\,y
c         = -\int_T  \phi (a w_x + b w_y) \ d\,x d\,y
c           +\int_T  w \phi ((a,b) \cdot {\bf n}) \ d\,x d\,y
c
c     Advection_matrix(n,m): contribution to n-th row of matrix by m-th col
c
      do n=1,3
         do m=1,3
            advection_matrix(n,m) =   
     >      (element_normal(m,1)*a_ave + element_normal(m,2)*b_ave)/6
         enddo
      enddo
c
c (5) (SUPG) least-squares stabilization. Let \lambda = (a,b)
c
c     \int_T  (\lambda \cdot \nabla u) \sigma (\lambda \cdot \nabla w) dA
c
c     set characteristic element length  h = sqrt(area)
c
c     \sigma = h / |\lambda|
c
      sigma = dsqrt(element_area/(a_ave**2+b_ave**2 + epsilon**2))
cc      sigma = dsqrt(element_area/(a_ave**2+b_ave**2 + epsilon**2))*
cc     >     1.11111
cc     >     2.d00
      do n=1,3
         do m=1,3
            advection_matrix(n,m) = advection_matrix(n,m) +
cc            advection_matrix(n,m) =-
     >      sigma/element_area/4*(
     >      (element_normal(n,1)*a_ave + element_normal(n,2)*b_ave)*
     >      (element_normal(m,1)*a_ave + element_normal(m,2)*b_ave))
c           write(*,*)
c     >      sigma/element_area/4*(
c     >      (element_normal(n,1)*a_ave + element_normal(n,2)*b_ave)*
c     >      (element_normal(m,1)*a_ave + element_normal(m,2)*b_ave))
         enddo
      enddo

      return
      end
C======================================================================
      real*8 function area(x1,y1,x2,y2,x3,y3)
C======================================================================
      implicit real*8 (a-h,o-z)
c----------------------------------------------------------------------     
c     This function determines the area of a right-handed triangle.
c----------------------------------------------------------------------
      area =  .5*(x2*(y3-y1) + x3*(y1-y2) + x1*(y2-y3))
      area =  dabs(area)
      return
      end
C======================================================================
      subroutine source_element_vector(x,y,source_vector)
C======================================================================
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3),source_vector(3)
C
      element_area = area(x(1),y(1),x(2),y(2),x(3),y(3))

      call coefficients(x(1),y(1),a1,b1,c1,d1,s1)
      call coefficients(x(2),y(2),a2,b2,c2,d2,s2)
      call coefficients(x(3),y(3),a3,b3,c3,d3,s3)

      source_vector(1) = element_area*s1/3.0
      source_vector(2) = element_area*s2/3.0
      source_vector(3) = element_area*s3/3.0

      return
      end
C======================================================================
      subroutine coefficients(x,y,a,b,c,d,s)
C======================================================================
      implicit real*8 (a-h,o-z)
c----------------------------------------------------------------------
c    model equation   a u_x + b u_y - (c u_x)_x - (d u_y)_y  = s
c
c    with   a(x,y),b(x,y),c(x,y),d(x,y),s(x,y) \in R^2
c----------------------------------------------------------------------

      xx = x - .5
      yy = y - .5
c
c     advection-diffusion-source
c     
      pi = 3.141592653589793d00
C     
      a =  beta(1,x,y)
      b =  beta(2,x,y)
      c =  eps(x,y)
      d = c
      s = rhs(x,y)
C
      return
      end
C======================================================================
