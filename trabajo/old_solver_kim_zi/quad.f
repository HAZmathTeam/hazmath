C=====================================================================
      subroutine quad_elt
C=====================================================================
      implicit real*8(a-h,o-z)
      common /quad/ wg(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  This subroutine supplies the quadrature data for numerical 
C...  integration on the reference triagle with vertices (0,0), (1,0), 
C...  and (0,1).
C...
C...  Input:
C...    LQUAD  - The choice of quadrature; 
C...             1 for 1 point, 3 for 3 points, and 7 for 7 points.
C...
C...  Parameter:
C...    WG     - the weights corresponding to the integation points.
C...    QUADPT - the integation points;
C...             QUADPT(1,j) = x-coordinate of j-th point,
C...             QUADPT(2,j) = y-coordinate of j-th point.
C...    PHI    - the function values of 3 basis functions at each of 
C...             integration points. The 3 basis functions are
C...             PHI(1,(x,y)) = 1 - x - y, 
C...             PHI(2,(x,y)) = x, and PHI(3,(x,y)) = y.
C...    PHIX   - the corresponding derivatives of PHI w.r.t. x.
C...    PHIY   - the corresponding derivatives of PHI w.r.t. y.
C---------------------------------------------------------------------
      lquad = ich(11)
      d12 = 0.5d0
      d13 = 1.d0/3
C
      if (lquad .ne. 1) go to 200
C
      wg(1) = d12
      quadpt(1,1) = d13
      quadpt(2,1) = d13
      phi(1,1) = d13
      phi(2,1) = d13
      phi(3,1) = d13
      phix(1,1) = - 1.d0
      phix(2,1) =   1.d0
      phix(3,1) =   0.d0
      phiy(1,1) = - 1.d0
      phiy(2,1) =   0.d0
      phiy(3,1) =   1.d0
      go to 500
C
 200  if (lquad .ne. 3) go to 300
C
      d16 = 1.d0/6
      wg(1) = d16
      wg(2) = d16
      wg(3) = d16
      quadpt(1,1) = d12
      quadpt(2,1) = 0.d0
      quadpt(1,2) = d12
      quadpt(2,2) = d12
      quadpt(1,3) = 0.d0
      quadpt(2,3) = d12
      phi(1,1) = d12
      phi(2,1) = d12
      phi(3,1) = 0.d0
      phi(1,2) = 0.d0
      phi(2,2) = d12
      phi(3,2) = d12
      phi(1,3) = d12
      phi(2,3) = 0.d0
      phi(3,3) = d12
      do 250 j = 1,3
         phix(1,j) = - 1.d0
         phix(2,j) =   1.d0
         phix(3,j) =   0.d0
         phiy(1,j) = - 1.d0
         phiy(2,j) =   0.d0
         phiy(3,j) =   1.d0
 250  continue
      go to 500
C
 300  if (lquad .ne. 7) go to 400
C
      d140 = 0.025d0
      d115 = 1.d0/15
      d940 = 0.225d0
      wg(1) = d140
      wg(2) = d140
      wg(3) = d140
      wg(4) = d115
      wg(5) = d115
      wg(6) = d115
      wg(7) = d940
      quadpt(1,1) = 0.d0
      quadpt(2,1) = 0.d0
      quadpt(1,2) = 1.d0
      quadpt(2,2) = 0.d0
      quadpt(1,3) = 0.d0
      quadpt(2,3) = 1.d0
      quadpt(1,4) = d12
      quadpt(2,4) = 0.d0
      quadpt(1,5) = d12
      quadpt(2,5) = d12
      quadpt(1,6) = 0.d0
      quadpt(2,6) = d12
      quadpt(1,7) = d13
      quadpt(2,7) = d13
      phi(1,1) = 1.d0
      phi(2,1) = 0.d0
      phi(3,1) = 0.d0
      phi(1,2) = 0.d0
      phi(2,2) = 1.d0
      phi(3,2) = 0.d0
      phi(1,3) = 0.d0
      phi(2,3) = 0.d0
      phi(3,3) = 1.d0
      phi(1,4) = d12
      phi(2,4) = d12
      phi(3,4) = 0.d0
      phi(1,5) = 0.d0
      phi(2,5) = d12
      phi(3,5) = d12
      phi(1,6) = d12
      phi(2,6) = 0.d0
      phi(3,6) = d12
      phi(1,7) = d13
      phi(2,7) = d13
      phi(3,7) = d13
      do 350 j = 1,7
         phix(1,j) = - 1.d0
         phix(2,j) =   1.d0
         phix(3,j) =   0.d0
         phiy(1,j) = - 1.d0
         phiy(2,j) =   0.d0
         phiy(3,j) =   1.d0
 350  continue
      go to 500
C     
 400  print*, 'Error: No such choice of lquad, lquad = ', lquad,'.'
      print*, 'Default value was used: lquad = 3.'
      lquad = 3
      go to 200
C
 500  return 
      end
C=====================================================================
      subroutine quad_data
C=====================================================================
      implicit real*8(a-h,o-z)
      common /quad_gs_ntrl/ z(7),zw(5,7)
      common /menu/ ich(200)
      dimension w(7)
C---------------------------------------------------------------------
C...  This subroutine supplies the Gauss quadrature data for numerical 
C...  integration on the reference edge I = [-1, 1]. The common data 
C...  set, QUAD_GS_NTRL, is used for the integration on the Natural
C...  (Neumann) boundary edges.
C...
C...  Input:
C...    LQUAD_GS - The choice of quadrature; it is N for N points.
C...
C...  Parameter:
C...    W(k) - the weights corresponding to the Gauss quadrature 
C...           points, G_k.
C...    Z(k) - the ratio between the length of the line segment
C...           [-1,G_k] and that of the line segment I, i.e., the 
C...           value, (G_k + 1)/2. Also, it can be interpreted in two 
C...           other ways: 
C...           1) If we think of the line segment J=[a,b] on which an
C...              integration will be performed, then the point
C...              a + Z(k)*(b-a) will the gauss point on the interval
C...              J corresponding to the gauss point G_k on I.
C...           2) If we think of the linear basis function PHI2(x) 
C...              which has the value 1 at the right end point of the
C...              edge I, then Z(k) = PHI2(G_k), the value of PHI2 at
C...              the k-th gauss point.
C...    ZW   - the multiplication of the values of basis functions at
C...           the Gauss points by the corresponding weights. Let PHI1
C...           and PHI2 are the basis functions which are 1 at the 
C...           right and at the left end points, respectively, of the
C...           interval I, then we define
C...              ZW(1,k) = PHI1(G_k)*w(k), 
C...              ZW(2,k) = PHI2(G_k)*w(k),
C...              ZW(3,k) = PHI1(G_k)*PHI1(G_k)*w(k),
C...              ZW(4,k) = PHI1(G_k)*PHI2(G_k)*w(k),
C...              ZW(5,k) = PHI2(G_k)*PHI2(G_k)*w(k).
C---------------------------------------------------------------------
C
      lquad_gs = ich(12)
 1    go to (10,20,30,40,50,60,70) lquad_gs
 5    print*, ' Invalid choice of LQUAD_GS in Gauss quadrature. '
      print*, ' Default -- two point. '
      lquad_gs = 2
      go to 1
C
 10   continue
      z(1) = 0.d0
      w(1) = 2.d0
      go to 100
C
 20   continue
      z(2) = 0.57735 02691 89625 76450 9149
      w(2) = 1.d0
      go to 100
C
 30   continue
      z(2) = 0.d0
      z(3) = 0.77459 66692 41483 37703 5835
      w(2) = 0.88888 88888 88888 88888 889
      w(3) = 0.55555 55555 55555 55555 556
      go to 100
C
 40   continue
      z(3) = 0.33998 10435 84856 26480 2666
      z(4) = 0.86113 63115 94052 57522 3946
      w(3) = 0.65214 51548 62546 14262 694
      w(4) = 0.34785 48451 37453 85737 306
      go to 100
C
 50   continue
      z(3) = 0.d0
      z(4) = 0.53846 93101 05683 09103 6314
      z(5) = 0.90617 98459 38663 99279 7627
      w(3) = 0.56888 88888 88888 88888 889
      w(4) = 0.47862 86704 99366 46804 129
      w(5) = 0.23692 68850 56189 08751 426
      go to 100
C
 60   continue
      z(4) = 0.23861 91860 83196 90863 0502
      z(5) = 0.66120 93864 66264 51366 1400
      z(6) = 0.93246 95142 03152 02781 2302
      w(4) = 0.46791 39345 72691 04738 987
      w(5) = 0.36076 15730 48138 60756 983
      w(6) = 0.17132 44923 79170 34504 030
      go to 100
C
 70   continue
      z(4) = 0.d0
      z(5) = 0.40584 51513 77397 16690 6607
      z(6) = 0.74153 11855 99394 43986 3865
      z(7) = 0.94910 79123 42758 52452 6190
      w(4) = 0.41795 91836 73469 38775 510
      w(5) = 0.38183 00505 05118 94495 037
      w(6) = 0.27970 53914 89276 66790 147
      w(7) = 0.12948 49661 68869 69327 061
      go to 100
C
 100  do 110 j = 1,int(lquad_gs/2)
         w(j) =  w(lquad_gs+1-j)
         z(j) = -z(lquad_gs+1-j)
 110  continue
C
      do 120 j = 1,lquad_gs
	 z(j) = (z(j) + 1.d0)*0.5d0
 120  continue
C
      do 130 j = 1,lquad_gs
	 zw(1,j) = z(lquad_gs+1-j)*w(j)
	 zw(2,j) = z(j)           *w(j)
	 zw(3,j) = z(lquad_gs+1-j)*zw(1,j)
	 zw(4,j) = z(j)           *zw(1,j)
	 zw(5,j) = z(j)           *zw(2,j)
 130  continue
C
      return
      end
C=====================================================================
      subroutine quad_data1
C=====================================================================
      implicit real*8(a-h,o-z)
      common /quad_gs/ w(7),z(7)
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  This subroutine supplies the Gauss quadrature data for numerical 
C...  integration on the reference edge I = [-1, 1].
C...
C...  Input:
C...    LQUAD_GS - The choice of quadrature; it is N for N points.
C...
C...  Parameter:
C...    W(k) - the weights corresponding to the Gauss quadrature 
C...           points, G_k.
C...    Z(k) - the ratio between the length of the line segment
C...           [-1,G_k] and that of the line segment I, i.e., the 
C...           value, (G_k + 1)/2. Also, it can be interpreted in 
C...           the another way; If we think of the line SEGMENT J=[A,B]
C...           on which an integration will be performed, then the 
C...           point a + z(k)*(b-a) will the gauss point on the 
C...           interval j corresponding to the gauss point g_k on i.
C---------------------------------------------------------------------
C
      lquad_gs = ich(13)
 1    go to (10,20,30,40,50,60,70) lquad_gs
 5    print*, ' Invalid choice of lquad_gs in gauss quadrature. '
      print*, ' Default -- 3 points. '
      lquad_gs = 3
      go to 1
C
 10   continue
      z(1) = 0.d0
      w(1) = 2.d0
      go to 100
C
 20   continue
      z(2) = 0.57735 02691 89625 76450 9149
      w(2) = 1.d0
      go to 100
C
 30   continue
      z(2) = 0.d0
      z(3) = 0.77459 66692 41483 37703 5835
      w(2) = 0.88888 88888 88888 88888 889
      w(3) = 0.55555 55555 55555 55555 556
      go to 100
C
 40   continue
      z(3) = 0.33998 10435 84856 26480 2666
      z(4) = 0.86113 63115 94052 57522 3946
      w(3) = 0.65214 51548 62546 14262 694
      w(4) = 0.34785 48451 37453 85737 306
      go to 100
C
 50   continue
      z(3) = 0.d0
      z(4) = 0.53846 93101 05683 09103 6314
      z(5) = 0.90617 98459 38663 99279 7627
      w(3) = 0.56888 88888 88888 88888 889
      w(4) = 0.47862 86704 99366 46804 129
      w(5) = 0.23692 68850 56189 08751 426
      go to 100
C
 60   continue
      z(4) = 0.23861 91860 83196 90863 0502
      z(5) = 0.66120 93864 66264 51366 1400
      z(6) = 0.93246 95142 03152 02781 2302
      w(4) = 0.46791 39345 72691 04738 987
      w(5) = 0.36076 15730 48138 60756 983
      w(6) = 0.17132 44923 79170 34504 030
      go to 100
C
 70   continue
      z(4) = 0.d0
      z(5) = 0.40584 51513 77397 16690 6607
      z(6) = 0.74153 11855 99394 43986 3865
      z(7) = 0.94910 79123 42758 52452 6190
      w(4) = 0.41795 91836 73469 38775 510
      w(5) = 0.38183 00505 05118 94495 037
      w(6) = 0.27970 53914 89276 66790 147
      w(7) = 0.12948 49661 68869 69327 061
      go to 100
C
 100  do 110 j = 1,int(lquad_gs/2)
         w(j) =  w(lquad_gs+1-j)
         z(j) = -z(lquad_gs+1-j)
 110  continue
C
      do 120 j = 1,lquad_gs
	 z(j) = (z(j) + 1.d0)*0.5d0
 120  continue
C
      return
      end
C=====================================================================














