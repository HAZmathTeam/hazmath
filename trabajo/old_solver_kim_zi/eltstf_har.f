*
* elt_stf_har.f
* This subroutine assembles the elementary stiffness matrix stiff(ndl,ndl)
* by harmonic average
* and the elementary right hand side rhs(ndl)
*
* xx,yy: the coordinates of the 3 nodes in the element

      subroutine elt_stf_har(xx,yy,ndl,stiff,rhs)
      implicit real*8(a-h,o-z)
      real*8 xx(3), yy(3), stiff(ndl,ndl), rhs(ndl)
      real*8 tr(2,3), cota(3,3)
      common /params/lbeta, lalpha, lmethod

*     compute the lengths of the edges
      e12 = sqrt((xx(2)-xx(1))**2
     $     + (yy(2)-yy(1))**2)
      e23 = sqrt((xx(3)-xx(2))**2
     $     + (yy(3)-yy(2))**2)
      e31 = sqrt((xx(1)-xx(3))**2
     $     + (yy(1)-yy(3))**2)

 3    continue

*     use Heron's formula to compute the area of the triangle
      s = (e12 + e23 + e31) * 0.5d00
      area = sqrt(s * (s-e12) * (s-e23) * (s-e31))

*     compute the transformation matrix from reference triangle to this one
      tr(1,1) = xx(2) - xx(1)
      tr(1,2) = xx(3) - xx(1)
      tr(1,3) = xx(1)
      tr(2,1) = yy(2) - yy(1)
      tr(2,2) = yy(3) - yy(1)
      tr(2,3) = yy(1)

*     compute the cot of each angle: cot A = (b^2+c^2-a^2)/(4*area)
      cot12 = (e23*e23 + e31*e31 - e12*e12) / (4. * area)
      cot23 = (e31*e31 + e12*e12 - e23*e23) / (4. * area)
      cot31 = (e12*e12 + e23*e23 - e31*e31) / (4. * area)

      cota(1,2) = cot12
      cota(2,1) = cot12
      cota(2,3) = cot23
      cota(3,2) = cot23
      cota(3,1) = cot31
      cota(1,3) = cot31

c      print '("elt_stf_har: area,cot: ", 4g12.6)',
c     $     area, cot12, cot23, cot31


*      print *, "check 3"
      do ii=1,ndl
         do jj=1,ndl
            stiff(ii,jj) = 0.0
         enddo
      enddo

      do ii=1,ndl


         do 10 jj=1,ndl

            if (jj .eq. ii) goto 10

C            print *, 'elt_stf_har: check b,ii,jj=', ii,jj

            xii = xx(ii)
            yii = yy(ii)
            xjj = xx(jj)
            yjj = yy(jj)

            xxx = (xii + xjj) * 0.5d00
            yyy = (yii + yjj) * 0.5d00
C
C            print *, 'elt_stf_har: xii,yii,xjj,yjj=',xii,yii,xjj,yjj
C
            a0e = (alpha(xii,yii) + alpha(xjj,yjj)) * 0.5d00
            btau = beta_x(xxx,yyy) * (xjj-xii)
     $           + beta_y(xxx,yyy) * (yjj-yii)

C            print *,'elt_stf_har: a0e, btau=',a0e,btau

            stiff(ii,ii) = stiff(ii,ii) +
     &           bernoulli(-btau/a0e) * a0e * cota(ii,jj) *0.5d00

            stiff(ii,jj) = stiff(ii,jj) -
     &           bernoulli(btau/a0e) * a0e * cota(ii,jj)* 0.5d00

C            print *, ' elt_stf_har: ', ii,jj,stiff(ii,jj)

 10      enddo

         nt = 3
         rhs(ii) = dbl_int(nt, ii, jj, tr, area,xx,yy)
C
c         print *, ' elt_stf_har: rhs ', rhs(ii)
C

 60   enddo

      return
      end

