
C====================================================================
      subroutine onpnam(name1)
C====================================================================
      character*100 name1
      l1 = len('//home/euler2/ltz/grids/')
      l2 = index(name1,' ')
      name1(l1+1:l1+l2)=name1(1:l2-1)
      name1(1:l1)='/home/euler2/ltz/grids/'
      write(*,*) ' INPUT file name: ', 
     >     name1(1:index(name1,' ')-1)
      return
      end
C====================================================================
      subroutine inpnam(name1,name2)
C====================================================================
      character*100 name1,name2
      l1 = len('/home/euler2/ltz/grids/')
      ll1 = len('/home/euler2/ltz/adjncy/')
      l2 = index(name1,' ')
      name1(l1+1:l1+l2-1)=name1(1:l2-1)
      name2(ll1+1:ll1+l2-1)=name1(1:l2-1)
      name2(ll1+l2:ll1+l2+3)='.adj'
      name1(1:l1)='/home/euler2/ltz/grids/'
      name2(1:ll1)='/home/euler2/ltz/adjncy/'
      do i = 1 , 100
c         write(*,*) ichar(name1(i:i)),ichar(' '),ichar(name2(i:i))
         if(ichar(name2(i:i)) .eq. 0) name2(i:i) = ' '
c         write(*,*) ichar(name1(i:i)),ichar(' '),ichar(name2(i:i))
c         read(*,*)
      end do
      m10 =     index(name1,' ')
      m10 =     index(name2,' ')
      write(*,*) ' INPUT and ADJ file names: ', 
     >     name1(1:index(name1,' ')-1),
     >     name2(1:index(name2,' ')-1)
      return
      end
C====================================================================
      subroutine vanek_data(name1)
C====================================================================
      character*100 name1,name2
      l1 = len('/home/euler2/ltz/vvanek/matrix.dat/')
      l2 = index(name1,' ')
      name1(l1+1:l1+l2-1)=name1(1:l2-1)
      name2(ll1+1:ll1+l2-1)=name1(1:l2-1)
      name2(ll1+l2:ll1+l2+3)='.adj'
      name1(1:l1)='/home/euler2/ltz/grids/'
      name2(1:ll1)='/home/euler2/ltz/adjncy/'
      write(*,*) ' INPUT and ADJ file names: ', 
     >     name1(1:index(name1,' ')-1),
     >     name2(1:index(name2,' ')-1)
      return
      end
C====================================================================
      subroutine outnam(name1)
C====================================================================
      character*100 name1
      l1 = len('/home/euler2/ltz/results/')
      l2 = index(name1,' ')
      name1(l1+1:l1+l2)=name1(1:l2-1)
      name1(1:l1)='/home/euler2/ltz/results/'
      write(*,*) ' OUTPUT file name: ', 
     >     name1(1:index(name1,' ')-1)
      return
      end
C
C============================================================
      subroutine aat_nokept (ia,ja,an,n,m)
C============================================================
      implicit real*8 (a-h,o-z)
      pointer (miat,iat),(mjat,jat),(mant,ant)
      dimension ia(1),ja(1),iat(*),jat(*),an(1),ant(*)
C--------------------------------------------------------------------
C...  TRANSPOSITION OF THE SPARSE MATRIX.=== matrix is not kept.
C--------------------------------------------------------------------
C....
C....
      mh = m + 1
      nh = n + 1
C
      call mem_allocate_integer(mh,miat)
C
      iab = ia(nh) - 1
C
      call mem_allocate_integer(iab,mjat)
C
      call mem_allocate_real(iab,mant)
C
      do i = 2,mh
         iat(i) = 0
      end do
C
      iab = ia(nh) - 1
C
      do i = 1 , iab
         j = ja(i) + 2
         if(j .le. mh) iat(j) = iat(j) + 1
      end do
C
      iat(1) = 1
      iat(2) = 1
      if(m .ne. 1) then
         do i = 3 , mh
            iat(i) = iat(i) + iat(i-1)
         end do
      end if
C
      do i = 1 , n
C
         iaa = ia(i)
         iab = ia(i+1) - 1
c
         if(iab .ge. iaa) then
c
            do jp = iaa , iab
               j = ja(jp) + 1
               k = iat(j) 
               jat(k) = i
               ant(k) = an(jp)
               iat(j) = k + 1
            end do
c
         end if
c
      end do
c

      if( m .gt. n) then 
         write(*,*)    
     >        ' some error might be occured here: n = ', n, 
     >        ' is less than m = ', m
      else
         call icopyv(iat,ia,mh)
         nnz = ia(mh)-1
         call icopyv(jat,ja,nnz)
         call copyv(ant,an,nnz)
      end if
      call free(miat)
      call free(mjat)
      call free(mant)
      return
      end

C====================================================================
      subroutine iit_nokept(ia,ja,n,m)
C====================================================================
      pointer (miat,iat), (mjat,jat)
      integer*4 ia(1),ja(1),iat(*),jat(*),n,m
C--------------------------------------------------------------------
C... TRANSPOSITION OF GRAPH (OR THE MATRIX --- SYMBOLICALLY)
C--------------------------------------------------------------------
c
c   n+1 is the dimension of ia
c
      mh = m + 1
C
      call mem_allocate_integer(mh,miat)
C
      nh = n + 1
c
      do i = 2,mh
         iat(i) = 0
      end do
c
      iab = ia(nh) - 1
C
      call mem_allocate_integer(iab,mjat)
C
      iab = ia(nh) - 1
c
      do i = 1 , iab
         j = ja(i) + 2
         if(j .le. mh) iat(j) = iat(j) + 1
      end do
c
      iat(1) = 1
      iat(2) = 1
      if(m .ne. 1) then
         do i = 3 , mh
            iat(i) = iat(i) + iat(i-1)
         end do
      end if
c
      do i = 1 , n
c
         iaa = ia(i)
         iab = ia(i+1) - 1
c
         if(iab .ge. iaa) then
c
            do jp = iaa , iab
               j = ja(jp) + 1
               k = iat(j) 
               jat(k) = i
               iat(j) = k + 1
            end do
c
         end if
c
      end do
c

      if( m .gt. n) then 
         write(*,*)    
     >        ' some error might be occured here: n = ', n, 
     >        ' is less than m = ', m
      else
         call icopyv(iat,ia,mh)
         nnz = ia(mh+1)-1
         call icopyv(jat,ja,nnz)
      end if
      call free(miat)
      call free(mjat)
      return
      end
C====================================================================
      subroutine abyvg_sym(ia,ja,an,b,n,c)
C====================================================================
      implicit real*8 (a-h,o-z)
      pointer(miaw,iaw),(mjaw,jaw),(manw,anw),(mdw,dw)
      dimension ia(1),ja(1),an(1),b(n),c(n)
      dimension iaw(*),jaw(*),anw(*),dw(*)
C--------------------------------------------------------------------
      call mem_allocate_integer(n+1,miaw)
C
      call mem_allocate_real(n,mdw)
C
      nwka = ia(n+1)-1
C
      nwkaw = (nwka - n) / 2
C
      call mem_allocate_integer(nwkaw,mjaw)
      call mem_allocate_real(nwkaw,manw)
C
      jp = 1
      iab = ia(1)-1
      iaw(1) = 1
      do  i = 1 , n
         dw(i) = 0.0d00
         iaa = iab + 1
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do ip = iaa,iab
               j = ja(ip)
               if(i-j .eq. 0) then
                  dw(i) = an(ip)
                  go to 100
               end if
               if(j-i .lt. 0) go to 100
               jaw(jp) = j
               anw(jp) = an(ip)
               jp = jp + 1
 100           continue
            end do
         end if
         iaw(i+1) = jp
      end do
C
      call sabyv(iaw,jaw,anw,dw,b,n,c)
C
      call my_free(4*n+4,miaw)
      call my_free(4*nwkaw,mjaw)
      call my_free(8*n,mdw)
      call my_free(8*nwkaw,manw)
C     
      return
      end
C====================================================================
      subroutine abyva0(ia,ja,an,b,n,c)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),c(n)
C--------------------------------------------------------------------
C... PRODUCT--- GENERAL block, SPARSE MATRIX BY VECTOR-Vector
C...c = a - A * b.
C--------------------------------------------------------------------
      do i = 1 , n
         u = 0.0d00
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do k = iaa,iab
               u = u - an(k) * b(ja(k))
            end do
         end if
         c(i) = u
      end do
      return
      end
C====================================================================
      subroutine rfabc_non_symmetric(f,ia,ja,an,b,n,c)
c====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(n),f(n)
c--------------------------------------------------------------------
c...    c = f-A*b.
c--------------------------------------------------------------------
      call abyvg(ia,ja,an,b,n,c)
ccccccc      call abyvg_sym(ia,ja,an,b,n,c)
C
      call vuminv(f,c,n)
C
      return
      end
C====================================================================
      real*8 function omsor(n)
C====================================================================
      implicit real*8(a-h,o-z)
      ddn = dsqrt(dble(n))+1.d00
      om = dsin(3.141592653589793/ddn)
      omsor = 2/(1.d00+om)
      return
      end
c
C====================================================================
      subroutine symass(ie,je,iet,jet,ia,ja,n,nwk)
C====================================================================
      dimension ia(1),ie(1),je(1),iet(1),jet(1),ja(1)
C--------------------------------------------------------------------
C...  ASSEMBLING OF A SYMMETRIC SPARSE MATRIX--- FEM
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C...  SYMBOLIC PART
C--------------------------------------------------------------------
      jp = 1
      nm = n - 1
      do  i = 1 , nm
         jpi = jp
         if(ia(i) .ne. n) then
            ieta = iet(i)
            ietb = iet(i+1) - 1
            do  ip = ieta,ietb
               j = jet(ip)
               iea = ie(j)
               ieb = ie(j+1) - 1
               do  kp = iea,ieb
                  k = je(kp)
                  if(k .gt. i) then
                     if(ia(k) .lt. i) then
                        ja(jp) = k
                        jp = jp + 1
                        ia(k) = i
                     end if
                  end if
               end do
            end do
         end if 
         ia(i) = jpi
      end do
      ia(n) = jp
      ia(n+1) = jp
      nwk = ia(n+1)
      return
      end
C====================================================================
      subroutine rcuruu(ia,ja,an,ad,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension an(1),ad(1),ia(1),ja(1)
C--------------------------------------------------------------------
C...  TRANSFORMATION FROM GENERAL UNORDERED MATRIX FORMAT
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C...  INTO UPPER TRIANGULAR ROW-WISE FORMAT.
C--------------------------------------------------------------------
      jp = 1
      iab = ia(1)-1
      do  i = 1 , n
         ad(i) = 0.0d00
         iaa = iab + 1
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do ip = iaa,iab
               j = ja(ip)
               if(i-j .eq. 0) then
                  ad(i) = an(ip)
                  go to 100
               end if
               if(j-i .lt. 0) go to 100
               ja(jp) = j
               an(jp) = an(ip)
               jp = jp + 1
 100           continue
            end do
         end if
         ia(i+1) = jp
      end do
      return
      end
C====================================================================
      subroutine rcuruu_envelope_size(ia,ja,an,ad,n,isize)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension an(1),ad(1),ia(1),ja(1)
C--------------------------------------------------------------------
C...  TRANSFORMATION FROM GENERAL UNORDERED MATRIX FORMAT
C...    the envelope size is computed ...
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C...  INTO UPPER TRIANGULAR ROW-WISE FORMAT.
C--------------------------------------------------------------------
      jp = 1
      isize = 2
      iab = ia(1)-1
      do  i = 1 , n
         ad(i) = 0.0d00
         iaa = iab + 1
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            izz = 0
            do ip = iaa,iab
               j = ja(ip)
               if(i-j .eq. 0) then
                  ad(i) = an(ip)
                  go to 200
               end if
               if(j-i .lt. 0) go to 100
               ja(jp) = j
               an(jp) = an(ip)
               jp = jp + 1
               go to 200
 100           continue
               izz =  max0(i-j,izz)
 200           continue
            end do
            isize = isize + izz + 1
         end if
         ia(i+1) = jp
      end do
      return
      end
C====================================================================
      subroutine shrink(nxadj,nadj,n)
C====================================================================
      integer nxadj(1),nadj(1),n
      kend = nxadj(n+1) - 1
         istrt = nxadj(1)
      do k = 1 , n
         iend = nxadj(k+1) - 1
         klngth= iend - istrt + 1
         do j = istrt,iend
            if(nadj(j) .eq. 0) klngth = klngth - 1
         end do
         nxadj(k+1) = nxadj(k) + klngth
         istrt = iend + 1
      end do
      l = 0
      do k = 1 , kend
         if(nadj(k) .ne. 0) then
            l = l + 1
            nadj(l) = nadj(k)
         end if
      end do
      return
      end
C====================================================================
      subroutine shreal(nxadj,nadj,a,n)
      implicit real*8(a-h,o-z)
C====================================================================
      dimension nxadj(1),nadj(1),a(1)
      kend = nxadj(n+1) - 1
         istrt = nxadj(1)
      do k = 1 , n
         iend = nxadj(k+1) - 1
         klngth= iend - istrt + 1
         do j = istrt,iend
            if(nadj(j) .lt. 1 .or. nadj(j) .gt. n) 
     >           klngth = klngth - 1
         end do
         nxadj(k+1) = nxadj(k) + klngth
         istrt = iend + 1
      end do
c
      l = 0

      do k = 1 , kend
         if(nadj(k) .gt. 0 .and. nadj(k) .le. n) then
            l = l + 1
            nadj(l) = nadj(k)
            a(l) = a(k)
         end if
      end do
      return
      end
C====================================================================
      subroutine shr_matrix(nxadj,nadj,a,n)
      implicit real*8(a-h,o-z)
C====================================================================
      dimension nxadj(1),nadj(1),a(1)
      do k = 1 , n
         do jk = nxadj(k),nxadj(k+1)-1
c
            if(dabs(a(jk)) .lt. 1.d-32) nadj(jk) = 0
         end do
      end do
      kend = nxadj(n+1) - 1
         istrt = nxadj(1)
      do k = 1 , n
         iend = nxadj(k+1) - 1
         klngth= iend - istrt + 1
         do j = istrt,iend
            if(nadj(j) .lt. 1 .or. nadj(j) .gt. n) 
     >           klngth = klngth - 1
         end do
         nxadj(k+1) = nxadj(k) + klngth
         istrt = iend + 1
      end do
c
      l = 0

      do k = 1 , kend
         if(nadj(k) .gt. 0 .and. nadj(k) .le. n) then
            l = l + 1
            nadj(l) = nadj(k)
            a(l) = a(k)
         end if
      end do
      return
      end
C====================================================================
      subroutine shreal_new(nxadj,nadj,a,n,n_max)
      implicit real*8(a-h,o-z)
C====================================================================
      dimension nxadj(1),nadj(1),a(1)
      kend = nxadj(n+1) - 1
         istrt = nxadj(1)
      do k = 1 , n
         iend = nxadj(k+1) - 1
         klngth= iend - istrt + 1
         do j = istrt,iend
            if(nadj(j) .lt. 1 .or. nadj(j) .gt. n_max) 
     >           klngth = klngth - 1
         end do
         nxadj(k+1) = nxadj(k) + klngth
         istrt = iend + 1
      end do
c
      l = 0

      do k = 1 , kend
         if(nadj(k) .gt. 0 .and. nadj(k) .le. n_max) then
            l = l + 1
            nadj(l) = nadj(k)
            a(l) = a(k)
         end if
      end do
      return
      end
C====================================================================
      subroutine dirsim(xa,row,nbb,n,ad)
C====================================================================
      integer*4 xa(1),row(1),nbb(1)
      real*8 ad(1)
c
      do i1 = 1 , n
         if(nbb(i1) .ne. 0) then
            ad(i1) = 1.D00
            call inullv(row(xa(i1)),xa(i1+1)-xa(i1))
         else
            do k = xa(i1),xa(i1+1)-1
               if(nbb(row(k)) .ne. 0) row(k) = 0
            end do
         end if
      end do
      return
      end
C====================================================================
      subroutine ircuru(ia,ja,n)
C--------------------------------------------------------------------
C...  TRANSFORMATION FROM GENERAL UNORDERED MATRIX FORMAT
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C...  INTO UPPER TRIANGULAR ROW-WISE FORMAT.
C...   SYMBOLIC SECTION, I.E. ONLY GRAPH STRUCTURE.
C--------------------------------------------------------------------
C====================================================================
      integer ia(1),ja(1)
      jp = 1
      iab = ia(1)-1
      do  i = 1 , n
         iaa = iab + 1
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do ip = iaa,iab
               j = ja(ip)
               if(j-i .le. 0) go to 100
               ja(jp) = j
               jp = jp + 1
 100           continue
            end do
         end if
         ia(i+1) = jp
      end do
      return
      end
C====================================================================
      subroutine dbya(ia,an,d,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),an(1),d(1)
C--------------------------------------------------------------------
C...  MULTIPLICATION DIAGONAL MATRIX BY GENERAL ONE
C--------------------------------------------------------------------
      do i = 1 , n
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .ge. iaa) then
            dd = d(i)
            do j = iaa,iab
               an(j) = an(j) * dd
            end do
         end if
      end do
      return
      end
      double precision function ddot_kz(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot_kz = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot_kz = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot_kz = dtemp
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
C====================================================================
