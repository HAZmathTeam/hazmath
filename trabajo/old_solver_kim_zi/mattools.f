C=====================================================================
      subroutine aplbs(ia,ja,ib,jb,n,m,ic,jc,ix)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),ib(1),jb(1),ic(1),jc(1),ix(m)
C--------------------------------------------------------------------
C...  SYMBOLLIC ADDITION OF TWO GENERAL SPARSE MATRICIES: A + B
C--------------------------------------------------------------------
      ip = 1
      do i = 1, m
         ix(i) = 0
      end do
C
      do i = 1, n
         ic(i) = ip
         iaa = ia(i)
         iab = ia(i+1) - 1
C
         if(iab .ge. iaa) then
            do jp = iaa, iab
               j = ja(jp)
               jc(ip) = j
               ip = ip + 1
               ix(j) = i
            end do
         end if
         iba = ib(i)
         ibb = ib(i+1) - 1
C
         if(ibb .ge. iba) then
            do jp = iba, ibb
               j = jb(jp)
               if(ix(j) .ne. i) then
                  jc(ip) = j
                  ip = ip + 1
               end if
            end do
         end if
      end do
C
      ic(n+1) = ip
      return
      end
C=====================================================================
      subroutine aplb(ia,ja,ib,jb,n,m,ic,jc,an,bn,cn,x)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),ib(1),jb(1),ic(1),jc(1)
      dimension an(1),bn(1),cn(1),x(m)
C--------------------------------------------------------------------
C...  ADDITION OF TWO SPARSE MATRICIES: A + B
C--------------------------------------------------------------------
      do i = 1, n
         ih = i + 1
         ica = ic(i)
         icb = ic(ih) - 1
         if(icb .ge. ica) then
            do ip = ica, icb
               x(jc(ip)) = 0.0d00
            end do
            iaa = ia(i)
            iab = ia(ih) - 1
            if(iab .ge. iaa) then
               do ip = iaa, iab
                  x(ja(ip)) = an(ip)
               end do
            end if
            iba = ib(i)
            ibb = ib(ih)-1
            if(ibb .ge. iba) then
               do ip = iba, ibb
                  j = jb(ip)
                  x(j) = x(j) + bn(ip)
               end do
            end if
            do ip = ica, icb
               cn(ip) = x(jc(ip))
            end do
         end if
      end do
      return
      end
C=====================================================================
      subroutine abyvg(ia,ja,an,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(n)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL SPARSE MATRIX BY VECTOR: c = A*b.
C--------------------------------------------------------------------
      do i = 1, n
         u = 0.0d00
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do k = iaa, iab
               u = u + an(k) * b(ja(k))
            end do
         end if
         c(i) = u
      end do
      return
      end
C=====================================================================
      subroutine abyv_unif(nx,ny,n,x,y)
C=====================================================================
      implicit none
      integer nx,ny,n,i,j,i0,i1,i2,i3,i4
      real*8  x(1),y(1)
      integer*4 nomxy
C--------------------------------------------------------------------
C...  Performs the action of Y = A*X on the structured triangular mesh.
C...  Note that nodes are numbered y-direction first from bottom to
C...  top and increased in x-direction from left to right. 
C...
C...  Parameters:
C...    NX,NY   - grid points in x- and y-directions
C...    N       - N = NX*NY, the number of unknowns
C...    X,Y     - input and output, Y = A*X.
C--------------------------------------------------------------------
C...  Initialize Y to be zero.
C
      call nullv(y,n)
C      
      do 50 i = 2, nx - 1
         do 40 j = 2, ny - 1
            i0 = nomxy(i,  j,    nx,ny)
            i1 = nomxy(i+1,j,    nx,ny)  
            i2 = nomxy(i,  j+1,  nx,ny)  
            i3 = nomxy(i,  j-1,  nx,ny)  
            i4 = nomxy(i-1,j  ,  nx,ny)  
            y(i0) = 4.0d0*x(i0) - x(i4) - x(i1) - x(i3) - x(i2) 
 40      continue
 50   continue
C
      return
      end
C======================================================================
      subroutine abyvcs(ia,ja,an,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(n)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, Vector + SPARSE MATRIX BY VECTOR
C...  c = c + A * b.
C--------------------------------------------------------------------
      do i = 1, n
         u = c(i)
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do k = iaa, iab
               u = u + an(k) * b(ja(k))
            end do
         end if
         c(i) = u
      end do
      return
      end
C=====================================================================
      subroutine abyvcm(ia,ja,an,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(n)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, Vector - SPARSE MATRIX BY VECTOR
C...  c = c - A * b.
C--------------------------------------------------------------------
      do i = 1, n
         u = c(i)
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do k = iaa, iab
               u = u - an(k) * b(ja(k))
            end do
         end if
         c(i) = u
      end do
      return
      end
C=====================================================================
      subroutine abyvas(a,ia,ja,an,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),c(n),a(n)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, Vector + SPARSE MATRIX BY VECTOR
C...  c = a + A * b.
C--------------------------------------------------------------------
      do i = 1, n
         u = a(i)
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do k = iaa, iab
               u = u + an(k) * b(ja(k))
            end do
         end if
         c(i) = u
      end do
      return
      end
C=====================================================================
      subroutine abyvam(a,ia,ja,an,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),c(n),a(n)
C---------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, Vector - SPARSE MATRIX BY VECTOR
C...  c = a - A * b.  It is the same as the subroutine RFABC_NS if 
C...  every entry of the matrix A is stored in the array AN.
C---------------------------------------------------------------------
      do i = 1, n
         u = a(i)
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do k = iaa, iab
               u = u - an(k) * b(ja(k))
            end do
         end if
         c(i) = u
      end do
      return
      end
C=====================================================================
      subroutine abyva0(ia,ja,an,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),c(n)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, 0 - SPARSE MATRIX BY VECTOR
C...  c = 0 - A * b.
C--------------------------------------------------------------------
      do i = 1, n
         u = 0.0d00
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do k = iaa, iab
               u = u - an(k) * b(ja(k))
            end do
         end if
         c(i) = u
      end do
      return
      end
C=====================================================================
      subroutine vbya(ia,ja,an,b,n,m,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(m)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, Vector^t by SPARSE MATRIX
C...  c = b^t*A.
C--------------------------------------------------------------------
      do i = 1, m
         c(i) = 0.0d00
      end do
      do i = 1, n
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            z = b(i)
            do k = iaa, iab
               j = ja(k)
               c(j) = c(j) + an(k)*z
            end do
         end if
      end do
      return
      end
C=====================================================================
      subroutine vbyaas(a,ia,ja,an,b,n,m,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(m),a(1)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, Vector + Vector^t by SPARSE MATRIX
C...  c = a + b^t*A.
C--------------------------------------------------------------------
      do i = 1, m
         c(i) = a(i)
      end do
      do i = 1, n
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            z = b(i)
            do k = iaa, iab
               j = ja(k)
               c(j) = c(j) + an(k)*z
            end do
         end if
      end do
      return
      end
C=====================================================================
      subroutine vbyacs(ia,ja,an,b,n,m,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(m),c(m)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, Vector + Vector^t by SPARSE MATRIX
C...  c = c + b^t*A.
C--------------------------------------------------------------------
      do i = 1, n
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            z = b(i)
            do k = iaa, iab
               j = ja(k)
               c(j) = c(j) + an(k)*z
            end do
         end if
      end do
      return
      end
C=====================================================================
      subroutine vbyaam(a,ia,ja,an,b,n,m,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(m),a(1)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, Vector - Vector^t by SPARSE MATRIX
C...  c = a - b^t*A.
C--------------------------------------------------------------------
      do i = 1,  m
         c(i) = a(i)
      end do
      do i = 1, n
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            z = b(i)
            do k = iaa, iab
               j = ja(k)
               c(j) = c(j) - an(k)*z
            end do
         end if
      end do
      return
      end
C=====================================================================
      subroutine vbyaa0(ia,ja,an,b,n,m,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(m)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, 0 - Vector^t by SPARSE MATRIX
C...  c = 0 - b^t*A.
C--------------------------------------------------------------------
      do i = 1, m
         c(i) = 0.0d00
      end do
      do i = 1, n
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            z = b(i)
            do k = iaa, iab
               j = ja(k)
               c(j) = c(j) - an(k)*z
            end do
         end if
      end do
      return
      end
C=====================================================================
      subroutine vbyacm(ia,ja,an,b,n,m,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(m),c(n)
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL block, Vector - Vector^t by SPARSE MATRIX
C...  c = c - b^t*A.
C--------------------------------------------------------------------
      do i = 1, n
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            z = b(i)
            do k = iaa, iab
               j = ja(k)
               c(j) = c(j) - an(k)*z
            end do
         end if
      end do
      return
      end
C=====================================================================
      subroutine dbyv(ad,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension b(n),c(n),ad(n)
C--------------------------------------------------------------------
C...  c = D*b, where D is a diagonal matrix.
C--------------------------------------------------------------------
      do i = 1, n
         c(i) = ad(i) * b(i)
      end do
      return
      end
C=====================================================================
      subroutine dbyvsc(ad,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension b(n),c(n),ad(n)
C--------------------------------------------------------------------
C...  c = c + D*b, where D is a diagonal matrix.
C--------------------------------------------------------------------
      do i = 1, n
         c(i) = c(i) + ad(i) * b(i)
      end do
      return
      end
C=====================================================================
      subroutine dbyvsa(a,ad,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension b(n),c(n),ad(n),a(n)
C--------------------------------------------------------------------
C...  c = a + D*b, where D is a diagonal matrix.
C--------------------------------------------------------------------
      do i = 1, n
         c(i) = a(i) + ad(i) * b(i)
      end do
      return
      end
C=====================================================================
      subroutine dbyvcm(ad,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension b(n),c(n),ad(n)
C--------------------------------------------------------------------
C...  c = c - D*b, where D is a diagonal matrix.
C--------------------------------------------------------------------
      do i = 1, n
         c(i) = c(i) - ad(i) * b(i)
      end do
      return
      end
C=====================================================================
      subroutine dbyvam(a,ad,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension b(n),c(n),ad(n),a(n)
C--------------------------------------------------------------------
C...  c = a - D*b, where D is a diagonal matrix.
C--------------------------------------------------------------------
      do i = 1, n
         c(i) = a(i) - ad(i) * b(i)
      end do
      return
      end
C=====================================================================
      subroutine dbyvm0(ad,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension b(n),c(n),ad(n)
C--------------------------------------------------------------------
C...  c = 0 - D*b, where D is a diagonal matrix.
C--------------------------------------------------------------------
      do i = 1, n
         c(i) = - ad(i) * b(i)
      end do
      return
      end
C=====================================================================
      subroutine sabyv(ia,ja,an,ad,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(n),ad(n)
C--------------------------------------------------------------------
C...  PRODUCT--- SYMMETRIC SPARSE MATRIX BY VECTOR. c = A*b.
C--------------------------------------------------------------------
      call dbyv(ad,b,n,c)
C
      call vbyacs(ia,ja,an,b,n,n,c)
C 
      call abyvcs(ia,ja,an,b,n,c)
C 
      return
      end
C=====================================================================
      subroutine rfabc(f,ia,ja,an,ad,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(n),ad(n),f(n)
C--------------------------------------------------------------------
C...  c = f - A*b.
C--------------------------------------------------------------------
      call sabyv(ia,ja,an,ad,b,n,c)
C
      call vuminv(f,c,n)
C
      return
      end
C=====================================================================
      subroutine rfabc_ns(f,ia,ja,an,b,n,c)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),b(n),c(n),f(n)
C--------------------------------------------------------------------
C...  c = f - A*b. It is the same as the subroutine ABYVAM if every
C...  entry of the matrix A is stored in the array AN.
C--------------------------------------------------------------------
      call abyvg(ia,ja,an,b,n,c)
C
      call vuminv(f,c,n)
C
      return
      end
C=====================================================================
      subroutine abybs(ia,ja,ib,jb,np,nq,nr,ic,jc,ix)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),ib(1),jb(1),ic(1),jc(1),ix(nr)
C--------------------------------------------------------------------
C...  SYMBOLIC MULTIPLICATION OF TWO GENERAL SPARSE MATRICIES
C...    np = number of rows of a
C...    nq = number of columns of a
C...    nr = number of columns of b
C--------------------------------------------------------------------
      ip = 1
      do i = 1, nr
         ix(i) = 0
      end do
C
      do i = 1, np
         ic(i) = ip
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge. iaa) then
            do jp = iaa,iab
               j = ja(jp)
               iba = ib(j)
               ibb = ib(j+1) - 1
               if(ibb .ge. iba) then
                  do kp = iba,ibb
                     k = jb(kp) 
                     if(ix(k) .ne. i) then
                        jc(ip) = k
                        ip = ip + 1
                        ix(k) = i
                     end if
                  end do
               end if
            end do
         end if
      end do
      ic(np+1) = ip
      return
      end
C=====================================================================
      subroutine abybs_new(ia,ja,ib,jb,np,nq,nr,ic,jc,ix,mem_chk)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),ib(1),jb(1),ic(1),jc(1),ix(nr)
      logical nflag
C--------------------------------------------------------------------
C...  SYMBOLIC MULTIPLICATION OF TWO GENERAL SPARSE MATRICIES
C...    np = number of rows of a
C...    nq = number of columns of a
C...    nr = number of columns of b
C--------------------------------------------------------------------
      ip = 1
      do i = 1 , nr
         ix(i) = 0
      end do
C     
      nflag = .true.
C
      do i = 1 , np
         ic(i) = ip
C...     check memory
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .ge.iaa) then
            do jp = iaa,iab
               j = ja(jp)
               iba = ib(j)
               ibb = ib(j+1) - 1
               if(ibb .ge. iba) then
                  do kp = iba,ibb
                     k = jb(kp) 
                     if(ix(k) .ne. i) then
                        if(ip .le. mem_chk) then
                           jc(ip) = k
                           ipp = ip
                        else
                           if(nflag) nflag = .false.
                        end if
                        ip = ip + 1
                        ix(k) = i
                     end if
                  end do
               end if
            end do
         end if
      end do
      ic(np+1) = ip
cc    if(nflag) then
cc       return
cc    else
cc       write(*,*) ' going out urgently ',ic(np+1),jc(ipp)
cc       read(*,*)
cc    end if
      return
      end
C=====================================================================
      subroutine abyb(ia,ja,ib,jb,np,ic,jc,an,bn,cn,x,nq)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),ib(1),jb(1),ic(1),jc(1)
      dimension an(1),bn(1),cn(1),x(nq)
C--------------------------------------------------------------------
C...  MULTIPLICATION OF TWO GENERAL SPARSE MATRICIES: C = A*B
C--------------------------------------------------------------------
      do i  = 1 , np
         ica = ic(i) 
         icb = ic(i+1) - 1
         if(icb .ge. ica) then
            do j = ica,icb
               x(jc(j)) = 0.0d00
            end do
            iaa = ia(i)
            iab = ia(i+1) - 1
            do jp = iaa,iab
               j = ja(jp)
               a = an(jp)
               iba = ib(j)
               ibb = ib(j+1) - 1
               if(ibb .ge. iba) then
                  do kp = iba,ibb
                     k = jb(kp)
                     x(k) = x(k) + a*bn(kp)
                  end do
               end if
            end do
            do j = ica,icb
               cn(j) = x(jc(j))
            end do
         end if
      end do
      return
      end
C=====================================================================
      subroutine sfactr(ia,ja,n,iu,ju,ip,nwku)
C=====================================================================
      dimension ia(1),ja(1),iu(1),ju(1),ip(n)
C--------------------------------------------------------------------
C...  SYMBOLIC FACTORIZATION OF A SYMMETRIC SPARSE MATRIX.
C--------------------------------------------------------------------
      nm = n - 1
      nh = n + 1
      do i = 1 , n
         iu(i) = 0
         ip(i) = 0
      end do
C
      jp = 1
      do i = 1 , nm
         jpi = jp
         jpp = n + jp - i
         min = nh
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .ge. iaa) then
            do j = iaa,iab
               jj = ja(j)
               ju(jp) = jj
               jp = jp + 1
               if(jj .lt. min) min = jj
               iu(jj) = i
            end do
         end if
          last = ip(i)
         if(last .eq. 0) go to 60
         l = last
 40      l = ip(l)
         lh = l + 1
         iua = iu(l)
         iub = iu(lh) - 1
         if(lh .eq. i) iub = jpi - 1
         iu(i) = i
         do j = iua,iub
            jj = ju(j)
            if(iu(jj) .ne. i) then
               ju(jp) = jj
               jp = jp + 1
               iu(jj) = i
               if(jj .lt. min)min = jj
            end if
         end do
         if(jp .eq. jpp) go to 70
         if(l.ne. last) go to 40
 60      if (min .eq. nh) go to 90
 70      l = ip(min)
         if(l .eq. 0) go to 80
         ip(i) = ip(l)
         ip(l) = i
         go to 90
 80      ip(min) = i
         ip(i) = i
 90      continue
         iu(i) = jpi
      end do
      iu(n) = jp
      iu(nh) = jp
      nwku = iu(n+1)
      return
      end
C=====================================================================
      subroutine sfactr_new(ia,ja,n,iu,ju,ip,nwku,mem_chk)
C=====================================================================
      dimension ia(1),ja(1),iu(1),ju(1),ip(n)
C--------------------------------------------------------------------
C... SYMBOLIC FACTORIZATION OF A SYMMETRIC SPARSE MATRIX.
C--------------------------------------------------------------------
      nm = n - 1
      nh = n + 1
      iu(nh) = 0
      do i = 1 , n
         iu(i) = 0
         ip(i) = 0
      end do
C
      jp = 1
      do i = 1 , nm
cc         if(mod(i,10) .eq. 1.and. nm .lt. 7000) write(*,*) i
         jpi = jp
         jpp = n + jp - i
         min = nh
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .ge. iaa) then
            do j = iaa,iab
               jj = ja(j)
C...           memory check, 
               if(jp .le. mem_chk) then
                  ju(jp) = jj
               else
                  iu(n+1) = mem_chk+3
                  return
               end if
               jp = jp + 1
               if(jj .lt. min) min = jj
               iu(jj) = i
            end do
         end if
         last = ip(i)
         if(last .eq. 0) go to 60
         l = last
 40      l = ip(l)
         lh = l + 1
         iua = iu(l)
         iub = iu(lh) - 1
         if(lh .eq. i) iub = jpi - 1
         iu(i) = i
         do j = iua,iub
C...        memory check, 
            if(j .le. mem_chk) then
               jj = ju(j)
            else
               iu(n+1) = mem_chk+3
               return
            end if
            if(iu(jj) .ne. i) then
C...           memory check, 
               if(jp .le. mem_chk) then
                  ju(jp) = jj
               else
                  iu(n+1) = mem_chk+3
                  return
               end if
               jp = jp + 1
               iu(jj) = i
               if(jj .lt. min) min = jj
            end if
         end do
         if(jp .eq. jpp) go to 70
         if(l.ne. last) go to 40
 60      if (min .eq. nh) go to 90
 70      l = ip(min)
         if(l .eq. 0) go to 80
         ip(i) = ip(l)
         ip(l) = i
         go to 90
 80      ip(min) = i
         ip(i) = i
 90      continue
         iu(i) = jpi
      end do
      iu(n) = jp
      iu(nh) = jp
      nwku = iu(n+1)
      return
      end
C=====================================================================
      subroutine factor(ia,ja,n,iu,ju,ip,iup,an,ad,un,di)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),iu(1),ju(1),ip(1),iup(1)
      dimension an(1),ad(1),di(1),un(1)
C--------------------------------------------------------------------
C...  Factorization of symmetric sparse matrix.
C--------------------------------------------------------------------
      do 10 j = 1 , n
         iup(j) = 0
         ip(j) = 0
 10   continue
C
      do 130 i = 1 , n
cc       if(mod(i,10).eq.1 .and. n.lt.7000) write(*,*) ' factoring',i
         ih = i + 1
         iua = iu(i)
         iub = iu(ih) - 1
         if(iub .lt. iua) go to 40
         do 20 j = iua,iub
            di(ju(j)) = 0.0d00
 20      continue
         iaa = ia(i)
         iab = ia(ih) - 1
         if(iab .lt. iaa) go to 40
         do 30 j = iaa,iab
            di(ja(j)) = an(j)
 30      continue
 40      di(i) = ad(i)
         last = ip(i)
         if(last .eq. 0) go to 90
         ln = ip(last)
 50      l = ln
         ln = ip(l)
         iuc = iup(l)
         iud = iu(l + 1) -1
         um = un(iuc) * di(l)
         do 60 j = iuc,iud
            jj = ju(j)
            di(jj) = di(jj) - un(j) * um
 60      continue
         un(iuc) = um
         iup(l) = iuc + 1
         if(iuc .eq. iud) go to 80
         j = ju(iuc+1)
         jj = ip(j)
         if(jj .eq. 0) go to 70
         ip(l) = ip(jj)
         ip(jj) = l
         go to 80
 70      ip(j) = l
         ip(l) = l
 80      if(l .ne. last) go to 50
C
 90      di(i) = 1.d00/di(i)
         if(iub .lt. iua) go to 120
         do 100 j = iua,iub
            un(j) = di(ju(j))
 100     continue
         j = ju(iua)
         jj = ip(j)
         if(jj .eq. 0) go to 110
         ip(i) = ip(jj)
         ip(jj) = i
         go to 120
 110     ip(j) = i
         ip(i) = i
 120     iup(i) = iua
 130  continue
C
      return
      end
C=====================================================================
      subroutine forbac(iu,ju,un,di,n,x)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension iu(1),ju(1),un(1),x(1),di(1)
C--------------------------------------------------------------------
C...  REDUCTION AND BACK SUBSTITUTION FOR SOLVING THE
C...  SYMMETRIC SYSTEM WITH DECOMPOSED (FACTORIZED MATRIX)
C--------------------------------------------------------------------
      nm = n-1
c     do i = 1 , n
c        x(i) = b(i)
c     end do
      do k = 1 , nm
         iua = iu(k)
         iub = iu(k+1) - 1
         xx = x(k)
         if(iub .ge. iua) then
            do i = iua,iub
               x(ju(i)) = x(ju(i))-un(i)*xx
            end do
         end if
         x(k) = xx*di(k)
      end do
      x(n) = x(n)*di(n)
      k = nm
 50   iua = iu(k)
      iub = iu(k+1) - 1
      if(iub .lt. iua) go to 70
      xx = x(k)
      do i = iua,iub
         xx = xx - un(i) *x(ju(i))
      end do
      x(k) = xx
 70   k = k - 1
      if(k .gt. 0) go to 50
C
      return
      end
C===================================================================== 
