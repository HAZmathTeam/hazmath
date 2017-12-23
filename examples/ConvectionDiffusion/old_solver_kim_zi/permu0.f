C=====================================================================
      subroutine cut_off(n,ia,ja,a,iaw,jaw,
     >     iord,isub,nblk,lwork1,lwork2,lwork3)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension iord(n),ia(*),ja(*),iaw(*),jaw(*),isub(n),
     >     lwork1(n),lwork2(n),lwork3(n)
      dimension a(*)
C---------------------------------------------------------------------
C...  INPUT is ia and ja; output is "iord()" which is the permutation"
C---------------------------------------------------------------------
      tol = 1.d-12
      nnz = ia(n+1)-1
C
      call icopyv(ia,iaw,n+1)
      call icopyv(ja,jaw,nnz)
C
      do 30 i = 2, n
         iaa = ia(i)
         iab = ia(i+1)-1
         do 20 ii = iaa,iab
            j = ja(ii)
            if (i .le. j) go to 20
            iac = ia(j)
            iad = ia(j+1)-1
            do 10 iik = iac, iad
               iii = iik
               if (ja(iii) .eq. i) go to 15
 10         continue
            write(*,*) ' Error j has no i, (i,j) = (',i,j,')'
            stop 
 15         continue
            aa = a(ii)
            bb = a(iii)
            if(i .ne. j .and. aa .gt. tol) aa = 0.0d00
            if(i .ne. j .and. bb .gt. tol) bb = 0.0d00
            if(dabs(aa) .lt. tol .and. dabs(bb) .lt. tol) then
               jaw(ii) = 0
               jaw(iii) = 0
               go to 20
            end if
            call chsize(aa,bb,tol,imin)
            if (imin .eq. 0) go to 20
            if (imin .eq. 1) then
               jaw(ii) = 0
            else
               jaw(iii) = 0
            end if
 20      continue
 30   continue
C
      call shift(iaw,jaw,n)
C
      call depth_first_search_block(n,iaw,jaw,
     >     iord,isub,nblk,lwork1,lwork2,lwork3)
cc      write(*,*) ' Ordering:::', nblk, ' blocks.'
C
      return
      end
C======================================================================
      subroutine chsize(a,b,tol,imin)
C======================================================================
      implicit real*8(a-h,o-z)
C---------------------------------------------------------------------
C...
C---------------------------------------------------------------------
      imin = 0
      if (dabs(a) .gt. 1.d-13 .and. dabs(b) .gt. 1.d-13) then
         ra = dabs(abs(a)/dabs(b)-1.d00)
         if(ra .lt. tol) then
            return
         else
            go to 10
         end if
      end if
C
      if(dabs(a) .gt. 1.d-13 .or. dabs(b) .gt. 1.d-13) go to 10
      return
C
 10   if(dabs(a) .gt. dabs(b)) then
         imin = 2
      else
         imin = 1
      end if
C
      return
      end
C====================================================================
      subroutine shift(nxadj,nadj,n)
C====================================================================
      integer nxadj(*),nadj(*),n
C---------------------------------------------------------------------
C...
C---------------------------------------------------------------------
      kend = nxadj(n+1) - 1
      istrt = nxadj(1)
      do k = 1, n
         iend = nxadj(k+1) - 1
         klngth = iend - istrt + 1
         do j = istrt,iend
            if(nadj(j) .eq. 0) klngth = klngth - 1
         end do
         nxadj(k+1) = nxadj(k) + klngth
         istrt = iend + 1
      end do
      l = 0
      do k = 1, kend
         if (nadj(k) .ne. 0) then
            l = l + 1
            nadj(l) = nadj(k)
         end if
      end do
      return
      end
C======================================================================
      subroutine depth_first_search_block(n,xadj,adjncy,
     >     p_vstack,lowlink,nblk,subtree_blk,edgeptr,lstack_number)
C======================================================================
      integer*4 xadj(*),adjncy(*),p_vstack(*),subtree_blk(*)
      integer*4 n,np1,nb,nblk,cnt,lfp1,sp,vp,v,wp,w,qvp1,i
      integer*4 lowlink(*),edgeptr(*),lstack_number(*)
cc      integer*4 q(*)
C---------------------------------------------------------------------
C...  This subroutine performs the Tarjan's algorithm for finding the
C...  strong components of a diggraph. 
C...  This is the fred gustavson's implementation of the algorithm. 
C...  This is an implementation without q(). q is commented out.
C---------------------------------------------------------------------
C...  Initialization.
C
      nblk = 0
      nb = 0
C
      cnt  = nb+1
      lfp1 = cnt
      np1 = n+1
      vp = np1
      sp = np1
C
      do i = 1 , n
         edgeptr(i) = xadj(i)
         lstack_number(i) = 0
         lowlink(i) = 0
      end do
C
      lstack_number(n+1) = 0
C
C...  Get out when the output renumbering is full;
C...  Otherwise begin the search at vertex LFP1 for another 
C...  tree in the forest
C...
C
 10   if(cnt .eq. np1) then
C
C...     WE DO NOT REVERSE THE ORDERING NOW!!!
C...     Reverse the ordering and exit.
cc         do i = 1,n/2
cc            ipp = p_vstack(i)
cc            p_vstack(i) = p_vstack(n-i+1)
cc            p_vstack(n-i+1) = ipp
cc         enddo
C
C...     Now the starting addresses for the blocks are in lowlink.
cc         subtree_blk(nblk+1) = n+1
cc         lowlink(1) = 1
cc         do i = 1, nblk
cc            lowlink(i+1) = lowlink(i)+
cc     >           subtree_blk(nblk-i+2)-subtree_blk(nblk-i+1)
cc         enddo
         call icopyv(subtree_blk,lowlink,nblk+1)
         return
      end if
C
      do i = lfp1,n
         if(lstack_number(i) .eq. 0) go to 30
      end do
      write(*,*) ' There is an error in DEPTH FIRST SEARCH.'
      stop 1
C
 30   continue
      v = i
      lfp1 = v + 1
      go to 50
C     
C...  Recursive call of sico or whatever:::
 40   continue
      vp = vp - 1
      subtree_blk(vp) = v
      v=w
C
C...  Add tree branch to the current tree::::
 50   continue
      nb = nb + 1
      lstack_number(v) = nb
      lowlink(v) = lstack_number(v)
      sp = sp - 1
      p_vstack(sp) = v
      qvp1 = v+1
cc      qvp1 = q(v)+1
C
C...  Examine all the edges leading out of v.
 60   continue
      wp = edgeptr(v)
      w = adjncy(wp)
      edgeptr(v) = wp+1
      if(lstack_number(w) .ge. lstack_number(v)) go to 70
      if(lstack_number(w) .eq. 0) go to 40
      lowlink(v) = min0(lowlink(v),lstack_number(w))
 70   continue
      if(edgeptr(v) .lt. xadj(qvp1)) go to 60
C
C...   Check if v is a strong component root:::
C...   Gather the component if it is:::
C  
      if(lowlink(v) .lt. lstack_number(v)) go to 90
C
      nblk = nblk + 1
      subtree_blk(nblk) = cnt
C
 80   continue
      w = p_vstack(sp)
      lstack_number(w) = np1
      sp = sp + 1
      p_vstack(cnt) = w
      cnt = cnt + 1
      if(v .ne. w) go to 80
C
C...  Go to the begining. The present tree is finished...
C
      if(sp .eq. np1) go to 10
 90   continue
      w = v
      v = subtree_blk(vp)
      vp = vp + 1
      qvp1 = v + 1
cc    qvp1 = q(v) + 1
      lowlink(v) = min0(lowlink(v),lowlink(w))
      go to 70
C
      end
C====================================================================
      subroutine permat(iord,ia,ja,an,n,m,iat,jat,ant)
C====================================================================
      implicit real*8(a-h,o-z)
      integer*4 iord(*),ia(*),ja(*),iat(*),jat(*),n,m
      dimension an(*),ant(*)
C---------------------------------------------------------------------
C...  Permutes the matrix (P^t*A*P)
C---------------------------------------------------------------------
cc      call outmat(ia,ja,an,n)
cc      write(*,*) ' iord :' 
cc      write(*,*) (iord(k),k=1,n)
cc      read(*,*)
C
      call perm0(iord,ia,ja,an,n,m,iat,jat,ant)
C
cc      call aat(iat,jat,ant,n,m,ia,ja,an)
cc      call outmat(ia,ja,an,n)
cc      read(*,*)
C
      call perm0(iord,iat,jat,ant,m,n,ia,ja,an)
C
      return
      end
C====================================================================
      subroutine pervec(iord,u1,u2,n)
C====================================================================
      implicit real*8(a-h,o-z)
      integer*4 iord(*)
      dimension u1(*),u2(*)
C---------------------------------------------------------------------
C...  Permutes the vector (u1 = P^t*u1)
C---------------------------------------------------------------------
      do k = 1, n
         u2(k) = u1(k)
      end do
      do k = 1, n
         u1(k) = u2(iord(k))
      end do
      return
      end
C====================================================================
      subroutine perback(iord,u1,u2,n)
C====================================================================
      implicit real*8(a-h,o-z)
      integer*4 iord(*)
      dimension u1(*),u2(*)
C---------------------------------------------------------------------
C...  Permutes to the original vector (u1 = P*u1)
C---------------------------------------------------------------------
      do k = 1, n
         u2(k) = u1(k)
      end do
      do k = 1 , n
         u1(iord(k)) = u2(k)
      end do
      return
      end
C====================================================================
      subroutine perm0(iord,ia,ja,an,n,m,iat,jat,ant)
C====================================================================
      implicit real*8(a-h,o-z)
      integer*4 iord(*),ia(*),ja(*),iat(*),jat(*),n,m
      dimension an(*),ant(*)
C--------------------------------------------------------------------
C...  Permutation of the columns of a general sparse matrix A.
C...
C...  Input:
C...    IA, JA, AN    - given matrix A in RRCU.
C...    IORD          - the permutation. 
C...    N             - number of rows of the matrix.
C...    M             - number of columns of the matrix.
C...
C...  Output:
C...    IAT, JAT, ANT - transposed matrix A^t in RRCO.
C...
C...  Note:
C...    N+1 is the dimension of IA.
C...    M+1 is the dimension of IAT.
C--------------------------------------------------------------------
      mh = m + 1
      nh = n + 1
      do 10 i = 2, mh
         iat(i) = 0
 10   continue
      iab = ia(nh) - 1
      do 20 i = 1, iab
         j = ja(i) + 2
         if(j .le. mh) iat(j) = iat(j) + 1
 20   continue
      iat(1) = 1
      iat(2) = 1
      if (m .ne. 1) then
         do 30 i = 3, mh
            iat(i) = iat(i) + iat(i-1)
 30      continue
      end if
      do 50 i = 1, n
         iaa = ia(iord(i))
         iab = ia(iord(i)+1) - 1
         if(iab .lt. iaa) go to 50
         do 40 jp = iaa, iab
            j = ja(jp) + 1
            k = iat(j) 
            jat(k) = (i)
	    ant(k) = an(jp)
            iat(j) = k + 1
 40      continue
 50   continue
C
      return
      end
C=====================================================================
      subroutine random(iperm,n)
C=====================================================================
      parameter (imax = 1 000 000 000)
      implicit real*8(a-h,o-z)
      integer*4 iperm(*)
C---------------------------------------------------------------------
C...  Generating a permutaion (quasi-random order) of N integers 
C...  {1,2,...,N}.
C...  
C...  Parameter
C...    IPERM - a permutation of N integers {1,2,...,N}.
C---------------------------------------------------------------------
      k  = 0
      n1 = n + 1
C
      do 10 i = 1, n
         iperm(i) = 0
 10   continue
C
      do 40 i = 1, imax
         x  = dble(i)
         x  = dsin(x)
         x  = dabs(x)
         k1 = n1*x
         k1 = mod(k1, n) + 1
         if (iperm(k1) .eq. 0) then
            k         = k + 1 
            iperm(k1) = k
            if (k .eq. n) go to 50
         end if
 40   continue
C 
 50   continue
cc      write(*,*) (i,'*',iperm(i), i=1,n)
C
      return
      end
C=====================================================================
