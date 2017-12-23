C====================================================================
      subroutine iit(ia,ja,n,m,iat,jat)
C====================================================================
      integer*4 ia(1),ja(1),iat(1),jat(1),n,m
C--------------------------------------------------------------------
C...  Transposition of a graph (or the matrix) symbolically.
C...
C...  Input:
C...    IA, JA   - given graph (or matrix).
C...    N        - number of rows of the matrix.
C...    M        - number of columns of the matrix.
C...
C...  Output:
C...    IAT, JAT - transposed graph (or matrix).
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
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .lt. iaa) go to 50
         do 40 jp = iaa, iab
            j = ja(jp) + 1
            k = iat(j) 
            jat(k) = i
            iat(j) = k + 1
 40      continue
 50   continue
C
      return
      end
C====================================================================
      subroutine aat(ia,ja,an,n,m,iat,jat,ant)
C====================================================================
      implicit real*8(a-h,o-z)
      integer*4 ia(1),ja(1),iat(1),jat(1),n,m
      dimension an(1),ant(1)
C--------------------------------------------------------------------
C...  Transposition of a general sparse matrix A to A^t.
C...
C...  Input:
C...    IA, JA, AN    - given matrix A in RRCU.
C...    N             - number of rows of the matrix.
C...    M             - number of columns of the matrix.
C...
C...  Output:
C...    IAT, JAT, ANT - transposed matrix A^t in RRCO.
C...
C...  Note:
C...    N+1 is the dimension of IA.
C...    M+1 is the dimension of IA.
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
         iaa = ia(i)
         iab = ia(i+1) - 1
         if(iab .lt. iaa) go to 50
         do 40 jp = iaa, iab
            j = ja(jp) + 1
            k = iat(j) 
            jat(k) = i
	    ant(k) = an(jp)
            iat(j) = k + 1
 40      continue
 50   continue
C
      return
      end
C====================================================================
