C=====================================================================
      subroutine ilu0(n,ia,ja,a,a_lu,lu_d,iw,icode)
C=====================================================================
      integer n,ia(n+1),ja(*),lu_d(n),iw(n)
      real*8  a(*),a_lu(*)
C---------------------------------------------------------------------
C...  Set-up routine for ILU(0) preconditioner. This routine computes 
C...  the L and U factors of the ILU(0) factorization of a general 
C...  sparse matrix A stored in CSR format. Since L is unit triangular, 
C...  the L and U factors can be stored as a single matrix which 
C...  occupies the same storage as A. The ia and ja arrays are not 
C...  needed for the LU matrix since the pattern of the LU matrix is 
C...  identical with that of A.
C...
C...  INPUT :
C...    n       - dimension of matrix A
C...    ia,ja,a - sparse matrix in CSR format
C...    iw      - integer work array of length n
C...
C...  Output :
C...    a_lu    - L/U matrices stored together. On return, ia, ja, 
C...              a_lu are the combined CSR data structure for the
C...              LU factors
C...    lu_d    - pointer to the diagonal elements in the CSR data
C...              structure ia, ja, a_lu
C...    icode   - integer indicating error code on return
C...              0 : normal return
C...              k : encountered a zero pivot at step k
C---------------------------------------------------------------------
C
      zero = 1.0d-14
C
C...  Initialize work array iw to zero and a_lu array to a.
C
CCCC      inullv(iw,n)
CCCC      copyv(a,a_lu,ia(n+1)-1)
C
      do 100 i = 1, n
         iw(i) = 0
 100  continue
C
      do 200 i = 1, ia(n+1) - 1
         a_lu(i) = a(i)
 200  continue
C
C...  Main loop
C
      do 800 k = 1, n
         j1 = ia(k)
         j2 = ia(k+1) - 1
         do 300 j = j1, j2
            iw(ja(j)) = j
 300     continue
C
         j = j1
 400     continue
         jrow = ja(j)
C
C...     Exit if diagonal element is reached.
         if (jrow .ge. k) goto 600
C
C...     Compute the multiplier for jrow.
         tl = a_lu(j)*a_lu(lu_d(jrow))
         a_lu(j) = tl
C
C...     Perform linear combination.
         do 500 jj = lu_d(jrow) + 1, ja(jrow+1) - 1
            jw = iw(ja(jj))
            if (jw .ne. 0) a_lu(jw) = a_lu(jw) - tl*a_lu(jj)
 500     continue
C
         j = j + 1
         if (j .le. j2) goto 400
C
C...     Store pointer to diagonal element.
 600     lu_d(k) = j
         if (jrow .ne. k .or. a_lu(j) .eq. zero) goto 900
         a_lu(j) = 1.0d0/a_lu(j)
C
C...     Refresh all entries of iw to zero.
         do 700 i = j1, j2
            iw(ja(i)) = 0
 700     continue
 800  continue
C
      icode = 0
      return
C
 900  icode = k
      write(*,'(10x,a,i7)') ' Error: zero pivot at Step : ', icode
C
      return
      end
C=====================================================================
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  Algorithm ILU(0)
cc  ----------------
cc  for i = 2, n 
cc      for k = 1,i-1 and for a_ik .ne. 0
cc          compute a_ik = a_ik / a_kk
cc          for j = k+1, n and for a_ij .ne. 0
cc              compute a_ij = a_ij - a_ik * a_kj
cc          end
cc      end
cc  end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

