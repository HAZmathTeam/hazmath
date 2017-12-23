C====================================================================
      subroutine mat_conv(ia,ja,a,n,nnz,ir,ic,aij)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),ir(1),ic(1)
      dimension a(1),aij(1)
C--------------------------------------------------------------------
C...  This subroutine converts the structures of a given matrix from
C...  (IR,IC,AIJ) to (IA,JA,A).
C...
C...  Input:
C...    IR, IC - the row and column indices of a given matrix. 
C...    AIJ    - numerical values of the matrix corresponding to IR 
C...             and IC. 
C...  Output: 
C...    IA, JA - structure of matrix A in RRCU. The order of A is N.
C...    A      - numerical values of nonzeros of A in RRCU.
C...    N      - the size of the matrix is N by N.
C...    NNZ    - number of nonzeros of A.
C--------------------------------------------------------------------
      open(20,file='mat_input',status='unknown',form='formatted')
cc      open(20,file='mat_data',status='unknown',form='formatted')
C
      ni   = 0
      nj   = 0
      nnz  = 0
      zero = 1.d-12
C
      k = 0
 10   read(20,*,END=20,ERR=20) irt, ict, aijt
      if (dabs(aijt) .gt. zero) then 
         k      = k + 1
         ir(k)  = irt
         ic(k)  = ict
         aij(k) = aijt
         ni     = max(ni,ir(k))
         nj     = max(nj,ic(k))
      end if
      go to 10
C
 20   continue
C
      if (ni .eq. nj) then
         n = ni
      else
         n = max(ni,nj)
         write(*,*), 'No. of rows and columns do not match:', ni, nj
         write(*,*)
      end if
C
      nnz = k
      do k = 1, n+1
         ia(k) = 0
      end do
cc      call inullv(ia,n+1)
C
      do 30 k = 1, nnz
         irk     = ir(k) + 1
         ia(irk) = ia(irk) + 1
 30   continue
C
      do 40 k = 2, n+1
         if (ia(k) .ne. 0) then
            nzk = k - 1
            go to 50
         end if
 40   continue
C
 50   continue
C
      do 60 k = 1, nzk
         ia(k) = 1
 60   continue
C
      do 70 k = nzk+1, n+1
         ia(k) = ia(k-1) + ia(k)
 70   continue
C
      do 80 k = 1, n
         ica = ia(k)
         icb = ia(k+1)
         if (icb .gt. ica)  ja(ia(k+1)-1) = ia(k)
 80   continue
C
      do 90 k = 1, nnz
         ik     = ir(k)
         iend   = ia(ik+1) - 1
         jp     = ja(iend)
         ja(jp) = ic(k)
         a(jp)  = aij(k)
         if (iend .ne. jp)  ja(iend) = ja(iend) + 1
 90   continue
C     
      return
      end
C====================================================================
