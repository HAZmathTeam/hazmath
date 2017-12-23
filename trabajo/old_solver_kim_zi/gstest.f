C====================================================================
      program gstest
C====================================================================
      implicit real*8(a-h,o-z)
      parameter(nsubmax = 10, max_dditer = 6)
      parameter (nlast = 5 000 000, nrlast = 2 500 000) !level = 8
      common /top/ m(nlast)              !To use the machine witt
      common /top1/ r(nrlast)            !To use the machine witt
cc      dimension m(nlast), r(nrlast)
      dimension atemp(3,4),soln(3)

      n = 10
      kk = 3
      ia = 1
      ja = 11
      kdinv = 21
      do ii = 1 , kk
         do jj = 1 , kk
            atemp(ii,jj) = dble(ii)/dble(jj)*(ii+jj)-jj
         end do
      end do  
      
      do 20 k = 1, kk
         do 10 k1 = 1, kk
            if (k1 .eq. k) then
               atemp(k1,kk+1) = 1.d0
            else 
               atemp(k1,kk+1) = 0.d0
            end if
 10      continue
         do ii =1,3
            write(*,*) ii, (atemp(ii,jj), jj=1,4)
         end do
         call gauspp(kk,kk+1,kk,atemp,soln)
         write(*,*) k, (soln(jj), jj=1,3)         
 20   continue
         
      
      write(*,*) n,kk

      
cc      call inv_blk_gauss(ia,ja,r(1),r(kdinv),n,kk)
      stop
      end

C=====================================================================
      subroutine inv_blk_gauss(ia,ja,an,dinv,n,m) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(m,m,1),atemp(3,4),soln(3)
      dimension dinv(m,m,1)
C---------------------------------------------------------------------
C...  This subroutine finds the inverse of a block diagonal matrix of
C...  given block matrix AN. Gaussian elimination with partial pivoting 
C...  will be used to invert the each diagonal block of AN.
C...
C...  Output
C...    DINV - inverse of a block diagonal matrix of a block matrix AN.
C---------------------------------------------------------------------
      do 100 i = 1, n
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 100
ccccccccccccccccc         write(*,*) " node ", i,iaa,iab 
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then  
               do ii = 1 , m
                  do jj = 1 , m
                     atemp(ii,jj) = ii*jj
cc                    write(*,*) ii,jj,an(ii,jj,j)
                    write(*,*) n,m
                  end do
               end do
               do 20 k = 1, m
                  do 10 kk = 1, m
                     if (kk .eq. k) then
                        atemp(kk,m+1) = 1.d0
                     else 
                        atemp(kk,m+1) = 0.d0
                     end if
 10               continue
cc                  call nullv(soln,3)
                  call gauspp(m,m+1,m,atemp,soln)
                  do ii = 1 , m
                     dinv(ii,k,i) = soln(ii)
                  end do
 20            continue
               do l1 = 1, 3
                  do l2 = 1,4
                     write(*,*) atemp(l1,l2)
                  end do
               end do

               do l1 = 1, 3
                  do l2 = 1,3
                     write(*,*) dinv(l1,l2,i)
                  end do
               end do

            end if
 30      continue
 40      continue
 100  continue
C
55      return 
      end
C=====================================================================
      subroutine gauspp(nrows, ncols, nd, cm, soln)
C=====================================================================
      implicit real*8(a-h, o-z)
      dimension cm(nrows, ncols), soln(nrows)
      logical error 
C---------------------------------------------------------------------
C...  This subroutine uses the Gaussian elimination with partial pivoting 
C...  to solve a linear system. All coefficients and the load are stored 
C...  in the augmented matrix cm.
C...  In general, NCOLS = NROWS + 1 and ND = NCOLS.
C---------------------------------------------------------------------
      npivot = 1
      error  = .false.
C
 10   if (npivot.lt.nd .and. .not.error) then
C...     !To reorder the equations so that the pivot position in the pivot 
C...     !equation has the maximum absolute value.
         maxrow = npivot
         do 20 irow = npivot + 1, nd
            if (dabs(cm(irow, npivot)) .gt. dabs(cm(maxrow, npivot))) 
     >           maxrow = irow
 20      continue
         if (dabs(cm(maxrow, npivot)) .lt. 1.d-15) then
            error = .true.
         else
            if (maxrow .ne. npivot) then
               do 30 k = 1, nd + 1
                  temp = cm(maxrow, k)
                  cm(maxrow, k) = cm(npivot, k)
                  cm(npivot, k) = temp
 30            continue
            end if 
         end if
C
         if (.not.error) then
C...        !To eliminate the element in the pivot position from rows following
C...        !the pivot equation.
            do 50 irow = npivot + 1, nd
               factor = cm(irow, npivot)/ cm(npivot, npivot)
               cm(irow, npivot) = 0.d0
               do 40 icol = npivot + 1, nd +1
                  cm(irow,icol) = cm(irow,icol) - cm(npivot,icol)*factor
 40            continue
 50         continue
            npivot = npivot + 1
         end if
         go to 10
      end if
C
      if (error) then
         print*, ' No unique solution exists.'
      else
C...     !To perform the back-substitution to determine the solution to the 
C...     !system of equations.
         do 70 irow = nd, 1, -1
            do 60 icol = nd, irow + 1, -1
               cm(irow,nd+1) = cm(irow,nd+1) - soln(icol)*cm(irow,icol)
 60         continue
            soln(irow) = cm(irow, nd+1)/cm(irow, irow)
 70      continue
      end if
C
      return
      end

