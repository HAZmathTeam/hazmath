C=====================================================================
      subroutine fwd_gs_ns(x,ia,ja,an,b,n,max_sweeps) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by backward
C...  Gauss-Seidel method.
C---------------------------------------------------------------------
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
      do 40 i = 1, n
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 40
         u = b(i)    
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
               u = u - an(j) * x(ja(j))
            end if
 30      continue
         u    = u/an(id)
         ud   = u - x(i)
         umax = dmax1(umax,dabs(ud)/dabs(an(id)))
         x(i) = u
 40   continue
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      return
      end
C=====================================================================
      subroutine fwd_gs_ns_blk(x,ia,ja,an,b,dinv,n,m,max_sweeps) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(m,m,1),b(m,1),x(m,1),u(10),v(10)
      dimension dinv(m,m,1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by forward
C...  block Gauss-Seidel method.
C---------------------------------------------------------------------
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
      do 40 i = 1, n
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 40
         call copyv(b(1,i),u,m)    
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
C...           U and X are vectors, and an(j) is an m by m matrix.
C...           U = U - A*X.
               call uuminav(u,an(1,1,j),x(1,ja(j)),m)
            end if
 30      continue
C
C...     Compute U = D^(-1)*U.
         call copyv(u,v,m)
         call uabyv(u,dinv(1,1,i),v,m)
C
         do ii = 1, m    
            ud = u(ii) - x(ii,i)
            umax = dmax1(umax,dabs(ud))
         end do
         call copyv(u,x(1,i),m)
 40   continue
C
      if (iteration .lt. max_sweeps .and. 
     >     umax .gt. 0.5d-13) go to 20
C
      return
      end
C=====================================================================
      subroutine bwd_gs_ns(x,ia,ja,an,b,n,max_sweeps) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by backward
C...  Gauss-Seidel method.
C---------------------------------------------------------------------
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
      do 40 i = n, 1, -1
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 40
         u = b(i)    
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
               u = u - an(j) * x(ja(j))
            end if
 30      continue
         u    = u/an(id)
         ud   = u - x(i)
         umax = dmax1(umax,dabs(ud)/dabs(an(id)))
         x(i) = u
 40   continue
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      return
      end
C=====================================================================
      subroutine bwd_gs_ns_blk(x,ia,ja,an,b,dinv,n,m,max_sweeps) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(m,m,1),b(m,1),x(m,1),u(10),v(10)
      dimension dinv(m,m,1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by backward
C...  block Gauss-Seidel method.
C---------------------------------------------------------------------
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
      do 40 i = n, 1, -1
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 40
         call copyv(b(1,i),u,m)    
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
C...           U and X are vectors, and an(j) is an m by m matrix.
C...           U = U - A*X.
               call uuminav(u,an(1,1,j),x(1,ja(j)),m)
            end if
 30      continue
C
C...     Compute U = D^(-1)*U.
         call copyv(u,v,m)
         call uabyv(u,dinv(1,1,i),v,m)
C
         do ii = 1, m    
            ud = u(ii) - x(ii,i)
            umax = dmax1(umax,dabs(ud))
         end do
         call copyv(u,x(1,i),m)
 40   continue
C
      if (iteration .lt. max_sweeps .and. 
     >     umax .gt. 0.5d-13) go to 20
C
      return
      end
C=====================================================================
      subroutine sym_gs_ns(x,ia,ja,an,b,n,max_sweeps) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by symmetric
C...  Gauss-Seidel method.
C---------------------------------------------------------------------
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
C
C...  Forward Gauss-Seidel
C  
      do 40 i = 1, n
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 40
         u = b(i)    
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
               u = u - an(j) * x(ja(j))
            end if
 30      continue
         x(i) = u/an(id)
 40   continue
C
C...  Backward Gauss-Seidel
C  
      do 60 i = n, 1, -1
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 60
         u = b(i)    
         do 50 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
               u = u - an(j) * x(ja(j))
            end if
 50      continue
         u    = u/an(id)
         ud   = u - x(i)
         umax = dmax1(umax,dabs(ud)/dabs(an(id)))
         x(i) = u
 60   continue
C
c      write(*,*) (x(i),i=1,n)
c      read(*,*)
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      return
      end
C=====================================================================
      subroutine sym_gs_ns_blk(x,ia,ja,an,b,dinv,n,m,max_sweeps) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(m,m,1),b(m,1),x(m,1),u(10),v(10)
      dimension dinv(m,m,1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by symmetricd
C...  block Gauss-Seidel method.
C---------------------------------------------------------------------
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
C
C...  Forward Gauss-Seidel
C 
      do 40 i = 1, n
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 40
         call copyv(b(1,i),u,m)    
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
C...           U and X are vectors, and an(j) is an m by m matrix.
C...           U = U - A*X.
               call uuminav(u,an(1,1,j),x(1,ja(j)),m)
            end if
 30      continue
C
C...     Compute X = D^(-1)*U.
         call uabyv(x(1,i),dinv(1,1,i),u,m)
 40   continue
C
C...  Backward Gauss-Seidel
C 
      do 60 i = n, 1, -1
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 40
         call copyv(b(1,i),u,m)    
         do 50 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
C...           U and X are vectors, and an(j) is an m by m matrix.
C...           U = U - A*X.
               call uuminav(u,an(1,1,j),x(1,ja(j)),m)
            end if
 50      continue
C
C...     Compute U = D^(-1)*U.
         call copyv(u,v,m)
         call uabyv(u,dinv(1,1,i),v,m)
C
         do ii = 1, m    
            ud = u(ii) - x(ii,i)
            umax = dmax1(umax,dabs(ud))
         end do
         call copyv(u,x(1,i),m)
 60   continue
C
      if (iteration .lt. max_sweeps .and. 
     >     umax .gt. 0.5d-13) go to 20
C
      return
      end
C=====================================================================
      subroutine jacobi(x,ia,ja,an,b,n,max_sweeps,xold,omega) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1),xold(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by Jacobi
C...  method.
C---------------------------------------------------------------------
      call copyv(x,xold,n)
C
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
      do 40 i = 1,n
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 40
         u = b(i)    
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
               u = u - an(j) * xold(ja(j))
            end if
 30      continue
         u    = u/an(id)
         ud   = u - x(i)
         umax = dmax1(umax,dabs(ud)/dabs(an(id)))
         x(i) = omega*u + (1-omega)*xold(i)
 40   continue
C
      call copyv(x,xold,n)
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      return
      end
C=====================================================================
      subroutine jacobi_blk(x,ia,ja,an,b,dinv,n,m,
     >     max_sweeps,xold,omega) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(m,m,1),b(m,1),x(m,1),u(10),v(10)
      dimension dinv(m,m,1),xold(m,1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by block Jacobi
C...  method.
C---------------------------------------------------------------------
      call copyv(x,xold,n*m)
C
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
      do 40 i = 1,n
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .lt. iaa) go to 40
         call copyv(b(1,i),u,m)
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then
               id = j
            else
C...           U and X are vectors, and an(j) is an m by m matrix.
C...           U = U - A*X.
               call uuminav(u,an(1,1,j),xold(1,ja(j)),m) 
            end if
 30      continue
C
C...     Compute U = D^(-1)*U.
         call copyv(u,v,m)
         call uabyv(u,dinv(1,1,i),v,m)
         do ii = 1, m    
            ud = u(ii) - x(ii,i)
            umax = dmax1(umax,dabs(ud))
         end do
         do k = 1, m
            x(k,i) = omega*u(k) + (1-omega)*xold(k,i)
         end do
 40   continue
C
      call copyv(x,xold,n*m)
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      return
      end
C=====================================================================
      subroutine uuminav(u,a,v,n)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension a(n,n),u(1),v(1)
C---------------------------------------------------------------------
C...  Compute  U = U - A*V.
C---------------------------------------------------------------------
      do i = 1 , n
         hold = 0d0
         do j = 1 , n
            hold = hold + a(i,j)*v(j)
         end do
         u(i) = u(i) - hold
      end do
      return
      end
C====================================================================
      subroutine uuplsav(u,a,v,n)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension a(n,n),u(1),v(1)
C---------------------------------------------------------------------
C...  Compute  U = U + A*V.
C---------------------------------------------------------------------
      do i = 1 , n
         hold = 0d0
         do j = 1 , n
            hold = hold + a(i,j)*v(j)
         end do
         u(i) = u(i) + hold
      end do
      return
      end
C====================================================================
      subroutine uabyv(u,a,v,n)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension a(n,n),u(1),v(1)
C---------------------------------------------------------------------
C...  Compute  U = A*V.
C---------------------------------------------------------------------
      do i = 1 , n
         hold = 0d0
         do j = 1 , n
            hold = hold + a(i,j)*v(j)
         end do
         u(i) = hold
      end do
      return
      end
C====================================================================

C=====================================================================
      subroutine inv__blk_gauss(ia,ja,an,dinv,n,m) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(m,m,1),atemp(3.4),soln(3)
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
         do 30 j = iaa,iab
            if(ja(j) .eq. i) then  
               do 20 k = 1, m
                  copyv(an(1,1,j),atemp(1,1),m*m)
                  do 10 kk = 1, m
                     if (kk .eq. k) then
                        atemp(kk,m+1) = 1.d0
                     else 
                        atemp(kk,m+1) = 0.d0
                     end if
 10               continue
                  call gauspp(m,m+1,m,atemp,soln)
                  copyv(soln,dinv(1,k,i),m)
 20            continue
            end if
 30      continue
 100  continue
C
      return 
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
C=====================================================================
