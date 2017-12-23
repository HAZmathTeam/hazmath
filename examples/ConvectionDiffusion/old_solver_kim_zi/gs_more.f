C=====================================================================
      subroutine gauss_seidel(ia,ja,mwk,ex_sol,sol,rhs,ra,wk,iexact,
     >     n,nnz,igs,max_sweeps,norder,tol,iwr)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),mwk(1)
      dimension sol(1),rhs(1),ra(1),ex_sol(1),wk(1)
C---------------------------------------------------------------------
C...  Gauss-seidel iterative methods with Tarjan's ordering.
C...  Initial guess should be an input in the vector sol.
C
C...  Parameters:
C...    MAX_SWEEPS - maximum number of sweeps
C...    IGS        - 1: forward G-S, 2: backward G-S, 3: sym G-S.
C...    NORDER     - 1: Tarjan's ordering, 0: no ordering
C...    TOL        - stopping tolerance for the iterations
C...    IWR        - write the convergence history or not: 1 or 0
C...    IEXACT     - compute the relative error with exact soln U, 
C...                 ||U-U_k||_0 / ||U||_0, if it is 2.
C---------------------------------------------------------------------
      if (norder .eq. 1) then
         iord   = 1
         iaw    = iord + n + 1
         jaw    = iaw + n + 1 + 1
         isub   = jaw + nnz 
         iwork1 = isub + n + 1 + 1
         iwork2 = iwork1 + n + 1
         iwork3 = iwork2 + n + 1
         lasti  = iwork3 + n + 1 + 1
         call cut_off(
     >        n,ia,ja,ra,mwk(iaw),mwk(jaw),mwk(iord),
     >        mwk(isub),nblk,mwk(iwork1),mwk(iwork2),mwk(iwork3))
      else
         iord   = 1
         mwk(iord) = 0
      end if
C
C...  Zero initial guess.
cc      call init_g0(ia,ja,ra,rhs,sol,n)
C
      zero = 1.0d-13
      iter = 1
      iwr_in = 0
C
      if (iexact .eq. 2) then
         call l2norm(ex_sol,soln_norm,n,1.0d0)
cc         write(*,*) ' soln_norm: ', soln_norm
         if (soln_norm .le. zero) then
            soln_norm = 1.d0
            write(*,*) 'Zero exact solution is detected!!!!.'
         end if
      end if
C
C...  To compute the initial residual: wk = rhs - A*sol.
      call abyvam(rhs,ia,ja,ra,sol,n,wk)
C
      call scpro(wk,wk,err_res0,n)
      err_res0 = dsqrt(err_res0)
C
      if (err_res0 .le. zero) then
         err_res0 = 1.d0
         write(*,*) 'Zero right hand side is detected!!!!.'
      end if
C
      err_relres = err_res0
      err_rel = 0.d0
C
      if (irw .eq. 1)
      write(*,'(2X,i5,2X,3(3X,e11.4))') 
     >     0,1.d0,err_res0,1.d0
C
C...  Iteration starts here.
C
      iteration = 0
 20   iteration = iteration + 1
C
      if (igs. eq. 1) then
         call fwd_gs_ord(sol,ia,ja,ra,rhs,mwk(iord),n,iter,iwr_in)
      else if (igs. eq. 2) then
         call bwd_gs_ord(sol,ia,ja,ra,rhs,mwk(iord),n,iter,iwr_in)
      else if (igs. eq. 3) then
         call sym_gs_ord(sol,ia,ja,ra,rhs,mwk(iord),n,iter,iwr_in)
      end if
C
      if (iteration .gt. 100) then
         nmod = mod(iteration,1000)
      else
         nmod = 1
      end if
C
C...  To compute the residual: wk = rhs - A*sol.
      call abyvam(rhs,ia,ja,ra,sol,n,wk)
C
C...  To compute the residual error.
      call scpro(wk,wk,err_res,n)
      err_res = dsqrt(err_res) 
      err_relres = err_res/err_res0
C
C...  To compute the relative error if the exact solution is known.
      if (iexact .eq. 2) then
         call wuminv(wk,ex_sol,sol,n)
         call l2norm(wk,xdiffn,n,1.0d0)
         err_rel = xdiffn/soln_norm  
      end if
C
      if (nmod .eq. 1 .and. iwr .eq. 1) then
         write(*,'(2X,i5,2X,3(3X,e11.4))')
     >        iteration,err_relres,err_res,err_rel
      end if
C
      if (iteration .gt. max_sweeps) then
         write(*,*), ' Iteration limit exceeded:', max_sweeps
         go to 100
      end if
C
      if  (iexact .eq. 2) then
         if (err_rel .gt. tol) go to 20
      else
         if (err_relres .gt. tol) go to 20
      end if
C
 100  continue
C
      if (iwr .eq. 1) then
         write(*,'(a,a)') ' =============================',
     >        '=================================================='
C
         arfac  = (err_relres)**(1./dble(kk))
         arfac1 = (err_rel)**(1./dble(kk))
         write(*,*) 
         write(*,'(a,2f12.5)'),' Average Reduction Factor:', arfac,arfac1
C
         write(*,*)
      end if
C
      return
      end
C=====================================================================
      subroutine init_g0(ia,ja,a,b,x,n)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),a(1),x(1),b(1)
C--------------------------------------------------------------------
C     Initial guess to be zero except on the dirichlet boundary.
C--------------------------------------------------------------------
      call nullv(x,n)
C
      do k = 1 , n
         i1 = ia(k)
         i2 = ia(k+1)-1
         if(i2-i1 .lt. 0) stop 10
         if(i2-i1 .eq. 0) then
            x(k) = b(k)
            a(i1) = 1.d0
         else
cc           isx = ifix(sngl(1d0/rng(k*k*k)))
cc           x(k)= dble(k)/dble(n)
cc           x(k) = rng(isx)*(1-2*mod(k,3))
           x(k) = 0.d0
         end if
      end do
      return
      end
C====================================================================
      subroutine init_g_random(ia,ja,a,b,x,n)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),a(1),x(1),b(1)
C--------------------------------------------------------------------
C     Initial guess randomly except on the dirichlet boundary.
C--------------------------------------------------------------------
      call nullv(x,n)
C
      do k = 1 , n
         i1 = ia(k)
         i2 = ia(k+1)-1
         if(i2-i1 .lt. 0) stop 10
         if(i2-i1 .eq. 0) then
            x(k) = b(k)
            a(i1) = 1.d0
         else
           isx = ifix(sngl(1d0/rng(k*k*k)))
           x(k)= dble(k)/dble(n)
           x(k) = rng(isx)*(1-2*mod(k,3))
cc           x(k) = 0.d0
         end if
      end do
      return
      end
C====================================================================
      REAL*8 FUNCTION RNG(ISEED)
C=====================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C--------------------------------------------------------------------
C     RANDOM NUMBER GENERATOR
C--------------------------------------------------------------------
      ISEED=MOD(ISEED*24+23,1361)
      RNG=ISEED/1359.*.99923+.00038
      RETURN
      END
C=====================================================================
      subroutine zero_dir(x,idir,n)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension idir(1),x(1)
C---------------------------------------------------------------------
C...  Assign components of x to be zero if they are Dirichlet nodes.
C---------------------------------------------------------------------
      do 10 k = 1,n
         if(idir(k) .gt. n) x(k) = 0.0d00
 10   continue
      return
      end
C=====================================================================
      subroutine fwd_gs_ord(x,ia,ja,an,b,iord,n,max_sweeps,iw) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1),iord(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by forward
C...  Gauss-Seidel method with Tarjan's ordering.
C...
C...  Parameter:
C...    IW   - To determine whether write the convergence history:
C...           IW = 0 : do not write it; IW = 1 : write it.
C...    IORD - a different ordering of unknowns will be used if the
C...           array contains some permutation of the original ordering.
C...           Please note that if no permutation is considered, then 
C...           simply assign iord(1) = 0. 
C---------------------------------------------------------------------
      if (iw .eq. 1) then
         nmodd = n/7
         if(nmodd .eq. 0) nmodd = 2
      end if
C
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
C
      if (iord(1) .eq. 0) then
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
 30         continue
            u    = u/an(id)
            ud   = u - x(i)
            umax = dmax1(umax,dabs(ud)/dabs(an(id)))
cc          umax = dmax1(umax,dabs(ud))
            x(i) = u
 40      continue
      else 
         do 60 i = 1, n
            io = iord(i)
            iaa = ia(io)
            iab = ia(io+1)-1
            if(iab .lt. iaa) go to 40
            u = b(io)    
            do 50 j = iaa,iab
               if(ja(j) .eq. io) then
                  id = j
               else
                  u = u - an(j) * x(ja(j))
               end if
 50         continue
            u    = u/an(id)
            ud   = u - x(io)
            umax = dmax1(umax,dabs(ud)/dabs(an(id)))
cc          umax = dmax1(umax,dabs(ud))
            x(io) = u
 60      continue
      end if
C
cc      nmodd = mod(iteration,nmodd)
      nmodd = 1
C
      if (nmodd .eq. 1 .and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' FGS : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in FGS: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine bwd_gs_ord(x,ia,ja,an,b,iord,n,max_sweeps,iw) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1),iord(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by backward
C...  Gauss-Seidel method with Tarjan's ordering.
C...
C...  Parameter:
C...    IW   - To determine whether write the convergence history:
C...           IW = 0 : do not write it; IW = 1 : write it.
C...    IORD - a different ordering of unknowns will be used if the
C...           array contains some permutation of the original ordering.
C...           Please note that if no permutation is considered, then 
C...           simply assign iord(1) = 0. 
C---------------------------------------------------------------------
      if (iw .eq. 1) then
         nmodd = n/7
         if(nmodd .eq. 0) nmodd = 2
      end if
C
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
C
      if (iord(1) .eq. 0) then
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
 30         continue
            u    = u/an(id)
            ud   = u - x(i)
            umax = dmax1(umax,dabs(ud)/dabs(an(id)))
cc          umax = dmax1(umax,dabs(ud))
            x(i) = u
 40      continue
      else
         do 60 i = n, 1, -1
            io = iord(i)
            iaa = ia(io)
            iab = ia(io+1)-1
            if(iab .lt. iaa) go to 40
            u = b(io)    
            do 50 j = iaa,iab
               if(ja(j) .eq. io) then
                  id = j
               else
                  u = u - an(j) * x(ja(j))
               end if
 50         continue
            u    = u/an(id)
            ud   = u - x(io)
            umax = dmax1(umax,dabs(ud)/dabs(an(id)))
cc          umax = dmax1(umax,dabs(ud))
            x(io) = u
 60      continue
      end if
C
cc      nmodd = mod(iteration,nmodd)
      nmodd = 1
C
      if (nmodd .eq. 1 .and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' BGS : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in BGS: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine sym_gs_ord(x,ia,ja,an,b,iord,n,max_sweeps,iw) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1),iord(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by symmetric
C...  Gauss-Seidel method with Tarjan's ordering.
C...
C...  Parameter:
C...    IW   - To determine whether write the convergence history:
C...           IW = 0 : do not write it; IW = 1 : write it.
C...    IORD - a different ordering of unknowns will be used if the
C...           array contains some permutation of the original ordering.
C...           Please note that if no permutation is considered, then 
C...           simply assign iord(1) = 0. 
C---------------------------------------------------------------------
      if (iw .eq. 1) then
         nmodd = n/7
         if(nmodd .eq. 0) nmodd = 2
      end if
C
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
C
      if (iord(1) .eq. 0) then
C
C...     Forward Gauss-Seidel
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
 30         continue
            x(i) = u/an(id)
 40      continue
C
C...     Backward Gauss-Seidel
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
 50         continue
            u    = u/an(id)
            ud   = u - x(i)
            umax = dmax1(umax,dabs(ud)/dabs(an(id)))
cc          umax = dmax1(umax,dabs(ud))
            x(i) = u
 60      continue
C
      else
C
C...     Forward Gauss-Seidel
C     
         do 80 i = 1, n
            io = iord(i)
            iaa = ia(io)
            iab = ia(io+1)-1
            if(iab .lt. iaa) go to 40
            u = b(io)    
            do 70 j = iaa,iab
               if(ja(j) .eq. io) then
                  id = j
               else
                  u = u - an(j) * x(ja(j))
               end if
 70         continue
            x(io) = u/an(id)
 80      continue
C
C...     Backward Gauss-Seidel
C  
         do 100 i = n, 1, -1
            io = iord(i)
            iaa = ia(io)
            iab = ia(io+1)-1
            if(iab .lt. iaa) go to 60
            u = b(io)    
            do 90 j = iaa,iab
               if(ja(j) .eq. io) then
                  id = j
               else
                  u = u - an(j) * x(ja(j))
               end if
 90         continue
            u    = u/an(id)
            ud   = u - x(io)
            umax = dmax1(umax,dabs(ud)/dabs(an(id)))
cc          umax = dmax1(umax,dabs(ud))
            x(io) = u
 100     continue
      end if
C
cc      nmodd = mod(iteration,nmodd)
      nmodd = 1
C
      if (nmodd .eq. 1 .and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' SGS : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in SGS: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine jacobi(x,ia,ja,an,b,n,xold,omega,max_sweeps,iw) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1),xold(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by Jacobi
C...  method.
C...
C...  Parameter:
C...    IW - To determine whether write the convergence history or not,
C...         IW = 0 : do not write it; IW = 1 : write it.
C----------------------------------------------------------------------
      if (iw .eq. 1) then
         nmodd = n/7
         if(nmodd .eq. 0) nmodd = 2
      end if
C
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
cc         umax = dmax1(umax,dabs(ud))
         x(i) = omega*u + (1-omega)*xold(i)
 40   continue
C
      call copyv(x,xold,n)
C
cc      nmodd = mod(iteration,nmodd)
      nmodd = 1
C
      if (nmodd .eq. 1 .and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' JACOBI : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in JACOBI: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
C=====================================================================
      subroutine fwd_gs_ord_blk(x,ia,ja,an,b,iord,n,kmatch,
     >     max_sweeps) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1),iord(1)
C---------------------------------------------------------------------
C...  This subroutine does MAX_SWEEPS steps on the linear system A*x = b 
C...  by forward  Gauss-Seidel method. IORD is the ordering if any. 
C...  If no ordering then iord(i) = i. This is a special smoother
C...  done for the matching coarsening. 
C---------------------------------------------------------------------
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
      do 40 i = kmatch+1,n
         io  = iord(i)
         iaa = ia(io)
         iab = ia(io+1)-1
         if(iab .lt. iaa) go to 40
         u = b(io)    
         do 30 j = iaa,iab
            if(ja(j) .eq. io) then
               id = j
            else
               u = u - an(j) * x(ja(j))
            end if
 30      continue
         u    = u/an(id)
         ud   = u - x(io)
         umax = dmax1(umax,dabs(ud)/dabs(an(id)))
cc         umax = dmax1(umax,dabs(ud))
         x(io) = u
 40   continue
C
C...  Solve directly 2 x 2 systems by inverting the matrices.
C
      do 140 i = 1, kmatch, 2
         io1  = iord(i)
         io2  = iord(i+1)
         iaa1 = ia(io1)
         iab1 = ia(io1+1)-1
         u1   = b(io1)    
         iaa2 = ia(io2)
         iab2 = ia(io2+1)-1
         u2   = b(io2)    
         do j = iaa1,iab1
            if(ja(j) .eq. io1) then
               a22 = an(j)
            else
               if(ja(j) .eq. io2) then
                  a12 = -an(j)
               else
                  u1 = u1 - an(j) * x(ja(j))
               end if
            end if
         end do
         do j = iaa2,iab2
            if(ja(j) .eq. io2) then
               a11 = an(j)
            else
               if(ja(j) .eq. io1) then
                  a21 = -an(j)
               else
                  u2 = u2 - an(j) * x(ja(j))
               end if
            end if
         end do
         delta = 1d0/(a11*a22-a12*a21)
         u    = delta*(a11*u1+a12*u2)
         u2   = delta*(a21*u1+a22*u2)
         u1   = u
         ud1  = u1 - x(io1)
         ud2  = u2 - x(io2)
         umax = dmax1(umax,dabs(ud1))
         umax = dmax1(umax,dabs(ud2))
cc         umax = dmax1(umax,dabs(ud))
         x(io1) = u1
         x(io2) = u2
 140  continue
C
cc      if (mod(iteration,10) .eq. 1) write(*,*) iteration,umax 
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      return
      end
C=====================================================================
      subroutine bwd_gs_ord_blk(x,ia,ja,an,b,iord,n,kmatch,
     >     max_sweeps) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1),iord(1)
C---------------------------------------------------------------------
C...  This subroutine does MAX_SWEEPS steps on the linear system A*x = b 
C...  by backward  Gauss-Seidel method. IORD is the ordering if any. 
C...  If no ordering then iord(i) = i. This is a special smoother
C...  done for the matching coarsening. 
C---------------------------------------------------------------------
      iteration = 0
 20   iteration = iteration + 1
      umax = -1.d20
cd      delta0 = -1d20
      do 140 i = kmatch-1,1,-2
         io1  = iord(i)
         io2  = iord(i+1)
cd         write(*,*) io1,io2
dc         read(*,*)
         iaa1 = ia(io1)
         iab1 = ia(io1+1)-1
         u1 = b(io1)    
         iaa2 = ia(io2)
         iab2 = ia(io2+1)-1
         u2 = b(io2)    
         do j = iaa1,iab1
            if(ja(j) .eq. io1) then
               a22 = an(j)
            else
               if(ja(j) .eq. io2) then
                  a12 = -an(j)
               else
                  u1 = u1 - an(j) * x(ja(j))
               end if
            end if
         end do
         do j = iaa2,iab2
            if(ja(j) .eq. io2) then
               a11 = an(j)
            else
               if(ja(j) .eq. io1) then
                  a21 = -an(j)
               else
                  u2 = u2 - an(j) * x(ja(j))
               end if
            end if
         end do
         delta = 1d0/(a11*a22-a12*a21)
cc         delta0 = dmax1(dabs(delta),delta0)
         u    = delta*(a11*u1+a12*u2)
         u2    = delta*(a21*u1+a22*u2)
         u1 = u
         ud1   = u1 - x(io1)
         ud2   = u2 - x(io2)
         umax = dmax1(umax,dabs(ud1))
         umax = dmax1(umax,dabs(ud2))
         x(io1) = u1
         x(io2) = u2
 140  continue
C
      do 40 i = n,kmatch+1,-1
         io  = iord(i)
         iaa = ia(io)
         iab = ia(io+1)-1
         if(iab .lt. iaa) go to 40
         u = b(io)    
         do 30 j = iaa,iab
            if(ja(j) .eq. io) then
               id = j
            else
               u = u - an(j) * x(ja(j))
            end if
 30      continue
         u    = u/an(id)
         ud   = u - x(io)
         umax = dmax1(umax,dabs(ud)/dabs(an(id)))
cc         umax = dmax1(umax,dabs(ud))
         x(io) = u
 40   continue
cc      if (mod(iteration,10) .eq. 1) write(*,*) iteration,umax,delta0 
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      return
      end
C====================================================================
      subroutine gauss_seidel_per(ia,ja,idir,mwk,sol,rhs,ra,rat,n,nnz)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),idir(1),mwk(1)
      dimension sol(1),rhs(1),ra(1),rat(1)
C---------------------------------------------------------------------
C...  Gauss-seidel iterative methods with Tarjan's ordering.
C...  IORDER = 1: ordering, IORDER = 0: no ordering.
C---------------------------------------------------------------------
      max_sweeps = 3*n
      iorder     = 1 
C
      if (iorder .eq. 1) then
         iaw    = 1
         jaw    = iaw + n + 1 + 1
         isub   = jaw + nnz 
         iord   = isub + n + 1 + 1
         iwork1 = iord + n + 1
         iwork2 = iwork1 + n + 1
         iwork3 = iwork2 + n + 1
         lasti  = iwork3 + n + 1 + 1
         call cut_off(
     >        n,ia,ja,ra,mwk(iaw),mwk(jaw),mwk(iord),
     >        mwk(isub),nblk,mwk(iwork1),mwk(iwork2),mwk(iwork3))
      else
         iord = 1
         call iseqv(mwk(iord),n)
      end if
C
      call permat(mwk(iord),ia,ja,ra,n,n,mwk(iaw),mwk(jaw),rat)
      call pervec(mwk(iord),rhs,sol,n)
C
      call init_guess(sol,rhs,idir,n)
C
      call outmat(ia,ja,n,202,ra,rhs)
      call fwd_gs_ord(sol,ia,ja,ra,rhs,n,max_sweeps,1)
C     
      call perback(mwk(iord),sol,rhs,n)
C
      return
      end
C=====================================================================
      subroutine init_guess(x,b,idir,n)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension idir(1),b(1),x(1)
C---------------------------------------------------------------------
C...  Initial guess for iterative methods. If there are Dirichlet 
C...  boundary nodes, then assign the initial solution to be the exact
C...  Dirichlet boundary data at those nodes, but to be zero at all 
C...  other nodes.
C---------------------------------------------------------------------
      do 10 k = 1,n
         if(idir(k) .gt. n) then
            x(k) = b(k) 
         else
            x(k) = 0.0d00
         end if
 10   continue
      return
      end
C====================================================================
      subroutine init_g_random0(ia,ja,a,b,x,n)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),a(1),x(1),b(1)
C--------------------------------------------------------------------
C     Initial guess is a sequence except on the dirichlet boundary.
C--------------------------------------------------------------------
      call nullv(x,n)
C
      do k = 1 , n
         i1 = ia(k)
         i2 = ia(k+1)-1
         if(i2-i1 .lt. 0) stop 10
         if(i2-i1 .eq. 0) then
            x(k) = b(k)
            a(i1) = 1.d0
         else
            x(k) = dble(k)
         end if
      end do
      return
      end
C=====================================================================
      subroutine fwd_gs(x,ia,ja,an,b,n,max_sweeps,iw) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by forward
C...  Gauss-Seidel method.
C...
C...  Parameter:
C...    IW - To determine whether write the convergence history or not,
C...         IW = 0 : do not write it; IW = 1 : write it.
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
cc         write(*,*) 'i,id, an(id), umax',i,id, an(id), umax
cc         umax = dmax1(umax,dabs(ud))
         x(i) = u
 40   continue
C
cc      nmodd = mod(iteration,nmodd)
      nmodd = 1
C
      if (nmodd .eq. 1 .and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' FGS : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in FGS: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine bwd_gs(x,ia,ja,an,b,n,max_sweeps,iw) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by backward
C...  Gauss-Seidel method.
C...
C...  Parameter:
C...    IW - To determine whether write the convergence history or not,
C...         IW = 0 : do not write it; IW = 1 : write it.
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
cc         write(*,*) 'i,id, an(id), umax',i,id, an(id), umax
cc         umax = dmax1(umax,dabs(ud))
         x(i) = u
 40   continue
C
cc      nmodd = mod(iteration,nmodd)
      nmodd = 1
C
      if (nmodd .eq. 1 .and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' BGS : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in BGS: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine sym_gs(x,ia,ja,an,b,n,max_sweeps,iw) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),an(1),b(1),x(1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by symmetric
C...  Gauss-Seidel method.
C...
C...  Parameter:
C...    IW - To determine whether write the convergence history or not,
C...         IW = 0 : do not write it; IW = 1 : write it.
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
cc         write(*,*) 'i,id, an(id), umax',i,id, an(id), umax
cc         umax = dmax1(umax,dabs(ud))
        x(i) = u
 60   continue
C
cc      nmodd = mod(iteration,100)
      nmodd = 1
C
      if (nmodd .eq. 1.and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' SGS : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-10) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in SGS: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
