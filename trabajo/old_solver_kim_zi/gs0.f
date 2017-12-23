C=====================================================================
      subroutine gauss_seidel(ia,ja,mwk,sol,rhs,a,wk,
     >     n,nnz,igs,max_it,ich_ord,tol,iwr)
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer ia(1),ja(1),mwk(1)
      integer iexact,n,nnz,igs,max_it,ich_ord,iwr
      double precision sol(1),rhs(1),a(1),wk(1),tol
C
      integer iord,iaw,jaw,isub,lasti,nblk,iter,iwr_in
      integer iwork1,iwork2,iwork3,nmod,iteration,kk
      double precision zero,soln_norm,xdiffn,arfac,arfac1
      double precision err_res0,err_relres,err_rel,err_res
C
C
C---------------------------------------------------------------------
C...  Gauss-seidel iterative methods. It uses a reordering of unknowns
C...  by Tarjan's ordering algorithm if ich_ord is 1. Forward, backward,
C...  and symmetric iterations are considered. Initial guess should be
C...  an input in the vector SOL.
C
C...  Parameters:
C...    IA,JA,A - contain the matrix data structure for A
C...    N       - the number of unknowns or the order of the Matrix A
C...    RHS     - the right-hand side vector
C...    SOL     - the initial quess as an input and solution upon return
C...    NNZ     - the number of non-zeros in the IA, JA, A storage 
C...              for the matrix A.
C...    IGS     - 1: forward G-S, 2: backward G-S, 3: sym G-S.
C...    ICH_ORD - a different ordering of unknowns is used for G-S
C...              smoothing steps if it is 1. Please note that if it 
C...              is not 1, then simply assign m(iord(k)) = 0 with
C...              array size being 1 for each level k before calling
C...              this subroutine.
C...    MAX_IT  - maximum number of sweeps
C...    TOL     - stopping criterion for iteration procedure.
C...    IEXACT  - to determine whether reading or writing the computed
C...              exact solution is considered or not;
C...              0 - no reading and writing the computed (exact) soln;
C...              1 - writing the computed soln;
C...              2 - reading the computed exact soln;
C...              3 - writing the computed exact soln;
C...              If IEXACT = 0,1,3, iteration stops when 
C...                         ||r||/||rhs|| < TOL,
C...              where r = rhs - A*sol.
C...              IEXACT = 2 is often useful for checking and comparing 
C...              different routines. For this case, the user must supply
C...              the "exact" solution or a very accurate approximation 
C...              (one with an error much less than TOL) through a common 
C...              block,  
C...                         COMMON /DSLBLK/ SOLN( )
C...              If IEXACT = 2, iteration stops when 
C...                         ||sol-soln||/||soln|| < tol.
C...              Note that this requires the user to set up the 
C...              "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.
C...              The routine with this declaration should be loaded before
C...              the stop test so that the correct length is used by the 
C...              loader.
C...    IWR     - write the convergence history or not; 1 : write it and 
C...              0 : do not write it
C...    MWK, WK - integer and real*8 working spaces, respectively
C---------------------------------------------------------------------
      if (ich_ord .eq. 1) then
         iord   = 1
         iaw    = iord + n + 1
         jaw    = iaw + n + 1 + 1
         isub   = jaw + nnz 
         iwork1 = isub + n + 1 + 1
         iwork2 = iwork1 + n + 1
         iwork3 = iwork2 + n + 1
         lasti  = iwork3 + n + 1 + 1
         call cut_off(
     >        n,ia,ja,a,mwk(iaw),mwk(jaw),mwk(iord),
     >        mwk(isub),nblk,mwk(iwork1),mwk(iwork2),mwk(iwork3))
      else
         iord   = 1
         mwk(iord) = 0
      end if
C
C...  Zero initial guess.
cc      call init_g0(ia,ja,a,rhs,sol,n)
C
      zero = 1.0d-13
      iter = 1
      iwr_in = 0
C
C...  To compute the initial residual: wk = rhs - A*sol.
      call abyvam(rhs,ia,ja,a,sol,n,wk)
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
      if (iwr .eq. 1)
     >     write(*,'(2X,i5,2X,3(3X,e11.4))') 0,1.d0,err_res0,1.d0
C
C...  Iteration starts here.
C
      iteration = 0
 20   iteration = iteration + 1
C
      if (igs. eq. 1) then
         call fwd_gs_ord(sol,ia,ja,a,rhs,mwk(iord),n,iter,iwr_in)
      else if (igs. eq. 2) then
         call bwd_gs_ord(sol,ia,ja,a,rhs,mwk(iord),n,iter,iwr_in)
      else if (igs. eq. 3) then
         call sym_gs_ord(sol,ia,ja,a,rhs,mwk(iord),n,iter,iwr_in)
      end if
C
      if (iteration .gt. 100) then
cc         nmod = mod(iteration,1000)
      else
         nmod = 1
      end if
C
C...  To compute the residual: wk = rhs - A*sol.
      call abyvam(rhs,ia,ja,a,sol,n,wk)
C
C...  To compute the residual error.
      call scpro(wk,wk,err_res,n)
      err_res = dsqrt(err_res) 
      err_relres = err_res/err_res0
C
C...  To compute the relative error if the exact solution is known.
      if (nmod .eq. 1 .and. iwr .eq. 1) then
         write(*,'(2X,i5,2X,3(3X,e11.4))')
     >        iteration,err_relres,err_res,err_rel
      end if
C
      if (iteration .gt. max_it) then
         write(*,*), ' Iteration limit exceeded:', max_it
         go to 100
      end if
C
      if (err_relres .gt. tol) go to 20
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
         write(*,'(a,2f12.5)'),'Average Reduction Factor:', arfac,arfac1
C
         write(*,*)
      end if
C
      return
      end
C=====================================================================
      subroutine init_g0(ia,ja,a,b,x,n)
C=====================================================================
cc      implicit real*8 (a-h,o-z)
      implicit none
C
      integer ia(1),ja(1),n
      double precision a(1),x(1),b(1)
C
      integer k,i1,i2
C--------------------------------------------------------------------
C...  Initial guess to be zero except on the dirichlet boundary where
C...  it is assigned to the exact value which is stored also in the
C...  right hand side load vector B.
C...  
C...  Parameters:
C...    IA,JA,A - contain the matrix data structure for A
C...    N       - the number of unknowns or the order of the Matrix A
C...    B       - the right-hand side vector
C...    X       - the initial quess
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
           x(k) = 0.d0
         end if
      end do
      return
      end
C====================================================================
      subroutine init_g_random(ia,ja,a,b,x,n)
C=====================================================================
cc      implicit real*8 (a-h,o-z)
      implicit none
C
      integer ia(1),ja(1),n
      double precision a(1),x(1),b(1)
C
      integer k,i1,i2,isx
C--------------------------------------------------------------------
C...  Initial guess randomly except on the dirichlet boundary where
C...  it is assigned to the exact value which is stored also in the
C...  right hand side load vector B.
C...  
C...  Parameters:
C...    IA,JA,A - contain the matrix data structure for A
C...    N       - the number of unknowns or the order of the Matrix A
C...    B       - the right-hand side vector
C...    X       - the initial quess
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
           x(k)= dble(k)/dble(n)
cc           x(k) = 0.d0
         end if
      end do
      return
      end
C=====================================================================
      subroutine zero_dir(x,idir,n)
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer idir(1),n
      double precision x(1)
C
      integer k
C---------------------------------------------------------------------
C...  Assign components of x to be zero if they are Dirichlet nodes.
C...  It performs same operation as the subroutine ZERO_BDRY. It is 
C...  recommended to use ZERO_BDRY rather than this subroutine.
C... 
C...  Parameters:
C...    IDIR - array which identifies Dirichlet nodes; it contains
C...           n+1 for a Dirichlet node, 0 otherwise.
C...    X    - it has zero components on Dirichlet nodes as an output
C...    N    - the number of unknowns or the order of the Matrix A
C---------------------------------------------------------------------
      do 10 k = 1,n
         if(idir(k) .gt. n) x(k) = 0.0d00
 10   continue
      return
      end
C=====================================================================
      subroutine zero_bdry(x,ia,n)
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer ia(1),n
      double precision x(1)
C
      integer k, i1, i2
C--------------------------------------------------------------------
C     Zero on the dirichlet boundary of the vector x. It performs 
C...  same operation as the subroutine ZERO_DIR. It is recommended to
C...  use this subroutine rather than ZERO_DIR because the current
C...  subroutine does not use the vector IDIR as in the ZERO_DIR.
C...  Instead it does get the information about the Dirichlet nodes
C...  from the data structure "IA" of the matrix A.
C... 
C...  Parameters:
C...    IA   - a part of the data structure of the matrix A
C...    X    - it has zero components on Dirichlet nodes as an output
C...    N    - the number of unknowns or the order of the Matrix A
CC--------------------------------------------------------------------
      do k = 1 , n
         i1 = ia(k)
         i2 = ia(k+1)-1
         if(i2-i1 .lt. 0) stop 10
         if(i2-i1 .eq. 0) x(k) = 0.0d0
      end do
      return
      end
C=====================================================================
      subroutine fwd_gs_ord(x,ia,ja,a,b,iord,n,max_it,iw) 
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer ia(1),ja(1),iord(1),n,max_it,iw
      real*8  a(1),b(1),x(1)
C
      integer nmodd,iteration,i,iaa,iab,j,id,io
      real*8  umax,u,ud
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by forward
C...  Gauss-Seidel method. It will use a rerodering of unknowns if
C...  IORD is a nonzero vector.
C...
C...  Parameter:
C...    IA,JA,A - contain the matrix data structure for A
C...    N       - the number of unknowns or the order of the Matrix A
C...    X       - initial solution as an input and updated solution as
C...              an output
C...    B       - right hand side
C...    IORD    - a different ordering of unknowns will be used if the
C...              array contains some permutation of the original ordering.
C...              Please note that if no permutation is considered, then 
C...              it should be assigned iord(1) = 0 as an input.
C...    MAX_IT  - maximum number of sweeps
C...    IW      - To determine whether write the convergence history:
C...              IW = 0 : do not write it; IW = 1 : write it.
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
                  u = u - a(j) * x(ja(j))
               end if
 30         continue
            u    = u/a(id)
            ud   = u - x(i)
            umax = dmax1(umax,dabs(ud)/dabs(a(id)))
cc          umax = dmax1(umax,dabs(ud))
            x(i) = u
 40      continue
      else 
         do 60 i = 1, n
            io = iord(i)
            iaa = ia(io)
            iab = ia(io+1)-1
            if(iab .lt. iaa) go to 60
            u = b(io)    
            do 50 j = iaa,iab
               if(ja(j) .eq. io) then
                  id = j
               else
                  u = u - a(j) * x(ja(j))
               end if
 50         continue
            u    = u/a(id)
            ud   = u - x(io)
            umax = dmax1(umax,dabs(ud)/dabs(a(id)))
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
      if (iteration .lt. max_it .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in FGS: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine bwd_gs_ord(x,ia,ja,a,b,iord,n,max_it,iw) 
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer ia(1),ja(1),iord(1),n,max_it,iw
      real*8  a(1),b(1),x(1)
C
      integer nmodd,iteration,i,iaa,iab,j,id,io
      real*8  umax,u,ud
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by backward
C...  Gauss-Seidel method. It will use a rerodering of unknowns if
C...  IORD is a nonzero vector.
C...
C...  Parameter:
C...    IA,JA,A - contain the matrix data structure for A
C...    N       - the number of unknowns or the order of the Matrix A
C...    X       - initial solution as an input and updated solution as
C...              an output
C...    B       - right hand side
C...    IORD    - a different ordering of unknowns will be used if the
C...              array contains some permutation of the original ordering.
C...              Please note that if no permutation is considered, then 
C...              it should be assigned iord(1) = 0 as an input.
C...    MAX_IT  - maximum number of sweeps
C...    IW      - To determine whether write the convergence history:
C...              IW = 0 : do not write it; IW = 1 : write it.
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
                  u = u - a(j) * x(ja(j))
               end if
 30         continue
            u    = u/a(id)
            ud   = u - x(i)
            umax = dmax1(umax,dabs(ud)/dabs(a(id)))
cc          umax = dmax1(umax,dabs(ud))
            x(i) = u
 40      continue
      else
         do 60 i = n, 1, -1
            io = iord(i)
            iaa = ia(io)
            iab = ia(io+1)-1
            if(iab .lt. iaa) go to 60
            u = b(io)    
            do 50 j = iaa,iab
               if(ja(j) .eq. io) then
                  id = j
               else
                  u = u - a(j) * x(ja(j))
               end if
 50         continue
            u    = u/a(id)
            ud   = u - x(io)
            umax = dmax1(umax,dabs(ud)/dabs(a(id)))
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
      if (iteration .lt. max_it .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in BGS: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine sym_gs_ord(x,ia,ja,a,b,iord,n,max_it,iw) 
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer ia(1),ja(1),iord(1),n,max_it,iw
      real*8  a(1),b(1),x(1)
C
      integer nmodd,iteration,i,iaa,iab,j,id,io
      real*8  umax,u,ud
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by symmetric
C...  Gauss-Seidel method. It will use a rerodering of unknowns if
C...  IORD is a nonzero vector.
C...
C...  Parameter:
C...    IA,JA,A - contain the matrix data structure for A
C...    N       - the number of unknowns or the order of the Matrix A
C...    X       - initial solution as an input and updated solution as
C...              an output
C...    B       - right hand side
C...    IORD    - a different ordering of unknowns will be used if the
C...              array contains some permutation of the original ordering.
C...              Please note that if no permutation is considered, then 
C...              it should be assigned iord(1) = 0 as an input.
C...    MAX_IT  - maximum number of sweeps
C...    IW      - To determine whether write the convergence history:
C...              IW = 0 : do not write it; IW = 1 : write it.
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
                  u = u - a(j) * x(ja(j))
               end if
 30         continue
            x(i) = u/a(id)
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
                  u = u - a(j) * x(ja(j))
               end if
 50         continue
            u    = u/a(id)
            ud   = u - x(i)
            umax = dmax1(umax,dabs(ud)/dabs(a(id)))
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
            if(iab .lt. iaa) go to 80
            u = b(io)    
            do 70 j = iaa,iab
               if(ja(j) .eq. io) then
                  id = j
               else
                  u = u - a(j) * x(ja(j))
               end if
 70         continue
            x(io) = u/a(id)
 80      continue
C
C...     Backward Gauss-Seidel
C  
         do 100 i = n, 1, -1
            io = iord(i)
            iaa = ia(io)
            iab = ia(io+1)-1
            if(iab .lt. iaa) go to 100
            u = b(io)    
            do 90 j = iaa,iab
               if(ja(j) .eq. io) then
                  id = j
               else
                  u = u - a(j) * x(ja(j))
               end if
 90         continue
            u    = u/a(id)
            ud   = u - x(io)
            umax = dmax1(umax,dabs(ud)/dabs(a(id)))
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
      if (iteration .lt. max_it .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in SGS: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine jacobi(x,ia,ja,a,b,n,xold,omega,max_it,iw) 
C=====================================================================
cc      implicit real*8(a-h,o-z)
      implicit none
C
      integer ia(1),ja(1),n,max_it,iw
      real*8  a(1),b(1),x(1),xold(1),omega
C
      integer nmodd,iteration,i,iaa,iab,j,id
      real*8  umax,u,ud
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by Jacobi method.
C...
C...  Parameter:
C...    IA,JA,A - contain the matrix data structure for A
C...    N       - the number of unknowns or the order of the Matrix A
C...    X       - initial solution as an input and updated solution as
C...              an output
C...    XOLD    - initial solution X will be copied to this vector
C...    B       - right hand side
C...    MAX_IT  - maximum number of sweeps
C...    IW      - To determine whether write the convergence history:
C...              IW = 0 : do not write it; IW = 1 : write it.
C---------------------------------------------------------------------
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
               u = u - a(j) * xold(ja(j))
            end if
 30      continue
         u    = u/a(id)
         ud   = u - x(i)
         umax = dmax1(umax,dabs(ud)/dabs(a(id)))
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
      if (iteration .lt. max_it .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST in JACOBI: ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
C====================================================================
      subroutine scpro(u,v,scpr,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n)
C--------------------------------------------------------------------
C...  INNER PRODUCT
C--------------------------------------------------------------------
cc      scpr = ddot_kz(n,u,1,v,1)
      scpr = 0.0d00
      do i = 1 , n
         scpr = scpr + u(i)*v(i)
      end do
C
      return
      end
C====================================================================
      subroutine uuminv(u,v,n)
C====================================================================
      implicit real*8(a-h,o-z)
C--------------------------------------------------------------------
C...  U <--- U - V
C--------------------------------------------------------------------
      dimension u(n),v(n)
      do i = 1 , n
         u(i) = u(i) - v(i)
      end do
      return
      end
C====================================================================
      subroutine vuminv(u,v,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n)
C--------------------------------------------------------------------
C...  V <--- U - V
C--------------------------------------------------------------------
      do i = 1 , n
         v(i) = u(i) - v(i)
      end do
      return
      end
C====================================================================
      subroutine uupluv(u,v,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n)
C--------------------------------------------------------------------
C...  U <--- U + V
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = u(i) + v(i)
      end do
      return
      end
C====================================================================
      subroutine uuplmv(u,v,n,smult)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n)
C--------------------------------------------------------------------
C...  U <--- U + SMULT*V
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = u(i) + v(i)*smult
      end do
      return
      end
C====================================================================
      subroutine umuplv(u,v,n,smult)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n)
C--------------------------------------------------------------------
C...  U <--- SMULT*U  + V
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = smult*u(i) + v(i)
      end do
      return
      end
C====================================================================
      subroutine usmultv(u,v,n,smult)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n)
C--------------------------------------------------------------------
C...  U <--- SMULT*V
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = v(i)*smult
      end do
      return
      end
C====================================================================
      subroutine usmultu(u,n,smult)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n)
C--------------------------------------------------------------------
C...  U <--- SMULT*U  
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = smult*u(i)
      end do
      return
      end
C====================================================================
      subroutine wuminv(w,u,v,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n),w(n)
C--------------------------------------------------------------------
C...  W <--- U - V
C--------------------------------------------------------------------
      do i = 1 , n
         w(i) = u(i) - v(i)
      end do
      return
      end
C====================================================================
      subroutine wupluv(w,u,v,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n),w(n)
C--------------------------------------------------------------------
C...  W <--- U + V
C--------------------------------------------------------------------
      do i = 1 , n
         w(i) = u(i) + v(i)
      end do
      return
      end
C====================================================================
      subroutine wuplmv(w,u,v,n,smult)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n),w(n)
C--------------------------------------------------------------------
C...  W <--- U + SMULT*V
C--------------------------------------------------------------------
      do i = 1 , n
         w(i) = u(i) + v(i)*smult
      end do
      return
      end
C====================================================================
      subroutine l2norm(u,l2nr,n,h)
C====================================================================
      implicit real*8(a-h,l,o-z)
      dimension u(n)
C--------------------------------------------------------------------
C...  ||U||_{0} (FOR TWO DIMENSIONS)
C--------------------------------------------------------------------
      call scpro(u,u,l2nr,n)
      l2nr = dsqrt(l2nr)*h
      return
      end
C====================================================================
      subroutine c0norm(u,c0nr,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n)
C--------------------------------------------------------------------
C...  ||U||_{\infty}
C--------------------------------------------------------------------
      c0nr = -1.d00
      do i = 1 , n
         c0nr = dmax1(c0nr,dabs(u(i)))
      end do   
      return
      end
C====================================================================
      subroutine c0norm_number(u,c0nr,n,i1)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n)
C--------------------------------------------------------------------
C...  ||U||_{\infty}
C--------------------------------------------------------------------
      c0nr = -1.d00
      i1 = 0
      do i = 1 , n
         if(c0nr .lt. dabs(u(i))) then
            c0nr = dabs(u(i))
            i1 = i
         end if
      end do   
      return
      end
C====================================================================
      subroutine copyv(u,v,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n)
C--------------------------------------------------------------------
C...  V <--- U 
C--------------------------------------------------------------------
      do i = 1 , n
         v(i) = u(i)
      end do
      return
      end
C====================================================================
      subroutine icopyv(iu,iv,n)
C====================================================================
      dimension iu(n),iv(n)
C--------------------------------------------------------------------
C...  IV <--- IU
C--------------------------------------------------------------------
      do i = 1 , n
         iv(i) = iu(i)
      end do
      return
      end
C====================================================================
      subroutine nullv(u,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n)
C--------------------------------------------------------------------
C...  U = 0.
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = 0.0d00
      end do
      return
      end
C====================================================================
      subroutine onesv(u,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n)
C--------------------------------------------------------------------
C...  U = [1,1,...,1]
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = 1.0d00
      end do
      return
      end
C====================================================================
      subroutine inullv(iu,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension iu(n)
C--------------------------------------------------------------------
C...  IU = 0
C--------------------------------------------------------------------
      do i = 1 , n
         iu(i) = 0
      end do
      return
      end
C====================================================================
      subroutine ionev(iu,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension iu(n)
C--------------------------------------------------------------------
C...  IU = [1,1,...,1]
C--------------------------------------------------------------------
      do i = 1 , n
         iu(i) = 1
      end do
      return
      end
C====================================================================
      subroutine innv(iu,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension iu(n)
C--------------------------------------------------------------------
C...  IU = [N,N,...,N]
C--------------------------------------------------------------------
      do i = 1 , n
         iu(i) = n
      end do
      return
      end
C====================================================================
      subroutine dbnv(u,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n)
C--------------------------------------------------------------------
C...  U = [N,N,...,N]
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = dble(n)
      end do
      return
      end
C====================================================================
      subroutine iseqv(iu,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension iu(n)
C--------------------------------------------------------------------
C...  IU = [1,2,...,N]
C--------------------------------------------------------------------
      do i = 1 , n
         iu(i) = i
      end do
      return
      end
C====================================================================
      subroutine iseqrv(iu,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension iu(n)
C--------------------------------------------------------------------
C...  IU = [N,N-1,...,1]
C--------------------------------------------------------------------
      n1 = n - 1
      do i = 0, n1
         iu(i) = n - i
      end do
      return
      end
C====================================================================
      subroutine seqv(u,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n)
C--------------------------------------------------------------------
C...  U = [1,2,...,N]
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = dble(i)
      end do
      return
      end
C====================================================================
      subroutine seqrv(u,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n)
C--------------------------------------------------------------------
C...  U = [N,N-1,...,1]
C--------------------------------------------------------------------
      n1 = n - 1
      do i = 0, n1
         u(i) = dble(n - i)
      end do
      return
      end
C=====================================================================
      subroutine ireverse(iu,n)
C====================================================================
      implicit real*8(a-h,o-z)
      integer iu(1),n,k,iu0
C--------------------------------------------------------------------
C...  Reverse the numbering of the vector IU.
C--------------------------------------------------------------------
      do k =  1 , n/2
         iu0 =  iu(n-k+1)
         iu(n-k+1) = iu(k)
         iu(k) = iu0
      end do
      return
      end
C====================================================================
      subroutine reverse(u,n)
C====================================================================
      implicit real*8(a-h,o-z)
      integer n,k
      dimension u(n)
C--------------------------------------------------------------------
C...  Reverse the numbering of the vector U.
C--------------------------------------------------------------------
      do k =  1 , n/2
         u0 =  u(n-k+1)
         u(n-k+1) = u(k)
         u(k) = u0
      end do
      return
      end
C====================================================================
      real*8 function ddot_kz(n,dx,incx,dy,incy)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension dx(1),dy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
C---------------------------------------------------------------------
C...  Forms the dot product of two vectors.
C...  Uses unrolled loops for increments equal to one.
C---------------------------------------------------------------------
      ddot_kz = 0.0d0
      dtemp = 0.0d0
      if(n.le.0) return
      if(incx.eq.1.and.incy.eq.1) go to 20
C
C...  Code for unequal increments or equal increments not equal to 1.
C
      ix = 1
      iy = 1
      if(incx.lt.0) ix = (-n+1)*incx + 1
      if(incy.lt.0) iy = (-n+1)*incy + 1
      do 10 i = 1,n
         dtemp = dtemp + dx(ix)*dy(iy)
         ix = ix + incx
         iy = iy + incy
 10   continue
      ddot_kz = dtemp
      return
C
C...  Code for both increments equal to 1. Clean-up loop
C
 20   m = mod(n,5)
      if(m .eq. 0) go to 40
      do 30 i = 1,m
         dtemp = dtemp + dx(i)*dy(i)
 30   continue
      if(n .lt. 5) go to 60
 40   mp1 = m + 1
      do 50 i = mp1,n,5
         dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1)
     >         + dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
 50   continue
 60   ddot_kz = dtemp
C
      return
      end
C=====================================================================
      subroutine interp_func(funct0,u,n,x,y)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension u(1),x(1),y(1)
      external funct0
C---------------------------------------------------------------------
C...  The linear interpolation of the given function FUNCT0.
C---------------------------------------------------------------------
      do k = 1 , n
         u(k) = funct0(x(k),y(k))
      end do
      return
      end
C====================================================================
      subroutine wumultv(w,u,v,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension w(n),u(n),v(n)
C--------------------------------------------------------------------
C...  Vector PRODUCT componentwise. w = u:v
C--------------------------------------------------------------------
      do i = 1 , n
         w(i) = u(i)*v(i)
      end do
C
      return
      end
C====================================================================
      subroutine uumultv(u,v,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension u(n),v(n)
C--------------------------------------------------------------------
C...  Vector PRODUCT componentwise. u = u:v
C--------------------------------------------------------------------
      do i = 1 , n
         u(i) = u(i)*v(i)
      end do
C
      return
      end
C=====================================================================
      integer*4 function nomxy(i,j,nx,ny)
C======================================================================
      implicit none
      integer*4 i,j,nx,ny
C--------------------------------------------------------------------
C...  NOMXY gives the global number of the node (i,j) in an NX x NY
C...  structured grid. Note that global numbers are enumerated 
C...  y-direction first from bottom to top and increased in 
C...  x-direction from left to right.
C--------------------------------------------------------------------
      if(i .lt. 1 .or. j .lt. 1 .or. 
     >     i .gt. nx .or. j .gt. ny  ) then
         nomxy = 0
      else
cc       nomxy = (j-1)*nx + i
         nomxy = (i-1)*ny + j
      end if
C
      RETURN
      END
C=====================================================================
      subroutine nomxy_inv(i,j,n,ny)
C======================================================================
      implicit none
      integer*4 i,j,n,ny
C--------------------------------------------------------------------
C...  NOMXY_INV gives the node (i,j) for the global number n in an
C...  NX x NY structured grid. Note that global numbers are enumerated
C...  y-direction first from bottom to top and increased in x-direction
C...  from left to right.
C--------------------------------------------------------------------
      j = mod(n,ny)
      if (j .eq. 0) then
         j = ny
         i = n/ny
      else
         i = n/ny + 1
      end if
C
      RETURN
      END
C======================================================================
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
