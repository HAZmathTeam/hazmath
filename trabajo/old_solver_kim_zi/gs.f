C=====================================================================
      subroutine gauss_seidel(ia,ja,mwk,sol,rhs,a,wk,iexact,
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
      double precision soln(1)
      COMMON /DSLBLK/ SOLN
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
      if (iexact .eq. 2) then
         call l2norm(soln,soln_norm,n,1.0d0)
cc         write(*,*) ' soln_norm: ', soln_norm
         if (soln_norm .le. zero) then
            soln_norm = 1.d0
            write(*,*) 'Zero exact solution is detected!!!!.'
         end if
      end if
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
      if (iexact .eq. 2) then
         call wuminv(wk,soln,sol,n)
         call l2norm(wk,xdiffn,n,1.0d0)
         err_rel = xdiffn/soln_norm  
      end if
C
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
