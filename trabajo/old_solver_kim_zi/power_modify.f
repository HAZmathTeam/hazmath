C=====================================================================
      subroutine power_new(m,iao,jao,idiro,r,sol,rhs,ao,
     >     n,nnz,maxit,tol,iexact,nel,jdf,ie,je,x,y,
     >     iao_p,jao_p,idiro_p,kal,ksol,krhs,kra) 
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1)
      dimension r(1)
C
      integer iao(1),jao(1),idiro(1)
      double precision sol(1),rhs(1),ao(1)
C
      integer lmax
      parameter (lmax = 20)
      integer nd, n1, n2, nz, idir, 
     >     ia, ja, ka, kag,
     >     ipp, jpp, nzp, 
     >     ku, kb, iperm, iord, 
     >     ic, jc, kc, ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ic, jc, kc, ifree, kfree, lf, lc
C
      dimension ie(1),je(1),x(1),y(1)
      external msolve
C
      common /prenormal/ ipre_normal
C---------------------------------------------------------------------
C...  Compute the norm of the operator T = I - BA, i.e., |||T|||.
C...  Here A is the EAFE scheme and B is the mg V-cycle on EAFE scheme.
C...  Let C = A_0 \bar{A}^(-t) T^t \bar{A}^t A_0^(-1) \bar{A} T \bar{A}^(-1)
C...  For this purpose, we use a power method described below:
C...    1.  x = a random vector for initial guess with ||x|| = 1 
C...    2.  v = a random vector such that (x, A_0v) is not zero, where
C...            A_0 is the Laplacian
C...    3.  w = A_0v
C...    4.  xw = (x, w)
C...    5.  Iterate until convergence
C...    6.     y = Cx
C...    7.     yw = (y, w)
C...    8.     r = yw / xw
C...    9.     r = sqrt(r)
C...   10.     yn = ||y||
C...   11.     x = y / yn
C...   12.     xw = yw / yn
C...   13.  End iterate 
C...   14.  r is the desired norm
C...
C...  Parameters:
C...    TOL    - stopping iteration procedure.
C---------------------------------------------------------------------
      kv    = kfree
      kw    = kv + n
      kasave = kw + n
      kfree = kasave + nnz
C
      iasave = ifree
      jasave = iasave + n + 1
      ifree = jasave + nnz
C     
C...  Temporary save for iao, jao, ao
C
      call icopyv(iao,m(iasave),n+1)
      call icopyv(jao,m(jasave),nnz)
      call copyv(ao,r(kasave),nnz)
C
      ifree_old = ifree
      kfree_old = kfree
C
C...  Initial guess randomly for eigenvector.
C
      call init_g_random(iao,jao,ao,rhs,sol,n)
cc      write(*,*) (sol(ii), ii = 1,n)
cc      write(*,*)
C
C...  Normalize sol, i.e., sol = sol / || sol ||
C
      call scpro(sol,sol,sol_norm,n)
      sol_norm = 1.0d0 / dsqrt(sol_norm)
      call usmultu(sol,n,sol_norm)
C
C...  Choose v, a random vector such that (sol, A_0v) is not zero, where
C...  A_0 is the Laplacian. 
C
      call copyv(sol,r(kv),n)
C
cc      call init_g_random00(iao,jao,ao,rhs,r(kv),n)
cc      write(*,*) (r(kv+ii-1), ii = 1,n)
C
C...  Compute the inner product, (sol, A_0v), to get x_dot_av.
      call abyvg(iao,jao,r(kal),r(kv),n,r(kw))
      call scpro(sol,r(kw),x_dot_av,n)
C
C...  One iteration of MG V-cycle to compute w = T*v. (This belongs
C...  to an old version.)
C
      ipre_normal = 0 !NO smoothings by G-S sweeps for the normal eq.
      mg_format = 1
      maxit_1 = 100
      tol_1   = 1.0d-11
      maxit_mg = 1
      tol_mg = -1.0d0
C
      r_new = 1.0d0
C
C...  Iterate until desired tolerance is achieved.
C
      do 100 kk = 1, maxit
C
         call copyv(sol,r(kv),n)
C
C======= Beginning of the action sol = C*v. =============================
C
C...     Step 1. Solve sol = \bar{A}^(-1)*v.
C
         call nullv(sol,n)
C
C...     Solve by MG or MG_normal.
C...     The following two subroutines are identical.
C
         if (mg_format .eq. 0) then
            call mg(m(1),iao,jao,idiro,
     >           r(1),sol,r(kv),ao,
     >           n,nnz,maxit_1,tol_1,iexact)
         else
            call premg_normal(m(1),iao,jao,idiro,
     >           r(1),sol,r(kv),ao,
     >           n,nnz,maxit_1,tol_1,iexact,msolve)
         end if
C
cc         write(100,*) 'sol 1',(sol(ii),ii=1,n)
C
         ifree = ifree_old
         kfree = kfree_old
C
C...     Step 2. One iteration of MG V-cycle to compute sol = T*sol.
C...     The following two subroutines are identical.
C
         if (mg_format .eq. 0) then
            call mg(m,iao,jao,idiro,
     >           r,sol,rhs,ao,
     >           n,nnz,maxit_mg,tol_mg,iexact)
         else
            call premg_normal(m,iao,jao,idiro,
     >           r,sol,rhs,ao,
     >           n,nnz,maxit_mg,tol_mg,iexact,msolve)
         end if
C
         ifree = ifree_old
         kfree = kfree_old
C
C...     Step 3. Compute v = \bar{A}*sol.
C
         call abyvg(iao,jao,ao,sol,n,r(kv))
C
C...     Step 4. Compute sol = A_0^(-1)*v by MG V-cycle. 
C
         call mg_s(m,iao_p,jao_p,idiro_p,ifree_old,
     >        r,ksol,kv,kal,kfree_old,
     >        n,n1(lf),n2(lf),nnz,nel,lf,1,tol_1) 
C
C...     Step 5. Compute v = \bar{A}^t*sol.
C
         call vbya(iao,jao,ao,sol,n,n,r(kv))
C
C...     Step 6. Compute the action sol = (T^*)v with zero initial quess.
C
         max_itstar = 1
         do k = 1,  maxit_mg
            call nullv(sol,n)
            call mg_star(m,iao,jao,idiro,
     >           r,sol,r(kv),ao,
     >           n,nnz,max_itstar,tol_mg,iexact,1)
C
            ifree = ifree_old
            kfree = kfree_old
C     
C...        Copy sol --> v.
            call copyv(sol,r(kv),n)
C     
            call icopyv(m(iasave),iao,n+1)
            call icopyv(m(jasave),jao,nnz)
            call copyv(r(kasave),ao,nnz)
         end do
C
C...     Step 7. Solve sol = \bar{A}^(-t)*v.
C           
C...     Compute at = A^t.
         call aat(m(iasave),m(jasave),r(kasave),n,n,iao,jao,ao)
C
         call nullv(sol,n)
C
C...     Solve by MG or MG_normal.
C
         call mg_tr1(m,iao,jao,idiro,
cc         call mg_tr(m,iao,jao,idiro,
     >        r,sol,r(kv),ao,
     >        n,nnz,maxit_1,tol_1,iexact)
C
         call icopyv(m(iasave),iao,n+1)
         call icopyv(m(jasave),jao,nnz)
         call copyv(r(kasave),ao,nnz)
C
         ifree = ifree_old
         kfree = kfree_old
C
C...     Step 8. Compute v = A_0*sol.
C
         call abyvg(iao,jao,r(kal),sol,n,r(kv))
         call copyv(r(kv),sol,n)
C
C======= End of the action sol = C*v. ==================================
C
C...     Compute the inner product, (y, A_0v), to get y_dot_av.
C
         call scpro(sol,r(kw),y_dot_av,n)
C
         write(*,*) 'x dot av : ', kk, '  :::  ', x_dot_av
         write(*,*) 'y dot av : ', kk, '  :::  ', y_dot_av
C
         r_old = r_new 
C
C...     Compute the ratio: r_new = y_dot_av / x_dot_av
C
         r_new = y_dot_av / x_dot_av
cc         r_new = dsqrt(r_new)
         write(*,*) '   Norm :  ', kk, '  :::  ', r_new
C
         if (dabs(r_new - r_old) / dabs(r_new) .lt. tol) go to 200
C
C...     Normalize sol, i.e., sol = sol / || sol ||
C
         call scpro(sol,sol,x_dot_av,n)
         x_dot_av = 1.0d0 / dsqrt(x_dot_av)
         call usmultv(sol,sol,n,x_dot_av)
C
C...     Here is the new x_dot_av.
         x_dot_av = y_dot_av * x_dot_av
C
 100  continue
C
cc      write(*,*) 'Maximum iteration count occurs at Power_Modified.f.'
C
 200  continue
C
      r_new = dsqrt(r_new)
      write(*,*) 'Final Norm ', kk, '  :::  ', r_new
C
      return
      end
C=====================================================================
      subroutine power_modify(m,iao,jao,idiro,r,sol,rhs,ao,
     >     n,nnz,maxit,tol,iexact,nel,jdf,ie,je,x,y,
     >     iao_p,jao_p,idiro_p,kal,ksol,krhs,kra)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1)
      dimension r(1)
C
      integer iao(1),jao(1),idiro(1)
      double precision sol(1),rhs(1),ao(1)
C
      integer lmax
      parameter (lmax = 20)
      integer nd, n1, n2, nz, idir, 
     >     ia, ja, ka, kag,
     >     ipp, jpp, nzp, 
     >     ku, kb, iperm, iord, 
     >     ic, jc, kc, ifree, kfree, lf, lc
      common /pointer/ 
     >     nd(lmax), n1(lmax), n2(lmax), nz(lmax), idir(lmax), 
     >     ia(lmax), ja(lmax), ka(lmax), kag(lmax), 
     >     ipp(lmax),jpp(lmax), nzp(lmax),
     >     ku(lmax), kb(lmax), iperm(lmax), iord(lmax),
     >     ic, jc, kc, ifree, kfree, lf, lc
C
      dimension ie(1),je(1),x(1),y(1)
      external msolve
C
      common /prenormal/ ipre_normal
C---------------------------------------------------------------------
C...  Compute the H^1 norm of the operator T = I - BA, i.e., ||T||_1.
C...  Here A is the EAFE scheme and B is the mg V-cycle on EAFE scheme.
C...  For this purpose, we use a power method described below:
C...    1.  x = a random vector for initial guess with ||x|| = 1 
C...    2.  v = a random vector such that (x, A_0v) is not zero, where
C...            A_0 is the Laplacian
C...    3.  w = A_0v
C...    4.  xw = (x, w)
C...    5.  Iterate until convergence
C...    6.     y = Cx = [(A_0^(-1))(T^t)(A_0)(T)]x
C...    7.     yw = (y, w)
C...    8.     r = yw / xw
C...    9.     r = sqrt(r)
C...   10.     yn = ||y||
C...   11.     x = y / yn
C...   12.     xw = yw / yn
C...   13.  End iterate 
C...   14.  r is the desired norm
C...
C...  Parameters:
C...    TOL    - stopping iteration procedure.
C---------------------------------------------------------------------
      kv    = kfree
      kw    = kv + n
      kasave = kw + n
      kfree = kasave + nnz
C
      iasave = ifree
      jasave = iasave + n + 1
      ifree = jasave + nnz
C
      ifree_old = ifree
      kfree_old = kfree
C
C...  Initial guess randomly for eigenvector.
C
      call init_g_random(iao,jao,ao,rhs,sol,n)
cc      write(*,*) (sol(ii), ii = 1,n)
cc      write(*,*)
C
C...  Normalize sol, i.e., sol = sol / || sol ||
C
      call scpro(sol,sol,sol_norm,n)
      sol_norm = 1.0d0 / dsqrt(sol_norm)
      call usmultu(sol,n,sol_norm)
C
C...  Choose v, a random vector such that (sol, A_0v) is not zero, where
C...  A_0 is the Laplacian. 
C
      call copyv(sol,r(kv),n)

cc      call init_g_random00(iao,jao,ao,rhs,r(kv),n)
cc      write(*,*) (r(kv+ii-1), ii = 1,n)

ccC...  Compute the inner product, (sol, A_0v), to get x_dot_av.
cc      call l2_h1s_product(1,sol,r(kv),n,x_dot_av,nel,jdf,ie,je,x,y)
cc      write(*,*) 'initial x dot av : ', x_dot_av
C
C...  Compute the inner product, (sol, A_0v), to get x_dot_av.
      call abyvg(iao,jao,r(kal),r(kv),n,r(kw))
      call scpro(sol,r(kw),x_dot_av,n)
C
C...  One iteration of MG V-cycle to compute w = T*v. (This belongs
C...  to an old version.)
C
      ipre_normal = 0 !NO smoothings by G-S sweeps for the normal eq.
      mg_format = 1
      max_it = 1
      tol_mg = -1.0d0
C
C...  Iterate until desired tolerance is achieved.
C
      r_new = 1.0d0
C
      do 100 kk = 1, maxit
C
C...     One iteration of MG V-cycle to compute sol = T*sol.
C...     The following two subroutines are identical.
C
         if (mg_format .eq. 0) then
            call mg(m,iao,jao,idiro,
     >           r,sol,rhs,ao,
     >           n,nnz,max_it,tol_mg,iexact)
         else
            call premg_normal(m,iao,jao,idiro,
     >           r,sol,rhs,ao,
     >           n,nnz,max_it,tol_mg,iexact,msolve)
         end if
C
         ifree = ifree_old
         kfree = kfree_old
C
ccC...     Compute the inner product, (sol, A_0Tv), to get y_dot_av.
cc         call l2_h1s_product(1,sol,r(kw),n,y_dot_av,nel,jdf,ie,je,x,y)
C
C...     Compute v = A_0*sol
C
         call abyvg(iao,jao,r(kal),sol,n,r(kv))
C
C...     Temporary save for iao, jao, ao
C
         call icopyv(iao,m(iasave),n+1)
         call icopyv(jao,m(jasave),nnz)
         call copyv(ao,r(kasave),nnz)
C
C...     Compute the action (T^*)v here with zero initial quess.
C
         max_itstar = 1
         do k = 1,  max_it
            call nullv(sol,n)
            call mg_star(m,iao,jao,idiro,
     >           r,sol,r(kv),ao,
     >           n,nnz,max_itstar,tol_mg,iexact,1)
            ifree = ifree_old
            kfree = kfree_old
C     
C...        Copy sol --> v.
C     
            call copyv(sol,r(kv),n)
C     
            call icopyv(m(iasave),iao,n+1)
            call icopyv(m(jasave),jao,nnz)
            call copyv(r(kasave),ao,nnz)
C     
         end do
C
C...     Compute y (= sol) = A_0^(-1)*v by MG V-cycle. 
C
         call mg_s(m,iao_p,jao_p,idiro_p,ifree_old,
     >        r,ksol,kv,kal,kfree_old,
     >        n,n1(lf),n2(lf),nnz,nel,lf,1,1.0d-11) 
C
C...     Compute the inner product, (y, A_0v), to get y_dot_av.
C
         call scpro(sol,r(kw),y_dot_av,n)
C
         write(*,*) 'x dot av : ', kk, '  :::  ', x_dot_av
         write(*,*) 'y dot av : ', kk, '  :::  ', y_dot_av
C
         r_old = r_new 
C
C...     Compute the ratio: r_new = y_dot_av / x_dot_av
C
         r_new = y_dot_av / x_dot_av
cc         r_new = dsqrt(r_new)
         write(*,*) '   Norm :  ', kk, '  :::  ', r_new
C
         if (dabs(r_new - r_old) / dabs(r_new) .lt. tol) go to 200
C
C...     Normalize sol, i.e., sol = sol / || sol ||
C
         call scpro(sol,sol,x_dot_av,n)
         x_dot_av = 1.0d0 / dsqrt(x_dot_av)
         call usmultv(sol,sol,n,x_dot_av)
C
C...     Here is the new x_dot_av.
         x_dot_av = y_dot_av * x_dot_av
C
ccC...     Compute the inner product, (sol, A_0v), to get x_dot_av.
cc         call l2_h1s_product(1,sol,r(kv),n,x_dot_av,nel,jdf,ie,je,x,y)
C
 100  continue
C
cc      write(*,*) 'Maximum iteration count occurs at Power_Modified.f.'
C
 200  continue
C
      r_new = dsqrt(r_new)
      write(*,*) 'Final Norm ', kk, '  :::  ', r_new
C
      return
      end
C=====================================================================
