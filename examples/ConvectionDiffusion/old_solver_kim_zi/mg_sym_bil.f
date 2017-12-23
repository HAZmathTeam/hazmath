C=====================================================================
      subroutine mg_s_bil(m,iao,jao,idiro,ifree,r,ksol,krhs,kra,kfree,
     >     n,nx,ny,nnz,nel,lf,lc,tol) 
C=====================================================================
      parameter (lmx = 20)
      implicit real*8(a-h,o-z)
      dimension m(1),ia(lmx),ja(lmx),idir(lmx),nd(lmx),nz(lmx)
      dimension ipp(lmx),jpp(lmx),nzp(lmx)
      dimension r(1),ku(lmx),kb(lmx),ka(lmx),n1(lmx),n2(lmx)
      character*100 meshf(0:lmx)
      common /menu/ ich(200)
      common /fname/ meshf
C---------------------------------------------------------------------
C...  Standard multigrid cycle for the SPD problems and bilinear element
C...  space is used.
C...
C...  Parameter:
C...    N,NX,NY - dimension of the mesh, N = NX x NY
C...    ICYCLE  - choice of MG:  1 = V-cycle, 2 = \-cycle, 
C...              3 = Variable V-cycle,(4 = W-cylce, 5 = Full MG)
C...    LF      - level of the finest mesh: 2,3,4, etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C...    MAXIT   - maximum number of MG iterations
C...    TOL     - stopping criterion for iteration procedure.
C---------------------------------------------------------------------
      mt     = 15
      icycle = ich(31)
      maxit  = ich(39)
C
cc      tol    = 0.5d-11
      zero   = 1.d-14
C
C...  Generating the stiffness matrix A and the load vector B at each level.
C
      lasti = ifree
      lastr = kfree
      nelnew = nel
      nd(lf) = n
      n1(lf) = nx
      n2(lf) = ny
C
      do 50 i = lf, lc, -1
         if (i .eq. lf) then
            nz(i)   = nnz
            ia(i)   = iao
            ja(i)   = jao
            idir(i) = idiro
            ipp(i)  = lasti
         else
            nelnew  = nelnew/4
            ia(i)   = lasti
            ja(i)   = ia(i) + nd(i) + 1
            idir(i) = ja(i) + 2*(nelnew + nd(i) + 5) + nd(i)
            ipp(i)  = idir(i) + nd(i) + 5
         end if
C
         jpp(i)  = ipp(i) + nd(i) + 1
         lasti   = jpp(i) + 2*nd(i)
C
         ku(i)   = lastr
         kb(i)   = ku(i) + nd(i)
         if (i .eq. lf) then
            ka(i)   = kra
            lastr   = kb(i) + nd(i)
         else
            ka(i)   = kb(i) + nd(i)
            lastr   = ka(i) + 2*(nelnew + nd(i) + 5) + nd(i)
         end if
C
C...     To compute the symbolic prolongation matrix P from level i-1
C...     to level i.
C
         if (i .eq. lc) go to 40
         call prolng_symb_bil(n1(i),n2(i),nd(i),n1(i-1),n2(i-1),
     >        nd(i-1),m(ipp(i)),m(jpp(i)),nzp(i))
 40      continue
C
C...     To compute stiffness matrices A at each level.
C
         if (i .eq. lf) then
            ia(i) = iao
            ja(i) = jao
            ka(i) = kra
            nz(i) = nnz
         else
            call formdir(m(idir(i+1)),m(ipp(i+1)),m(jpp(i+1)),
     >           nd(i+1),m(idir(i)),nd(i))
C
            ic = lasti
            ix = ic + nd(i+1) + 1
            jc = ix + nd(i) + 1

            call abybs(m(ia(i+1)),m(ja(i+1)),m(ipp(i+1)),m(jpp(i+1)),
     >           nd(i+1),nd(i+1),nd(i),m(ic),m(jc),m(ix))
C
            jce = m(ic+nd(i+1)) + 1
            kc  = lastr
            kx  = kc + jce
            call abyp_bil(m(ia(i+1)),m(ja(i+1)),m(ipp(i+1)),
     >           m(jpp(i+1)),
     >           nd(i+1),m(ic),m(jc),r(ka(i+1)),r(kc),r(kx),nd(i+1))
C
            ict = jc + jce + 1
            jct = ict + nd(i+1) + 1
            kct = kx + nd(i+1) + 1
            call aat(m(ic),m(jc),r(kc),nd(i+1),nd(i),m(ict),m(jct),
     >           r(kct))
c
            call abybs(m(ict),m(jct),m(ipp(i+1)),m(jpp(i+1)),nd(i),
     >           nd(i+1),nd(i),m(ia(i)),m(ja(i)),m(ix))
C     
            call abyp_bil(m(ict),m(jct),m(ipp(i+1)),m(jpp(i+1)),
     >           nd(i),
     >           m(ia(i)),m(ja(i)),r(kct),r(ka(i)),r(kx),nd(i+1))
C     
            nz(i) = m(ia(i)+nd(i))
C
            call adir_new(m(ia(i)),m(ja(i)),r(ka(i)),m(idir(i)),nd(i))
C
cc            call lprr(m(ia(i)),m(ja(i)),r(ka(i)),nd(i))
cc            call outmat1(m(ia(i)),m(ja(i)),r(ka(i)),nd(i),nz(i))
         end if
 50   continue
C
C...  ------------------- MG starts here. ---------------------------
C
cc      write(*,*)
cc      write(*,*) '=============================', 
cc     >           '=================================================='
cc      if (icycle .eq. 1) then
cc         write(*,*) '     METHOD :          MG : V-cycle '
cc      else if (icycle .eq. 2) then
cc         write(*,*) '     METHOD :          MG : Backslash cycle '
cc      else if (icycle .eq. 3) then
cc         write(*,*) '     METHOD :          MG : Variable V-cycle '
cc      end if
ccC
cc      write(*,*) '     Residual :        r := b-Au; '
cc      write(*,'(a,a,e12.4)')
cc     >     '      Convergence :     ||r||_0/||r_0||_0  &',
cc     >     '  ||Br||_0/||u||_0 <', tol
cc      write(*,'(a,2(5X,i3))') '      Coarse and fine levels : ',lc,lf
cc      write(*,*) '-----------------------------',
cc     >           '--------------------------------------------------'
cc      write(*,*) 'Iteration  ||r||/||r_0||    ||r||  ',
cc     >                  '   ||Br||/||u||      ||Br||   ||Br||/||r_0||'
cc      write(*,*) '-----------------------------',
cc     >           '--------------------------------------------------'
C
      kbf   = kb(lf)
      kuf   = ku(lf)
C
cc      write(*,*)
cc      write(*,*) (r(krhs+ii-1), ii = 1,n)
C
      call init_g0(m(iao),m(jao),r(kra),r(krhs),r(ksol),n)
C
cc      write(*,*) (r(ksol+ii-1), ii = 1,n)
cc      write(*,*)
cc      write(*,*) (r(krhs+ii-1), ii = 1,n)
C
C...  To compute the initial residual: b = rhs - A*sol.
      call abyvam(r(krhs),m(iao),m(jao),r(kra),r(ksol),n,r(kbf))
C
      call scpro(r(kbf),r(kbf),err_res0,n)
      err_res0 = dsqrt(err_res0)
C
      if (err_res0 .lt. zero)  err_res0 = 1.d0
C
      kk = 0
cc      write(*,'(2X,i5,2X,5(3X,e11.4))') 
cc     >     kk,1.d0,err_res0,1.d0,1.d0,1.d0
C
C...  Iterate until convergence:  sol = sol + B*b = sol + B(rhs - A*sol).
C
      do 200 kk = 1, maxit
C
C...     Iterator B using an MG-cycle, i.e., u = B*b = B(rhs - A*sol).
C
         call premg_bil(m(1),ia(1),ja(1),idir(1),
     >        ipp(1),jpp(1),nzp(1),nz(1),nd(1),
     >        r(1),ku(1),kb(1),ka(1),lmx,lf,lc,lasti,lastr)
C
C...     New solution:  sol = sol + u = sol + B(rhs - A*sol).
         call uupluv(r(ksol),r(kuf),n)
C
C...     To compute the relative error.
cc 	   call scpro(r(kuf),r(kuf),xdiffn,n)
cc         call scpro(r(ksol),r(ksol),xnewn,n)
cc         xdiffn  = dsqrt(xdiffn)
cc         xnewn   = dsqrt(xnewn)
C
cc         if (xnewn .gt. xdiffn) then
cc            err_rel = xdiffn/xnewn
cc         else
cc            err_rel = xdiffn  
cc         end if
cc         err_brr = xdiffn/err_res0
C
C...     To compute the residual: b = rhs - A*sol.
         call abyvam(r(krhs),m(iao),m(jao),r(kra),r(ksol),n,r(kbf))
C
C...     To compute the L_2 norm of the residual.
         call scpro(r(kbf),r(kbf),err_res,n)
         err_res = dsqrt(err_res) 
         err_relres = err_res/err_res0
C
         nmod = mod(kk,3)
C
         if (nmod .eq. 1) then
            write(*,'(2X,i5,2X,5(3X,e11.4))')
     >           kk,err_relres,err_res,err_rel,xdiffn,err_brr
         end if
         if (err_res .lt. tol) go to 333
 200  continue
C
cc      write(*,*), ' Iteration limit exceeded:', maxit
 333  continue
      niter = min(kk,maxit)
      write(*,'(2X,i5,2X,5(3X,e11.4))')
     >     niter,err_relres,err_res,err_rel,xdiffn,err_brr
C
cc      write(*,*) '=============================', 
cc     >           '=================================================='
C
 500  return
      end
C=====================================================================
      subroutine premg_bil(m,ia,ja,idir,ipp,jpp,nzp,nz,nd,
     >     r,ku,kb,ka,lmx,lf,lc,lasti,lastr)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1),ia(lmx),ja(lmx),idir(lmx),nd(lmx),nz(lmx)
      dimension r(1),ku(lmx),kb(lmx),ka(lmx)
      dimension ipp(lmx),jpp(lmx),nzp(lmx)
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  PREMG stands for an multigrid algorithm as a preconditioner, in 
C...  other words, it is an iterator B in an MG-cycle. For the solver
C...  of linear systems various multigrid method will be used: 
C...  V-cycle, \-cycle, W-cycle, etc.
C...
C...  Parameter:
C...    ICYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle,
C...              3 = Variable V-cycle, (4 = W-cycle, 5 = full MG)
C...    IPRSM   - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPSSM   - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRSM   - number of pre-smoothings: 1,2,3, etc
C...    NPSSM   - number of post-smoothings: 1,2,3, etc
C...    LF      - level of the finest mesh: 2,3,4, etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C...    ISOLV   - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C---------------------------------------------------------------------
      icycle = ich(31)
      iprsm = ich(32)
      ipssm = ich(33)
      nprsm = ich(34)
      npssm = ich(35)
      lsolv = ich(38)
C
      call nullv(r(ku(lf)),nd(lf))
      do 10 i = lc, lf - 1
         call nullv(r(ku(i)),nd(i))
         call nullv(r(kb(i)),nd(i))
 10   continue
C
C...  V(\)-cycle of the MG, i.e., the action of the iterator B.
C...  Going downward in the V-cycle.
C
      icount = -1
C
      do 100 k = lf, lc+1, -1
         iak = ia(k)
         jak = ja(k)
         idc = idir(k-1)
         kuk = ku(k)
         kbc = kb(k-1)
         kbk = kb(k)
         kak = ka(k)
         ndc = nd(k-1)
         ndk = nd(k)
         nzk = nz(k)
         ipk = ipp(k)
         jpk = jpp(k)
         nzpk  = nzp(k)
         iend  = lasti
         kwk1  = lastr
         kend  = kwk1 + ndk*7

cc         write(*,*)'iend, kend in down-cycle', iend, kend
C
C...     Presmoothing by Gauss-seidel smoother nprsm times: u=u+R(b-Au)
C
         if (icycle .eq. 3) then !For variable V-cyle
            icount = icount + 1
            noofsm = nprsm*2**icount
         else
            noofsm = nprsm      !For V or \-cyle
         end if
C
	 if (iprsm .eq. 1) then 
            call fwd_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
	 else if (iprsm .eq. 2) then 
            call bwd_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
	 else if (iprsm .eq. 3) then 
            call sym_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
         end if
C
C...     To compute the residual wk = b - Au.
C
         call abyvam(r(kbk),m(iak),m(jak),r(kak),r(kuk),ndk,r(kwk1))
C
C...     Restriction to the lower level: b = P^t*wk = wk^t*P.
C
         call rvbyp_bil(m(ipk),m(jpk),r(kwk1),ndk,r(kbc),ndc)
         call zero_dir(r(kbc),m(idc),ndc)
 100  continue
C
C...  Solving exactly at the coarsest level, i.e., u = A^(-1)*b.
C
      itmax = 1000
      kuc   = ku(lc)
      iac   = ia(lc)
      jac   = ja(lc)
      kac   = ka(lc)
C
      if (lsolv .eq. 1) then
         call fwd_gs_ns_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        ndc,itmax,0)
      else if (lsolv .eq. 2) then
         call bwd_gs_ns_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        ndc,itmax,0)
      else if (lsolv .eq. 3) then
         call sym_gs_ns_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        ndc,itmax,0)
cc      else if (lsolv .eq. 4) then
cc         maxa = iend
cc         kau  = kwk1
cc         call sgauss(m(iac),m(jac),r(kac),r(kbc),m(maxa),r(kau),ndc)
cc         call copyv(r(kbc),r(kuc),ndc)
      end if
C
C...  Going upward in the V-cycle.
C
      do 200 k = lc+1, lf
         iak = ia(k)
         jak = ja(k)
         idc = idir(k-1)
         kuc = ku(k-1)
         kuk = ku(k)
         kbk = kb(k)
         kak = ka(k)
         ndc = nd(k-1)
         ndk = nd(k)
         nzk = nz(k)
         ipk = ipp(k)
         jpk = jpp(k)
         nzpk  = nzp(k)
         iend  = lasti
         kend  = lastr
C
C...     Correction with prolongation from the lower level: 
C...     i.e., u_k = u_k + P*u_{k-1}.
C
         call pbyvcs_bil(m(ipk),m(jpk),r(kuc),ndc,r(kuk),ndk)
C
C...     Postsmoothing by Gauss-Seidel smoother npssm times: u=u+R(b-Au)
C
         if (icycle .eq. 2)  go to 200        !No post smoothing in \-cycle
C
         if (icycle .eq. 3) then !For variable V-cyle
            noofsm = npssm*2**icount
            icount = icount - 1
         else
            noofsm = npssm      !For V or \-cyle
         end if
C
	 if (ipssm .eq. 1) then 
            call fwd_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
	 else if (ipssm .eq. 2) then 
            call bwd_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
	 else if (ipssm .eq. 3) then 
            call sym_gs_ns_wr(r(kuk),m(iak),m(jak),r(kak),r(kbk),
     >           ndk,noofsm,0) 
         end if
 200  continue 
C
      return
      end
C=====================================================================
C=====================================================================
      subroutine fwd_gs_ns_wr(x,ia,ja,an,b,n,max_sweeps,iw) 
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
      if (iw .eq. 1) then
         nmodd = n/7
         if(nmodd .eq. 0) nmodd = 2
      end if
C
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
cc         umax = dmax1(umax,dabs(ud))
         x(i) = u
 40   continue
C
cc      nmodd = mod(iteration,nmodd)
      nmodd = 1
      if (nmodd .eq. 1 .and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' FGS : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST : ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine bwd_gs_ns(x,ia,ja,an,b,n,max_sweeps,iw) 
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
      if (iw .eq. 1) then
         nmodd = n/7
         if(nmodd .eq. 0) nmodd = 2
      end if
C
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
cc      nmodd = mod(iteration,nmodd)
      nmodd = 1
      if (nmodd .eq. 1 .and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' BGS : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-13) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST : ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
C
      return
      end
C=====================================================================
      subroutine sym_gs_ns(x,ia,ja,an,b,n,max_sweeps,iw) 
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
      if (iw .eq. 1) then
         nmodd = n/7
         if(nmodd .eq. 0) nmodd = 2
      end if
C
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
      nmodd = mod(iteration,100)
cc      nmodd = 1
      if (nmodd .eq. 1.and. iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' SGS : ',iteration, ' |x_n-x_{n+1}| = ',umax  
      end if
C
      if (iteration .lt. max_sweeps .and. umax .gt. 0.5d-10) go to 20
C
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST : ',iteration, ' |x_n-x_{n+1}| = ',umax 
      end if
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
