C=======================================================================
      program Lap_mg
C=======================================================================
      implicit real*8(a-h,o-z)
      parameter (nlast = 10 000 000, nrlast = 50 000 000)  !lvl = 11
      common /top/ m(nlast)
      common /top1/ r(nrlast)
cc      dimension m(nlast), r(nrlast)
      common /quad/ wg(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /quad_gs_ntrl/ z_ntrl(7),zw(5,7)
      common /quad_gs/ w(7),z(7)
      common /menu/ ich(200)
C---------------------------------------------------------------------
C...  To solve the Laplace equation:
C...
C...     - Laplace(U) = f(x,y)     on the unit square,
C...               U  = g(x,y)     on the boundary of the unit square.
C...
C...  Various MG cycles are used to solve the linear system:
C...       V-cycle, \-cycle, variable V-cycle.
C...  Since this program is working on the structured mesh grids, no 
C...  matrix and mesh information are stored. Instead, only actions
C...  are performed.
C---------------------------------------------------------------------
      call input_data
      call quad_elt
      call quad_data
      call quad_data1
C
      icycle = ich(31)
      iprsm = ich(32)
      ipssm = ich(33)
      nprsm = ich(34)
      npssm = ich(35)
      lsolv = ich(38)
      lf    = ich(36)
      lc    = ich(37)
      maxit = ich(39)
      tol   = 1.d-6 
      n     = (2**lf + 1)**2
C
      ksol  = 1
      krhs  = ksol + n
      kfree = krhs + n
C
C...  Get the RHS vector.
C
      call loadb_tri(lf,n,r(krhs))
C
C...  To take care of the Dirichlet boundary condition.
      call ess_bdry_actn(lf,n,r(krhs),1)
C
C...  Solve by a MG cycle.
C
      call mg_actn(r(1),r(ksol),r(krhs),kfree,n,lf,lc,
     >     icycle,iprsm,ipssm,nprsm,npssm,lsolv,maxit,tol)
C
C********************************************************************
C...  To print out the computed solution SOL.
C********************************************************************
      iexact = 0
      if (iexact .eq. 1) then
         lvl = 0
         call write_soln(lvl,r(ksol),n,iformatted)
      end if
C
      if (n .le. 65*65) then
         call wrsol01(r(ksol),n,100)
         endfile(100)
         close(100)
      end if
C
      smin = r(ksol)
      smax = smin
C
      do i = 2 , n
cc         if( i.le. 65*65) write(100,'(e15.7$)') (r(ksol+i-1))
         smin = dmin1(smin,r(ksol+i-1))
         smax = dmax1(smax,r(ksol+i-1))
cc         if(dabs(r(ksol+i-1)) .gt. 1d-05) then
cc             write(*,*) i,r(ksol+i-1)
cc             read(*,*)
cc         end if
      end do
C
      write(*,*) 
     >'=============================================================',
     >'================='
      write(*,*)
      write(*,'(a,f20.12)') ' MINIMUM of the solution: ', smin
      write(*,'(a,f20.12)') ' MAXIMUM of the solution: ', smax
      write(*,*)
      write(*,*) 
     >'=============================================================',
     >'================='

      write(*,'(10x,a)') ' CPU time:  '
      write(*,'(10x,a)') ' ======================== '
      write(*,'(10x,a,f10.3,a)') '   SETUP: ',t_setup, ' s '
      write(*,'(10x,a)') ' ------------------------'
      write(*,'(10x,a,f10.3,a)') '   SOLVE: ',t_iter, ' s '
      write(*,'(10x,a)') ' ------------------------'
      write(*,'(10x,a,f10.3,a)') '   TOTAL: ',t_setup+t_iter, ' s '
      write(*,'(10x,a)') ' ======================== '
C
 999  continue
C
      stop
      end
C=======================================================================
      subroutine loadb_tri(lvl,n,b)
C=======================================================================
      implicit real*8(a-h,o-z)
      dimension b(1)
cc      common /menu/ ich(200)
C-----------------------------------------------------------------------
C...  This subroutine computes the vector b = (b(j)) in the structured
C...  triangular mesh, where b(j) = integral_omega(f*phi_j).
C...
C...  Note that no Dirichlet boundary condition is considered at this
C...  subroutine. But, it will be considered by calling the subroutine
C...  ESS_BDRY_ACTN.
C-----------------------------------------------------------------------
      ldm = 2**lvl
      ld  = ldm + 1
      h   = 1.d00/ldm
      isource = 2
C
      call nullv(b,n)
C
C...  Integration in Omega
C
      do 20 i = 1, ldm
         kt = (i-1)*ld
         do 10 j = 1, ldm
            k = kt + j
C
C...        The first triangle having K as the smallest node numbering.
C
            k1 = k
            k2 = k + ld
            k3 = k2 + 1
            ai11 = (i-1)*h
            ai21 = (j-1)*h
            ai12 = i*h
            ai22 = (j-1)*h
            ai13 = i*h
            ai23 = j*h
            x2m1 = ai12 - ai11
            x3m1 = ai13 - ai11
            y2m1 = ai22 - ai21
            y3m1 = ai23 - ai21
            det  = dabs(x2m1*y3m1 - x3m1*y2m1)
            element_area = det*0.5
C
            if (isource .eq. 1) then
               b123 = element_area*
     >              rhs(ai11+(x2m1+x3m1)/3,ai21+(y2m1+y3m1)/3)/3
               b(k1) = b(k1) + b123
               b(k2) = b(k2) + b123
               b(k3) = b(k3) + b123
            else if(isource .eq. 2) then
               b1 = element_area*rhs(ai11,ai21)/3
               b2 = element_area*rhs(ai12,ai22)/3
               b3 = element_area*rhs(ai13,ai23)/3
               b(k1) = b(k1) + b1
               b(k2) = b(k2) + b2
               b(k3) = b(k3) + b3
            end if
C
C...        The second triangle having K as the smallest node numbering.
C
            k1 = k
            k2 = k + ld + 1
            k3 = k + 1
            ai11 = (i-1)*h
            ai21 = (j-1)*h
            ai12 = i*h
            ai22 = j*h
            ai13 = (i-1)*h
            ai23 = j*h
            x2m1 = ai12 - ai11
            x3m1 = ai13 - ai11
            y2m1 = ai22 - ai21
            y3m1 = ai23 - ai21
            det  = dabs(x2m1*y3m1 - x3m1*y2m1)
            element_area = det*0.5
C
            if (isource .eq. 1) then
               b123 = element_area*
     >              rhs(ai11+(x2m1+x3m1)/3,ai21+(y2m1+y3m1)/3)/3
               b(k1) = b(k1) + b123
               b(k2) = b(k2) + b123
               b(k3) = b(k3) + b123
            else if(isource .eq. 2) then
               b1 = element_area*rhs(ai11,ai21)/3
               b2 = element_area*rhs(ai12,ai22)/3
               b3 = element_area*rhs(ai13,ai23)/3
               b(k1) = b(k1) + b1
               b(k2) = b(k2) + b2
               b(k3) = b(k3) + b3
            end if
 10      continue
 20   continue
C     
      return
      end
C=======================================================================
      subroutine form_dir(lvl,n,idir,idir_opt)
C=======================================================================
      implicit real*8(a-h,o-z)
      dimension idir(1)
C-----------------------------------------------------------------------
C... Dirichlet nodes are formed on the boundary of the unit square domain.
C... If i is a Dirichlet boundary node, then set IDIR(i) = n+1, otherwise 
C... set 0.
C...
C... Parameter:
C...   IDIR_OPT = 0 : No node on the boudary is Dirichlet node.
C...   IDIR_OPT = 1 : All nodes on the boudary are Dirichlet nodes.
C...   IDIR_OPT = 2 : The bottom portion is Dirichlet boundary.
C...   IDIR_OPT = 3 : The top portion is Dirichlet boundary.
C...   IDIR_OPT = 4 : The left side is Dirichlet boundary.
C...   IDIR_OPT = 5 : The right side is Dirichlet boundary.
C...   IDIR_OPT = 6 : The bottom and the left side are Dirichlet boundary.
C...   IDIR_OPT = 7 : The bottom and the right side are Dirichlet boundary.
C...   IDIR_OPT = 8 : The bottom and the top side are Dirichlet boundary.
C...   IDIR_OPT = 9 : The top and the left side are Dirichlet boundary.
C...   IDIR_OPT = 10 : The top and the right side are Dirichlet boundary.
C...   IDIR_OPT = 11 : The left and the right side are Dirichlet boundary.
C...   IDIR_OPT = 12 : Only the bottom side is Neumann boundary.
C...   IDIR_OPT = 13 : Only the top side is Neumann boundary.
C...   IDIR_OPT = 14 : Only the left side is Neumann boundary.
C...   IDIR_OPT = 15 : Only the right side is Neumann boundary.
C-----------------------------------------------------------------------
      ld = 2**lvl + 1
      call inullv(idir,n)
C
      if (idir_opt.eq.1 .or. idir_opt.eq.2 .or. idir_opt.eq.6 .or.
     >     idir_opt.eq.7 .or. idir_opt.eq.8 .or. idir_opt.eq.13 .or. 
     >     idir_opt.eq.14 .or. idir_opt.eq.15) then
C...     The bottom portion is Dirichlet boundary.
         do 10 i = 1, ld
            k = (i-1)*ld + 1
            idir(k) = n + 1
 10      continue
      else if (idir_opt.eq.1 .or. idir_opt.eq.3 .or. idir_opt.eq.8 .or.
     >     idir_opt.eq.9 .or. idir_opt.eq.10 .or. idir_opt.eq.12 .or. 
     >     idir_opt.eq.14 .or. idir_opt.eq.15) then
C...     The top portion is Dirichlet boundary.
         do 20 i = 1, ld
            k = i*ld 
            idir(k) = n + 1
 20      continue
      else if (idir_opt.eq.1 .or. idir_opt.eq.4 .or. idir_opt.eq.6 .or.
     >     idir_opt.eq.9 .or. idir_opt.eq.11 .or. idir_opt.eq.12 .or. 
     >     idir_opt.eq.13 .or. idir_opt.eq.15) then
C...     The left side is Dirichlet boundary.
         do 30 i = 1, ld
            k = i
            idir(k) = n + 1
 30      continue
      else if (idir_opt.eq.1 .or. idir_opt.eq.5 .or. idir_opt.eq.7 .or.
     >     idir_opt.eq.10 .or. idir_opt.eq.11 .or. idir_opt.eq.12 .or. 
     >     idir_opt.eq.13 .or. idir_opt.eq.14) then
C...     The right side is Dirichlet boundary.
         do 40 i = 1, ld
            k = (ld-1)*ld + i
            idir(k) = n + 1
 40      continue
      end if
C
      return
      end
C=======================================================================
C=====================================================================
      subroutine mg_actn(r,sol,rhs,kfree,n,lf,lc,
     >     icycle,iprsm,ipssm,nprsm,npssm,lsolv,maxit,tol)
C=====================================================================
      parameter (lmx = 15)
      implicit real*8(a-h,o-z)
      dimension nd(lmx),r(1),ku(lmx),kb(lmx),sol(1),rhs(1)
C---------------------------------------------------------------------
C...  Standard multigrid cycles for the SPD problems on the structured
C...  grid of the unit square.
C...
C...  Parameter:
C...    N       - dimension of the mesh, N = (2**LF + 1)**2
C...    ICYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle,
C...              3 = Variable V-cycle
C...    IPRSM   - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPSSM   - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRSM   - number of pre-smoothings: 1,2,3, etc
C...    NPSSM   - number of post-smoothings: 1,2,3, etc
C...    LF      - level of the finest mesh: 2,3,4, etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C...    MAXIT   - maximum number of MG iterations
C...    TOL     - stopping criterion for iteration procedure.
C---------------------------------------------------------------------
cc      maxit  = 100
cc      tol    = 0.5d-6
      zero   = 1.d-13
C
      lastr = kfree
C
      do 50 i = lf, lc, -1
         nd(i) = (2**i + 1)**2
         ku(i) = lastr
         kb(i) = ku(i) + nd(i)
         lastr = kb(i) + nd(i)
 50   continue
C
      kwk = lastr
C
C...  ------------------- MG starts here. ---------------------------
C
      write(*,*)
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
      call init_guess_actn(lf,n,rhs,sol)
C
C...  To compute the initial residual: b = rhs - A*sol:
C...  Compute b = A*sol.
      call abyvg_actn(lf,n,sol,r(kwk),r(kbf))
C
C...  Compute b = rhs - b =  rhs - A*sol.
      call vuminv(rhs,r(kbf),n)
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
         call premg_actn(nd(1),r(1),ku(1),kb(1),lastr,lmx,lf,lc,
     >        icycle,iprsm,ipssm,nprsm,npssm,lsolv)
C
C...     New solution:  sol = sol + u = sol + B(rhs - A*sol).
         call uupluv(sol,r(kuf),n)
C
C...     To compute the relative error.
	 call scpro(r(kuf),r(kuf),xdiffn,n)
         call scpro(sol,sol,xnewn,n)
C
         xdiffn  = dsqrt(xdiffn)
         xnewn   = dsqrt(xnewn)
C
         if (xnewn .gt. xdiffn) then
            err_rel = xdiffn/xnewn
         else
            err_rel = xdiffn  
         end if
         err_brr = xdiffn/err_res0
C
C...     To compute the residual: b = rhs - A*sol.
C...     Compute b = A*sol.
         call abyvg_actn(lf,n,sol,r(kwk),r(kbf))

C...     Compute b = rhs - b =  rhs - A*sol.
         call vuminv(rhs,r(kbf),n)
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
cc         if (err_relres .lt. tol .and. err_rel .lt. tol) go to 333
         if (err_rel .lt. tol) go to 333
 200  continue
C
 333  continue
      niter = min(kk,maxit)
      write(*,'(2X,i5,2X,5(3X,e11.4))')
     >     niter,err_relres,err_res,err_rel,xdiffn,err_brr
C
cc      write(*,*) '=============================', 
cc     >           '=================================================='
C
cc      arfac  = (err_relres)**(1./dble(kk))
cc      arfac1 = (xdiffn/err_res0)**(1./dble(kk))
cc      write(*,*) 
cc      write(*,'(a,2f12.5)'),' Average Reduction Factor:', arfac,arfac1
cc      write(*,*)
C
 500  return
      end
C=====================================================================
      subroutine premg_actn(nd,r,ku,kb,lastr,lmx,lf,lc,
     >     icycle,iprsm,ipssm,nprsm,npssm,lsolv)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension nd(lmx),r(1),ku(lmx),kb(lmx)
C---------------------------------------------------------------------
C...  PREMG stands for an multigrid algorithm as a preconditioner, in 
C...  other words, it is an iterator B in an MG-cycle. Various multigrid 
C...  cycles will be used: V-cycle, \-cycle, variable V-cycle
C...
C...  Parameter:
C...    ICYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle,
C...              3 = Variable V-cycle
C...    IPRSM   - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPSSM   - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRSM   - number of pre-smoothings: 1,2,3, etc
C...    NPSSM   - number of post-smoothings: 1,2,3, etc
C...    LF      - level of the finest mesh: 2,3,4, etc
C...    LC      - level of the coarsest mesh: 1,2,3, etc
C...    ISOLV   - choice of the coarsest grid solver:
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C---------------------------------------------------------------------
      kwk1  = lastr
      kwk   = kwk1 + nd(lf)
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
         kuk = ku(k)
         kbc = kb(k-1)
         kbk = kb(k)
         ndc = nd(k-1)
         ndk = nd(k)
C
C...     Presmoothing by Gauss-seidel smoother noofsm times: u=u+R(b-Au)
C
         if (icycle .eq. 3) then !For variable V-cyle
            icount = icount + 1
            noofsm = nprsm*2**icount
         else
            noofsm = nprsm      !For V or \-cyle
         end if
C
         call gs_actn(k,ndk,r(kuk),r(kbk),r(kwk),iprsm,noofsm,0) 
C
C...     To compute the residual wk1 = b - Au:
C...     Compute wk1 = A*u.
         call abyvg_actn(k,ndk,r(kuk),r(kwk),r(kwk1))

C...     Compute wk1 = b - wk1 =  b - A*u.
         call vuminv(r(kbk),r(kwk1),ndk)
C
C...     Restriction to the lower level: bc = P^t*wk1 = wk1^t*P.
C
         call restrict_actn(k,ndk,r(kwk1),r(kwk),r(kbc))
 100  continue
C
C...  Solving exactly at the coarsest level by GS iteration, 
C...  i.e., u = A^(-1)*b.
C
      itmax = 1000
      kuc   = ku(lc)
C
      call gs_actn(lc,ndc,r(kuc),r(kbc),r(kwk),lsolv,itmax,0) 
C
C...  Going upward in the V-cycle.
C
      do 200 k = lc+1, lf
         kuc = ku(k-1)
         kuk = ku(k)
         kbk = kb(k)
         ndc = nd(k-1)
         ndk = nd(k)
C
C...     Correction with prolongation from the lower level: 
C...     i.e., u_k = u_k + P*u_{k-1}.
C
         call prolong_actn(k,ndk,r(kuc),r(kwk),r(kwk1))
         call uupluv(r(kuk),r(kwk1),ndk)
C
C...     Postsmoothing by Gauss-Seidel smoother noofsm times: u=u+R(b-Au)
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
         call gs_actn(k,ndk,r(kuk),r(kbk),r(kwk),ipssm,noofsm,0)
 200  continue 
C
      return
      end
C=======================================================================
      subroutine init_guess_actn(lvl,n,b,x)
C=======================================================================
      implicit real*8(a-h,o-z)
      dimension  x(1),b(1)
C-----------------------------------------------------------------------
C...  Initial guess.
C-----------------------------------------------------------------------
      ld = 2**lvl + 1
C
      do 20 i = 1, ld
         kt = (i-1)*ld
         do 10 j = 1, ld
            k = kt + j
            if (j.eq.1 .or. j.eq.ld .or. i.eq.1 .or. i.eq.ld) then
               x(k) = b(k)
            else
               x(k) = 0.d0
            end if
 10      continue
 20   continue
C
      return
      end
C=======================================================================
      subroutine ess_bdry_actn(lvl,n,q,ibdry)
C=======================================================================
      implicit real*8(a-h,o-z)
      dimension  q(1)
cc      common /menu/ ich(200)
C-----------------------------------------------------------------------
C...  This subroutine adjusts the LVL-th level vector q
C...  so that it satisfies the Dirichlet boundary condition.
C...
C...  Parameters:
C...    IBDRY = 0 : Zero boundary condition
C...    IBDRY = 1 : Nonzero boundary condition with the function g(x,y)
C-----------------------------------------------------------------------
      ld = 2**lvl + 1
      h  = 1.d0/(ld-1)
C
      if (ibdry .eq. 0) then
         do 20 i = 1, ld
            kt = (i-1)*ld
            do 10 j = 1, ld
C...           On the boundary, q(k) = 0.
	       if (j.eq.1 .or. j.eq.ld .or. i.eq.1 .or. i.eq.ld) then
                  k = kt + j
	          q(k) = 0.d0
	       end if
 10         continue
 20      continue
      else if (ibdry .eq. 1) then
         do 40 i = 1, ld
            kt = (i-1)*ld
            do 30 j = 1, ld
C...           On the boundary, q(k) = g(x_k,y_k).
	       if (j.eq.1 .or. j.eq.ld .or. i.eq.1 .or. i.eq.ld) then
                  k = kt + j
                  x = (i-1)*h
                  y = (j-1)*h
	          q(k) = g(x,y)
	       end if
 30         continue
 40      continue         
      end if
C
      return
      end
C=======================================================================
      subroutine restrict_actn(lvl,n,p,pm,q)
C=======================================================================
      implicit real*8(a-h,o-z)
      dimension p(1),q(1),pm(2**lvl+1, 2**lvl+1)
C-----------------------------------------------------------------------
C...                                  l   t               l   t
C...  This subroutine computes q = ( I   )  * p, where ( I   )  is the
C...                                  l-1                 l-1
C...  l-th level restriction operator.
C-----------------------------------------------------------------------
      ld  = 2**lvl+1
      ldc = 2**(lvl-1)+1
C
      call nullv(q,ldc*ldc)
C
C...  Matrix representation(pm) of the vector p.
C
      do 20 i = 1,ld
         kt = (i-1)*ld
	 do 10 j = 1,ld
	    k = kt + j
            pm(i,j) = p(k)
 10      continue
 20   continue
C
      do 50 ic = 1, ldc
         ib = 2*ic - 1
	 kt = (ic-1)*ldc
	 do 40 jc = 1, ldc
            jb = 2*jc - 1
	    k  = kt + jc
            if (ic.eq.1 .or. ic.eq.ldc .or. jc.eq.1 .or. jc.eq.ldc) then
               q(k) = pm(ib,jb)      !For the Dirichlet boundary nodes.
            else
               q(k) = pm(ib,jb) + (pm(ib,jb-1) + pm(ib,jb+1) 
     >              + pm(ib-1,jb) + pm(ib+1,jb) 
     >              + pm(ib-1,jb-1) + pm(ib+1,jb+1))/2 
            end if
 40      continue
 50   continue
C
      return
      end
C=======================================================================
      subroutine prolong_actn(lvl,n,p,pm,q)
C=======================================================================
      implicit real*8(a-h,o-z)
      dimension p(1),q(1),pm(2**(lvl-1)+1, 2**(lvl-1)+1)
C-----------------------------------------------------------------------
C...                                  l                   l
C...  This subroutine computes q = ( I   )  * p, where ( I   ) is the
C...                                  l-1                 l-1
C...  l-th level prolongation operator.
C-----------------------------------------------------------------------
      ldcm = 2**(lvl-1)
      ldc = ldcm + 1
      ld = 2**lvl + 1
C     
      call nullv(q,n)
C
C...  Matrix representation(pm) of the vector p.
C  
      do 20 i = 1,ldc
         kt = (i-1)*ldc
	 do 10 j = 1,ldc
	    k = kt + j
            pm(i,j) = p(k)
 10      continue
 20   continue
C
C...  Consider the odd columns first.  
C
      do 40 i = 1, ldc
CCC         jb = 2*i - 1, kt = (jb-1)*ld
         kt = (2*i-2)*ld
	 do 30 j = 1, ldcm
	    k      = kt + 2*j
	    q(k-1) = pm(i,j)
	    q(k)   = (pm(i,j) + pm(i,j+1))/2      
 30      continue
	 q(kt+ld)  = pm(i,ldc)
 40   continue
C
C...  Now consider the even columns.  
C
      do 60 i = 1, ldcm
CCC         jb = 2*i, kt = (jb-1)*ld
         kt = (2*i-1)*ld
	 do 50 j = 1, ldcm
	    k      = kt + 2*j
	    q(k-1) = (pm(i,j) + pm(i+1,j))/2
	    q(k)   = (pm(i,j) + pm(i+1,j+1))/2      
 50      continue
         q(kt+ld)  = (pm(i,ldc) + pm(i+1,ldc))/2
 60   continue
C
      return
      end
C=======================================================================
      subroutine gs_actn(lvl,n,x,b,xm,ifbs,max_sweeps,iw) 
C=======================================================================
      implicit real*8(a-h,o-z)
      dimension x(1),b(1),xm(2**lvl + 1, 2**lvl + 1)
C---------------------------------------------------------------------
C...  This subroutine solves the linear system A*x = b by Gauss-Seidel 
C...  method in the structured mesh.
C...
C...  Parameters:
C...    IFBS = 1, 2, 3 : Forward, Backward, Symmetric GS
C...    IW             : Whether write the convergence history or not,
C...                     IW = 0 : do not write it; IW = 1 : write it.
C---------------------------------------------------------------------
      tol = 0.5d-13         !tol - stopping criteria for iterations
      ld  = 2**lvl + 1
      ldm = ld - 1
C
C...  Matrix representation(xm) of the vector x.
C  
      do 20 i = 1, ld
	 kt = (i-1)*ld
	 do 10 j = 1, ld
	    k = kt + j
            if (j.eq.1 .or. j.eq.ld .or. i.eq.1 .or. i.eq.ld) then
               xm(i,j) = 0.d0
            else
               xm(i,j) = x(k)
            end if
 10      continue
 20   continue
C
      do 200 kk = 1, max_sweeps
ccc      xdiffn = 0.d0
         xmax = -1.d20
	 if (ifbs .eq. 2)  go to 60
C
         do 50 i = 2, ldm          !Forward GS
            kt = (i-1)*ld
            do 40 j = 2, ldm
               k = kt + j
               xt = (b(k) + xm(i-1,j) + xm(i,j+1) 
     .              + xm(i+1,j) + xm(i,j-1))/4.d0
               if (ifbs .eq. 3)  go to 30
               xdiff = xt - x(k)
               xmax = dmax1(xmax,dabs(xdiff))
ccc            xdiffn = xdiffn + xdiff*xdiff
 30            xm(i,j) = xt
 40         continue
 50      continue
C
	 if (ifbs .eq. 1)  go to 90
 60      continue
C
         do 80 i = ldm, 2, -1           !Backward GS
            kt = (i-1)*ld
            do 70 j = ldm, 2, -1
               k = kt + j
               xt = (b(k) + xm(i-1,j) + xm(i,j+1) 
     >              + xm(i+1,j) + xm(i,j-1))/4.d0
               xdiff = xt - x(k)
               xmax = dmax1(xmax,dabs(xdiff))
ccc            xdiffn = xdiffn + xdiff*xdiff
               xm(i,j) = xt
 70         continue
 80      continue
C
 90      continue
C
         do 110 i = 2, ldm
            kt = (i-1)*ld
            do 100 j = 2, ldm
               k = kt + j
               x(k) = xm(i,j)
 100        continue
 110     continue
C
ccc      call scpro(x,x,xnn,n)
ccc	 if (xnn .lt. 1.d-14) then
ccc         stopl = dsqrt(xdiffn)
ccc	 else
ccc	    stopl  = dsqrt(xdiffn/xnn)
ccc	 end if
C
cc         nmodd = mod(kk,100)
         nmodd = 1
         if (nmodd .eq. 1 .and. iw .eq. 1) then
            write(*,'(a,i7,a,e20.12)') 
     >           ' GS : ',kk, ' |x_n-x_{n+1}| = ',xmax  
         end if
         if (xmax .lt. tol)  go to 333
 200  continue
C
 333  if (kk .gt. max_sweeps) kk = max_sweeps
      if (iw .eq. 1) then
         write(*,'(a,i7,a,e20.12)') 
     >        ' LAST : ',kk, ' |x_n-x_{n+1}| = ',xmax 
      end if
C
      return
      end
C=======================================================================
      subroutine abyvg_actn(lvl,n,x,xm,w)
C=======================================================================
      implicit real*8(a-h,o-z)
      dimension x(1),w(1),xm(2**lvl + 1, 2**lvl + 1)
C-----------------------------------------------------------------------
C...  This subroutine computes matrix_vector multiplication(W=B*x) in the
C...  structured grid, where B is the stiffness matrix.
C...
C...  Parameters:     
C...    LVL, N   : N = (2**lvl + 1)^2
C...    XM       : matrix representation of X
C-----------------------------------------------------------------------
      ld  = 2**lvl + 1
      ldm = ld - 1
C
      call nullv(w,n)
C
C...  Matrix representation xm of the vector x.
C
      do 30 i = 1, ld
	 kt = (i-1)*ld
	 do 20 j = 1, ld
	    k = kt + j
            if (j.eq.1 .or. j.eq.ld .or. i.eq.1 .or. i.eq.ld) then
               xm(i,j) = 0.d0
               w(k)    = x(k)
            else
               xm(i,j) = x(k)
            end if
 20      continue
 30   continue
C      
      do 60 i = 2, ldm
	 kt = (i-1)*ld
	 do 50 j = 2, ldm
	    k = kt + j
            w(k) = 4*xm(i,j) - xm(i-1,j) - xm(i+1,j)
     >           - xm(i,j-1) - xm(i,j+1) 
 50      continue
 60   continue
C     
      return
      end
C=======================================================================
      subroutine abyvg_actn_more(lvl,n,x,xm,w,imat,c,d,f)
C=======================================================================
      implicit real*8(a-h,o-z)
      dimension x(1),w(1),xm(2**lvl + 1, 2**lvl + 1)
C-----------------------------------------------------------------------
C...  This subroutine computes matrix_vector multiplication(W=B*x) in the
C...  structured grid, where B is defined below.
C...
C...  Parameters:     
C...    imat = 1 : B = A;
C...    imat = 2 : B = A + cA_x + dA_y;
C...    imat = 3 : B = A + cA_x + dA_y + fM;
C...    imat = 4 : B = cA_x + dA_y;
C...    imat = 5 : B = cA_x + dA_y + fM;
C...    imat = 6 : B = fM;
C...    imat = 7 : B = A + fM,
C...    A        : the stiffness matrix 
C...    A_x, A_y : lower order matrices 
C...    M        : the mass matix
C...    C, D, F  : constants
C...    LVL, N   : N = (2**lvl + 1)^2
C...    XM       : matrix representation of X
C-----------------------------------------------------------------------
      ld = 2**lvl + 1
      ldm = ld - 1
      ldm2 = ldm*ldm
C
      call nullv(w,n)
C
C...  Matrix representation xm of the vector x.
C
      do 30 i = 1, ld
	 kt = (i-1)*ld
	 do 20 j = 1, ld
	    k = kt + j
            if (j.eq.1 .or. j.eq.ld .or. i.eq.1 .or. i.eq.ld) then
               xm(i,j) = 0.d0
               w(k)    = x(k)
            else
               xm(i,j) = x(k)
            end if
 20      continue
 30   continue
C      
      do 60 i = 2, ldm
	 kt = (i-1)*ld
	 do 50 j = 2, ldm
	    k = kt + j
            if (imat.ne.4 .and. imat.ne.5 .and. imat.ne.6) then
	       w(k) = 4*xm(i,j) - xm(i-1,j) - xm(i+1,j)
     >              - xm(i,j-1) - xm(i,j+1) 
            end if
	    if (imat.ne.1 .and. imat.ne.6 .and. imat.ne.7) then
	       w(k) = w(k) + c*((-2*xm(i-1,j) + 2*xm(i+1,j)-xm(i,j-1)
     >              + xm(i+1,j-1) - xm(i-1,j+1) + xm(i,j+1))/(6*ldm))
     >              + d*((-xm(i-1,j) + xm(i+1,j) - 2*xm(i,j-1) 
     >              - xm(i+1,j-1) + xm(i-1,j+1) + 2*xm(i,j+1))/(6*ldm))
	    end if
	    if (imat.ne.1 .and. imat.ne.2 .and. imat.ne.4) then
	       w(k) = w(k) + f*((6*xm(i,j) + xm(i-1,j) + xm(i+1,j)
     >              + xm(i,j-1) + xm(i+1,j-1)
     >              + xm(i-1,j+1) + xm(i,j+1))/(12*ldm2))
	    end if
 50      continue
 60   continue
C      
      return
      end
C=======================================================================
C     END of PROGRAM                                                   C
C=======================================================================
C=======================================================================
C=======================================================================
