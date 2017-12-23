C====================================================================
      subroutine match10(
     >     n,ia,ja,a,nnz,jp,jpt,iord,kmatch,
     >     nc,nnzc,iac,jac,kac,
     >     m,r,ifree,kfree)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1),ia(1),ja(1),r(1),a(1),iord(1)
      external heavedg,lighedg
C
      jp =  1
      jpt = jp + n + 1
      imask = jpt + 2*n
C
cc      call wrxy(x,y,n,100)
C
      maxadj = 100
C
cc      call wradj_r(ia,ja,a,n,200)
C
      call inullv(m(imask),n)
      call inullv(m(jp),n)
C
      kfree0 = 1
      nc = 0
C
cc      call lprr(ia,ja,a,n)
C
      iheav = 1
      if(iheav .ne. 0) then
         call match00(heavedg,
     >        ia,ja,a,m(jp),m(jpt),n,nc,lenpp,
     >        r(kfree0),maxadj,m(imask),
     >        iord,kmatch,kiso,kdir)
      else
         call match00(lighedg,
     >        ia,ja,a,m(jp),m(jpt),n,nc,lenpp,
     >        r(kfree0),maxadj,m(imask),
     >        iord,kmatch,kiso,kdir)
      end if
C     
      kac = 1
      iac = jpt + 2*nc
      jac = iac + nc + 1
      imask = jac + nnz
C
      call ptap(n,nc,ia,ja,a,m(jp),m(jpt),
     >     m(iac),m(jac),r(kac),m(imask))
C
cc      call lprr(m(iac),m(jac),r(kac),ndc)
cc      read(*,*)
C
      nnzc = m(iac+nc)-1
      jp   = jp  + ifree-1
      jpt  = jpt + ifree-1
      iac  = iac + ifree-1
      jac  = jac + ifree-1
      ifree = jac + nnzc
C
      kac = kfree 
      kfree = kac + nnzc
C
      return 
      end
C====================================================================
      subroutine match00(getedg,ia,ja,a,jp,jpt,
     >     n,nc,lpro,work,maxadj,mask,
     >     iord,kmatch,kiso,kdir)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),a(1),jp(1),jpt(1),mask(1),iord(1)
      dimension iwork(100),work(maxadj)
      external getedg
C
cc      call lpri(ia,ja,n)
      if(nc .ne. 0) stop 77
C
      kdir = 0
      do k = 1,n
         if(ia(k+1)-ia(k) .le. 1) then
            mask(k) = 2*n+2
cc            iord(n-kdir) = k
            kdir = kdir + 1
         end if
      end do
C
      kmatch = 0
      kiso = 0
      kbeg = n-kdir
cc      call lpri(ia,ja,n)
      do k = n,1,-1
         if(mask(k) .lt. n) then
            iz = 0
            do jk = ia(k),ia(k+1)-1
               j = ja(jk)
               ajk = a(jk)
               if(j .eq. k) then
                  ak = ajk
                  dk = 1.d00/ak
                  go to 10
               end if
            end do
            write(*,*) ' no diagonal' , k
            stop 55
 10         do jk = ia(k),ia(k+1)-1
               j = ja(jk)
               ajk = a(jk)
               if(mask(j) .eq. 0) then
                  iz = iz + 1
                  if(j .ne. k) then
                     do i123 = ia(j),ia(j+1)-1
                        if(ja(i123) .eq. j) then
                           aj = a(i123)
                           dj = 1.d00/aj
                           go to 20
                        end if
                     end do
                     write(*,*) ' no symmetric pattern ',j,k
                     stop 66
 20                  work(iz) = ajk*ajk*dj*dk
                     iwork(iz) = j
                  else
                     iwork(iz) = 0
                     work(iz) = -1d15
                  end if
               end if
            end do
            nc = nc + 1
            if(iz .eq. 1) then
C...           Isolated points.
               jp(k) = nc
               jpt(2*nc-1) = k
               jpt(2*nc) = 0
               mask(k) = n + mask(k)
cc               iord(kbeg-kiso) = k
               kiso = kiso + 1
            else
C...           Matched edges.
               call getedg(work,iwork,ipick,iz)
               jp(k) = nc
               jp(ipick) = nc
               jpt(2*nc-1) = k
               jpt(2*nc) = ipick
               mask(k) = n + mask(k)
               mask(ipick) = n+mask(ipick)
               kmatch = kmatch + 1
cc               iord(kmatch) = k
               kmatch = kmatch + 1
cc               iord(kmatch) = ipick
            end if
         end if
      end do
C     
      lpro = 2*nc
C
cc      write(*,*) " edges in the matching ", nc
cc      write(*,*) " prolongation length ", lpro
C
      return
      end
C====================================================================
      subroutine heavedg(wei,numb,iheav,n)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension wei(1),numb(1)
C--------------------------------------------------------------------
C...  Pick the heaviest WEI.
C--------------------------------------------------------------------
cc      write(*,*) (numb(k),k=1,n)
cc      read(*,*)
      j = 1
      iheav = numb(j)
      temp = wei(j)
      do while (numb(j) .eq. 0)
         j = j + 1
         iheav = numb(j)
         temp = wei(j)
      end do
      do k = j+1 , n
         if(wei(k) .gt. temp .and. numb(k) .ne. 0) then
            temp = wei(k)
            iheav = numb(k)
         end if
      end do
      return
      end
C====================================================================
      subroutine lighedg(wei,numb,iligh,n)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension wei(1),numb(1)
C--------------------------------------------------------------------
C...  Pick the lightest WEI.
C--------------------------------------------------------------------
c      write(*,*) (numb(k),k=1,n)
c      read(*,*)
      j = 1
      iligh = numb(j)
      temp = wei(j)
      do while (numb(j) .eq. 0)
         j = j + 1
         iligh = numb(j)
         temp = wei(j)
      end do
      do k = j+1 , n
         if(wei(k) .lt. temp .and. numb(k) .ne. 0) then
            temp = wei(k)
            iligh = numb(k)
         end if
      end do
      return
      end
C====================================================================
      subroutine ptap(n,nc,ia,ja,a,jp,jpt,iac,jac,ac,pntr)
C====================================================================
      real*8 a(1),ac(1),hold
      integer*4 ia(1),ja(1),jp(1),jpt(1),pntr(1),iac(1),jac(1)
C--------------------------------------------------------------------
C...  Special C = (P^t)*A*P for the matching
C--------------------------------------------------------------------
      do i = 1 , nc
         pntr(i) = 0
         iac(i) = 0
      end do
      iac(nc+1) = 0
C
      iacp = 1
      do i = 1 , nc
         iac(i) = iacp
         do ik = 1 , 2
            i1 = jpt(2*i+ik-2) 
cc            write(*,*)i,i1, " i,i1"
            if(i1 .eq. 0) go to 10
            do ji = ia(i1),ia(i1+1)-1
               j = ja(ji)
               k = jp(j)
               if(k .ne. 0) then
                  if(pntr(k) .ne. i) then
                     jac(iacp) = k
cc                     write(*,*)j,k, "j,k"
                     iacp = iacp + 1
                     pntr(k) = i
                  end if
               end if
            end do
 10         continue
         end do
cc         read(*,*)
      end do
      iac(nc+1) = iacp
C
cc      call lpri(iac,jac,nc)
C...  End of symbolic part...
C
      do i = 1 , nc
         ibegin = iac(i)
         iend = iac(i+1)-1
         i1 = jpt(2*i-1)
         i2 = jpt(2*i)
         do ji = ibegin,iend
            j = jac(ji)
            j1 = jpt(2*j-1)
            j2 = jpt(2*j)
C
            hold = 0d0
            ik = ia(i1)
            ic = 0
            do while (ik .le. ia(i1+1)-1 .and. ic .lt. 2)
               if(ja(ik) .eq. j1 .or. ja(ik) .eq. j2) then
                  hold = hold + a(ik)
                  ic = ic + 1
               end if
               ik = ik + 1
            end do
C
            if(i2 .eq. 0) go to 20
            ik = ia(i2)
            ic = 0
            do while (ik .le. ia(i2+1)-1 .and. ic .lt. 2)
               if(ja(ik) .eq. j1 .or. ja(ik) .eq. j2) then
                  hold = hold + a(ik)
                  ic = ic + 1
               end if
               ik = ik + 1
            end do
 20         continue
            ac(ji) = hold
         end do
      end do
      return
      end
C====================================================================
      subroutine prolo(iord,jp,jpt,n,nc,kmatch,kiso,kdir)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension iord(1),jp(1),jpt(1)
C
      do k = 1,n
         jp(k) = 0
         jpt(k) = 0
      end do
      ncc = 0
      do k = 1 , kmatch , 2
         ncc = ncc + 1
         jpt(2*ncc-1) = iord(k)
         jpt(2*ncc) = iord(k+1)
         jp(iord(k)) = ncc
         jp(iord(k+1)) = ncc
      end do
      do k = kmatch + 1 , kmatch + kiso
         ncc = ncc + 1
         jpt(2*ncc-1) = iord(k)
         jpt(2*ncc) = 0
         jp(iord(k)) = ncc
      end do
      write(*,*) kmatch,kiso,ncc,nc
      read(*,*) iii
      write(*,*)  " jp ", (jp(kk),kk=1,n)
      read(*,*) iii
      write(*,*)  " jpt ", (jpt(kk),kk=1,ncc)
      return
      end
C====================================================================
      subroutine prolong(jp,x,y,nc,n)
C====================================================================
      implicit  real*8 (a-h,o-z)
      dimension  jp(1),x(1),y(1)
C--------------------------------------------------------------------
C...  Prolongation y = P*x.
C...
C...  Parameters:
C...    n  - no. of rows of P
C...    nc - no. of columns of P
C--------------------------------------------------------------------
      do i=1,n
         if (jp(i) .le. 0) then
            y(i)=0d0
         else
            y(i)=x(jp(i))
         end if
      end do
C     
      return 
      end 
C====================================================================
      subroutine prolong_plus(jp,x,y,nc,n)
C====================================================================
      implicit  real*8 (a-h,o-z)
      dimension  jp(1),x(1),y(1)
C--------------------------------------------------------------------
C...  Adding prolongation y = y + P*x.
C...
C...  Parameters:
C...    n  - no. of rows of P
C...    nc - no. of columns of P
C--------------------------------------------------------------------
      do i=1,n
         if (jp(i) .gt. 0) then
            y(i)=y(i) + x(jp(i))
         end if
      end do
C     
      return 
      end 
C====================================================================
      subroutine restr(jp,x,y,nc,n)
C====================================================================
      implicit  real*8 (a-h,o-z)
      dimension  jp(1),x(1),y(1)
C--------------------------------------------------------------------
C...  Restriction x = P^t*y.
C...
C...  Parameters:
C...    n  - no. of rows of P
C...    nc - no. of columns of P
C--------------------------------------------------------------------
      call nullv(x,nc)
C
      do i=1,n
         if(jp(i) .gt. 0) then
            x(jp(i))=x(jp(i))+y(i)
         end if
      end do
      return 
      end
C=====================================================================
      subroutine mg_match(sol,rhs,n,max_iter,m,r,
     >     tol,ex_sol,iexact)
C=====================================================================
      implicit none
C
      include 'paramsmg.h'
C
      integer maximum_levels
      integer nd,iord, 
     > ia, ja, ka, 
     > ipp,jpp, kpp,
     > ku , kb, 
     > ifree , kfree,lf,lc
      parameter (maximum_levels = 100)
C
      common /point/ 
     > nd(maximum_levels), iord(maximum_levels), 
     > ia(maximum_levels), ja(maximum_levels), ka(maximum_levels), 
     > ipp(maximum_levels),jpp(maximum_levels), kpp(maximum_levels),
     > ku(maximum_levels) , kb(maximum_levels), 
     > ifree , kfree ,lf,lc

      double precision r(1),sol(1),rhs(1),ex_sol(1),zero,one,nr0
      double precision xdiffn,xnewn,err_rel,err_res,err_relres,tol,h
      double precision ex_sol_norm
      integer m(1)
      integer*4 kmatch(100),nmod,maxa,kau,noofsm,kk,k2
      integer i,k,iak,jak,kuk,kbc,kbk,kak,ndc,ndk,jpk,
     >     iordk,kwk1,n,max_iter,itmax,iexact,
     >     iac,jac,kuc,kac,jiter,jptk,nnzk,nnzc,kfrold,ifrold
      integer iaw,jaw,isub,iwork1,iwork2,iwork3,nblk,nh
      common /onezero/ zero,one
C---------------------------------------------------------------------
C...  MG_MATCH stands for an multigrid algorithm (see also MG_MATCH_MODIFY
C...  which is identical to MG_MATCH), and it performs
C...  u^{k+1} = u^k + B(f - A*u^k), where B an iterator B in an MG-cycle. 
C...  Various multigrid method will be used: V-cycle, \-cycle, variable
C...  V-cycle, etc.
C...
C...  Parameter:
C...    ICYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle, 
C...              3 = variable V-cycle
C...    IPRSM   - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPSSM   - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRSM   - number of pre-smoothings: 1,2,3, etc
C...    NPSSM   - number of post-smoothings: 1,2,3, etc
C...    LF      - level of the finest mesh, which is backward in order: 1
C...    LC      - level of the coarsest mesh: 2,3, etc
C...    ISOLV_COARSE  - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C...              5 = CGM
C...    ICH_ORD - To determine whether Tarjan's ordering is considered 
C...              or not. If ordering is considered, then it is 1.
C...    IEXACT  - To determine the stoppping criterion;
C...              0,1,3 - compute the relative error using 
C...                      ||sol_k - sol_{k-1}|| / ||sol_k||;
C...              2     - compute the relative error using the computed 
C...                      exact soln, i.e., ||ex_sol - sol_k|| / ||ex_sol||;
C--------------------------------------------------------------------
cc      call wrsol01(sol,n,150)
C
      lf = 1
      k  = 1
      nd(1) = n
C
      do while (nd(k) .gt. 10)
         ndk = nd(k)
         iak = ia(k)
         jak = ja(k)
         kak = ka(k)
         nnzk = m(iak+ndk) - 1
C
C...     To determine whether Tarjan's ordering is considered or not.
C...     If ordering is considered, then set ICH_ORD = 1.
C
         if (ich_ord .eq. 1) then
            iordk  = ifree
            iaw    = iordk + ndk + 1
            jaw    = iaw + ndk + 1 + 1
            isub   = jaw + nnzk
            iwork1 = isub + ndk + 1
            iwork2 = iwork1 + ndk + 1
            iwork3 = iwork2 + ndk + 1
            ifree  = iaw
C     
            call cut_off(ndk,m(iak),m(jak),r(kak),m(iaw),m(jaw),
     >           m(iordk),m(isub),nblk,m(iwork1),m(iwork2),m(iwork3))
         else
            iordk  = ifree
            ifree  = iordk + ndk + 1
C     
            call iseqv(m(iordk),ndk)
         end if
C
         call match10(ndk,m(iak),m(jak),r(kak),
     >        nnzk,jpk,jptk,m(iordk),kmatch(k),
     >        ndc,nnzc,iac,jac,kac,
     >        m(ifree),r(kfree),ifree,kfree)
C
         iord(k) = iordk
         jpp(k)  = jpk
         nd(k+1) = ndc
         ia(k+1) = iac
         ja(k+1) = jac
         ka(k+1) = kac
C
         k = k + 1
cc         write(*,*) 'ndk, ndc, k', ndk, ndc, k 
      end do
C
      lc = k
C
      kwk1 = kfree
      kfree = kwk1 + n
C
      jiter = 0
      kfrold = kfree
      ifrold = ifree
C
C...  Compute the L^2 norm of the computed exact solution.
C
      if (iexact .eq. 2) then
         call l2norm(ex_sol,ex_sol_norm,n,1d0)
         if (ex_sol_norm .le. 1.d-13) then
            ex_sol_norm = 1.d0
            write(*,*) 'Zero exact solution is detected!!!!.'
         end if
      end if
C
      err_relres = 1.d0
      err_rel    = 1.d0
C
cc      nh  = dsqrt(dble(n))
cc      h   = 1/dble(nh)
cc      tol = 1.d-3/n
cc      write(*,*) 'n,sqrt(n),h,1.d-3*h,tol:',n,nh,h,1.d-3/nh,tol
      if (iexact .eq. 1)  tol = 1.d-13
      do while (err_rel .gt. tol .and. err_rel .lt.1.d20)
cc      do while ((err_relres .gt. tol .or. err_rel .gt. tol) 
cc     >     .and. err_relres .lt. 10.d0 .and. jiter .lt.5000)
cc      do while (err_rel .gt. tol)
C
         jiter = jiter + 1
C
         kuk = ku(1)
         kbk = kb(1)
         iak = ia(1)
         jak = ja(1)
         kak = ka(1)
C
         if (jiter .eq. 1) then
C...        To compute the residual: b = rhs - A*sol.
            call abyvam(rhs,m(iak),m(jak),r(kak),sol,n,r(kbk))
C
            if (iexact .eq. 2) then 
               write(*,'(2X,i5,2X,1(3X,e11.4))')  0,1.d0
            else
C...           To compute the L_2 norm of the initial residual.
               call l2norm(r(kbk),nr0,n,1d0)
               if (nr0 .le. 1.d-13) nr0 = 1.d0
               err_relres = 1.d0
               write(*,'(2X,i5,2X,4(3X,e11.4))')
     >              0,err_relres,nr0,1.d0,1.d0
            end if
         end if
C
C************************************************************************
C...     MG-cycle, i.e., the action of the iterator B: B*b.
C************************************************************************
C
C----------------------------------------------------------------------
C...     Going downward in the V-cycle, or (\)-cycle.
C----------------------------------------------------------------------
C
         do 50 k = 1, lc-1
            ndk = nd(k)
            iak = ia(k)
            jak = ja(k)
            kuk = ku(k)
            kbk = kb(k)
            kak = ka(k)
            iordk = iord(k)
C     
cc            call lpri(m(iac),m(jac),ndc)
C
C...        Zero initial guess.
            call nullv(r(kuk),ndk)
C
C...        Presmoothing by Gauss-seidel smoother nprsm times: u=u+R(b-Au)
C
            if (icycle .eq. 3) then                    !For variable V-cyle
               k2 = k/2
               noofsm = 2**k2*nprsm
            else
               noofsm = nprsm                              !For V or \-cyle
            end if
C
cc            write(*,*) 'No. of pre-smoothings:', k,noofsm
C
            if (iprsm .eq. 1) then 
               call fwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
            else if (iprsm .eq. 2) then 
               call bwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
            else if (iprsm .eq. 3) then 
               call sym_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
cc            else if (iprsm .eq. 4) then 
cc               call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc            else if (iprsm .eq. 5) then 
cc               call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc	      else if (iprsm .eq. 6) then 
cc               call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc               call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
            end if
C
cc            if (k .eq. 1) then
cc               call wruv_plus(sol,r(kuk),n,151)
cc            end if
C
C...        To compute the residual wk = b - Au.
            call abyvam(r(kbk),m(iak),m(jak),r(kak),r(kuk),ndk,r(kwk1))
C     
            ndc = nd(k+1)
            jpk = jpp(k)
C     
cc            write(*,*) kfree,kac+nnzc
cc            write(*,*) ifree,jac+nnzc
C
            kbc     = kfree
            kb(k+1) = kbc
            ku(k+1) = kbc + nd(k+1)
            kfree   = ku(k+1) + nd(k+1)
C
C...        Restriction to the lower level: b = P^t*wk = wk^t*P.
            call restr(m(jpk),r(kbc),r(kwk1),ndc,ndk)
 50      continue
C
C----------------------------------------------------------------------
C...     Solving exactly at the coarsest level, i.e., u = A^(-1)*b.
C----------------------------------------------------------------------
C     
         itmax = 300*ndc+1
         kuc   = ku(lc)
         iac   = ia(lc)
         jac   = ja(lc)
         kac   = ka(lc)
C     
cc         call lprr(m(iac),m(jac),r(kac),ndc)
C
C...     Zero initial guess.
         call nullv(r(kuc),ndk)
C
         if (isolv_coarse .eq. 1) then
            call fwd_gs_ns(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >           ndc,itmax)
cc            call fwd_gs_ns_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
cc     >           ndc,itmax,1)
         else if (isolv_coarse .eq. 2) then
            call bwd_gs_ns(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >           ndc,itmax)
         else if (isolv_coarse .eq. 3) then
            call sym_gs_ns(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >           ndc,itmax)
         else if (isolv_coarse .eq. 4) then
            maxa = ifree
            kau  = kfree
            call sgauss(m(iac),m(jac),r(kac),r(kbc),m(maxa),r(kau),ndc)
            call copyv(r(kbc),r(kuc),ndc)
cc         else if (isolv_coarse .eq. 5) then
cc            call cg_silent(m(iac),m(jac),r(kac),r(kuc),r(kbc),
cc     >           r(kfree),ndc,zero)
         end if
C
cc         read(*,*)

C----------------------------------------------------------------------
C...     Going upward in the V-cycle, or /-cycle.
C----------------------------------------------------------------------
C
         do 100 k = lc-1, 1, -1
            iak = ia(k)
            jak = ja(k)
            kuc = ku(k+1)
            kuk = ku(k)
            kbk = kb(k)
            kak = ka(k)
            ndc = nd(k+1)
            ndk = nd(k)
            jpk = jpp(k)
            iordk = iord(k)
C
cc            write(*,*) "k, jpp ", k,jpp(k)
cc            write(*,*) " kak ", kak
cc            read(*,*)
C
C...        Correction with prolongation from the lower level: 
C...        i.e., u_k = u_k + P*u_{k-1}.
C
            call prolong_plus(m(jpk),r(kuc),r(kuk),ndc,ndk)
C
cc            write(*,*) (r(kuk+i-1),i=1,ndk)
cc            read(*,*)
cc            write(*,*) kmatch
C
            if (icycle .eq. 2)  go to 70      !No post smoothing in \-cycle
C
C...        Postsmoothing by Gauss-Seidel smoother npssm times: u=u+R(b-Au)
C
            if (icycle .eq. 3) then                 !For variable V-cyle
               k2 = k/2
               noofsm = 2**k2*npssm
            else
               noofsm = npssm                       !For V or \-cyle
            end if 
C
cc            write(*,*) 'No. of post-smoothings:', noofsm
C
            if (ipssm .eq. 1) then 
               call fwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
            else if (ipssm .eq. 2) then 
               call bwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
            else if (ipssm .eq. 3) then 
               call sym_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
cc            else if (ipssm .eq. 4) then 
cc               call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc            else if (ipssm .eq. 5) then 
cc               call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc            else if (ipssm .eq. 6) then 
cc               call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc               call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
            end if
C
 70         continue
 100     continue
C
C************************************************************************
C...     End of MG-cycle, i.e., end of the action of the iterator B: B*b.
C************************************************************************
C
cc         if(jiter .eq. max_iter-1) then
cc            call wrsol01(r(ku(1)),n,152)
cc         end if
cc         if(jiter .eq. 4)stop
C
C...     Updating the solution: sol = sol + B*b.
         call uupluv(sol,r(ku(1)),n)
C
cc         call wrsol01(sol,n,153)
C
C...     To compute the relative error.
         if (iexact .eq. 2) then
            call wuminv(r(ku(1)),ex_sol,sol,n)
            call l2norm(r(ku(1)),xdiffn,n,1d0)
            err_rel = xdiffn/ex_sol_norm  
         else
            call l2norm(r(ku(1)),xdiffn,n,1d0)
            call l2norm(sol,xnewn,n,1d0)
            if (xnewn .lt. xdiffn .or. xnewn .lt. 1.d-4) then
               err_rel = xdiffn/nr0  
            else
               err_rel = xdiffn/xnewn
            end if
         end if
C
C...     To compute the residual: b = rhs - A*sol.
         call abyvam(rhs,m(iak),m(jak),r(kak),sol,n,r(kbk))
C
         if (jiter .gt. 100) then
            nmod = mod(jiter,100)
         else
            nmod = 1
         end if
C
         if (iexact .eq. 2) then
            if (nmod .eq. 1) then
               write(*,'(2X,i7,2X,1(3X,e11.4))')
     >              jiter,err_rel
            end if
         else
C...        To compute the L_2 norm of the residual.
            call l2norm(r(kbk),err_res,n,1d0)
            err_relres = err_res/nr0
            if (nmod .eq. 1) then
               write(*,'(2X,i7,2X,4(3X,e11.4))')
     >              jiter,err_relres,err_res,err_rel,xdiffn
            end if
         end if
C
         kfree = kfrold
         ifree = ifrold
C
      end do
C
C...  Print the error at the last iteration. 
C
      if (iexact .eq. 2) then
         write(*,'(a,i7,a,e11.4)') 
     >        ' LAST MG : ',jiter, ' Rel_err = ',err_rel 
      else
         write(*,'(2X,i7,2X,4(3X,e11.4))')
     >              jiter,err_relres,err_res,err_rel,xdiffn 
      end if
C     
      return
      end
C====================================================================
      subroutine mg_match_power(sol,rhs,n,max_iter,m,r,tol)
C=====================================================================
      implicit none
C
      include 'paramsmg.h'
C
      integer maximum_levels
      integer nd,iord, 
     > ia, ja, ka, 
     > ipp,jpp, kpp,
     > ku , kb, 
     > ifree , kfree,lf,lc
      parameter (maximum_levels = 100)
C
      common /point/ 
     > nd(maximum_levels),iord(maximum_levels), 
     > ia(maximum_levels), ja(maximum_levels), ka(maximum_levels), 
     > ipp(maximum_levels),jpp(maximum_levels), kpp(maximum_levels),
     > ku(maximum_levels) , kb(maximum_levels), 
     > ifree , kfree ,lf,lc

      double precision r(1),sol(1),rhs(1),zero,one,tol
      double precision eig_val_rel,sol_norm,eig_val_old,eig_val_new
      integer m(1)
      integer*4 kmatch(100),nmod,maxa,kau,kk,noofsm,k2
      integer i,k,iak,jak,kuk,kbc,kbk,kak,ndc,ndk,jpk,
     >     iordk,iend ,kwk1,kend,n,max_iter,itmax,
     >     iac,jac,kuc,kac,jiter,jptk,nnzk,nnzc,kfrold,ifrold,
     >     keig,l2_or_max,kkk
      integer iaw,jaw,isub,iwork1,iwork2,iwork3,nblk
      common /onezero/ zero,one
C---------------------------------------------------------------------
C...  MG_MATCH_POWER finds the spectral radius (the largest eigenvalue)
C...  of the operator I - BA by power method, where B is the iterator
C...  in the MG-cycle. Various multigrid method will be used: 
C...  V-cycle, \-cycle, variable V-cycle, etc.
C...
C...  Parameter:
C...    ICYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle
C...    IPRSM   - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPSSM   - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRSM   - number of pre-smoothings: 1,2,3, etc
C...    NPSSM   - number of post-smoothings: 1,2,3, etc
C...    LF      - level of the finest mesh, which is backward in order: 1
C...    LC      - level of the coarsest mesh: 2,3, etc
C...    ISOLV_COARSE  - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C...              5 = CGM
C...    ICH_ORD - To determine whether Tarjan's ordering is considered 
C...              or not. If ordering is considered, then it is 1.
C---------------------------------------------------------------------
cc      call wrsol01(sol,n,150)
C
      lf = 1
      k  = 1
      nd(1) = n
C
      do while (nd(k) .gt. 10)
         ndk = nd(k)
         iak = ia(k)
         jak = ja(k)
         kak = ka(k)
         nnzk = m(iak+ndk)-1
C
C...     To determine whether Tarjan's ordering is considered or not.
C...     If ordering is considered, then set ICH_ORD = 1.
C
         if (ich_ord .eq. 1) then
            iordk  = ifree
            iaw    = iordk + ndk + 1
            jaw    = iaw + ndk + 1 + 1
            isub   = jaw + nnzk
            iwork1 = isub + ndk + 1
            iwork2 = iwork1 + ndk + 1
            iwork3 = iwork2 + ndk + 1
            ifree  = iaw
C     
            call cut_off(ndk,m(iak),m(jak),r(kak),m(iaw),m(jaw),
     >           m(iordk),m(isub),nblk,m(iwork1),m(iwork2),m(iwork3))
         else
            iordk  = ifree
            ifree  = iordk + ndk + 1
C     
            call iseqv(m(iordk),ndk)
         end if
C
         call match10(ndk,m(iak),m(jak),r(kak),
     >        nnzk,jpk,jptk,m(iordk),kmatch(k),
     >        ndc,nnzc,iac,jac,kac,
     >        m(ifree),r(kfree),ifree,kfree)
C
         iord(k) = iordk
         jpp(k)  = jpk
         nd(k+1) = ndc
         ia(k+1) = iac
         ja(k+1) = jac
         ka(k+1) = kac
C
         k = k + 1
cc         write(*,*) 'ndk, ndc, k', ndk, ndc, k 
      end do
C
      lc = k
C
      jiter = 0
      keig = kfree
      kwk1 = keig + n 
      kfree = kwk1 + n
      kfrold = kfree
      ifrold = ifree
C
      eig_val_old = 1.d0
      eig_val_new = 1.d0
      eig_val_rel = 1.d0
C
C...  Choose the norm used: 0 - L^2 norm,  1 - Maximum norm.
      l2_or_max = 0              !The choice 1 doesn't work well.
C
      do while (eig_val_rel .gt. tol)
CC      do while (jiter. lt.100)
         kfree = kfrold
         ifree = ifrold
C
         jiter = jiter + 1
C
         kuk = ku(1)
         kbk = kb(1)
         iak = ia(1)
         jak = ja(1)
         kak = ka(1)
C
         if (jiter .eq. 1) then
            if (l2_or_max .eq. 0) then
C...           To compute the L^2 norm of the initial guess.
               call l2norm(sol,sol_norm,n,1d0)
            else
C...           To compute the L^{infty} norm of the initial guess.
               call c0norm(sol,sol_norm,n)
            end if
C
C...        Compute eig_vec = sol/sol_norm, i.e., normalization of sol.
            call usmultv(r(keig),sol,n,1.d0/sol_norm)
         end if
C
C...     To compute the residual: b = rhs - A*eig_vec.
         call abyvam(rhs,m(iak),m(jak),r(kak),r(keig),n,r(kbk))
C
C************************************************************************
C...     MG-cycle, i.e., the action of the iterator B: B*b.
C************************************************************************
C
C----------------------------------------------------------------------
C...     Going downward in the V-cycle, or (\)-cycle.
C----------------------------------------------------------------------
C
         do 50 k = 1, lc-1
            ndk = nd(k)
            iak = ia(k)
            jak = ja(k)
            kuk = ku(k)
            kbk = kb(k)
            kak = ka(k)
            iordk = iord(k)
            iend  = ifree
            kend  = kfree
C     
cc            call lpri(m(iac),m(jac),ndc)
C
C...        Zero initial guess.
            call nullv(r(kuk),ndk)
C
C...        Presmoothing by Gauss-seidel smoother nprsm times: u=u+R(b-Au)
C
            if (icycle .eq. 3) then                    !For variable V-cyle
               k2 = k/2
               noofsm = 2**k2*nprsm
            else
               noofsm = nprsm                              !For V or \-cyle
            end if
C
cc            write(*,*) 'No. of pre-smoothings:', k,noofsm
C
            if (iprsm .eq. 1) then 
               call fwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
            else if (iprsm .eq. 2) then 
               call bwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
            else if (iprsm .eq. 3) then 
               call sym_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
cc            else if (iprsm .eq. 4) then 
cc               call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc            else if (iprsm .eq. 5) then 
cc               call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc	      else if (iprsm .eq. 6) then 
cc               call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc               call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
            end if
C
cc            if (k .eq. 1) then
cc               call wruv_plus(sol,r(kuk),n,151)
cc            end if
C
C...        To compute the residual wk = b - Au.
            call abyvam(r(kbk),m(iak),m(jak),r(kak),r(kuk),ndk,r(kwk1))
C     
            ndc = nd(k+1)
            jpk = jpp(k)
C     
cc            write(*,*) kfree,kac+nnzc
cc            write(*,*) ifree,jac+nnzc
C
            kbc     = kfree
            kb(k+1) = kbc
            ku(k+1) = kbc + nd(k+1)
            kfree   = ku(k+1) + nd(k+1)
C
C...        Restriction to the lower level: b = P^t*wk = wk^t*P.
            call restr(m(jpk),r(kbc),r(kwk1),ndc,ndk)
 50      continue
C
C----------------------------------------------------------------------
C...     Solving exactly at the coarsest level, i.e., u = A^(-1)*b.
C----------------------------------------------------------------------
C     
         itmax = 2*ndc+1
         kuc   = ku(lc)
         iac   = ia(lc)
         jac   = ja(lc)
         kac   = ka(lc)
C     
cc         call lprr(m(iac),m(jac),r(kac),ndc)
C
C...     Zero initial guess.
         call nullv(r(kuc),ndk)
C
         if (isolv_coarse .eq. 1) then
            call fwd_gs_ns(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >           ndc,itmax)
         else if (isolv_coarse .eq. 2) then
            call bwd_gs_ns(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >           ndc,itmax)
         else if (isolv_coarse .eq. 3) then
            call sym_gs_ns(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >           ndc,itmax)
         else if (isolv_coarse .eq. 4) then
            maxa = iend
            kau  = kwk1
            call sgauss(m(iac),m(jac),r(kac),r(kbc),m(maxa),r(kau),ndc)
            call copyv(r(kbc),r(kuc),ndc)
cc         else if (isolv_coarse .eq. 5) then
cc            call cg_silent(m(iac),m(jac),r(kac),r(kuc),r(kbc),
cc     >           r(kfree),ndc,zero)
         end if
C     
C----------------------------------------------------------------------
C...     Going upward in the V-cycle, or /-cycle.
C----------------------------------------------------------------------
C
         do 100 k = lc-1, 1, -1
            iak = ia(k)
            jak = ja(k)
            kuc = ku(k+1)
            kuk = ku(k)
            kbk = kb(k)
            kak = ka(k)
            ndc = nd(k+1)
            ndk = nd(k)
            jpk = jpp(k)
            iordk = iord(k)
            iend  = ifree
            kend  = kfree
C
cc            write(*,*) "k, jpp ", k,jpp(k)
cc            write(*,*) " kak ", kak
cc            read(*,*)
C
C...        Correction with prolongation from the lower level: 
C...        i.e., u_k = u_k + P*u_{k-1}.
C
            call prolong_plus(m(jpk),r(kuc),r(kuk),ndc,ndk)
C
cc            write(*,*) (r(kuk+i-1),i=1,ndk)
cc            read(*,*)
cc            write(*,*) kmatch
C
            if (icycle .eq. 2)  go to 70       !No post smoothing in \-cycle
C
C...        Postsmoothing by Gauss-Seidel smoother npssm times: u=u+R(b-Au)
C
            if (icycle .eq. 3) then                 !For variable V-cyle
               k2 = k/2
               noofsm = 2**k2*npssm
            else
               noofsm = npssm                       !For V or \-cyle
            end if 
C
cc            write(*,*) 'No. of post-smoothings:', noofsm
C
            if (ipssm .eq. 1) then 
               call fwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
            else if (ipssm .eq. 2) then 
               call bwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
            else if (ipssm .eq. 3) then 
               call sym_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >              r(kbk),m(iordk),ndk,noofsm) 
cc            else if (ipssm .eq. 4) then 
cc               call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc            else if (ipssm .eq. 5) then 
cc               call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc            else if (ipssm .eq. 6) then 
cc               call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc               call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >              r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
            end if
C     
 70         continue
 100     continue
C
C************************************************************************
C...     End of MG-cycle, i.e., end of the action of the iterator B: B*b.
C************************************************************************
C
cc         if(jiter .eq. max_iter-1) then
cc            call wrsol01(r(ku(1)),n,152)
cc         end if
cc         if(jiter .eq. 4)stop
C
C...     Updating the solution: sol = eig_vec + B*b.
         call wupluv(sol,r(keig),r(ku(1)),n)
C
cc         call wrsol01(sol,n,153)
C
         if (l2_or_max .eq. 0) then
C...        To compute the L^2 norm of the updated solution.
            call l2norm(sol,sol_norm,n,1d0)
         else
C...        To compute the L^{infty} norm of the updated solution.
            call c0norm(sol,sol_norm,n)
         end if
C
         eig_val_old = eig_val_new
C
C...     Compute the scalar product eig_val_new = (r(keig),sol)
         call scpro(r(keig),sol,eig_val_new,n)
C
         eig_val_rel = dabs((eig_val_old - eig_val_new)/eig_val_new)
C
         if (mod(jiter,2) .eq. 0) then
            write(*,'(2X,i5,2X,2(3X,e17.10))') 
     >           jiter,eig_val_new,eig_val_rel
         end if
C
C...     Compute eig_vec = sol/sol_norm, i.e., normalization of sol.
         call usmultv(r(keig),sol,n,1.d0/sol_norm)
      end do
C
      write(*,'(2X,i5,2X,2(3X,e17.10))') 
     >     jiter,eig_val_new,eig_val_rel 

      return
      end
C====================================================================
      subroutine cg_silent(ia,ja,a,yn,b,ra,n,tol)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),a(1),b(1),yn(1),ra(1)
C--------------------------------------------------------------------
C...  Conjugate gradient method without write statements.
C--------------------------------------------------------------------
      if(n .eq. 1) then
         yn(1) = b(1)/a(1)
         return
      end if
      maxitr = n*n
      iter = 0
C
C...  Initial guess
C
      irn = 1
      izn = irn + n
      iapn = izn + n
      ipn = iapn + n
C      
      call abyvam(b,ia,ja,a,yn,n,ra(irn))
C
      call l2norm(ra(irn),res1,n,1d0)
      resss = res1
C
C.... z0 = B*r0
C
      call copyv(ra(irn),ra(izn),n)
C
C...  p0 = z0
      call copyv(ra(izn),ra(ipn),n)
C
      call scpro(ra(izn),ra(irn),resin,n)
      resin = dsqrt(resin)
C
      if(resin .lt. tol) then
         resss = 1
         res = resin
         go to 10
      end if
C
 1    iter = iter + 1
C
C...  apn = A*pn
C
      call abyvg(ia,ja,a,ra(ipn),n,ra(iapn))
C
C.... alphau = <ra(irn),ra(izn)>    
C     
      call scpro(ra(irn),ra(izn),alphau,n)
C
C.... alphad = <Apn,pn>    
C     
      call scpro(ra(iapn),ra(ipn),alphad,n)
C
      alpha = alphau/alphad
C
      call uuplmv(yn,ra(ipn),n,alpha)
      call uuplmv(ra(irn),ra(iapn),n,-alpha)
C
      call copyv(ra(irn),ra(izn),n)
C
C...  betad = alphau
C...  betau = <zn,rn>
C
C...  check convergence
C
      call scpro(ra(izn),ra(irn),betau,n)
      res = dsqrt(betau)/resin
      call l2norm(ra(irn),res1,n,1d0)
      if(res .lt. tol) go to 10
C
      beta = -betau/alphau
C
      call umuplv(ra(ipn),ra(izn),n,-beta)
C     
      if(iter .eq. maxitr) go to 20
C----------------------END OF ITERATION LOOP----------------------
      go to 1
C-----------------------------------------------------------------
 10   continue
C
C        write(*,'(3X,i4,5X,4(2X, e14.6,1X))')
C    >     iter, res,res*resin,res1/resss,res1*n
C
      return
C
 20   write(*,*) ' !!!! ', iter, 
     >   '  residual = ',res
C
      return
      end
C=====================================================================
C=====================================================================
      subroutine mg_match_modify(sol,rhs,n,max_iter,m,r,
     >     tol,ex_sol,iexact)
C=====================================================================
      implicit none
C
      include 'paramsmg.h'
C
      integer maximum_levels
      integer nd,iord, 
     > ia, ja, ka, 
     > ipp,jpp, kpp,
     > ku , kb, 
     > ifree , kfree,lf,lc
      parameter (maximum_levels = 100)
C
      common /point/ 
     > nd(maximum_levels), iord(maximum_levels), 
     > ia(maximum_levels), ja(maximum_levels), ka(maximum_levels), 
     > ipp(maximum_levels),jpp(maximum_levels), kpp(maximum_levels),
     > ku(maximum_levels) , kb(maximum_levels), 
     > ifree , kfree ,lf,lc

      double precision r(1),sol(1),rhs(1),ex_sol(1),zero,one,nr0
      double precision xdiffn,xnewn,err_rel,err_res,err_relres,tol,h
      double precision ex_sol_norm
      integer m(1)
      integer*4 kmatch(100),nmod
      integer i,k,iak,jak,kuk,kbk,kak,
     >     n,max_iter,iexact,
     >     jiter
      integer nh
      common /onezero/ zero,one
C---------------------------------------------------------------------
C...  MG_MATCH_MODIFY is identical to MG_MATCH, and it performs
C...  u^{k+1} = u^k + B(f - A*u^k), where B an iterator B in an MG-cycle. 
C...  Various multigrid method will be used: V-cycle, \-cycle, variable
C...  V-cycle, etc.
C...
C...  Parameter:
C...    ICYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle, 
C...              3 = variable V-cycle
C...    IPRSM   - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPSSM   - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRSM   - number of pre-smoothings: 1,2,3, etc
C...    NPSSM   - number of post-smoothings: 1,2,3, etc
C...    LF      - level of the finest mesh, which is backward in order: 1
C...    LC      - level of the coarsest mesh: 2,3, etc
C...    ISOLV_COARSE  - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C...              5 = CGM
C...    ICH_ORD - To determine whether Tarjan's ordering is considered 
C...              or not. If ordering is considered, then it is 1.
C...    IEXACT  - To determine the stoppping criterion;
C...              0,1,3 - compute the relative error using 
C...                      ||sol_k - sol_{k-1}|| / ||sol_k||;
C...              2     - compute the relative error using the computed 
C...                      exact soln, i.e., ||ex_sol - sol_k|| / ||ex_sol||;
C--------------------------------------------------------------------
cc      call wrsol01(sol,n,150)
C
C...  Generate lower(coarse) level matrices using graph matching.
C
      call gen_coarse_mat(n,m(1),r(1))
C
      jiter = 0
C
C...  Compute the L^2 norm of the computed exact solution.
C
      if (iexact .eq. 2) then
         call l2norm(ex_sol,ex_sol_norm,n,1d0)
         if (ex_sol_norm .le. 1.d-13) then
            ex_sol_norm = 1.d0
            write(*,*) 'Zero exact solution is detected!!!!.'
         end if
      end if
C
      err_relres = 1.d0
      err_rel    = 1.d0
C
cc      nh  = dsqrt(dble(n))
cc      h   = 1/dble(nh)
cc      tol = 1.d-3/n
cc      write(*,*) 'n,sqrt(n),h,1.d-3*h,tol:',n,nh,h,1.d-3/nh,tol
      if (iexact .eq. 1)  tol = 1.d-13
C
      do while (err_rel .gt. tol .and. err_rel .lt.1.d20)
cc      do while ((err_relres .gt. tol .or. err_rel .gt. tol) 
cc     >     .and. err_relres .lt. 10.d0 .and. jiter .lt.5000)
cc      do while (err_rel .gt. tol)
C
         jiter = jiter + 1
C
         kuk = ku(1)
         kbk = kb(1)
         iak = ia(1)
         jak = ja(1)
         kak = ka(1)
C
         if (jiter .eq. 1) then
C...        To compute the residual: b = rhs - A*sol.
            call abyvam(rhs,m(iak),m(jak),r(kak),sol,n,r(kbk))
C
            if (iexact .eq. 2) then 
               write(*,'(2X,i5,2X,1(3X,e11.4))')  0,1.d0
            else
C...           To compute the L_2 norm of the initial residual.
               call l2norm(r(kbk),nr0,n,1d0)
               if (nr0 .le. 1.d-13) nr0 = 1.d0
               err_relres = 1.d0
               write(*,'(2X,i5,2X,4(3X,e11.4))')
     >              0,err_relres,nr0,1.d0,1.d0
            end if
         end if
C
C------------------------------------------------------------------------
C...     Call MG-cycle, i.e., the action of the iterator B: B*b.
C------------------------------------------------------------------------
         call premg_match(n,m(1),r(1))
C------------------------------------------------------------------------
C...     End of MG-cycle, i.e., end of the action of the iterator B: B*b.
C------------------------------------------------------------------------
C
cc         if(jiter .eq. max_iter-1) then
cc            call wrsol01(r(ku(1)),n,152)
cc         end if
cc         if(jiter .eq. 4)stop
C
C...     Updating the solution: sol = sol + B*b.
         call uupluv(sol,r(ku(1)),n)
C
cc         call wrsol01(sol,n,153)
C
C...     To compute the relative error.
         if (iexact .eq. 2) then
            call wuminv(r(ku(1)),ex_sol,sol,n)
            call l2norm(r(ku(1)),xdiffn,n,1d0)
            err_rel = xdiffn/ex_sol_norm  
         else
            call l2norm(r(ku(1)),xdiffn,n,1d0)
            call l2norm(sol,xnewn,n,1d0)
            if (xnewn .lt. xdiffn .or. xnewn .lt. 1.d-4) then
               err_rel = xdiffn/nr0  
            else
               err_rel = xdiffn/xnewn
            end if
         end if
C
C...     To compute the residual: b = rhs - A*sol.
         call abyvam(rhs,m(iak),m(jak),r(kak),sol,n,r(kbk))
C
         if (jiter .gt. 100) then
            nmod = mod(jiter,100)
         else
            nmod = 1
         end if
C
         if (iexact .eq. 2) then
            if (nmod .eq. 1) then
               write(*,'(2X,i7,2X,1(3X,e11.4))')
     >              jiter,err_rel
            end if
         else
C...        To compute the L_2 norm of the residual.
            call l2norm(r(kbk),err_res,n,1d0)
            err_relres = err_res/nr0
            if (nmod .eq. 1) then
               write(*,'(2X,i7,2X,4(3X,e11.4))')
     >              jiter,err_relres,err_res,err_rel,xdiffn
            end if
         end if
C
      end do
C
C...  Print the error at the last iteration. 
C
      if (iexact .eq. 2) then
         write(*,'(a,i7,a,e11.4)') 
     >        ' LAST MG : ',jiter, ' Rel_err = ',err_rel 
      else
         write(*,'(2X,i7,2X,4(3X,e11.4))')
     >              jiter,err_relres,err_res,err_rel,xdiffn 
      end if
C
      return
      end
C=====================================================================
      subroutine gen_coarse_mat(n,m,r)
C=====================================================================
      implicit none
C
      include 'paramsmg.h'
C
      integer maximum_levels
      integer nd,iord, 
     > ia, ja, ka, 
     > ipp, jpp, kpp,
     > ku, kb, 
     > ifree, kfree, lf, lc
      parameter (maximum_levels = 100)
C
      common /point/ 
     > nd(maximum_levels), iord(maximum_levels), 
     > ia(maximum_levels), ja(maximum_levels), ka(maximum_levels), 
     > ipp(maximum_levels), jpp(maximum_levels), kpp(maximum_levels),
     > ku(maximum_levels), kb(maximum_levels), 
     > ifree, kfree, lf, lc

      double precision r(1)
      integer m(1)
      integer*4 kmatch(100)
      integer k,iak,jak,kak,ndc,ndk,jpk,
     >     iordk,n,iac,jac,kac,jptk,nnzk,nnzc
      integer iaw,jaw,isub,iwork1,iwork2,iwork3,nblk
C---------------------------------------------------------------------
C...  Generate lower(coarse) level matrices using graph matching
C...  and other information. 
C...  The following finite element method can be applied: 
C...  Standard Galerkin, EAFE (Edge Average Finite Element).
C...  These will be used in the subroutine PREMG_MATCH
C...
C...  Parameter:
C...    LF      - level of the finest mesh, which is backward in order: 1
C...    LC      - level of the coarsest mesh: 2,3, etc
C...    ICH_ORD - To determine whether Tarjan's ordering is considered 
C...              or not. If ordering is considered, then it is 1.
C--------------------------------------------------------------------
      lf = 1
      k  = 1
      nd(1) = n
C
      do while (nd(k) .gt. 10)
         ndk = nd(k)
         iak = ia(k)
         jak = ja(k)
         kak = ka(k)
         nnzk = m(iak+ndk) - 1
C
C...     To determine whether Tarjan's ordering is considered or not.
C...     If ordering is considered, then set ICH_ORD = 1.
C
         if (ich_ord .eq. 1) then
            iordk  = ifree
            iaw    = iordk + ndk + 1
            jaw    = iaw + ndk + 1 + 1
            isub   = jaw + nnzk
            iwork1 = isub + ndk + 1
            iwork2 = iwork1 + ndk + 1
            iwork3 = iwork2 + ndk + 1
            ifree  = iaw
C     
            call cut_off(ndk,m(iak),m(jak),r(kak),m(iaw),m(jaw),
     >           m(iordk),m(isub),nblk,m(iwork1),m(iwork2),m(iwork3))
         else
            iordk  = ifree
            ifree  = iordk + ndk + 1
C     
            call iseqv(m(iordk),ndk)
         end if
C
         call match10(ndk,m(iak),m(jak),r(kak),
     >        nnzk,jpk,jptk,m(iordk),kmatch(k),
     >        ndc,nnzc,iac,jac,kac,
     >        m(ifree),r(kfree),ifree,kfree)
C
         iord(k) = iordk
         jpp(k)  = jpk
         nd(k+1) = ndc
         ia(k+1) = iac
         ja(k+1) = jac
         ka(k+1) = kac
C
         k = k + 1
cc         write(*,*) 'ndk, ndc, k', ndk, ndc, k 
      end do
C
      lc = k
C
      return
      end
C=====================================================================
      subroutine premg_match(n,m,r)
C=====================================================================
      implicit none
C
      include 'paramsmg.h'
C
      integer maximum_levels
      integer nd,iord, 
     > ia, ja, ka, 
     > ipp,jpp, kpp,
     > ku, kb, 
     > ifree, kfree, lf, lc
      parameter (maximum_levels = 100)
C
      common /point/ 
     > nd(maximum_levels), iord(maximum_levels), 
     > ia(maximum_levels), ja(maximum_levels), ka(maximum_levels), 
     > ipp(maximum_levels),jpp(maximum_levels), kpp(maximum_levels),
     > ku(maximum_levels), kb(maximum_levels), 
     > ifree, kfree, lf, lc
C
      double precision r(1)
      integer m(1)
      integer*4 kmatch(100),maxa,kau,noofsm,k2
      integer k,iak,jak,kuk,kbc,kbk,kak,ndc,ndk,jpk,
     >     iordk,kwk1,n,itmax,
     >     iac,jac,iordc,kuc,kac,nnzc,kfrold,ifrold
C---------------------------------------------------------------------
C...  PREMG_MATCH stands for a multigrid cycle as a preconditioner, 
C...  and it performs the action u = B*rhs, where B is an iterator B  
C...  in an MG-cycle. Various multigrid method will be used: V-cycle, 
C...  \-cycle, variable, V-cycle, etc.
C...
C...  Parameter:
C...    ICYCLE  - choice of MG: 1 = V-cycle, 2 = \-cycle, 
C...              3 = variable V-cycle
C...    IPRSM   - choice of pre-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    IPSSM   - choice of post-smoothing: 
C...              1 = fwd G-S, 2 = bwd G-S, 3 = sym G-S
C...    NPRSM   - number of pre-smoothings: 1,2,3, etc
C...    NPSSM   - number of post-smoothings: 1,2,3, etc
C...    LF      - level of the finest mesh, which is backward in order: 1
C...    LC      - level of the coarsest mesh: 2,3, etc
C...    ISOLV_COARSE  - choice of the coarsest grid solver: 1 = fwd G-S,
C...              2 = bwd G-S, 3 = sym G-S, 4 = Gaussian elimination
C...              5 = CGM
C...    ICH_ORD - To determine whether Tarjan's ordering is considered 
C...              or not. If ordering is considered, then it is 1.
C--------------------------------------------------------------------
cc      call wrsol01(sol,n,150)
C
      ifrold = ifree
      kfrold = kfree
C
      kwk1 = kfree
      kfree = kwk1 + n
C
C------------------------------------------------------------------------
C...  MG-cycle, i.e., the action of the iterator B: B*b.
C------------------------------------------------------------------------
C
C----------------------------------------------------------------------
C...  Going downward in the V-cycle, or (\)-cycle.
C----------------------------------------------------------------------

C
      do 50 k = 1, lc-1
         ndk = nd(k)
         iak = ia(k)
         jak = ja(k)
         kuk = ku(k)
         kbk = kb(k)
         kak = ka(k)
         iordk = iord(k)
C     
cc         call lpri(m(iac),m(jac),ndc)
C
C...     Zero initial guess.
         call nullv(r(kuk),ndk)
C
C...     Presmoothing by Gauss-seidel smoother nprsm times: u=u+R(b-Au)
C
         if (icycle .eq. 3) then !For variable V-cyle
            k2 = k/2
            noofsm = 2**k2*nprsm
         else
            noofsm = nprsm      !For V or \-cyle
         end if
C
cc         write(*,*) 'No. of pre-smoothings:', k,noofsm
C
         if (iprsm .eq. 1) then 
            call fwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >           r(kbk),m(iordk),ndk,noofsm) 
         else if (iprsm .eq. 2) then 
            call bwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >           r(kbk),m(iordk),ndk,noofsm) 
         else if (iprsm .eq. 3) then 
            call sym_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >           r(kbk),m(iordk),ndk,noofsm) 
cc         else if (iprsm .eq. 4) then 
cc            call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >           r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc         else if (iprsm .eq. 5) then 
cc            call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >           r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc	   else if (iprsm .eq. 6) then 
cc            call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >           r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc            call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >           r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
         end if
C
cc         if (k .eq. 1) then
cc            call wruv_plus(sol,r(kuk),n,151)
cc         end if
C
C...     To compute the residual wk = b - Au.
         call abyvam(r(kbk),m(iak),m(jak),r(kak),r(kuk),ndk,r(kwk1))
C     
         ndc = nd(k+1)
         jpk = jpp(k)
C     
cc         write(*,*) kfree,kac+nnzc
cc         write(*,*) ifree,jac+nnzc
C
         kbc     = kfree
         kb(k+1) = kbc
         ku(k+1) = kbc + nd(k+1)
         kfree   = ku(k+1) + nd(k+1)
C
C...     Restriction to the lower level: b = P^t*wk = wk^t*P.
         call restr(m(jpk),r(kbc),r(kwk1),ndc,ndk)
 50   continue
C
C----------------------------------------------------------------------
C...  Solving exactly at the coarsest level, i.e., u = A^(-1)*b.
C----------------------------------------------------------------------
C     
      itmax = 300*ndc+1
      kuc   = ku(lc)
      iac   = ia(lc)
      jac   = ja(lc)
      kac   = ka(lc)
C     
cc      call lprr(m(iac),m(jac),r(kac),ndc)
C
C...  Zero initial guess.
      call nullv(r(kuc),ndk)
C
      if (isolv_coarse .eq. 1) then
         call fwd_gs_ns(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        ndc,itmax)
cc         call fwd_gs_ns_wr(r(kuc),m(iac),m(jac),r(kac),r(kbc),
cc     >        ndc,itmax,1)
      else if (isolv_coarse .eq. 2) then
         call bwd_gs_ns(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        ndc,itmax)
      else if (isolv_coarse .eq. 3) then
         call sym_gs_ns(r(kuc),m(iac),m(jac),r(kac),r(kbc),
     >        ndc,itmax)
      else if (isolv_coarse .eq. 4) then
         maxa = ifree
         kau  = kfree
         call sgauss(m(iac),m(jac),r(kac),r(kbc),m(maxa),r(kau),ndc)
         call copyv(r(kbc),r(kuc),ndc)
cc      else if (isolv_coarse .eq. 5) then
cc         call cg_silent(m(iac),m(jac),r(kac),r(kuc),r(kbc),
cc     >        r(kfree),ndc,zero)
      end if
C
C----------------------------------------------------------------------
C...  Going upward in the V-cycle, or /-cycle.
C----------------------------------------------------------------------
C
      do 100 k = lc-1, 1, -1
         iak = ia(k)
         jak = ja(k)
         kuc = ku(k+1)
         kuk = ku(k)
         kbk = kb(k)
         kak = ka(k)
         ndc = nd(k+1)
         ndk = nd(k)
         jpk = jpp(k)
         iordk = iord(k)
C
C...     Correction with prolongation from the lower level: 
C...     i.e., u_k = u_k + P*u_{k-1}.
C
         call prolong_plus(m(jpk),r(kuc),r(kuk),ndc,ndk)
C
         if (icycle .eq. 2)  go to 70 !No post smoothing in \-cycle
C
C...     Postsmoothing by Gauss-Seidel smoother npssm times: u=u+R(b-Au)
C
         if (icycle .eq. 3) then                 !For variable V-cyle
            k2 = k/2
            noofsm = 2**k2*npssm
         else
            noofsm = npssm                       !For V or \-cyle
         end if 
C
cc         write(*,*) 'No. of post-smoothings:', noofsm
C
         if (ipssm .eq. 1) then 
            call fwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >           r(kbk),m(iordk),ndk,noofsm) 
         else if (ipssm .eq. 2) then 
            call bwd_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >           r(kbk),m(iordk),ndk,noofsm) 
         else if (ipssm .eq. 3) then 
            call sym_gs_ns_ord(r(kuk),m(iak),m(jak),r(kak),
     >           r(kbk),m(iordk),ndk,noofsm) 
cc         else if (ipssm .eq. 4) then 
cc            call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >           r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc         else if (ipssm .eq. 5) then 
cc            call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >           r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc         else if (ipssm .eq. 6) then 
cc            call fwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >           r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
cc            call bwd_gs_ns_ord_blk(r(kuk),m(iak),m(jak),r(kak),
cc     >           r(kbk),m(iordk),ndk,kmatch(k),noofsm) 
         end if
C
 70      continue
 100  continue
C
C------------------------------------------------------------------------
C...     End of MG-cycle, i.e., end of the action of the iterator B: B*b.
C------------------------------------------------------------------------
C
      ifree = ifrold
      kfree = kfrold
C
      return
      end
C=====================================================================
      subroutine msolve(n, rhs, sol, nnz, iao, jao, ao, isym, r, m)
C=====================================================================
      implicit real*8(a-h,o-z)
C
      integer n, nnz, iao(1), jao(1), isym, m(1)
      double precision rhs(1), sol(1), ao(1), r(1)
C
      integer maximum_levels
      integer nd,iord, 
     > ia, ja, ka, 
     > ipp,jpp, kpp,
     > ku, kb, 
     > ifree, kfree, lf, lc
      parameter (maximum_levels = 100)
C
      common /point/ 
     > nd(maximum_levels), iord(maximum_levels), 
     > ia(maximum_levels), ja(maximum_levels), ka(maximum_levels), 
     > ipp(maximum_levels),jpp(maximum_levels), kpp(maximum_levels),
     > ku(maximum_levels), kb(maximum_levels), 
     > ifree, kfree, lf, lc
C--------------------------------------------------------------------
C...  MSOLVE solves a linear system P*sol = rhs for sol given rhs with
C...  the preconditioning matrix P (P is supplied via R and M arrays). 
C...  (As an example, it is used as an external subroutine for DGMRES.)
C...
C...  Parameters:
C...    IAO,JAO,AO - contain the matrix data structure for A (It could 
C...                 take any form.)
C...    N       - the number of unknowns
C...    RHS     - the right-hand side vector
C...    SOL     - the solution upon return
C...    NNZ     - the number of non-zeros in the SLAP IAO, JAO, AO storage 
C...              for the matrix A. (It could take any form.)
C...    ISYM    - a flag which, if non-zero, denotes that A is symmetric
C...              and only the lower or upper triangle is stored. If it
C...              is zero, all non-zero entries of the matrix are stored.
C...    R       - can be used to pass necessary preconditioning 
C...              information and/or workspace 
C...    M       - the same purpose as R.
C--------------------------------------------------------------------
C
      call copyv(rhs,r(kb(lf)),n)
C     
C...  Iterator B using an MG-cycle, i.e., u = B*b.
C     
      call premg_match(n,m(1),r(1))
C     
C...  Copy the solution.
      call copyv(r(ku(lf)),sol,n)
C
      return
      end
C=====================================================================
      subroutine matvec(n,x,y,nelt,ia,ja,a,isym)
C=====================================================================
      integer n, nelt, ia(1), ja(1), isym
      double precision x(n), y(n), a(1)
cc      integer i, iaa, iab, k
cc      double precision u
C--------------------------------------------------------------------
C...  PRODUCT--- GENERAL SPARSE MATRIX BY VECTOR: y = A*x.
C...  MATVEC is used as an external subroutine for DGMRES subroutine.
C...
C...  Parameters:
C...    IA,JA,A - contain the matrix data structure for A (It could 
C...                take any form.)
C...    Y       - the product A*X
C...    X       - an input vector
C...    NELT    - the number of non-zeros in the SLAP IA, JA, A 
C...              storage for the matrix A. 
C...    ISYM    - a flag which, if non-zero, denotes that A is symmetric
C...              and only the lower or upper triangle is stored. If it
C...              is zero, all non-zero entries of the matrix are stored.
C--------------------------------------------------------------------
      call abyvg(ia,ja,a,x,n,y)
C
cc      do i = 1, n
cc        u = 0.0d00
cc         iaa = ia(i)
cc         iab = ia(i+1) - 1
cc         if(iab .ge. iaa) then
cc            do k = iaa, iab
cc               u = u + a(k) * x(ja(k))
cc            end do
cc         end if
cc         y(i) = u
cc      end do
C
      return
      end
C=====================================================================
