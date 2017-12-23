C=====================================================================
      program TGMG
C=====================================================================
      implicit real*8(a-h,o-z)
      parameter (iter_tg_max = 12)
cc      parameter (nlast = 2 000 000, nrlast = 6 000 000) !level = 8
cc      parameter (nlast = 10 000 000, nrlast =20 000 000) !level = 9
cc      parameter (nlast = 40 000 000, nrlast = 55 000 000) !level = 10
cc      parameter (nlast = 130 000 000, nrlast = 105 000 000) !level = 11
      parameter (nlast = 141 000 000, nrlast = 100 000 000) !level = 11
      common /top/ m(nlast)              !To use the machine witt
      common /top1/ r(nrlast)            !To use the machine witt
cc      dimension m(nlast), r(nrlast)
      character*100 meshf(0:20)
      common /quad/ wg(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /quad_gs_ntrl/ z_ntrl(7),zw(5,7)
      common /quad_gs/ w(7),z(7)
      common /menu/ ich(200)
      common /fname/ meshf
      parameter (lmx = 20)
      dimension multl(lmx),mcont(0:128)
      dimension n1(lmx),n2(lmx),nel(lmx),nd(lmx),ned(lmx)
      dimension ia(lmx),ja(lmx),idir(lmx),nz(lmx)
      dimension ku(lmx),kb(lmx),ka(lmx),kal(lmx),krhs(lmx),ksol(lmx)
C---------------------------------------------------------------------
C...  This program solves the 2nd order linear elliptic PDEs. A 
C...  multilevel version of two grid method is considered. It is similar
C...  to the standard MG W-cycle schem. SPD systems
C...  are solved by the standard MG V-cycle scheme on all the levels 
C...  except the coarsest level where a non-SPD system is solved by the
C...  Gaussian elimination. A Standard Galerkin scheme is used.
C...
C...  Parameters:
C...    ISTIFF - choice of stiffness matrices.
C...             1 = Standard Galerkin, (2 = EAFE.)
C...    LVL, KEND, JCONT, LMT   
C...           - to control the levels of W-cycle type scheme.
C---------------------------------------------------------------------
      call input_data
      call quad_elt
      call quad_data
      call quad_data1
C 
      lvl = ich(8)
      do 10 i = 1, lvl
         multl(i) = ich(40 + i)
 10   continue
C
C...  Call LVLCONTROL to control the levels of W-cycle type multigrid 
C...  version of two-grid method.
C
      call wcycle_lvl(lvl,kend,mcont)
C
      mt      = 15
      lpde    = ich(1)
      jdf     = ich(3)
      lmethod = ich(4)
      lread   = ich(5)
      lf      = multl(lvl)
      lcmg    = ich(37)
      tol     = 1.d-12
C
      write(*,1000) lread,lpde,lmethod,lf,lcmg
 1000 format(2x,' lread   = ', i7 / 2x,
     >          ' lpde    = ', i7 / 2x,
     >          ' lmethod = ', i7 / 2x,
     >          ' lf      = ', i7 / 2x,
     >          ' lc      = ', i7)
C
      write(*,*) 
      write(*,*) '================================================'
      write(*,'(2(a,i3))') ' Fine Level:',lf,'  Coarse Level:',multl(1)
      write(*,*) '================================================'
      write(*,*)
C
C...  Generating the stiffness matrices ^A, A, and the load vector B 
C...  at each level.
C
      ifree = 1
      kfree = 1
      lasti = ifree
      lastr = kfree
C
      do 50 i = lvl,1,-1
         if(lread .eq. 1) then
            open(mt,file=meshf(multl(i)),
     >           status='unknown',form='formatted')
            read(mt,*) n1(i),n2(i)
            read(mt,*) nel(i),nd(i),ned(i)
         else
            open(mt,file=meshf(multl(i)),
     >           status='unknown',form='unformatted')
            read(mt) n1(i),n2(i)
            read(mt) nel(i),nd(i),ned(i)
         end if
C     
         write(*,*)
     >        ' No. of elements, nodes, & bdry edges:',
     >        nel(i),nd(i),ned(i)
C     
C...     Compute the addresses of each array in the arrays m and r.
C
          if (i .eq. lvl) then
            ief      = lasti
            jef      = ief + nel(i) + 1
            lasti    = jef + nel(i)*3
C
            kxf      = lastr
            kyf      = kxf + nd(i)
            lastr    = kyf + nd(i)
         end if
C
         idir(i) = lasti
         ia(i)   = idir(i) + nd(i)
         ja(i)   = ia(i) + nd(i) + 1
         lasti   = ja(i) + 2*(nel(i) + nd(i) + 5) + nd(i)
C
         if (i .eq. lvl) then
            inedg   = lasti
         else
            iec     = lasti
            jec     = iec + nel(i) + 1
            inedg   = jec + nel(i)*3
         end if
C
         iet     = inedg + ned(i) * 2
         jet     = iet + nd(i) + 1
         jat     = jet + nel(i)*3
         iend    = jat + 2*(nel(i) + nd(i) + 5) + nd(i)
C
         ku(i)   = lastr
         kb(i)   = ku(i) + nd(i)
         ksol(i) = kb(i) + nd(i)
         krhs(i) = ksol(i) + nd(i)
         ka(i)   = krhs(i) + nd(i)
         if (i .gt. 1) then
            kal(i) = ka(i) + 2*(nel(i) + nd(i) + 5) + nd(i)
            lastr  = kal(i) + 2*(nel(i) + nd(i) + 5) + nd(i)
         else
            lastr  = ka(i) + 2*(nel(i) + nd(i) + 5) + nd(i)
         end if
C
         if (i .eq. lvl) then
            krat = lastr
         else
            kxc  = lastr
            kyc  = kxc + nd(i)
            krat = kyc + nd(i)
         end if
         kend    = krat + 2*(nel(i) + nd(i) + 5) + nd(i)
C
         write(*,*) 'lasti,lastr,iend,kend', lasti,lastr,iend,kend
C
         if (i .eq. lvl) then
            ie = ief
            je = jef
            kx = kxf
            ky = kyf
         else
            ie = iec
            je = jec
            kx = kxc
            ky = kyc
         end if
C
         call rdmesh(lread,mt,m(ie),m(je),ned(i),m(inedg),
     >        m(idir(i)),r(kx),r(ky),nel(i),nd(i),jdf)
         close(mt)
C
         call iit(m(ie),m(je),nel(i),nd(i),m(iet),m(jet))
C
         nz(i) = 0
         call smbasg(m(ie),m(je),m(iet),m(jet),m(ia(i)),m(ja(i)),
     >        nd(i),nz(i),m(idir(i)))
C
C...     To compute stiffness matrix ^A at each level.
C
         iip  = iet
         call pde_choice(lmethod,
     I        nel(i),nd(i),jdf,m(ie),m(je),r(kx),r(ky),
     O        m(ia(i)),m(ja(i)),r(ka(i)),nz(i),r(krhs(i)),
     W        ned(i),m(inedg),m(idir(i)),m(iip))
C
         lasti = ja(i) + nz(i)
C
         if (i .gt. 1) then         
            lastr = kal(i) + nz(i)
         else
            lastr = ka(i) + nz(i)
         end if
C
C...     To compute stiffness matrix (Laplacian) A at each level.
C
         if (i .gt. 1) then
            lm = 1             !1 - stiffness matrix,  2 -  mass matrix.
            call lapl_mass(lm,nel(i),nd(i),jdf,m(ie),m(je),r(kx),r(ky),
     >           m(ia(i)),m(ja(i)),r(kal(i)),nz(i),m(idir(i)),m(iip))
         end if
C
cc         write(*,*) '============== Discretization 11111111 ==========='
cc         call outmat1(m(ia(i)),m(ja(i)),r(kal(i)),nd(i),nz(i))
cc         call lprr(m(ia(i)),m(ja(i)),r(kal(i)),nd(i))
cc         call wradj_r(m(ia(i)),m(ja(i)),r(kal(i)),
cc     >        min0(nd(i),65*65),202)
C
cc         write(*,*) '============== Discretization 22222222 ==========='
cc         call outmat1(m(ia(i)),m(ja(i)),r(ka(i)),nd(i),nz(i))
cc         stop 111
 50   continue 
C
cc      call outmat1(m(ia(lvl)),m(ja(lvl)),r(ka(lvl)),nd(lvl),nz(lvl))
cc      call outmat1(m(ia(1)),m(ja(1)),r(ka(1)),nd(1),nz(1))
cc      stop 111
C
      ifree = lasti
      kfree = lastr
C
C********************** Multilevel Starts Here ******************************
C...  We consider zero boundary condition only for temporaray.              *
C****************************************************************************
C
C...  Initialize SOL to be zero.
      do 55 ii = 1, lvl
         call nullv(r(ksol(ii)),nd(ii))
 55   continue
C
      k = 0
      iter_tg = 1
 60   continue
      k = k + 1
      lasti = ifree
      lastr = kfree
C
C****************************************************************************
C...  STEP 1. Exact solver of the non-SPD system ^A*u_c = b at the coarsest *
C...          level. It is done by Gaussian elimination.                    *
C****************************************************************************
C...  Compute the RHS for each level first.
C
      mk = mcont(k-1)
C...  To compute the residual: b = rhs - ^A*sol.
      call abyvam(r(krhs(mk)),m(ia(mk)),m(ja(mk)),r(ka(mk)),
     >     r(ksol(mk)),nd(mk),r(kb(mk)))
      
      write(*,*) 'bbbb4'

cc      write(*,*) 'bbbb', (r(kb(mk)+i-1), i=1,nd(mk))

C
      ipp = lasti
      jpp = ipp + nd(lvl)  + 1
      idirtwo = jpp + 2*nd(lvl)
C     
      ktemp = lastr
C     
      do 80 ii = mk, 2, -1
         nf = nd(ii)
         n1f = n1(ii)
         n2f = n2(ii)
         lftwo = multl(ii)
         lctwo = multl(ii-1)
C     
         if (ii .lt. mk) then
            call copyv(r(krhs(ii)),r(kb(ii)),nf)
         end if
         call icopyv(m(idir(ii)),m(idirtwo),nf)
C     
         do 70 i = lftwo, lctwo+1, -1
            call prolng_symb(n1f,n2f,nf,n1c,n2c,nc,m(ipp),m(jpp),nzp)
            call rvbyp(m(ipp),m(jpp),r(kb(ii)),nf,r(ktemp),nc)
            call copyv(r(ktemp),r(kb(ii)),nc)
            call dir_two(r(kb(ii)),m(ipp),m(jpp),nf,nc,m(idirtwo))
C     
            n1f = n1c
            n2f = n2c
            nf = nc
 70      continue
C     
         call copyv(r(kb(ii)),r(krhs(ii-1)),nc)
 80   continue
C
C...  Solve the non-SPD system ^A*e_1 = b by Gaussian Elimination.
C
      call copyv(r(krhs(1)),r(kb(1)),nd(1))

      write(*,*) 'bbbb5'

      call sgauss(m(ia(1)),m(ja(1)),r(ka(1)),r(kb(1)),
     >     m(lasti),r(lastr),nd(1))

      write(*,*) 'bbbb6'
cc      write(*,*) 'bbbb', (r(kb(1)+i-1), i=1,nd(1))
C
C...  Interpolation (prolongation) to the next fine grid: sol = P*u_c.
C
      lctwo = multl(1)
      lftwo = multl(2)
      n1c   = n1(1)
      n2c   = n2(1)
      nc    = nd(1)
C
      call copyv(r(kb(1)),r(kb(2)),nc)
C
      do 90 i = lctwo+1,lftwo
         n1f = n1c*2 - 1
         n2f = n2c*2 - 1
         nf  = n1f*n2f
C     
         call prolng_symb(n1f,n2f,nf,n1c,n2c,nc,m(ipp),m(jpp),nzp)
         call pbyvg(m(ipp),m(jpp),r(kb(2)),nc,r(ktemp),nf)
         call copyv(r(ktemp),r(kb(2)),nf)
C     
         n1c = n1f
         n2c = n2f
         nc  = nf
 90   continue
C     
      if (lvl .eq. 2) then
C...     Update the solution: sol = sol + b.
         call uupluv(r(ksol(2)),r(kb(2)),nd(2))    
      else
         call copyv(r(kb(2)),r(ksol(2)),nd(2))
      end if
C
cc      write(*,*) 'bbbb', (r(ksol(2)+i-1), i=1,nd(2))
C
C****************************************************************************
C...  STEP 2. Exact solver of the SPD system A*e_2 = b at the second        *
C...          coarsest level. It is done by an MG V-cycle.                  *
C****************************************************************************
C
 95   continue
C
C...  To compute the residual: b = rhs - ^A*sol.
      call abyvam(r(krhs(2)),m(ia(2)),m(ja(2)),r(ka(2)),
     >     r(ksol(2)),nd(2),r(kb(2)))
C
C...  Solve of the SPD system A*e_2 = b by an MG V-cycle.
C
cc      tol  = 1.d-6
      lfmg = multl(2)
      lcmg = 1
C
      call mg_s(m(1),ia(2),ja(2),idir(2),lasti,
     >     r(1),ku(2),kb(2),kal(2),lastr,
     >     nd(2),n1(2),n2(2),nz(2),nel(2),lfmg,lcmg,tol) 
C
C...  Update the solution: sol = sol + u.
      call uupluv(r(ksol(2)),r(ku(2)),nd(2))      
C
C.... Compute the H^1 and L^2 norms of the error.
C
      if (lvl .eq. 2) then
         ksol_exact = lastr
         kerr       = ksol_exact + nd(lvl)
C
         ich_sol = 1
C
         if (ich_sol .eq. 1) then
            call read_soln(multl(lvl),r(ksol_exact),nd(lvl),lread)
         end if
C
         call wuminv(r(kerr),r(ksol_exact),r(ksol(lvl)),nd(lvl))
         call l2_h1s_norm(0,r(kerr),nd(lvl),sl2norm,nel(lvl),jdf,
     >        m(ief),m(jef),r(kxf),r(kyf))
         call l2_h1s_norm(1,r(kerr),nd(lvl),h1snorm,nel(lvl),jdf,
     >        m(ief),m(jef),r(kxf),r(kyf))
         h1norm = dsqrt(sl2norm*sl2norm + h1snorm*h1snorm)
C
         iter_tg = iter_tg + 1
C     
         write(*,'(a,i6,2(e14.4))') 
     >        ' No.iter     L_2      H_1: ',iter_tg,sl2norm,h1norm
C
cc         write(*,'(i6,2(e14.4))') iter_tg, sl2norm, h1norm
C
         if (iter_tg .eq. iter_tg_max) go to 300
      end if
C
C****************************************************************************
C...  STEP 3. Exact solver of the non-SPD system ^A*e_1 = b at the coarsest *
C...          level. It is done by Gaussian elimination.                    *
C****************************************************************************
C
C...  Compute b = ^A*sol.
      call abyvg(m(ia(2)),m(ja(2)),r(ka(2)),r(ksol(2)),nd(2),r(kb(2)))
C
C...  Restriction to the coarsest grid: b = P^t*b.
C
      nf = nd(2)
      n1f = n1(2)
      n2f = n2(2)
      lftwo = multl(2)
      lctwo = multl(1)
C     
      call icopyv(m(idir(2)),m(idirtwo),nf)
C     
      do 100 i = lftwo, lctwo+1, -1
         call prolng_symb(n1f,n2f,nf,n1c,n2c,nc,m(ipp),m(jpp),nzp)
         call rvbyp(m(ipp),m(jpp),r(kb(2)),nf,r(ktemp),nc)
         call copyv(r(ktemp),r(kb(2)),nc)
         call dir_two(r(kb(2)),m(ipp),m(jpp),nf,nc,m(idirtwo))
C     
         n1f = n1c
         n2f = n2c
         nf = nc
 100  continue
C     
      call copyv(r(kb(2)),r(kb(1)),nc)
C
C...  Compute the residual b = rhs - ^A*sol (= rhs - b).
      call vuminv(r(krhs(1)),r(kb(1)),nd(1))
C
C...  Solve the non-SPD system ^A*e_1 = b by Gaussian Elimination.
      call sgauss(m(ia(1)),m(ja(1)),r(ka(1)),r(kb(1)),
     >     m(lasti),r(lastr),nd(1))
C
C...  Interpolation (prolongation) to the next fine grid: b = P*u_c.
C
      lctwo = multl(1)
      lftwo = multl(2)
      n1c   = n1(1)
      n2c   = n2(1)
      nc    = nd(1)
C     
      call copyv(r(kb(1)),r(kb(2)),nc)
C     
      do 110 i = lctwo+1,lftwo
         n1f = n1c*2 - 1
         n2f = n2c*2 - 1
         nf  = n1f*n2f
C     
         call prolng_symb(n1f,n2f,nf,n1c,n2c,nc,m(ipp),m(jpp),nzp)
         call pbyvg(m(ipp),m(jpp),r(kb(2)),nc,r(ktemp),nf)
         call copyv(r(ktemp),r(kb(2)),nf)
C     
         n1c = n1f
         n2c = n2f
         nc  = nf
 110   continue
C
C...  Update the solution: sol = sol + b.
      call uupluv(r(ksol(2)),r(kb(2)),nd(2))
      
cc      write(*,*) 'bbbb', (r(ksol(2)+i-1), i=1,nd(2))
C
C****************************************************************************
C...  STEP 4. Finish the algorithm after updating the solution at the       *
C...          finest grid. Otherwise, solve the SPD system A*e_mk = b at   * 
C...          finer grid and go back to the Step 1.                         *
C****************************************************************************
C
      if (lvl .eq. 2) go to 200
C
      mk = mcont(k)
      if (mk .eq. 0) then
         mkt = lvl
      else
         mkt = mk - 1
      end if
C
      do 130 ii = 3, mkt
C
C...     Interpolation (prolongation) to the next fine grid: sol = P*sol.
C
         lctwo = multl(ii-1) 
         lftwo = multl(ii)
         n1c   = n1(ii-1)
         n2c   = n2(ii-1)
         nc    = nd(ii-1)
C     
         call copyv(r(ksol(ii-1)),r(kb(ii)),nc)
C     
         do 120 i = lctwo+1,lftwo
            n1f = n1c*2 - 1
            n2f = n2c*2 - 1
            nf  = n1f*n2f
C     
            call prolng_symb(n1f,n2f,nf,n1c,n2c,nc,m(ipp),m(jpp),nzp)
            call pbyvg(m(ipp),m(jpp),r(kb(ii)),nc,r(ktemp),nf)
            call copyv(r(ktemp),r(kb(ii)),nf)
C     
            n1c = n1f
            n2c = n2f
            nc  = nf
 120     continue
C     
C...     Update the solution: sol = sol + b.
         call uupluv(r(ksol(ii)),r(kb(ii)),nd(ii))
 130  continue
C
      if (mk .eq. 0) go to 200
C
C...  Interpolation (prolongation) to next finer level, MK : sol = P*sol.
C
      lctwo = multl(mk-1) 
      lftwo = multl(mk)
      n1c   = n1(mk-1)
      n2c   = n2(mk-1)
      nc    = nd(mk-1)
C     
      call copyv(r(ksol(mk-1)),r(kb(mk)),nc)
C     
      do 140 i = lctwo+1,lftwo
         n1f = n1c*2 - 1
         n2f = n2c*2 - 1
         nf  = n1f*n2f
C     
         call prolng_symb(n1f,n2f,nf,n1c,n2c,nc,m(ipp),m(jpp),nzp)
         call pbyvg(m(ipp),m(jpp),r(kb(mk)),nc,r(ktemp),nf)
         call copyv(r(ktemp),r(kb(mk)),nf)
C     
         n1c = n1f
         n2c = n2f
         nc  = nf
 140  continue
C     
      call copyv(r(kb(mk)),r(ksol(mk)),nd(mk))
C
C...  To compute the residual: b = rhs - ^A*sol.
      call abyvam(r(krhs(mk)),m(ia(mk)),m(ja(mk)),r(ka(mk)),
     >     r(ksol(mk)),nd(mk),r(kb(mk)))
C
C...  Solve of the SPD system A*e_mk = b by an MG V-cycle.
C
cc      tol  = 1.d-6
      lfmg = multl(mk)
      lcmg = 1
C
      call mg_s(m(1),ia(mk),ja(mk),idir(mk),lasti,
     >     r(1),ku(mk),kb(mk),kal(mk),lastr,
     >     nd(mk),n1(mk),n2(mk),nz(mk),nel(mk),lfmg,lcmg,tol) 
C
C...  Update the solution: sol = sol + u.
      call uupluv(r(ksol(mk)),r(ku(mk)),nd(mk))      
C
 200  continue
C
C.... Compute the H^1 and L^2 norms of the error.
C
      if (mk .eq. 0 .or. mk .eq. lvl .or. lvl .eq. 2) then
         ksol_exact = lastr
         kerr       = ksol_exact + nd(lvl)
C
         ich_sol = 1
C
         if (ich_sol .eq. 1) then
            call read_soln(multl(lvl),r(ksol_exact),nd(lvl),lread)
         end if
C
         call wuminv(r(kerr),r(ksol_exact),r(ksol(lvl)),nd(lvl))
         call l2_h1s_norm(0,r(kerr),nd(lvl),sl2norm,nel(lvl),jdf,
     >        m(ief),m(jef),r(kxf),r(kyf))
         call l2_h1s_norm(1,r(kerr),nd(lvl),h1snorm,nel(lvl),jdf,
     >        m(ief),m(jef),r(kxf),r(kyf))
         h1norm = dsqrt(sl2norm*sl2norm + h1snorm*h1snorm)
C
         if (lvl .eq. 2) k = iter_tg + 1
C
         write(*,'(a,i6,2(e14.4))') 
     >        ' No.iter     L_2      H_1: ',k,sl2norm,h1norm
C
cc         write(*,'(i6,2(e14.4))') k, sl2norm, h1norm
      end if
C
      if (lvl .eq. 2) then
         iter_tg = iter_tg + 1
         if (iter_tg .eq. iter_tg_max) then
            go to 300
         else
            go to 95
         end if
      else
         if (mk .eq. 0) then
            go to 300
         else
            go to 60
         end if
      end if
C
 300  continue
C
C****************************************************************************
C...  End of the algorithm.                                                 *
C****************************************************************************
C
C****************************************************************************
C...  To print out the computed solution SOL.                               *
C****************************************************************************
      iexact = 1
      if (iexact .eq. 1) then
         call write_soln(multl(lvl),r(ksol(lvl)),nd(lvl),lread)
      end if
C
      if (nd(lvl) .le. 65*65) then
         call wrsol01(r(ksol(lvl)),nd(lvl),100)
         endfile(100)
         close(100)
      end if
C
      smin = r(ksol(lvl))
      smax = smin
C
      do i = 2 , nd(lvl)
cc         if( i.le. 65*65) write(100,'(e15.7$)') (r(ksol(lvl)+i-1))
         smin = dmin1(smin,r(ksol(lvl)+i-1))
         smax = dmax1(smax,r(ksol(lvl)+i-1))
cc         if(dabs(r(ksol(lvl)+i-1)) .gt. 1d-05) then
cc             write(*,*) i,r(ksol(lvl)+i-1)
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
C
      write(*,'(10x,a)') ' CPU time:  '
      write(*,'(10x,a)') ' ======================== '
      write(*,'(10x,a,f10.3,a)') '   SETUP: ',t_setup, ' s '
      write(*,'(10x,a)') ' ------------------------'
      write(*,'(10x,a,f10.3,a)') '   SOLVE: ',t_iter, ' s '
      write(*,'(10x,a)') ' ------------------------'
      write(*,'(10x,a,f10.3,a)') '   TOTAL: ',t_setup+t_iter, ' s '
C
      stop
      end
C=====================================================================

