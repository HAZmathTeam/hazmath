
      program nstgb

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     NSTGB uses bilinear finite elements on uniform square mesh to solve
c     the nonsymmetric/indefinite elliptic Dirichlet boundary value problem by 
c     Two-Grid method:
c         {  -Laplacian(U) + cU_x + dU_y + fU = fxy   in the domain      }
c         {                                 U = hxy   on the boundary    },
c     where c, d, and f are constants.
c
c     First, the user should describe the functions f(x,y), h(x,y) in the
c     function subprograms 'fxy', 'hxy'.
c
c     During execution of the program the user specifies the level, lvl
c     (lvl = 4, 6, 8, or 10) of the  so that h = 1/(2^lvl), where
c     h is the length of each side of squares (parallel to x,y-axes) in the 
c     uniform mesh of the unit square.
c
c     leng = number of nodes in the finest structured grid
c              (if lvl = 8,  set leng =    66050)
c              (if lvl = 9,  set leng =   263170)
c              (if lvl = 10, set leng =  1050630)
c              (if lvl = 11, set leng =  4198410)
c     lvl = level of the structured space
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     icdomn : choice of domains
c              1 = the quarter circle.
c              3 = the unit square with zeros only on the boundary.
c              4 = the unit square with zeros near/on the boundary.
c     icalgm : choice of algorithms
c              1 = Algorithm 1
c              2 = Algorithm 2
c              3 = Algorithm 3
c        A_0 : stiffness matrix of the uniform mesh in the structured space
c              (bilinear elements on squares)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (leng = 1050630)
      implicit double precision (a-h, o-z)
      dimension b(leng), x(leng), exsoln(leng), wk(leng),
     &          bwk(2*leng), xwk(2*leng), wk1(2*leng), wk2(2*leng),
     &          itbl(15,2,20), dtbl(15,6,20), jd(14)
      real      time(2), etime, dtime
      logical   done1, done2

      open(unit = 20, file = 'nsck2', status = 'unknown')

 10   continue

      t11 = etime(time)
      t12 = dtime(time)

      c      =      -5000.d0
      d      =      -3000.d0
      f      =     0.d0

      icalgm = 2
      icsolv = 2              !Choice of an exact solver for the coarse grids.

cc      print*, ' Enter the choice of the domain, 1, 3, or 4:'
cc      read*, icdomn
      icdomn = 3

      print 20
 20   format (/' Enter the level (3 < lvl < 11) of the fine' 
     &        /' grid so that h =(approx.) 1/(2^lvl):')
      read*, lvl

      if (lvl .lt. 1 .or. lvl .gt. 11) then
         print*, ' Invalid choice. Reenter!!'
	 go to 10
      end if

      k  = 2**lvl 
      k1 = k - 1
      k2 = k + 1
      n  = k2*k2

      jd(1) = 1
      do 30 i = 1, 13
         jd(i+1) = jd(i) + (2**i+1)**2
 30   continue

      print*, ' Please wait a while if the level is large.'
      print*

      icb = 2                                       !Choice of load vector b.
      if (icb .eq. 1) then
c        To compute the load vector bwk of the linear system
c        (A_0 + cA_0x + dA_0y + fM_0)*u = bwk in the structured mesh directly.
         call loadbs(lvl, b, icdomn, c, d, f)
      else if (icb .eq. 2) then
c        To supply an exact FE solution and a correspoding load vector in order
c        to test the algorithms. (This is a temporary process.)
         call exctxb(lvl, n, exsoln, b, icdomn, c, d, f)
      end if 

c*************************New Approach******************************************

      ipb = 6                                   !ipb = 0, 1, 2, 3, 4, 5, 6
      if (ipb .eq. 0) go to 40

c     ipcnd=1 : bpx precond. 
c     ipcnd=2 : symmetric mg \-cycle as a precond. 
c     ipcnd=3 : symmetric mg V-cycle as a precond.
c     ipcnd=4 : multiplicative precond with mg \-cycle and a coarse grid solver.
c     ipcnd=5 : multiplicative precond with mg V-cycle and a coarse grid solver.
c     ipcnd=6 : nonsymmetric/indefinite mg \-cycle as a precond. 
c     ipcnd=7 : nonsymmetric/indefinite mg V-cycle as a precond.
      ipcnd = 7                !ipcnd = 1, 2, 3, 4, 5, 6, 7

	print*, 'Enter jcbegin'
	read*, jcb
      if (lvl .eq. 2) then
c	 jcb = 2
	 jce = 1      
      else if (lvl .eq. 3) then
c	 jcb = 2
         jce = 2
      else if (lvl .eq. 4) then
c	 jcb = 2
         jce = 3
      else if (lvl .eq. 5) then
c	 jcb = 2
         jce = 3
      else if (lvl .eq. 6) then
c	 jcb = 2
         jce = 3
      else if (lvl .eq. 7) then
c	 jcb = 2
         jce = 3
      else if (lvl .eq. 8) then
c	 jcb = 2
         jce = 3
      else if (lvl .eq. 9) then
c	 jcb = 2
         jce = 3
      else if (lvl .eq. 10) then
c	 jcb = 2
         jce = 3
      end if

      do 31 jc = jcb, jce                   !Choice of the coarse level
         nc  = jd(jc+1) - jd(jc)             !Dimension of the coarse level.
	 if (ipb .le. 2)    jc = lvl
	 if (ipb .ge. 3 .and. ipb .le. 5 .and. ipcnd .le. 3) jc = lvl

         t51 = etime(time)
         t52 = dtime(time)

         if (ipb .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Bi-CGSTAB method in the structured mesh.
            call bcgst(lvl, n, b, x, 1, nit, icdomn, 3, c, d, f)
         else if (ipb .eq. 2) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Bi-CGSTAB method in the structured mesh.
            call bcgsto(lvl, n, b, x, 1, nit, icdomn, 3, c, d, f)
         else if (ipb .eq. 3) then
	    if (lvl .le. 3 .and. ipcnd .ge. 4) go to 170
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Preconditioned Bi-CGSTAB method in the structured mesh.
            call bcgstp(lvl, n, b, x, 1, nit, icdomn, 3, c, d, f,
     &                  jd, bwk, xwk, wk, wk1, wk2, ipcnd, jc, nc)
         else if (ipb .eq. 4) then
	    if (lvl .le. 3 .and. ipcnd .ge. 4) go to 170
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Preconditioned Bi-CGSTAB method in the structured mesh.
            call bcgstq(lvl, n, b, x, 1, nit, icdomn, 3, c, d, f,
     &                  jd, bwk, xwk, wk, wk1, wk2, ipcnd, jc, nc)
         else if (ipb .eq. 5) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Preconditioned Bi-CGSTAB method in the structured mesh.
            call bcgsts(lvl, n, b, x, 1, nit, icdomn, 3, c, d, f,
     &                  jd, bwk, xwk, wk, wk1, wk2, ipcnd, jc, nc)
         else if (ipb .eq. 6) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by multigrid V- or \- cycle method in the structured mesh.
	    ipde   = 3
	    nit    = 100000
	    ivbs   = 1
cccc	    jc     = 1          !For an exact solver at the level 1.
	    nsmthg = 2
	    ifbs1  = 3
	    ifbs2  = 3
	    isolv  = 2
            call mgvbs(lvl, n, icdomn, ipde, nit, c, d, f, b, x, isolv,
     &               ivbs, jc, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
cccc	    jc = lvl            !For an exact solver at the level 1.
         end if

         t61 = etime(time)
         t62 = dtime(time)

c        Compute maximum, L^2, H^1 seminorm, and H^1 norm errors.
         call error0(lvl, n, jc, x, nit, t62, icdomn, c, d, f, exsoln, 
     &                  wk, xwk, 2, itbl, dtbl)
31    continue

      if (ipb .ne. 0) go to 170
c*******************************************************************************

40    do 50 kk = 2, lvl-1                    !Choice of the coarse level
	 if (lvl .eq. 10 .and. kk .eq. 9) go to 170
         nc  = jd(kk+1) - jd(kk)             !Dimension of the coarse level.

         t21 = etime(time)
         t22 = dtime(time)

         if (icalgm .eq. 1 .or. icalgm .eq. 2) then
            nit = 6
   	    call algm2(lvl, n, kk, nc, icdomn, icsolv, nit, c, d, f, 
     &           b, x, jd, bwk, xwk, wk, wk1, wk2, exsoln, itbl, dtbl)
         else if (icalgm .eq. 3) then
            nit = 1
	    call algm3(lvl, n, kk, nc, icdomn, icsolv, nit, c, d, f, 
     &                            b, x, jd, bwk, xwk, wk, wk1, wk2)
         end if

         t31 = etime(time)
         t32 = dtime(time)

	 if (icalgm .eq. 2) go to 50

c        Compute maximum, L^2, H^1 seminorm, and H^1 norm errors.
         call error0(lvl, n, kk, x, nit, t32, icdomn, c, d, f, exsoln, 
     &                  wk, xwk, 2, itbl, dtbl)
 50   continue

c      print*, 'CPU time = ', t11, t21, t31
c      print*, 'CPU time = ', t12, t22, t32
  
c     Presentation of result (computed solution).

 170  Print 180
 180  format(/' 2 Display the solution vector X with its coordinates.'
     &       /'   (various 3D views using MATLAB tool)'
     &       /' 3 Compute max.norm, L^2 norm, and H^1 (semi)norm errors'
     &	     /'   between computed solution(u_h) and the exact finite'
     &       /'   element solution(u_e) of the discretized problem.'
     &       /' 4 Compute max norm, L^2 norm, and H^1 (semi)norm errors'
     &	     /'   between computed solution(u_h) and interpolation(u_I)'
     &       /'   of continuous exact solution at grid points, if exact'
     &       /'   solution is known. (Look at subroutine EXTSOL!!)'
     &       /' 0 Quit.'
     &      //' Enter choice :')

      read*, iquery
      if (iquery .eq. 1) then
         continue
      else if (iquery .eq. 2) then
c        print*,' Solution at each node : ', (i,'*',x(i),i=1, nonode)
c        To display various 3D views of the computed solution using MATLAB tool.
         call grflab(lvl, n, x, icdomn)
      else if (iquery .eq. 3) then
         continue
      else if (iquery .eq. 4) then
         continue
      else if (iquery .eq. 0) then
  200    continue
	 print*, ' Do you want to run the program again with '
         print*, ' the different choice of levels lvlaux & lvluns?'
	 print*, ' If yes, enter 1:'
	 print*, ' If no, enter 99:'
         read*, kchoic
         if (kchoic .eq. 1) then
            go to 10
         else if (kchoic .eq. 99) then
            close(unit = 20)
            print*, ' GOOD BYE!!'
            stop
         else
            print*,' Invalid choice. Reenter!!'
            go to 200
	 end if
      else
         print*,' Invalid choice. Reenter!!'
         go to 170
      end if

      go to 170
      end
c=======================================================================

      subroutine loadbs(j, b, ichoic, c, d, f)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes the vector b = (b(j)) in the structured mesh,
c     where b(j) = integral_omega(f*phi_j).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension b((2**j + 1)**2)
      external psi, fxy

      m  = 2**j + 1
      m1 = m - 1
      h  = 1.d0/m1

      do 10 i = 1, m*m
	 b(i) = 0.d0
   10 continue

c*    Integration in Omega

      do 40 l = 1, m1
	 it = (l-1)*m
	 do 30 k = 1, m1
	    i  = it + k
            xa = (k-1)*h
            ya = (l-1)*h
	    xb = xa + h
	    yb = ya + h

c           Call INTGRL to compute integral on a rectangle.
            call intgrl(2, xa, xb, ya, yb,4,psi,fxy,result,ichoic,c,d,f)
            b(i)     = b(i)     + result

            call intgrl(2, xa, xb, ya, yb,2,psi,fxy,result,ichoic,c,d,f)
            b(i+1)   = b(i+1)   + result

            call intgrl(2, xa, xb, ya, yb,3,psi,fxy,result,ichoic,c,d,f)
            b(i+m)   = b(i+m)   + result

            call intgrl(2, xa, xb, ya, yb,1,psi,fxy,result,ichoic,c,d,f)
            b(i+m+1) = b(i+m+1) + result 
   30    continue
   40 continue

      call bpxin(j, b, ichoic)

      return
      end
c=======================================================================

      subroutine intgrl(n, xa, xb, ya,yb,jp,hxy,fxy,result,ichoic,c,d,f)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine calculates double integration of hxy*fxy on the rectangle
c     with vertices (xa,ya), (xa,yb), (xb,yb), (xb, ya), using Gauss qudrature.
c
c     n  -  the number of Gauss points in each side
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension w(4),z(4),t(4),s(4)

      xl = xb - xa
      yl = yb - ya

      if (n .ne. 1) go to 2
	 z(1) = 0.d0
	 w(1) = 2.d0
	 go to 15

    2 if (n .ne. 2) go to 3
	 z(2) = 1/dsqrt(3.d0)
	 w(2) = 1.d0
	 go to 5

    3 if (n .ne. 3) go to 4
	 z(2) = 0.d0
	 z(3) = 0.774596669241483d0
	 w(2) = 0.88888888888888888889d0
	 w(3) = 0.55555555555555555556d0
	 go to 5

    4 zt   = 4.d0*dsqrt(0.3d0)
      z(3) = dsqrt((3.d0-zt)/7.d0)
      z(4) = dsqrt((3.d0+zt)/7.d0)
      wt   = dsqrt(10.d0/3.d0)/12.d0
      w(3) = 0.5d0 + wt
      w(4) = 0.5d0 - wt

    5 do 10 j = 1,int(n/2)
         w(j) = w(n+1-j)
         z(j) = -z(n+1-j)
   10 continue

   15 do 20 j = 1,n
	 t(j) = (z(j)*xl+xa+xb)/2.d0
         s(j) = (z(j)*yl+ya+yb)/2.d0
   20 continue

      result = 0.d0
      do 40 j = 1,n
	 x    = t(j)
         temp = 0.d0
         do 30 i = 1,n
	    y = s(i)
c           Call functions FXY and HXY.
            gxy  = fxy(x,y,ichoic,c,d,f)*hxy(xa, xl, ya, yl, x, y, jp)
            temp = temp + gxy*w(i)
   30    continue
         temp   = temp*yl/2.d0
         result = result + temp*w(j)
   40 continue
      result = result*xl/2.d0

      return
      end
c=======================================================================

      subroutine algm1(j, n, jc, nc, icdomn, nit, c, d, f, b, x,isolv,
     &           ivbs, nsmthg, ifbs1, ifbs2, jd, bwk, xwk, wk, wk1, wk2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ALGM1 solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a Two-Grid method.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  x(n), b(n), wk(n), xwk(2*n), bwk(2*n)
      dimension  wk1(2*n), wk2(2*n), jd(14)

      noit = 1                       !No. of iterations of Algorithm 1.
c     isolv = 2                      !1 or 2

      do 20 i = 1, jd(j)-1
         xwk(i)  = 0.d0
         bwk(i)  = 0.d0
 20   continue

      do 30 i = 1, n
         x(i)           = 0.d0
	 xwk(jd(j)+i-1) = 0.d0
	 bwk(jd(j)+i-1) = b(i)
 30   continue

      do 200 k = 1, noit
c********Step 1. Solve the nonsymmetric problem exactly at coarser space.*******
         if (k .ne. 1) then
c           Compute wk = (A + cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 3, c, d, f)

c           Compute bwk = b - wk = b - (A + cA_x + dA_y + fM)*x.
            do 40 i = 1, n
               bwk(jd(j)+i-1) = b(i) - wk(i)
 40         continue
         end if

         do 60 i = j, jc+1, -1
            call restrn(i, bwk(jd(i)), bwk(jd(i-1))) !Restric. to lower level.
            call bpxin(i-1, bwk(jd(i-1)), icdomn)
 60      continue

	 print*

	 if (isolv .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Bi-CGSTAB method in the structured mesh.
            call bcgsto(jc, nc, bwk(jd(jc)), xwk(jd(jc)), 0, nit0,
     &                  icdomn, 3, c, d, f) 
	 else if (isolv .eq. 2) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by GSUPWD method in the structured mesh.
	    call gsupwd(jc, nc, bwk(jd(jc)), xwk(jd(jc)), icdomn,
     &                  3, 1, c, d, f, wk, 1500) 
	 end if

c********Step 2. Correction step. i.e., new half-step solution.*****************
         do 80 i = jc+1, j
            call prolng(i,xwk(jd(i-1)),xwk(jd(i))) !Prolongation to upper level.
            call bpxin(i, xwk(jd(i)), icdomn)
 80      continue

         do 100 i = 1, n
            x(i) = x(i) + xwk(jd(j)+i-1)
 100     continue

c********Step 3. Solve the symmetric problem at finer space.********************

         icit = 1                               !icit = 1 is better somehow.
         if (icit .eq. 1) then
c           To compute wk = (A_h + cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 3, c, d, f)
         else if (icit .eq. 2) then
c           To compute wk = (cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 5, c, d, f)
         end if

c        Compute bwk = b - wk = b - (A_h + cA_x + dA_y + fM)*x.
         do 120 i = 1, n
            bwk(jd(j)+i-1) = b(i) - wk(i)
 120     continue

         icmg = 2
c        To solve the SPD linear system A_h*xwkt = bwk.
         if (icmg .eq. 1) then
            nit1 = 1
            call mgvbs(j, n, icdomn, 1, nit1, 0.d0,0.d0,0.d0,bwk(jd(j)), 
     &      xwk(jd(j)),0, ivbs, 1, nsmthg, ifbs1, ifbs2, jd,wk,wk1,wk2)
         else if (icmg .eq. 2) then
            call premg(j, n, icdomn, 1, 0.d0, 0.d0, 0.d0, bwk(jd(j)), 
     &      xwk(jd(j)),0, ivbs,1,nsmthg,ifbs1,ifbs2,jd,wk,wk1,wk2)
         end if

c********Step 4. Correction step. i.e., new solution****************************
         if (icit .eq. 1) then
            do 140 i = 1, n
               x(i) = x(i) + xwk(jd(j)+i-1)
 140        continue
         else if (icit .eq. 2) then
            do 160 i = 1, n
               x(i) = xwk(jd(j)+i-1)
 160        continue
         end if

c********End of Step 4.*********************************************************
 200  continue

 300  return
      end
c=======================================================================

      subroutine algm2(j, n, jc, nc, icdomn, icsolv, nit, c, d, f, 
     &           b, x, jd, bwk, xwk, wk, wk1, wk2, exsoln, itbl, dtbl)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ALGM1 solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a Two-Grid method.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  x(n), b(n), wk(n), xwk(2*n), bwk(2*n), exsoln(n)
      dimension  wk1(2*n), wk2(2*n), jd(14), itbl(15,2,20),dtbl(15,6,20)

      t21 = etime(time)
      t22 = dtime(time)
      t32 = 0.d0

      nsmthg = 1
      ifbs1 = 1
      ifbs2 = 2

      ihalf = 0
      noit = nit                        !No. of iterations of Algorithm 2.

      do 20 i = 1, jd(j)-1
         xwk(i)  = 0.d0
         bwk(i)  = 0.d0
 20   continue

      do 30 i = 1, n
         x(i)           = 0.d0
	 xwk(jd(j)+i-1) = 0.d0
	 bwk(jd(j)+i-1) = b(i)
 30   continue

      do 200 k = 1, noit
c********Step 1. Solve the nonsymmetric problem exactly at coarser space.*******
         if (k .ne. 1) then
c           Compute wk = (A + cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 3, c, d, f)

c           Compute bwk = b - wk = b - (A + cA_x + dA_y + fM)*x.
            do 40 i = 1, n
               bwk(jd(j)+i-1) = b(i) - wk(i)
 40         continue
         end if

c	 if (k .eq. 3 .or. k .eq. 5) go to 130

         do 60 i = j, jc+1, -1
            call restrn(i, bwk(jd(i)), bwk(jd(i-1))) !Restric. to lower level.
            call bpxin(i-1, bwk(jd(i-1)), icdomn)
 60      continue

         if (icsolv .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c           by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
            call gssor(jc, nc, bwk(jd(jc)), xwk(jd(jc)),
     &                    icdomn, 3, 3, c, d, f, wk, nit0)
         else if (icsolv .eq. 2) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Bi-CGSTAB method in the structured mesh.
            call bcgsto(jc, nc, bwk(jd(jc)), xwk(jd(jc)), 0, nit0,
     &                  icdomn, 3, c, d, f) 
         end if

c********Step 2. Correction step. i.e., new half-step solution.*****************
         do 80 i = jc+1, j
            call prolng(i,xwk(jd(i-1)),xwk(jd(i))) !Prolongation to upper level.
            call bpxin(i, xwk(jd(i)), icdomn)
 80      continue

         do 100 i = 1, n
            x(i) = x(i) + xwk(jd(j)+i-1)
 100     continue

         t31 = etime(time)
         t32t = dtime(time)
	 t32  = t32 + t32t

c        Compute maximum, L^2, H^1 seminorm, and H^1 norm errors.
         call error0(j, n, jc, x, 10+k, t32, icdomn, c, d, f, exsoln, 
     &                  wk, xwk, 2, itbl, dtbl)

	 if (k .eq. noit .and. ihalf .eq. 1) go to 300

c********Step 3. Solve the symmetric problem at finer space.********************

         icit = 1                               !icit = 1 is better somehow.
         if (icit .eq. 1) then
c           To compute wk = (A_h + cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 3, c, d, f)
         else if (icit .eq. 2) then
c           To compute wk = (cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 5, c, d, f)
         end if

c        Compute bwk = b - wk = b - (A_h + cA_x + dA_y + fM)*x.
         do 120 i = 1, n
            bwk(jd(j)+i-1) = b(i) - wk(i)
 120     continue

 130     continue

c        To solve the SPD linear system A_h*xwkt = bwk.
         nit1 = 30
         call mgvbs(j, n, icdomn, 1, nit1, 0.d0, 0.d0,0.d0, bwk(jd(j)), 
     &      xwk(jd(j)), 1, 1, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)


c********Step 4. Correction step. i.e., new solution****************************
         if (icit .eq. 1) then
            do 140 i = 1, n
               x(i) = x(i) + xwk(jd(j)+i-1)
 140        continue
         else if (icit .eq. 2) then
            do 160 i = 1, n
               x(i) = xwk(jd(j)+i-1)
 160        continue
         end if

c********End of Step 4.*********************************************************

         t41 = etime(time)
         t42t = dtime(time)
	 t32  = t32 + t42t

c        Compute maximum, L^2, H^1 seminorm, and H^1 norm errors.
         call error0(j, n, jc, x, k, t32, icdomn, c, d, f, exsoln, 
     &                  wk, xwk, 2, itbl, dtbl)

 200  continue

 300  return
      end
c=======================================================================

      subroutine algm3(j, n, jc, nc, icdomn, icsolv, nit, c, d, f, 
     &                             b, x, jd, bwk, xwk, wk, wk1, wk2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ALGM1 solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a Two-Grid method.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  x(n), b(n), wk(n), xwk(2*n), bwk(2*n)
      dimension  wk1(2*n), wk2(2*n), jd(14)

      nsmthg = 1
      ifbs1  = 1
      ifbs2  = 2
      noit = nit                        !No. of iterations of Algorithm 3.

      do 20 i = 1, jd(j)-1
         xwk(i)  = 0.d0
         bwk(i)  = 0.d0
 20   continue

      do 30 i = 1, n
         x(i)           = 0.d0
	 xwk(jd(j)+i-1) = 0.d0
	 bwk(jd(j)+i-1) = b(i)
 30   continue

      do 200 k = 1, noit
c********Step 1. Solve the nonsymmetric problem exactly at coarser space.*******
         if (k .ne. 1) then
c           Compute wk = (A + cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 3, c, d, f)

c           Compute bwk = b - wk = b - (A + cA_x + dA_y + fM)*x.
            do 40 i = 1, n
               bwk(jd(j)+i-1) = b(i) - wk(i)
 40         continue
         end if

         do 60 i = j, jc+1, -1
            call restrn(i, bwk(jd(i)), bwk(jd(i-1))) !Restric. to lower level.
            call bpxin(i-1, bwk(jd(i-1)), icdomn)
 60      continue

         if (icsolv .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c           by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
            call gssor(jc, nc, bwk(jd(jc)), xwk(jd(jc)),
     &                    icdomn, 3, 3, c, d, f, wk, nit0)
         else if (icsolv .eq. 2) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Bi-CGSTAB method in the structured mesh.
            call bcgsto(jc, nc, bwk(jd(jc)), xwk(jd(jc)), 0, nit0,
     &                  icdomn, 3, c, d, f) 
         end if

c********Step 2. Solve the symmetric problem at finer space.********************
         do 80 i = jc+1, j
            call prolng(i,xwk(jd(i-1)),xwk(jd(i))) !Prolongation to upper level.
            call bpxin(i, xwk(jd(i)), icdomn)
 80      continue

         do 100 i = 1, n
            x(i) = x(i) + xwk(jd(j)+i-1)
 100     continue

         icit = 1                               !icit = 1 is better somehow.
         if (icit .eq. 1) then
c           To compute wk = (A_h + cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 3, c, d, f)
         else if (icit .eq. 2) then
c           To compute wk = (cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 5, c, d, f)
         end if

c        Compute bwk = b - wk = b - (A_h + cA_x + dA_y + fM)*x.
         do 120 i = 1, n
            bwk(jd(j)+i-1) = b(i) - wk(i)
 120     continue

c        To solve the SPD linear system A_h*xwkt = bwk.
         nit1 = 30
         call mgvbs(j, n, icdomn, 1, nit1, 0.d0, 0.d0,0.d0, bwk(jd(j)), 
     &      xwk(jd(j)), 1, 1, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)

c********Step 3. Solve the nonsymmetric problem at coarse space.****************

         if (icit .eq. 1) then
            do 130 i = 1, n
               x(i) = x(i) + xwk(jd(j)+i-1)
 130        continue

c           To compute wk = (A_h + cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, wk, icdomn, 3, c, d, f)

c           Compute bwk = b - wk = b - (A_h + cA_x + dA_y + fM)*x.
            do 135 i = 1, n
               bwk(jd(j)+i-1) = b(i) - wk(i)
 135        continue
	 else if (icit .eq. 2) then
            do 140 i = 1, n
               x(i) = x(i) - xwk(jd(j)+i-1)
 140        continue

c           To compute bwk = (cA_x + dA_y + fM)*x at level j.
            call matvc2(j, x, bwk, icdomn, 5, c, d, f)
	 end if

         do 150 i = j, jc+1, -1
            call restrn(i, bwk(jd(i)), bwk(jd(i-1))) !Restric. to lower level.
            call bpxin(i-1, bwk(jd(i-1)), icdomn)
 150      continue

         if (icsolv .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c           by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
            call gssor(jc, nc, bwk(jd(jc)), xwk(jd(jc)),
     &                    icdomn, 3, 3, c, d, f, wk, nit0)
         else if (icsolv .eq. 2) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Bi-CGSTAB method in the structured mesh.
            call bcgsto(jc, nc, bwk(jd(jc)), xwk(jd(jc)), 0, nit0,
     &                  icdomn, 3, c, d, f) 
         end if

c********Step 4. Correction step. i.e., new solution****************************

         do 160 i = jc+1, j
            call prolng(i,xwk(jd(i-1)),xwk(jd(i))) !Prolongation to upper level.
            call bpxin(i, xwk(jd(i)), icdomn)
 160      continue

         do 170 i = 1, n
            x(i) = x(i) + xwk(jd(j)+i-1)
 170     continue

c********End of Step 4.*********************************************************
 200  continue

      return
      end
c=======================================================================

      subroutine mgvbs(j, n, icdomn, ipde, niter, c, d, f, b, u,isolv,
     &            ivbs,lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, rhs, actnB)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MGVBS stands for multigrid V- or \-cycle  and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a multigrid method using V- (or\-) 
c     cycle algorithm in the structured mesh.
c     ivbs = 1 : MG V-cycle,  ivbs = 2 : MG \-cycle
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  u(n), b(n), wk(n), rhs(2*n), actnB(2*n), jd(14)

      ismth = 3                  !Choice of smoother, 1, 2, 3, 4.
      imatv = 3                  !Choice of matrix, FE or UPWD, 2, 3
      ifbs3 = 3                  !Choice of smoother in exact solver, 1, 2, 3.
c     isolv = 0, 1, 2, 3, 4      !Choice of solver at the coarest level.

c      if (ismth .eq. 3) then
c         tol   = 1.d-12
c      else 
         tol   = 1.d-10      !tol - stopping criteria for iterations
c      end if

      itmax = niter
      xnewn = 0.d0
      stopl = 1.d0

      do 30 i = 1, jd(j) - 1
         rhs(i)    = 0.d0
 30   continue

      do 40 i = 1, n
         u(i)  = 0.d0
         wk(i) = 0.d0
         rhs(jd(j)-1+i) = b(i)
 40   continue

c*    Iterate until convergence:  u = u + B(b - (A_0 + cA_x0 + dA_y0 + fM_0)u).

      do 180 kk = 1, itmax
c*       V-cycle of the MG, i.e., the action of the iterator B.
c        Going downward in the V-cycle.

         do 50 i = 1, jd(j+1) - 1
            actnB(i)  = 0.d0
 50      continue

         do 100 k = j, lvlcst+1, -1
            kd = jd(k+1) - jd(k)

c           Presmoothing by the Gauss-seidel smoother nsmthg times.
	    if (ismth.eq.1 .or. ismth.eq.4) then 
               if (ismth .eq. 1) then
                  call gausei(k,rhs(jd(k)),actnB(jd(k)),icdomn,ipde,
     &                        c,d,f,ifbs1)
               else if (ismth .eq. 4) then
                  call gsupwd1(k,kd,rhs(jd(k)),actnB(jd(k)),icdomn,ipde,
     &                         ifbs1,c,d,f,)
               end if
c              call bpxin(k, actnB(jd(k)), icdomn)
	       do 80 ii = 2, nsmthg
c                 Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0+fM_0)*actnB.
                  call matvc2(k, actnB(jd(k)), wk, icdomn, ipde,c,d,f)

c                 Residual: wk = rhs - wk = rhs-(A_0+cA_x0+dA_y0+fM_0)*actnB.
                  do 60 i = 1, kd 
                     wk(i) = rhs(jd(k)-1+i) - wk(i)
 60               continue

                  if (ismth .eq. 1) then
                     call gausei(k,wk,wk,icdomn,ipde,c,d,f,ifbs1)
                  else if (ismth .eq. 4) then
                     call gsupwd1(k,kd,wk,wk,icdomn,ipde,ifbs1,c,d,f)
                  end if
c                 call bpxin(k, wk, icdomn)

c                 Correction step.
                  do 70 i = 1, kd 
                     actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 70               continue
 80            continue
	    else if (ismth .eq. 2) then
               call gssor(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                    ipde, ifbs1, c, d, f, wk, nsmthg)
	    else if (ismth .eq. 3) then
               call gsupwd(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                    ipde, ifbs1, c, d, f, wk, nsmthg)
	    end if 

            if (imatv .eq. 2) then
c              Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
               call matvc2(k, actnB(jd(k)), wk, icdomn, ipde, c, d, f)
            else if (imatv .eq. 3) then
c              Call matvc3 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
               call matvc3(k, actnB(jd(k)), wk, icdomn, ipde, c, d, f)
            end if

c           Residual: wk = rhs - wk =  rhs - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
            do 90 i = 1, kd 
               wk(i) = rhs(jd(k)-1+i) - wk(i)
 90         continue

c           Restriction to the lower level.
            call restrn(k, wk, rhs(jd(k-1)))
            call bpxin(k-1, rhs(jd(k-1)), icdomn)
 100     continue

         if (isolv .eq. 0) then
c           Solving exactly at the coarsest level, i.e., 
c           actnB = (A_0 + cA_x0 + dA_y0 + fM_0)^(-1)*rhs.
            call solexc(rhs(jd(lvlcst)), actnB(jd(lvlcst)), icdomn,
     &                  lvlcst, c, d, f)
         else if (isolv .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x+dA_0y+fM_0)*actnB = rhs
c           by Bi-CGSTAB method in the structured mesh.
            nd = jd(lvlcst+1) - jd(lvlcst)
      	    call bcgsto(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &                  0, nit, icdomn, ipde, c, d, f) 
         else if (isolv .eq. 2) then
c           Solve the discretized problem (A_0 + cA_0x+dA_0y+fM_0)*actnB = rhs
c           by GSUPWD method in the structured mesh.
            nd = jd(lvlcst+1) - jd(lvlcst)
	    call gsupwd(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &                  icdomn, ipde, ifbs3, c, d, f, wk, 2500)
         else if (isolv .eq. 3) then
c           Solve the discretized problem (A_0 + cA_0x+dA_0y+fM_0)*actnB = rhs
c           by Gauss-Seidel method in the structured mesh.
            nd = jd(lvlcst+1) - jd(lvlcst)
            call gssor(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &              icdomn, ipde, ifbs3, c, d, f, wk, 2500)
         end if

c        Going upward in the V-cycle.
         do 150 k = lvlcst+1, j 
            kd = jd(k+1) - jd(k)

c           Prolongation to the upper level.
            call prolng(k, actnB(jd(k-1)), wk)
            call bpxin(k, wk, icdomn)

c           Correction step.
            do 110 i = 1, kd 
               actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 110        continue

            if (ivbs .eq. 2)  go to 150           !No post smoothing in \-cycle

	    if (ismth.eq.1 .or. ismth.eq.4) then
               do 140 ii = 1, nsmthg
c                 Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0+fM_0)*actnB.
                  call matvc2(k, actnB(jd(k)), wk, icdomn, ipde, c,d,f)

c                 Residual: wk = rhs - wk =  rhs - (A_0+cA_x0+dA_y0+fM_0)*u.
                  do 120 i = 1, kd 
                     wk(i) = rhs(jd(k)-1+i) - wk(i)
 120              continue

c                 Postsmoothing by the Gauss-seidel smoother.
                  if (ismth .eq. 1) then
                     call gausei(k,wk,wk,icdomn,ipde,c,d,f,ifbs2)
                  else if (ismth .eq. 4) then
                     call gsupwd1(k,kd,wk,wk,icdomn,ipde,ifbs2,c,d,f)
                  end if
c                 call bpxin(k, wk, icdomn)

c                 Final correction step.
                  do 130 i = 1, kd 
                     actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 130              continue
 140           continue
	    else if (ismth .eq. 2) then
               call gssor(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                    ipde, ifbs2, c, d, f, wk, nsmthg)
	    else if (ismth .eq. 3) then
               call gsupwd(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                    ipde, ifbs2, c, d, f, wk, nsmthg)
	    end if
 150     continue 

c*       New solution:  u = u + B(b - (A_0 + cA_x0 + dA_y0 + fM_0)u).

         do 160 i = 1, n
            u(i) = u(i) + actnB(jd(j)-1+i)
 160     continue

         xoldn = xnewn
         xnewn = dot(u, u, n)

	 xdiffn = dot(actnB(jd(j)), actnB(jd(j)), n)
	 stopl  = xdiffn/xnewn

c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         call matvc2(j, u, wk, icdomn, ipde, c, d, f)

c        Residual: rhs =  b - wk =  b - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         do 170 i = 1, n
            rhs(jd(j)-1+i) = b(i) - wk(i)
 170     continue

         rnorm = dot(rhs(jd(j)), rhs(jd(j)), n)

	 km = mod(kk, 100)

	 if (ivbs .eq. 1 .and. km .eq. 0) then
	    print 300, j, kk, stopl, rnorm
	 else if (ivbs .eq. 2 .and. km .eq. 0) then
	    print 310, j, kk, stopl, rnorm
	 end if

         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 180  continue

      print*, 'Iteration limit exceeded:', itmax

 300  format(' Level, no. of iters, relative, residual at MG V-cycle:',
     &       i2, i5, 2e10.3)
 310  format(' Level, no. of iters, relative, residual at MG BS-cycle:',
     &       i2, i5, 2e10.3)

 333  niter = min(kk, itmax)

      return
      end
c=======================================================================

      subroutine mgv1(j, n, icdomn, ipde, niter, c, d, f, b, u, ivbs,
     &                lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, rhs, actnB)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MGV1 stands for multigrid V-cycle and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a multigrid method using V-cycle
c     algorithm in the structured mesh. It is same as MGV.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim=1050630) 

      implicit double precision (a-h, o-z)
      dimension  u(n), b(n), wk(n), rhs(2*n), actnB(2*n), jd(14)
      dimension  un(lim), bn(lim)

      tol   = 1.d-12/(4**j)            !tol - stopping criteria for iterations
c      itmax = 50                        !max. iteration
      itmax = niter
      xnewn = 0.d0
      stopl = 1.d0

      do 10 i = 1, n
         u(i)  = 0.d0
         bn(i) = b(i)
 10   continue

c*    Iterate until convergence:  u = u + B(b - (A_0 + cA_x0 + dA_y0 + fM_0)u).

      do 180 kk = 1, itmax
c        Iterator using MG V-cycle, i.e., un = B*bn.
         call premg(j, n, icdomn, ipde, c, d, f, bn, un, 0,ivbs,
     &               lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, rhs, actnB)

c*       New solution:  u = u + B(b - (A_0 + cA_x0 + dA_y0 + fM_0)u).

         do 160 i = 1, n
            u(i) = u(i) + un(i)
 160     continue

         xoldn = xnewn
         xnewn = dot(u, u, n)

	 xdiffn = dot(un, un, n)
	 stopl  = xdiffn/xnewn

c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         call matvc2(j, u, wk, icdomn, ipde, c, d, f)

c        Residual: bn =  b - wk =  b - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         do 170 i = 1, n
            bn(i) = b(i) - wk(i)
 170     continue

         rnorm = dot(bn, bn, n)
	 if (ivbs .eq. 1) then
            print 300, j, kk, stopl, rnorm
	 else if (ivbs .eq. 2) then
            print 310, j, kk, stopl, rnorm
	 end if

         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 180  continue

      print*, 'Iteration limit exceeded:', itmax

 300  format(' Level, no. of iters, relative, residual at MG V-cycle:',
     &       i2, i5, 2e10.3)
 310  format(' Level, no. of iters, relative, residual at MG BS-cycle:',
     &       i2, i5, 2e10.3)

 333  niter = min(kk, itmax)
      print*

      return
      end
c=======================================================================

      subroutine mgvf(j, n, icdomn, ipde, niter, c, d, f, b, u, ivbs,
     &           lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, rhs, actnB, bwk)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MGVF stands for full multigrid V-cycle and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a full multigrid method using
c     V-cycle algorithm in the structured mesh.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim=1050630) 

      implicit double precision (a-h, o-z)
      dimension  u(n), b(n), wk(n), rhs(2*n), actnB(2*n), jd(14)
      dimension  un(lim), bn(lim), bwk(2*n)

      isolv = 0                  !Choice of solver at the coarest level.
      ifbs3 = 1                  !Choice of smoother in exact solver, 1, 2, 3.

      do 30 i = 1, jd(j) - 1
         bwk(i) = 0.d0
 30   continue

      do 40 i = 1, n
         u(i)           = 0.d0
         wk(i)          = 0.d0
         un(i)          = 0.d0
         bn(i)          = 0.d0
         bwk(jd(j)-1+i) = b(i)
 40   continue

      do 60 i = j, lvlcst+1, -1
c        Restriction of b to the lower levels.
         call restrn(i,bwk(jd(i)),bwk(jd(i-1)))
         call bpxin(i-1, bwk(jd(i-1)), icdomn)
 60   continue

      if (isolv .eq. 0) then
c        Solving exactly at the coarsest level, i.e., 
c        u = (A_0 + cA_x0 + dA_y0 + fM_0)^(-1)*bwk.
         call solexc(bwk(jd(lvlcst)), u, icdomn, lvlcst, c, d, f)
      else if (isolv .eq. 1) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*u = bwk
c        by Bi-CGSTAB method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
      	 call bcgsto(lvlcst, nd, bwk(jd(lvlcst)), u,
     &               0, nit, icdomn, ipde, c, d, f) 
      else if (isolv .eq. 2) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*u = bwk
c        by GSUPWD method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
	 call gsupwd(lvlcst, nd, bwk(jd(lvlcst)), u,
     &               icdomn, ipde, ifbs3, c, d, f, wk, 2500)
      else if (isolv .eq. 3) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*u = bwk
c        by Gauss-Seidel method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
         call gssor(lvlcst, nd, bwk(jd(lvlcst)), u,
     &              icdomn, ipde, ifbs3, c, d, f, wk, 2500)
      end if

      do 200 k = lvlcst+1, j                                !Full MG starts.
         tol   = 1.d-50/(4**k)      !tol - stopping criteria for iterations
c        itmax = 50                                         !max. iteration
         itmax = niter
         xnewn = 0.d0
         stopl = 1.d0
         kd    = jd(k+1) - jd(k)
         kdl   = jd(k) - jd(k-1)

         do 70 ii = 1, kdl
            un(ii) = u(ii) 
 70      continue

c        Prolongation to the upper level.
         call prolng(k, un, u)
         call bpxin(k, u, icdomn)

c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         call matvc2(k, u, wk, icdomn, ipde, c, d, f)

c        Residual: bn = bwk - wk =  bwk - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         do 80 ii = 1, kd 
            bn(ii) = bwk(jd(k)-1+ii) - wk(ii)
 80      continue

c*       Iterate until convergence: u = u+B(bwk -(A_0 + cA_x0 + dA_y0 + fM_0)u).

         do 180 kk = 1, itmax
c           Iterator using MG V-cycle, i.e., un = B*bn.
            call premg(k, kd, icdomn, ipde, c, d, f, bn, un,isolv,ivbs,
     &                lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, rhs, actnB)

c*          New solution (final correction step):  u = u + un.
            do 160 i = 1, kd
               u(i) = u(i) + un(i)
 160        continue

            xoldn = xnewn
            xnewn = dot(u, u, kd)

            xdiffn = dot(un, un , kd)
            stopl  = xdiffn/xnewn
 
c           Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*u.
            call matvc2(k, u, wk, icdomn, ipde, c, d, f)

c           Residual: bn =  bwk - wk =  b - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
            do 170 i = 1, kd
               bn(i) = bwk(jd(k)-1+i) - wk(i)
 170        continue

            rnorm = dot(bn, bn, kd)
c	    print*, ' rnorm at level', k,':', rnorm

ccc         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 180     continue

c        print*, ' Iteration limit exceeded:', itmax

 333     nit1 = min(kk, itmax)
c        print 190, k, nit1
c190     format(' No. of iterations of Full MG at level ',i2,':',i3)
c        print*
 200  continue

      return
      end
c=======================================================================

      subroutine premg(j, n, icdomn, ipde, c, d, f, b, u, isolv, ivbs,
     &               lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, rhs, actnB)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     PREMG stands for multigrid algorithms as preconditioners, i.e., u=Bb
c     ivbs = 1 : V-cycle is used;   ivbs = 2 : \-cycle is used.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  u(n), b(n), wk(n), rhs(2*n), actnB(2*n), jd(14)

      ismth = 3                  !Choice of smoother, 1, 2, 3, 4.
      imatv = 2                  !Choice of matrix, FE or UPWD, 2, 3, 4
      ifbs3 = 3                  !Choice of smoother in exact solver, 1, 2, 3.
c     isolv = 0, 1, 2, 3, 4      !Choice of solver at the coarest level.

      do 30 i = 1, jd(j+1) - 1
         rhs(i)    = 0.d0
         actnB(i)  = 0.d0
 30   continue

      do 40 i = 1, n
         u(i)  = 0.d0
         wk(i) = 0.d0
         rhs(jd(j)-1+i) = b(i)
 40   continue

c*    V-cycle of the MG, i.e., the action of the iterator B.
c     Going downward in the V-cycle.

      do 100 k = j, lvlcst+1, -1
         kd = jd(k+1) - jd(k)

c        Presmoothing by the (symmetric) Gauss-seidel smoother nsmthg times.
	 if (ismth .eq. 1) then 
            call gausei(k, rhs(jd(k)), actnB(jd(k)), icdomn, ipde,
     &                  c, d, f, ifbs1)
c           call bpxin(k, actnB(jd(k)), icdomn)

cc	    nt = nsmthg*(2**(j-k))             !Variable \- or V-cycle version.
cc	    do 80 ii = 2, nt                   !Variable \- or V-cycle version.
	    do 80 ii = 2, nsmthg
c              Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0+fM_0)*actnB.
               call matvc2(k, actnB(jd(k)), wk, icdomn, ipde,c,d,f)

c              Residual: wk = rhs - wk = rhs-(A_0+cA_x0+dA_y0+fM_0)*actnB.
               do 60 i = 1, kd 
                  wk(i) = rhs(jd(k)-1+i) - wk(i)
 60            continue

	       call gausei(k, wk, wk, icdomn, ipde, c, d, f, ifbs1)
c              call bpxin(k, wk, icdomn)

c              Correction step.
               do 70 i = 1, kd 
                  actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 70            continue
 80         continue
	 else if (ismth .eq. 2) then
            call gssor(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                 ipde, ifbs1, c, d, f, wk, nsmthg)
	 else if (ismth .eq. 3) then
            call gsupwd(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                 ipde, ifbs1, c, d, f, wk, nsmthg)
	 else if (ismth .eq. 4) then
	    if (k .eq. j) then
               call gssor(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                    ipde, ifbs1, c, d, f, wk, nsmthg)
	    else
               call gsupwd(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                     ipde, ifbs1, c, d, f, wk, nsmthg)
	    end if  
	 end if

         if (imatv .eq. 2) then
c           Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
            call matvc2(k, actnB(jd(k)), wk, icdomn, ipde, c, d, f)
         else if (imatv .eq. 3) then
c           Call matvc3 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
            call matvc3(k, actnB(jd(k)), wk, icdomn, ipde, c, d, f)
         else if (imatv .eq. 4) then
            if (k .eq. j) then
c              Call matvc2 to compute wk=(A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
               call matvc2(k, actnB(jd(k)), wk, icdomn, ipde,c,d,f)
            else
c              Call matvc3 to compute wk=(A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
               call matvc3(k, actnB(jd(k)), wk, icdomn, ipde,c,d,f)
            end if
         end if
            
c        Residual: wk = rhs - wk =  rhs - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         do 90 i = 1, kd 
            wk(i) = rhs(jd(k)-1+i) - wk(i)
 90      continue

c        Restriction to the lower level.
         call restrn(k, wk, rhs(jd(k-1)))
         call bpxin(k-1, rhs(jd(k-1)), icdomn)
 100  continue

      if (isolv .eq. 0) then
c        Solving exactly at the coarsest level, i.e., 
c        actnB = (A_0 + cA_x0 + dA_y0 + fM_0)^(-1)*rhs.
         call solexc(rhs(jd(lvlcst)), actnB(jd(lvlcst)), icdomn,
     &               lvlcst, c, d, f)
      else if (isolv .eq. 1) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*actnB = rhs
c        by Bi-CGSTAB method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
      	 call bcgsto(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &               0, nit, icdomn, ipde, c, d, f) 
      else if (isolv .eq. 2) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*actnB = rhs
c        by GSUPWD method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
	 call gsupwd(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &               icdomn, ipde, ifbs3, c, d, f, wk, 2500)
      else if (isolv .eq. 3) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*actnB = rhs
c        by Gauss-Seidel method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
         call gssor(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &              icdomn, ipde, ifbs3, c, d, f, wk, 2500)
      end if

c     Going upward in the V-cycle.
      do 150 k = lvlcst+1, j 
         kd = jd(k+1) - jd(k)

c        Prolongation to the upper level.
         call prolng(k, actnB(jd(k-1)), wk)
         call bpxin(k, wk, icdomn)

c        Correction step.
         do 110 i = 1, kd 
            actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 110     continue

         if (ivbs .eq. 2)  go to 150           !No post smoothing in \-cycle

	 if (ismth .eq. 1) then
cc  	    nt = nsmthg*(2**(j-k))             !Variable V-cycle version.
cc	    do 140 ii = 1, nt                  !Variable V-cycle version.
            do 140 ii = 1, nsmthg
c              Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0+fM_0)*actnB.
               call matvc2(k, actnB(jd(k)), wk, icdomn, ipde, c,d,f)

c              Residual: wk = rhs - wk =  rhs - (A_0+cA_x0+dA_y0+fM_0)*u.
               do 120 i = 1, kd 
                  wk(i) = rhs(jd(k)-1+i) - wk(i)
 120           continue

c              Postsmoothing by the symmetric Gauss-seidel smoother.
               call gausei(k, wk, wk, icdomn, ipde, c, d, f, ifbs2)
c              call bpxin(k, wk, icdomn)

c              Final correction step.
               do 130 i = 1, kd 
                  actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 130           continue
 140        continue
	 else if (ismth .eq. 2) then
            call gssor(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                 ipde, ifbs2, c, d, f, wk, nsmthg)
	 else if (ismth .eq. 3) then
            call gsupwd(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                 ipde, ifbs2, c, d, f, wk, nsmthg)
	 else if (ismth .eq. 4) then
	    if (k .eq. j) then
               call gssor(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                    ipde, ifbs2, c, d, f, wk, nsmthg)
	    else
               call gsupwd(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                     ipde, ifbs2, c, d, f, wk, nsmthg)
	    end if 
	 end if

 150  continue 

c*    Action of preconditioner:  u = B(b).
      do 160 i = 1, n
         u(i) = actnB(jd(j)-1+i)
 160  continue

      return
      end
c=======================================================================

      subroutine solexc(b, x, icdomn, iexct, c, d, f)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the system (A_h)*actnB = rhs exactly at the 
c     coarsest level by using Gaussian elimination with partial pivoting.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  b(25), x(25), acm(3,4), acm2(9,10), soln(9)

c     Solving exactly at the coarsest level, i.e.,
c     x = (A_0 + cA_x0 + dA_y0 + fM_0)^(-1)*b.

      if (icdomn .eq. 3) then
         if (iexct .eq. 1) then
c           Compute b*1/(8/3 + f/9)
            x(5) = b(5)*9/(24 + f)
         else if (iexct .eq. 2) then
c           Call crsmst to generate the coarsest level matrix. 
            call crsmst(9, 10, 9, acm2, b, c, d, f)

c           To solve a linear system by Gauss elimination with partial pivoting.
            call gauspp(9, 10, 9, acm2, soln)

            x( 7) = soln(1)
            x( 8) = soln(2)
            x( 9) = soln(3)
            x(12) = soln(4)
            x(13) = soln(5)
            x(14) = soln(6)
            x(17) = soln(7)
            x(18) = soln(8)
            x(19) = soln(9)
         end if
      else if (icdomn .eq. 4) then
         x(13) = b(13)*36/(96 + f)
      else if (icdomn .eq. 1) then
c        Call crsmun to generate the coarsest level matrix. 
         call crsmun(3, 4, 3, acm, b, c, d, f)
c        To solve a linear system by Gauss elimination with partial pivoting.
         call gauspp(3, 4, 3, acm, soln)
         x( 7) = soln(1)
         x( 8) = soln(2)
         x(12) = soln(3)
      end if

      return
      end
c=======================================================================

      subroutine crsmun(nrows, ncols, ndim, cm, b, c, d, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine generates the coarsest level matrix in multigrid 
c     algorithm.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension cm(nrows, ncols), b(12)

      do 2 i=1,nrows
	 do 1 j=1,ncols
	    cm(i,j) = 0.d0
  1      continue
  2   continue

      e1d3   = 1.d0/3
      e8d3   = 8*e1d3
      cd12   = c/12
      dd12   = d/12
      cmdd48 = (c-d)/48
      fd36   = f/36
      fd144  = f/144
      fd576  = f/576
      
      cm(1,1) = e8d3 + fd36
      cm(2,1) = - e1d3 - cd12 + fd144
      cm(3,1) = - e1d3 - dd12 + fd144
      cm(1,2) = - e1d3 + cd12 + fd144
      cm(2,2) = cm(1,1)
      cm(3,2) = - e1d3 + cmdd48 + fd576
      cm(1,3) = - e1d3 + dd12 + fd144
      cm(2,3) = - e1d3 - cmdd48 + fd576
      cm(3,3) = cm(1,1)
      cm(1,4) = b(7)
      cm(2,4) = b(8)
      cm(3,4) = b(12)

      return
      end
c=======================================================================

      subroutine crsmst(nrows, ncols, ndim, cm, b, c, d, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine generates the coarsest level matrix in multigrid 
c     algorithm.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension cm(nrows,ncols), b(25)

      do 2 i=1,nrows
	 do 1 j=1,ncols
	    cm(i,j) = 0.d0
  1      continue
  2   continue
	
      e1d3   = 1.d0/3
      e8d3   = 8*e1d3
      cd12   = c/12
      dd12   = d/12
      cmdd48 = (c-d)/48
      cpdd48 = (c+d)/48
      fd36   = f/36
      fd144  = f/144
      fd576  = f/576
      cmt    = e8d3 + fd36
      cmt0   = - e1d3 + fd144
      cmt1   = cmt0 + cd12
      cmt2   = cmt0 - cd12
      cmt3   = cmt0 + dd12
      cmt4   = cmt0 - dd12

      do 10 i = 1, 9
         cm(i, i) = cmt
 10   continue

      do 30 i = 1, 3
         kt = 3*(i-1)
         lt = 2*(i-1)
         do 20 j = 1, 2
            k = kt + j
            l = lt + j
            cm(k, k+1) = cmt1
            cm(k+1, k) = cmt2
            cm(l, l+3) = cmt3
            cm(l+3, l) = cmt4
 20      continue
 30   continue

      cmt0   = - e1d3 + fd576
      cmt1   = cmt0 + cpdd48
      cmt2   = cmt0 - cmdd48
      cmt3   = cmt0 + cmdd48
      cmt4   = cmt0 - cpdd48

      do 50 i = 1, 2
         kt = 3*(i-1)
         do 40 j = 1, 2
            k = kt + j
            cm(k, k+4)   = cmt1
            cm(k+1, k+3) = cmt2
            cm(k+3, k+1) = cmt3
            cm(k+4, k)   = cmt4
 40      continue
 50   continue

      cm(1,10) = b( 7)
      cm(2,10) = b( 8)
      cm(3,10) = b( 9)
      cm(4,10) = b(12)
      cm(5,10) = b(13)
      cm(6,10) = b(14)
      cm(7,10) = b(17)
      cm(8,10) = b(18)
      cm(9,10) = b(19)

      return
      end
c=======================================================================

      subroutine gauspp(nrows, ncols, nd, cm, soln)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine uses the Gaussian elimination with partial pivoting to 
c     solve a linear system. All coefficients and the load are stored in the 
c     augmented matrix cm.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension cm(nrows, ncols), soln(nrows)
      logical error 

      npivot = 1
      error  = .false.

 10   if (npivot.lt.nd .and. .not.error) then
c        !To reorder the equations so that the pivot position in the pivot 
c        !equation has the maximum absolute value.
         maxrow = npivot
         do 20 irow = npivot + 1, nd
            if (dabs(cm(irow, npivot)) .gt. dabs(cm(maxrow, npivot))) 
     &           maxrow = irow
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

         if (.not.error) then
c           !To eliminate the element in the pivot position from rows following
c           !the pivot equation.
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

      if (error) then
         print*, ' No unique solution exists.'
      else
c        !To perform the back-substitution to determine the solution to the 
c        !system of equations.
         do 70 irow = nd, 1, -1
            do 60 icol = nd, irow + 1, -1
               cm(irow,nd+1) = cm(irow,nd+1) - soln(icol)*cm(irow,icol)
 60         continue
            soln(irow) = cm(irow, nd+1)/cm(irow, irow)
 70      continue
      end if

      return
      end
c=======================================================================

      subroutine gs(l, n, b, x, ichoic, kchoic, c, d, f, wk, niter) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the linear system A*x = b by Gauss-Seidel
c     method in the structured mesh.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension b(n), x(n), wk(n), xm(513, 513)

      tol   = 1.d-50/(4**l)             !tol - stopping criteria for iterations
c      itmax = 200000                  !max. iteration
      itmax = 1
      ld = 2**l + 1
      m1 = ld - 1
      m2 = m1*m1
      xnn   = 0.d0

      do 30 j = 1,ld
         kt = (j-1)*ld
	 do 20 i = 1,ld
            k = kt + i
            xm(i,j) = x(k)
 20      continue
 30   continue

      do 70 kk = 1, itmax
         xdiffn = 0.d0
         do 50 j = 2, m1
            npost = j*j
            kt = (j-1)*ld
            do 40 i = 2, m1
               npos = npost + i*i
               k = kt + i
               if (ichoic .eq. 1 .and. npos .gt. m2) then
                  xm(i,j) = 0.d0
               else if (ichoic.eq.4 .and. 
     &              (j.eq.2 .or. j.eq.m1 .or. i.eq.2 .or. i.eq.m1)) then
                  xm(i,j) = 0.d0
               else
                  xm(i,j) = (3*b(k) + xm(i-1,j-1) + xm(i-1,j)
     &                    + xm(i-1,j+1) + xm(i,j+1) + xm(i+1,j+1) 
     &                    + xm(i+1,j) + xm(i+1,j-1) + xm(i,j-1))/8
                  xdiff = xm(i,j) - x(k)
                  xdiffn = xdiffn + xdiff*xdiff
               end if
               x(k) = xm(i,j)
 40         continue
 50      continue

         xon = xnn
	 xnn = dot(x, x, n)

	 stopl  = xdiffn/xnn

c        Call matvc2 to compute wk = (A_0)*x.
         call matvc2(l, x, wk, ichoic, kchoic, c, d, f)

c        Residual: wk =  b - wk =  b - (A_0)*x.
         do 60 i = 1, n
            wk(i) = b(i) - wk(i)
 60      continue

         rnorm = dot(wk, wk, n)
c	 print*, ' rnorm in GS', rnorm
cc         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 100
 70   continue

c      print*, ' Iteration limit exceeded:', itmax

 100  niter = min(kk, itmax)
c      print 110, l, niter
c 110  format(' No. of iterations in G-S Method at level ', i2, ':',i6)
c      print*

      return
      end
c=======================================================================

      subroutine gssor(l, n, b, x, icdomn, ipde,ifbs,c,d,f,wk,niter) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the linear system A*x = b by Gauss-Seidel or SOR
c     method in the structured mesh.
c     ifbs = 1, 2, 3 : Forward SOR, Backward SOR, SSOR
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension b(n), x(n), wk(n), xm(1025, 1025)

      tol   = 1.d-24           !tol - stopping criteria for iterations
c      itmax = 1500                     !max. iteration
      itmax = niter
      ld = 2**l + 1
      m1 = ld - 1
      m2 = m1*m1               !h**2 = 1/m2
      xnn   = 0.d0

c     c=-1000,d=1002,w=0.03811(3),0.06363(4),0.02634(3s),0.04326(4s),0.08938(5s)
c     c=-100, d=200, w=0.067(2s),w=0.1533(3s),w=0.2687(4s),w=0.5000(5s)
c     c=-30, d=40, w=0.067(2s),w=0.1533(3s),w=0.2687(4s),w=0.5000(5s)

cc      if (l .eq. 2) then
cc         w = 0.067d0                           !relaxation parameter
cc      else if (l .eq. 3) then
cc         w = 0.1533d0
cc      else if (l .eq. 4) then
cc         w = 0.2687d0
cc      else if (l .eq. 5) then
cc         w = 0.5000d0
cc      else if (l .eq. 8) then
cc         w = 1.900d0
cc      end if
cc      print*, w

      w = 1.d0
                
      do 30 j = 1,ld
         kt = (j-1)*ld
	 do 20 i = 1,ld
            k = kt + i
            xm(i,j) = x(k)
 20      continue
 30   continue

      e1d3   = 1.d0/3
      chd12  = c/(12*m1)
      dhd12  = d/(12*m1)
      fhhd36 = f/(36*m2)
      amm    = - e1d3 - chd12 - dhd12 + fhhd36
      a0m    = - e1d3 - 4*dhd12 + 4*fhhd36
      apm    = - e1d3 + chd12 - dhd12 + fhhd36
      am0    = - e1d3 - 4*chd12 + 4*fhhd36
      a00    = 8*e1d3 + 16*fhhd36
      ap0    = - e1d3 + 4*chd12 + 4*fhhd36
      amp    = - e1d3 - chd12 + dhd12 + fhhd36
      a0p    = - e1d3 + 4*dhd12 + 4*fhhd36
      app    = - e1d3 + chd12 + dhd12 + fhhd36

      do 200 kk = 1, itmax
         xdiffn = 0.d0
	 if (ifbs .eq. 2)  go to 60

         do 50 j = 2, m1          !Forward SOR
            npost = j*j
            kt = (j-1)*ld
            do 40 i = 2, m1
               npos = npost + i*i
               k = kt + i
               if (icdomn .eq. 1 .and. npos .gt. m2) then
                  xm(i,j) = 0.d0
               else if (icdomn.eq.4 .and. 
     &              (j.eq.2 .or. j.eq.m1 .or. i.eq.2 .or. i.eq.m1)) then
                  xm(i,j) = 0.d0
               else
                  xt = w*(b(k) - amm*xm(i-1,j-1) - am0*xm(i-1,j)
     &               - amp*xm(i-1,j+1) - a0p*xm(i,j+1)-app*xm(i+1,j+1) 
     &               - ap0*xm(i+1,j) - apm*xm(i+1,j-1) 
     &               - a0m*xm(i,j-1) - a00*xm(i,j))/a00
                  xm(i,j) = xm(i,j) + xt
                  xdiffn = xdiffn + xt*xt
               end if
	       x(k) = xm(i,j)
 40         continue
 50      continue

 	 if (ifbs .eq. 1)  go to 90

 60	 continue
         do 80 j = m1, 2, -1           !Backward SOR
            npost = j*j
            kt = (j-1)*ld
            do 70 i = m1, 2, -1
               npos = npost + i*i
               k = kt + i
               if (icdomn .eq. 1 .and. npos .gt. m2) then
                  xm(i,j) = 0.d0
               else if (icdomn.eq.4 .and. 
     &              (j.eq.2 .or. j.eq.m1 .or. i.eq.2 .or. i.eq.m1)) then
                  xm(i,j) = 0.d0
               else
                  xt = w*(b(k) - amm*xm(i-1,j-1) - am0*xm(i-1,j)
     &               - amp*xm(i-1,j+1) - a0p*xm(i,j+1)-app*xm(i+1,j+1) 
     &               - ap0*xm(i+1,j) - apm*xm(i+1,j-1) 
     &               - a0m*xm(i,j-1) - a00*xm(i,j))/a00
                  xm(i,j) = xm(i,j) + xt
                  xdiffn = xdiffn + xt*xt
               end if
	       x(k) = xm(i,j)
 70         continue
 80      continue

 90      continue

	 if (itmax .le. 16) go to 200

         xon = xnn
	 xnn = dot(x, x, n)
	 stopl  = xdiffn/xnn

         if (stopl .lt. tol .or. kk .eq. 1) then
c           Call matvc2 to compute wk = (A_0 + cA_0x + dA_0y + fM_0)*x.
            call matvc2(l, x, wk, icdomn, ipde, c, d, f)

c           Residual: wk =  b - wk =  b - (A_0)*x.
            do 120 i = 1, n
               wk(i) = b(i) - wk(i)
 120        continue

            rnorm = dot(wk, wk, n)
	    print 300, l, kk, stopl, rnorm
            if (rnorm .lt. tol)  go to 333
	 end if
 200  continue

 300  format(' Level, no. of iters, relative, residual at GSSOR: ',
     &       i2, i5, 2e10.3)

      if (itmax .le. 16) go to 555
      print*, ' Iteration limit exceeded:', itmax
 333  nit = min(kk, itmax)
      print*

 555  return
      end
c=======================================================================

      subroutine gsupwd(l, n, b, x, icdomn, ipde,ifbs,c,d,f,wk,niter) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the linear system A*x = b by Gauss-Seidel
c     method in the structured mesh, where A is formed by five point FDS and
c     upwind scheme.
c     ifbs = 1, 2, 3 : Forward SOR, Backward SOR, SSOR
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension b(n), x(n), wk(n), xm(1025, 1025)

      tol   = 1.d-26           !tol - stopping criteria for iterations
c      itmax = 1500                     !max. iteration
      itmax = niter
      ld = 2**l + 1
      m1 = ld - 1
      m2 = m1*m1               !h**2 = 1/m2
      xnn   = 0.d0

c      do 10 i = 1, n
c	 x(i) = 0.d0
c 10   continue

      do 30 j = 1,ld
         kt = (j-1)*ld
	 do 20 i = 1,ld
            k = kt + i
            xm(i,j) = x(k)
 20      continue
 30   continue

      ch  = c/m1
      dh  = d/m1
      fhh = f/m2

      if (c .gt. 0.d0) then
	 chm = - ch
	 ch0 =   ch
	 chp = 0.d0
      else
	 chm = 0.d0
	 ch0 = - ch
	 chp =   ch
      end if

      if (d .gt. 0.d0) then
	 dhm = - dh
	 dh0 =   dh
	 dhp = 0.d0
      else
	 dhm = 0.d0
	 dh0 = - dh
	 dhp =   dh
      end if

      a0m    = - 1.d0 + dhm
      am0    = - 1.d0 + chm
      a00    =   4.d0 + ch0 + dh0 + fhh
      ap0    = - 1.d0 + chp
      a0p    = - 1.d0 + dhp

      do 200 kk = 1, itmax
         xdiffn = 0.d0
	 if (ifbs .eq. 2)  go to 60

         do 50 j = 2, m1          !Forward SOR
            npost = j*j
            kt = (j-1)*ld
            do 40 i = 2, m1
               npos = npost + i*i
               k = kt + i
               if (icdomn .eq. 1 .and. npos .gt. m2) then
                  xm(i,j) = 0.d0
               else if (icdomn.eq.4 .and. 
     &              (j.eq.2 .or. j.eq.m1 .or. i.eq.2 .or. i.eq.m1)) then
                  xm(i,j) = 0.d0
               else
                  xt = (b(k) - am0*xm(i-1,j) - a0p*xm(i,j+1)
     &               - ap0*xm(i+1,j) - a0m*xm(i,j-1) - a00*xm(i,j))/a00
                  xm(i,j) = xm(i,j) + xt
                  xdiffn = xdiffn + xt*xt
               end if
	       x(k) = xm(i,j)
 40         continue
 50      continue

 	 if (ifbs .eq. 1)  go to 90

 60	 continue
         do 80 j = m1, 2, -1           !Backward SOR
            npost = j*j
            kt = (j-1)*ld
            do 70 i = m1, 2, -1
               npos = npost + i*i
               k = kt + i
               if (icdomn .eq. 1 .and. npos .gt. m2) then
                  xm(i,j) = 0.d0
               else if (icdomn.eq.4 .and. 
     &              (j.eq.2 .or. j.eq.m1 .or. i.eq.2 .or. i.eq.m1)) then
                  xm(i,j) = 0.d0
               else
                  xt = (b(k) - am0*xm(i-1,j) - a0p*xm(i,j+1)
     &               - ap0*xm(i+1,j) - a0m*xm(i,j-1) - a00*xm(i,j))/a00
                  xm(i,j) = xm(i,j) + xt
                  xdiffn = xdiffn + xt*xt
               end if
	       x(k) = xm(i,j)
 70         continue
 80      continue

 90      continue

	 if (itmax .le. 100) go to 200

         xon = xnn
	 xnn = dot(x, x, n)
	 stopl  = xdiffn/xnn

ccc	 if (stopl .lt. tol .or. l .eq. 1) then
ccc	    print 300, l, kk, stopl, rnorm
ccc	    go to 333
ccc	 end if
ccc	 go to 200

         if (stopl .lt. tol .or. l .eq. 1) then
c           Call matvc3 to compute wk = (A_0 + cA_0x + dA_0y + fM_0)*x.
            call matvc3(l, x, wk, icdomn, ipde, c, d, f)

c           Residual: wk =  b - wk =  b - (A_0)*x.
            do 120 i = 1, n
               wk(i) = b(i) - wk(i)
 120        continue

            rnorm = dot(wk, wk, n)
ccc	    print 300, l, kk, stopl, rnorm
	    go to 333
	 end if
 200  continue

 300  format(' Level, no. of iters, relative, residual at GSUPWD: ',
     &       i2, i5, 2e10.3)

      if (itmax .le. 100) go to 555
      print*, ' Iteration limit exceeded:', itmax
 333  nit = min(kk, itmax)

 555  return
      end
c=======================================================================

       subroutine bcgstr(j, n, b, x, itol, nit, icdomn, ipde,c,d,f) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the non-symmetric/indefinite linear system Ax = b 
c     by using Bi-CGSTAB (bi-conjugate gradient stabilized) method. It is a
c     restarted version of the subroutine BCGST.
c     lim = 1089(5), 4225(6), 16641(7), 66049(8), 263169(9), 1050625(10)
c       irst = 0 : start with initial x = 0 
c       irst = 1 : restart with updated x.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      parameter (lim = 263169) 
      parameter (lim = 1050625) 
      implicit double precision (a-h, o-z)
      dimension x(n), b(n), r(lim), p(lim), rs(lim), ap(lim)
      dimension t(lim), ax(lim), wk(lim)

      tol  = 1.d-16/(4**j)       !stopping criterion for iterations
      itmax = 1500               !max iterations
c      itmax = nit
      irst = 0
      kt = 0

 10   if (irst .eq. 0) then                         !Initial start
         do 20 i = 1, n
            x(i) = 0.d0
            r(i) = b(i)
 20      continue
      else if (irst .eq. 1) then                    !Restart
         call matvc2(j, x, ax, icdomn, ipde, c, d, f)
c        To compute r = b - ax.
	 call saxpy(-1.d0, ax, b, n, r)
      end if

      do 30 i = 1, n
	 p(i)  = 0.d0
         ap(i) = 0.d0
	 rs(i) = r(i)
 30   continue

      alpha = 1.d0
      omega = 1.d0
      rhoo  = 1.d0
      rhon  = dot(rs, r, n)
      rnorm = rhon
      xnn   = 0.d0

      do 200 k = 1, itmax
c         print*
c         if (k .ge. 200) print*, 'omega', omega

	 if (dabs(omega) .le. 1.d-22) then
	    print*, ' Bi-CGSTAB breakdown, omega is too small.'
            print*, ' Restart with updated solution.'
            irst = 1
	    kt = kt + k - 1
	    go to 10 
	 end if

         beta = (rhon/rhoo)*(alpha/omega)

c        To compute p = r + beta*(p-omega*ap).
	 call saxpy(-omega, ap, p, n, p)
	 call saxpy(beta, p, r, n, p)

c        To compute ap = A*p.
         call matvc2(j, p, ap, icdomn, ipde, c, d, f)

         tau = dot(rs, ap, n)
c	 if (k .ge. 200) print*, 'tau1', tau

	 if (dabs(tau) .le. 1.d-22) then
	    print*, ' Bi-CGSTAB breakdown, tau#1 is too small.'
            print*, ' Restart with updated solution.'
            irst = 1
	    kt = kt + k - 1
	    go to 10
	 end if
         alpha = rhon/tau

c        To compute r = r - alpha*ap.
	 call saxpy(-alpha, ap, r, n, r)

c        To compute t = A*r.
         call matvc2(j, r, t, icdomn, ipde, c, d, f)

         tau = dot(t, t, n)
c	 if (k .ge. 200) print*, 'tau2', tau

	 if (dabs(tau) .le. 1.d-22) then
	    print*, ' Bi-CGSTAB breakdown, tau#2 is too small.'
            print*, ' Restart with updated solution.'
            irst = 1
	    kt = kt + k - 1
	    go to 10
	 end if
         omega = dot(t, r, n)/tau

         rhoo = rhon
         rhon = -omega*dot(rs, t, n)

c        To compute x = x + alpha*p + omega*r
         xon = xnn
         xdn = 0.d0
         xnn = 0.d0
         do 40 i = 1, n
            xt   = alpha*p(i) + omega*r(i)
	    x(i) = x(i) + xt
            xdn  = xdn + xt*xt 
            xnn  = xnn + x(i)*x(i)
   40    continue

	 if (xnn .lt. 1.d-40) then
	    print*, 'Norm of the new solution x is to small.'
	 end if
	 stopl  = xdn/xnn

c	 To compute wk = Ax.
	 call matvc2(j, x, wk, icdomn, ipde, c, d, f)

c	 To compute the residual wk = b - wk = b - Ax.
	 call saxpy(-1.d0, wk, b, n, wk)

	 arnorm = dot(wk, wk, n)
	 print*, 'Actual residual norm at Bi-CGSTAB-R:', arnorm

c        To compute r = r - omega*t.
         call saxpy(-omega, t, r, n, r)

         rnorm = dot(r, r, n)
	 print*, 'rnorm at Bi-CGSTAB-R:', rnorm

c        print*, 'x = ', (jj, '*', x(jj), jj=1,n)
	 if (k .eq. 1) then
            if (arnorm .lt. tol)  go to 333
	 else 
            if (stopl .lt. tol .and. arnorm .lt. tol)  go to 333
	 end if
 200  continue
         
 333  if (irst .eq. 0)  nit1 = min(k, itmax)
      if (irst .eq. 1)  nit = kt + k - 1

      print 444, j, nit, arnorm
 444  format(' No. of iters in Bi-CGSTAB-R at level ',i2,':',i5, e10.3)
      print*

      return
      end
c=======================================================================

      subroutine bcgst(j, n, b, x, itol, nit, icdomn, ipde, c, d, f) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the non-symmetric/indefinite linear system Ax = b 
c     by using Bi-CGSTAB (bi-conjugate gradient stabilized) method. It is a
c     variant of the subroutine BCGSTO.
c     lim = 1089(5), 4225(6), 16641(7), 66049(8), 263169(9), 1050625(10)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      parameter (lim = 263169) 
      parameter (lim = 1050625) 
      implicit double precision (a-h, o-z)
      dimension x(n), b(n), r(lim), p(lim), rs(lim), ap(lim), t(lim)
      dimension wk(lim)

c      tol  = 1.d-17/(4**j)       !stopping criterion for iterations
      tol = 1.d-16
c      if (j .eq. 8)  tol = 1.d-20

      itmax = 2000               !max iterations
c      itmax = nit
	
      do 10 i = 1, n
	 x(i)  = 0.d0
	 p(i)  = 0.d0
         ap(i) = 0.d0
         r(i)  = b(i)
	 rs(i) = r(i)
   10 continue

      alpha = 1.d0
      omega = 1.d0
      rhoo  = 1.d0
      rhon  = dot(rs, r, n)
      rnorm = rhon
      xnn   = 0.d0

      do 200 k = 1, itmax
c	 print*
c	 print*, 'omega', omega

	 if (dabs(omega) .le. 1.d-50) then
	    print*, 'Bi-CGSTAB breakdown, omega is too small'
	    goto 333 
	 end if

         beta = (rhon/rhoo)*(alpha/omega)

c        To compute p = r + beta*(p-omega*ap).
	 call saxpy(-omega, ap, p, n, p)
	 call saxpy(beta, p, r, n, p)

c        To compute ap = A*p.
         call matvc2(j, p, ap, icdomn, ipde, c, d, f)

         tau = dot(rs, ap, n)
c	 print*, 'tau1', tau

	 if (dabs(tau) .le. 1.d-50) then
	    print*, 'Bi-CGSTAB breakdown, tau#1 is too small'
	    goto 333 
	 end if
         alpha = rhon/tau

c        To compute r = r - alpha*ap.
	 call saxpy(-alpha, ap, r, n, r)

c        To compute t = A*r.
         call matvc2(j, r, t, icdomn, ipde, c, d, f)

         tau = dot(t, t, n)
c	 print*, 'tau2', tau

	 if (dabs(tau) .le. 1.d-50) then
	    print*, 'Bi-CGSTAB breakdown, tau#2 is too small'
	    goto 333 
	 end if
         omega = dot(t, r, n)/tau

         rhoo = rhon
         rhon = -omega*dot(rs, t, n)

c        To compute x = x + alpha*p + omega*r
         xon = xnn
         xdn = 0.d0
         xnn = 0.d0
         do 40 i = 1, n
            xt   = alpha*p(i) + omega*r(i)
	    x(i) = x(i) + xt
            xdn  = xdn + xt*xt 
            xnn  = xnn + x(i)*x(i)
   40    continue

	 if (xnn .lt. 1.d-40) then
	    print*, 'Norm of the new solution x is to small.'
	 end if
	 stopl  = xdn/xnn

c	 To compute wk = Ax.
	 call matvc2(j, x, wk, icdomn, ipde, c, d, f)

c	 To compute the residual wk = b - wk = b - Ax.
	 do 50 i = 1, n
	    wk(i) = b(i) - wk(i)
 50	 continue

	 arnorm = dot(wk, wk, n)
cc	 print*, 'Actual residual norm at Bi-CGSTAB-M:', arnorm

c        To compute r = r - omega*t.
         call saxpy(-omega, t, r, n, r)

         rnorm = dot(r, r, n)
cc	 print*, 'rnorm at Bi-CGSTAB-M:', rnorm

	 if (stopl .lt. tol .and. arnorm .lt. tol)  go to 333
 200  continue
         
 333  nit = min(k, itmax)
      print 444, j, nit, arnorm, rnorm
 444  format(' No. of it. in Bi-CGSTAB-M at level ',i2,':',i5, 2e10.3)
      print*
         
      return
      end
c=======================================================================
 
      subroutine bcgsto(j, n, b, x, itol, nit, icdomn, ipde, c, d, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the non-symmetric/indefinite linear system Ax = b 
c     by using Bi-CGSTAB (bi-conjugate gradient stabilized) method.
c     lim = 1089(5), 4225(6), 16641(7), 66049(8), 263169(9), 1050625(10)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      parameter (lim = 263169) 
      parameter (lim = 1050625) 
      implicit double precision (a-h, o-z)
      dimension x(n), b(n), r(lim), p(lim), rs(lim), ap(lim), t(lim)
      dimension wk(lim)

      if (itol .eq. 1) then
c        tol  = 1.d-17/(4**j)   !stopping criterion for iterations
         tol = 1.d-18
      else if (itol .eq. 0) then 
c        tol  = 1.d-17/(4**j) !stopping criterion for iterations
         tol = 1.d-22
      end if 

      itmax = 2000               !max iterations
c      itmax = nit

      do 10 i = 1, n
	 x(i)  = 0.d0
	 p(i)  = 0.d0
         ap(i) = 0.d0
         r(i)  = b(i)
	 rs(i) = r(i)
   10 continue

      alpha = 1.d0
      omega = 1.d0
      rhon  = 1.d0
      rnorm = dot(r, r, n)
      xnn   = 0.d0

      do 200 k = 1, itmax
c	print*, 'omega', omega

	 if (dabs(omega) .le. 1.d-80) then
	    print*, 'Bi-CGSTAB-O breakdown, omega is too small'
	    goto 333 
	 end if

         rhoo = rhon
	 rhon = dot(rs, r, n)

         beta = (rhon/rhoo)*(alpha/omega)

c        To compute p = r + beta*(p-omega*ap).
	 call saxpy(-omega, ap, p, n, p)

	 call saxpy(beta, p, r, n, p)

c        To compute ap = A*p.
         call matvc2(j, p, ap, icdomn, ipde, c, d, f)

         tau = dot(rs, ap, n)
c	 print*, 'tau1', tau

	 if (dabs(tau) .le. 1.d-80) then
	    print*, 'Bi-CGSTAB-O breakdown, tau#1 is too small'
	    goto 333 
	 end if
         alpha = rhon/tau

c        To compute r = r - alpha*ap.
	 call saxpy(-alpha, ap, r, n, r)

c        To compute t = A*r.
         call matvc2(j, r, t, icdomn, ipde, c, d, f)

         tau = dot(t, t, n)
c	 print*, 'tau2', tau

	 if (dabs(tau) .le. 1.d-80) then
	    print*, 'Bi-CGSTAB-O breakdown, tau#2 is too small'
	    goto 333 
	 end if
         omega = dot(t, r, n)/tau

c        To compute x = x + alpha*p + omega*r
         xon = xnn
         xdn = 0.d0
         xnn = 0.d0
         do 40 i = 1, n
            xt   = alpha*p(i) + omega*r(i)
	    x(i) = x(i) + xt
            xdn  = xdn + xt*xt 
            xnn  = xnn + x(i)*x(i)
   40    continue

c        To compute r = r - omega*t.
         call saxpy(-omega, t, r, n, r)

cc	 if (xnn .lt. 1.d-60) then
cc	    print*, 'Norm of the new solution x is to small.'
cc	 end if

	 stopl  = xdn/xnn

ccc	 print*, 'stop loop', stopl
	 if (stopl .lt. tol) then
            rnorm = dot(r, r, n)
  	    print*, 'Norm of r at Bi-CGSTAB-O:', k, rnorm

	    if (rnorm .lt. tol) go to 333

            if (itol .eq. 0) go to 200

c	    To compute wk = Ax.
	    call matvc2(j, x, wk, icdomn, ipde, c, d, f)

c	    To compute the residual wk = b - wk = b - Ax.
	    do 50 i = 1, n
	       wk(i) = b(i) - wk(i)
 50	    continue

	    arnorm = dot(wk, wk, n)
	    print*, 'Actual residual norm at Bi-CGSTAB-O:', k,arnorm

	    if (arnorm .lt. tol) go to 333
	 end if
 200  continue

 333  nit = min(k, itmax)
      print 444, j, nit, rnorm
 444  format(' No. of it. in Bi-CGSTAB-O at level ',i2,':',i5, e10.3)
      print*

      return
      end
c=======================================================================
 
      subroutine bcgstp(j, n, b, x, itol, nit, icdomn, ipde, c, d, f,
     &                  jd, bwk, xwk, wk, wk1, wk2, ipcnd, jc, nc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the non-symmetric/indefinite linear system Ax = b 
c     by using preconditioned Bi-CGSTAB (bi-conjugate gradient stabilized)
c     method.
c     lim = 1089(5), 4225(6), 16641(7), 66049(8), 263169(9), 1050625(10)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      parameter (lim = 263169) 
      parameter (lim = 1050625) 
      implicit double precision (a-h, o-z)
      dimension x(n), b(n), r(lim), p(lim), rs(lim), ap(lim), t(lim)
      dimension br(lim), bp(lim), bt(lim), jd(14)
      dimension bwk(2*n), xwk(2*n), wk1(2*n), wk2(2*n), wk(n)

c      tol  = 1.d-17/(4**j)       !stopping criterion for iterations
      tol = 1.d-20
      tol1 = 1.d-18
      itmax = 1500               !max iterations
c      itmax = nit

      lvlcst = 1
      nsmthg = 1
      ifbs1  = 1
      ifbs2  = 1

c     Choice of mg-cycle as a precond; ivbs=1 : V-cycle, ivbs=2 : \-cycle.
      if (ipcnd .eq. 2 .or. ipcnd .eq. 4) then 
         ivbs   = 2                
      else if (ipcnd .eq. 3 .or. ipcnd .eq. 5) then 
         ivbs   = 1
      end if 
  
      do 10 i = 1, n
	 x(i)  = 0.d0
	 p(i)  = 0.d0
	 t(i)  = 0.d0
         ap(i) = 0.d0
         bp(i) = 0.d0
         br(i) = 0.d0
         bt(i) = 0.d0
         r(i)  = b(i)
	 rs(i) = r(i)
   10 continue

c      do 11 i = 1, n/10
c      do 11 i = 1, 100*j
c	 rs(i*10) = rs(i*10) + 10.0
c   11 continue

      alpha = 1.d0
      omega = 1.d0
      rhon  = 1.d0
      rnorm = dot(r, r, n)
      xnn   = 0.d0

      do 200 k = 1, itmax
c	print*, 'omega', omega

	 if (dabs(omega) .le. 1.d-50) then
	    print*, 'Bi-CGSTAB breakdown, omega is too small'
	    goto 333 
	 end if

         rhoo = rhon
	 rhon = dot(rs, r, n)

         beta = (rhon/rhoo)*(alpha/omega)

c        To compute p = r + beta*(p-omega*ap).
	 call saxpy(-omega, ap, p, n, p)

	 call saxpy(beta, p, r, n, p)

         if (ipcnd .eq. 1) then
c           To compute bp = B*p, where B is a BPX preconditioner.
            call bpxg(j, n, jd, p, bp, icdomn)
         else if (ipcnd .eq. 2 .or. ipcnd .eq. 3) then
c           To compute bp = B*p, where B is an MG V- or \-cycle preconditioner.
            call premg(j, n, icdomn, 1, 0.d0, 0.d0, 0.d0, p, bp, ivbs,
     &                  lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         else if (ipcnd .eq. 4 .or. ipcnd .eq. 5) then
c           To compute bp = B*p, where B is an MG V- or \-cycle preconditioner
c           with a coarse grid solver.
            call algm1(j, n, jc, nc, icdomn, 1, c, d, f, p, bp,
     &           ivbs, nsmthg, ifbs1, ifbs2, jd, bwk, xwk, wk, wk1, wk2)
         end if

c        To compute ap = A*bp.
         call matvc2(j, bp, ap, icdomn, ipde, c, d, f)

         tau = dot(rs, ap, n)
c	 print*, 'tau1', tau

	 if (dabs(tau) .le. 1.d-50) then
	    print*, 'Bi-CGSTAB breakdown, tau#1 is too small'
	    goto 333 
	 end if
         alpha = rhon/tau

c        To compute r = r - alpha*ap.
	 call saxpy(-alpha, ap, r, n, r)

         if (ipcnd .eq. 1) then
c           To compute br = B*r, where B is a BPX preconditioner.
            call bpxg(j, n, jd, r, br, icdomn)
         else if (ipcnd .eq. 2 .or. ipcnd .eq. 3) then
c           To compute br = B*r, where B is an MG V- or \-cycle preconditioner.
            call premg(j, n, icdomn, 1, 0.d0, 0.d0, 0.d0, r, br, ivbs, 
     &                  lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         else if (ipcnd .eq. 4 .or. ipcnd .eq. 5) then
c           To compute br = B*r, where B is an MG V- or \-cycle preconditioner
c           with a coarse grid solver.
            call algm1(j, n, jc, nc, icdomn, 1, c, d, f, r, br,
     &           ivbs, nsmthg, ifbs1, ifbs2, jd, bwk, xwk, wk, wk1, wk2)
         end if

c        To compute t = A*br.
         call matvc2(j, br, t, icdomn, ipde, c, d, f)

         if (ipcnd .eq. 1) then
c           To compute bt = B*t, where B is a BPX preconditioner.
            call bpxg(j, n, jd, t, bt, icdomn)
         else if (ipcnd .eq. 2 .or. ipcnd .eq. 3) then
c           To compute bt = B*t, where B is an MG V- or \-cycle preconditioner.
            call premg(j, n, icdomn, 1, 0.d0, 0.d0, 0.d0, t, bt, ivbs,
     &                  lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         else if (ipcnd .eq. 4 .or. ipcnd .eq. 5) then
c           To compute bt = B*t, where B is an MG V- or \-cycle preconditioner
c           with a coarse grid solver.
            call algm1(j, n, jc, nc, icdomn, 1, c, d, f, t, bt,
     &           ivbs, nsmthg, ifbs1, ifbs2, jd, bwk, xwk, wk, wk1, wk2)
         end if

         tau = dot(bt, bt, n)
c	 print*, 'tau2', tau

	 if (dabs(tau) .le. 1.d-50) then
	    print*, 'Bi-CGSTAB breakdown, tau#2 is too small'
	    goto 333 
	 end if
         omega = dot(bt, br, n)/tau

c        To compute x = x + alpha*bp + omega*br
         xon = xnn
         xdn = 0.d0
         xnn = 0.d0
         do 40 i = 1, n
            xt   = alpha*bp(i) + omega*br(i)
	    x(i) = x(i) + xt
            xdn  = xdn + xt*xt 
            xnn  = xnn + x(i)*x(i)
   40    continue

c        To compute r = r - omega*t.
         call saxpy(-omega, t, r, n, r)

cc	 if (xnn .lt. 1.d-40) then
cc	    print*, 'Norm of the new solution x is to small.'
cc	 end if

	 stopl  = xdn/xnn

	 if (stopl .lt. tol) then
c	    To compute wk = Ax.
	    call matvc2(j, x, wk, icdomn, ipde, c, d, f)

c	    To compute the residual wk = b - wk = b - Ax.
	    do 50 i = 1, n
	       wk(i) = b(i) - wk(i)
 50	    continue

	    arnorm = dot(wk, wk, n)
	    print*, 'Actual residual norm at Bi-CGSTAB-P:', arnorm

            rnorm = dot(r, r, n)
  	    print*, 'Norm of r at Bi-CGSTAB-P:', rnorm

  	    if (arnorm .lt. tol .and. rnorm .lt. tol1)  go to 333
	 end if
 200  continue
         
 333  nit = min(k, itmax)
      print 444, j, nit, arnorm, rnorm
 444  format(' No. of it. in Bi-CGSTAB-P at level ',i2,':',i5, 2e10.3)
      print*

      return
      end
c=======================================================================

      subroutine bcgstq(j, n, b, x, itol, nit, icdomn, ipde, c, d, f,
     &                  jd, bwk, xwk, wk, wk1, wk2, ipcnd, jc, nc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the non-symmetric/indefinite linear system Ax = b 
c     by using preconditioned Bi-CGSTAB (bi-conjugate gradient stabilized)
c     method. It is a variant of the subroutine BCGSTP.
c     lim = 1089(5), 4225(6), 16641(7), 66049(8), 263169(9), 1050625(10)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      parameter (lim = 263169)
      parameter (lim = 1050625)  
      implicit double precision (a-h, o-z)
      dimension x(n), b(n), r(lim), p(lim), rs(lim), ap(lim), t(lim)
      dimension br(lim), bp(lim), jd(14)
      dimension bwk(2*n), xwk(2*n), wk1(2*n), wk2(2*n), wk(n)

c      tol  = 1.d-17/(4**j)       !stopping criterion for iterations
      tol = 1.d-20	
      tol1 = 1.d-18
      itmax = 1500               !max iterations
c      itmax = nit

      lvlcst = 1
      nsmthg = 1
      ifbs1  = 1
      ifbs2  = 1

c     Choice of mg-cycle as a precond; ivbs=1 : V-cycle, ivbs=2 : \-cycle.
      if (ipcnd .eq. 2 .or. ipcnd .eq. 4) then 
         ivbs   = 2                
      else if (ipcnd .eq. 3 .or. ipcnd .eq. 5) then 
         ivbs   = 1
      end if 

      do 10 i = 1, n
	 x(i)  = 0.d0
	 p(i)  = 0.d0
	 t(i)  = 0.d0
         ap(i) = 0.d0
         bp(i) = 0.d0
         br(i) = 0.d0
         r(i)  = b(i)
	 rs(i) = r(i)
   10 continue

c      do 11 i = 1, n/10
c      do 11 i = 1, 100*j
c	 rs(i*10) = rs(i*10) + 10.0
c   11 continue

      alpha = 1.d0
      omega = 1.d0
      rhon  = 1.d0
      rnorm = dot(r, r, n)
      xnn   = 0.d0

      do 200 k = 1, itmax
c	print*, 'omega', omega

	 if (dabs(omega) .le. 1.d-50) then
	    print*, 'Bi-CGSTAB breakdown, omega is too small'
	    goto 333 
	 end if

         rhoo = rhon
	 rhon = dot(rs, r, n)

         beta = (rhon/rhoo)*(alpha/omega)

c        To compute p = r + beta*(p-omega*ap).
	 call saxpy(-omega, ap, p, n, p)

	 call saxpy(beta, p, r, n, p)

         if (ipcnd .eq. 1) then
c           To compute bp = B*p, where B is a BPX preconditioner.
            call bpxg(j, n, jd, p, bp, icdomn)
         else if (ipcnd .eq. 2 .or. ipcnd .eq. 3) then
c           To compute bp = B*p, where B is an MG V- or \-cycle preconditioner.
            call premg(j, n, icdomn, 1, 0.d0, 0.d0, 0.d0, p, bp, ivbs,
     &                  lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         else if (ipcnd .eq. 4 .or. ipcnd .eq. 5) then
c           To compute bp = B*p, where B is an MG V- or \-cycle preconditioner
c           with a coarse grid solver.
            call algm1(j, n, jc, nc, icdomn, 1, c, d, f, p, bp,
     &           ivbs, nsmthg, ifbs1, ifbs2, jd, bwk, xwk, wk, wk1, wk2)
         end if

c        To compute ap = A*bp.
         call matvc2(j, bp, ap, icdomn, ipde, c, d, f)

         tau = dot(rs, ap, n)
c	 print*, 'tau1', tau

	 if (dabs(tau) .le. 1.d-50) then
	    print*, 'Bi-CGSTAB breakdown, tau#1 is too small'
	    goto 333 
	 end if
         alpha = rhon/tau

c        To compute r = r - alpha*ap.
	 call saxpy(-alpha, ap, r, n, r)

         if (ipcnd .eq. 1) then
c           To compute br = B*r, where B is a BPX preconditioner.
            call bpxg(j, n, jd, r, br, icdomn)
         else if (ipcnd .eq. 2 .or. ipcnd .eq. 3) then
c           To compute br = B*r, where B is an MG V- or \-cycle preconditioner.
            call premg(j, n, icdomn, 1, 0.d0, 0.d0, 0.d0, r, br, ivbs,
     &                  lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         else if (ipcnd .eq. 4 .or. ipcnd .eq. 5) then
c           To compute br = B*r, where B is an MG V- or \-cycle preconditioner
c           with a coarse grid solver.
            call algm1(j, n, jc, nc, icdomn, 1, c, d, f, r, br,
     &           ivbs, nsmthg, ifbs1, ifbs2, jd, bwk, xwk, wk, wk1, wk2)
         end if

c        To compute t = A*br.
         call matvc2(j, br, t, icdomn, ipde, c, d, f)

         tau = dot(t, t, n)
c	 print*, 'tau2', tau

	 if (dabs(tau) .le. 1.d-50) then
	    print*, 'Bi-CGSTAB breakdown, tau#2 is too small'
	    goto 333 
	 end if
         omega = dot(t, r, n)/tau

c        To compute x = x + alpha*bp + omega*br
         xon = xnn
         xdn = 0.d0
         xnn = 0.d0
         do 40 i = 1, n
            xt   = alpha*bp(i) + omega*br(i)
	    x(i) = x(i) + xt
            xdn  = xdn + xt*xt 
            xnn  = xnn + x(i)*x(i)
   40    continue

c        To compute r = r - omega*t.
         call saxpy(-omega, t, r, n, r)

cc	 if (xnn .lt. 1.d-40) then
cc	    print*, 'Norm of the new solution x is to small.'
cc	 end if

	 stopl  = xdn/xnn

	 if (stopl .lt. tol) then
c	    To compute wk = Ax.
	    call matvc2(j, x, wk, icdomn, ipde, c, d, f)

c	    To compute the residual wk = b - wk = b - Ax.
	    do 50 i = 1, n
	       wk(i) = b(i) - wk(i)
 50	    continue

	    arnorm = dot(wk, wk, n)
	    print*, 'Actual residual norm at Bi-CGSTAB-Q:', arnorm

            rnorm = dot(r, r, n)
  	    print*, 'Norm of r at Bi-CGSTAB-Q:', rnorm

  	    if (arnorm .lt. tol .and. rnorm .lt. tol1)  go to 333
	 end if
 200  continue
         
 333  nit = min(k, itmax)
      print 444, j, nit, arnorm, rnorm
 444  format(' No. of it. in Bi-CGSTAB-Q at level ',i2,':',i5, 2e10.3)
      print*

      return
      end
c=======================================================================

      subroutine bcgsts(j, n, b, x, itol, nit, icdomn, ipde, c, d, f,
     &                  jd, bwk, xwk, wk, wk1, wk2, ipcnd, jc, nc)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the non-symmetric/indefinite linear system Ax = b' 
c     by using Bi-CGSTAB (bi-conjugate gradient stabilized) method. In this 
c     subroutine a preconditioned system is solved, i.e., BAx = b = Bb'.
c     lim = 1089(5), 4225(6), 16641(7), 66049(8), 263169(9), 1050625(10)
c     isolv : choice of coarse grid solver
c     isolv = 0 : solve directly at level 1 or at level 2 by Gaussian 
c                 elimination with partial pivoting
c     isolv = 1 : solve by Bi-CGSTAB
c     isolv = 2 : solve by GSUPWD(Upwind scheme + G-S iteration)
c     isolv = 3 :
c     isolv = 4 : 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      parameter (lim = 263169)
      parameter (lim = 1050625)  
      implicit double precision (a-h, o-z)
      dimension x(n), b(n), r(lim), p(lim), rs(lim), ap(lim), t(lim)
      dimension bp(lim), jd(14)
      dimension bwk(2*n), xwk(2*n), wk1(2*n), wk2(2*n), wk(n)

      tol  = 1.d-8/(4**j)       !stopping criterion for iterations
c      tol = 1.d-12
      itmax = 2500               !max iterations
c      itmax = nit

      nsmthg = 1
      ifbs1  = 3
      ifbs2  = 3
      isolv  = 2
cc      if (jc .eq. 1) isolv = 0

c     Choice of mg-cycle as a precond; ivbs=1 : V-cycle, ivbs=2 : \-cycle.
      if (ipcnd .eq. 2 .or. ipcnd .eq. 4 .or. ipcnd .eq. 6) then 
         ivbs   = 2                
      else if (ipcnd .eq. 3 .or. ipcnd .eq. 5 .or. ipcnd .eq. 7) then 
         ivbs   = 1
      end if 

      if (ipcnd .eq. 1) then
c        To compute r = B*b, where B is a BPX preconditioner.
         call bpxg(j, n, jd, b, r, icdomn)
      else if (ipcnd .eq. 2 .or. ipcnd .eq. 3) then
c        To compute r = B*b, where B is an MG V- or \-cycle preconditioner.
         call premg(j, n, icdomn, 1, 0.d0, 0.d0, 0.d0, b, r, 0,ivbs, 
     &              1, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
      else if (ipcnd .eq. 4 .or. ipcnd .eq. 5) then
c        To compute r = B*b, where B is an MG V- or \-cycle preconditioner
c        with a coarse grid solver.
         call algm1(j, n, jc, nc, icdomn, 1, c, d, f, b, r,isolv,
     &          ivbs, nsmthg, ifbs1, ifbs2, jd, bwk, xwk, wk, wk1, wk2)
      else if (ipcnd .eq. 6 .or. ipcnd .eq. 7) then
c        To compute r = B*b, where B is a nonsymmetric/indefinite MG V- or
c        \-cycle preconditioner.
         call premg(j, n, icdomn, 3, c, d, f, b, r, isolv, ivbs, 
     &              jc, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
      end if

      do 10 i = 1, n
	 x(i)  = 0.d0
	 p(i)  = 0.d0
         ap(i) = 0.d0
         bp(i) = 0.d0
c	 rs(i) = b(i)
	 rs(i) = r(i)
   10 continue

c      do 11 i = 1, n/10
c	 rs(i*10) = rs(i*10)/2
c   11 continue

      alpha = 1.d0
      omega = 1.d0
      rhon  = 1.d0
      rnorm = dot(r, r, n)
      xnn   = 0.d0

      do 200 k = 1, itmax
c	 print*, 'omega', omega

	 if (dabs(omega) .le. 1.d-80) then
	    print*, 'Bi-CGSTAB-S breakdown, omega is too small'
	    goto 333
	 end if

         rhoo = rhon
	 rhon = dot(rs, r, n)

         beta = (rhon/rhoo)*(alpha/omega)

c        To compute p = r + beta*(p-omega*ap).
	 call saxpy(-omega, ap, p, n, p)

	 call saxpy(beta, p, r, n, p)

c        To compute bp = A*p.
         call matvc2(j, p, bp, icdomn, ipde, c, d, f)

         if (ipcnd .eq. 1) then
c           To compute ap = B*bp, where B is a BPX preconditioner.
            call bpxg(j, n, jd, bp, ap, icdomn)
         else if (ipcnd .eq. 2 .or. ipcnd .eq. 3) then
c           To compute ap = B*bp, where B is an MG V- or \-cycle preconditioner.
            call premg(j, n, icdomn, 1, 0.d0, 0.d0, 0.d0, bp, ap,0,ivbs, 
     &                  1, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         else if (ipcnd .eq. 4 .or. ipcnd .eq. 5) then
c           To compute ap = B*bp, where B is an MG V- or \-cycle preconditioner
c           with a coarse grid solver.
            call algm1(j, n, jc, nc, icdomn, 1, c, d, f, bp, ap,isolv,
     &           ivbs, nsmthg, ifbs1, ifbs2, jd, bwk, xwk, wk, wk1, wk2)
         else if (ipcnd .eq. 6 .or. ipcnd .eq. 7) then
c           To compute ap = B*bp, where B is a nonsymmetric/indefinite MG V- or
c           \-cycle preconditioner.
            call premg(j, n, icdomn, 3, c, d, f, bp, ap, isolv,ivbs, 
     &                 jc, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         end if

         tau = dot(rs, ap, n)
c	 print*, 'tau1', tau

	 if (dabs(tau) .le. 1.d-80) then
	    print*, 'Bi-CGSTAB-S breakdown, tau#1 is too small'
	    goto 333 
	 end if
         alpha = rhon/tau

c        To compute r = r - alpha*ap.
	 call saxpy(-alpha, ap, r, n, r)

c        To compute bp = A*r.
         call matvc2(j, r, bp, icdomn, ipde, c, d, f)

         if (ipcnd .eq. 1) then
c           To compute t = B*bp, where B is a BPX preconditioner.
            call bpxg(j, n, jd, bp, t, icdomn)
         else if (ipcnd .eq. 2 .or. ipcnd .eq. 3) then
c           To compute t = B*bp, where B is an MG V- or \-cycle preconditioner.
            call premg(j, n, icdomn, 1, 0.d0, 0.d0, 0.d0, bp, t,0,ivbs,
     &                 1, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         else if (ipcnd .eq. 4 .or. ipcnd .eq. 5) then
c           To compute t = B*bp, where B is an MG V- or \-cycle preconditioner
c           with a coarse grid solver.
            call algm1(j, n, jc, nc, icdomn, 1, c, d, f, bp, t,isolv,
     &           ivbs, nsmthg, ifbs1, ifbs2, jd, bwk, xwk, wk, wk1, wk2)
         else if (ipcnd .eq. 6 .or. ipcnd .eq. 7) then
c           To compute t = B*bp, where B is a nonsymmetric/indefinite MG V- or
c           \-cycle preconditioner.
            call premg(j, n, icdomn, 3, c, d, f, bp, t, isolv,ivbs, 
     &                 jc, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         end if

         tau = dot(t, t, n)
c	 print*, 'tau2', tau

	 if (dabs(tau) .le. 1.d-80) then
	    print*, 'Bi-CGSTAB-S breakdown, tau#2 is too small'
	    goto 333 
	 end if
         omega = dot(t, r, n)/tau

c        To compute x = x + alpha*p + omega*r
         xon = xnn
         xdn = 0.d0
         xnn = 0.d0
         do 40 i = 1, n
            xt   = alpha*p(i) + omega*r(i)
	    x(i) = x(i) + xt
            xdn  = xdn + xt*xt 
            xnn  = xnn + x(i)*x(i)
   40    continue

c        To compute r = r - omega*t.
         call saxpy(-omega, t, r, n, r)

	 stopl  = xdn/xnn

	 if (stopl .lt. tol .or. k .eq. 1) then
            rnorm = dot(r, r, n)

ccc	    if (rnorm .lt. tol) go to 333

c	    To compute wk = Ax.
	    call matvc2(j, x, wk, icdomn, ipde, c, d, f)

c	    To compute the residual wk = b - wk = b - Ax.
	    do 50 i = 1, n
	       wk(i) = b(i) - wk(i)
 50	    continue

	    arnorm = dot(wk, wk, n)
	    print 300, j, k, stopl, rnorm, arnorm   

	    if (arnorm .lt. tol) go to 333
	 end if
 200  continue

      print*, 'Iteration limit exceeded:', itmax

 300  format(' Level, no. of iters, relative, temporary residual,
     & actual residual at Bi-CGSTAB-S: ', i2, i5, 3e10.3)

 333  nit = min(k, itmax)
      print*

      return
      end
c=======================================================================

      subroutine bpxg(j, n, jd, x, y, ichoic)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes y = B*x, where B is a multilevel nodal basis
c     preconditioner (BPX) with a symmetric Gauss-seidel smoother.
c         nsmthg = the number of smoothings, 1, 2, or 3.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	

      parameter( nsmthg = 1, ifbs = 1, lim = 1050630)
      implicit double precision (a-h, o-z)
      dimension x(n), y(n), wk(2*lim), jd(14), bwk(lim), wk1(lim)

      igs = 1                      !igs = 1 or 2 : choice of G-S iteration.

      do 20 i = 1, jd(j)
        wk(i) = 0.d0
   20 continue

      do 30 i = 1, n
         y(i) = 0.d0
         wk(jd(j)-1+i) = x(i)
   30 continue

      do 40 l = j-1, 1, -1
	 call restrn(l+1, wk(jd(l+1)), wk(jd(l)))
c*       Call BPXIN to apply the BPX preconditioner inside the domain only.
	 call bpxin(l, wk(jd(l)), ichoic)
   40 continue

      jl  = jd(2) - jd(1)
      if (igs .eq. 1) then
         do 50 i = 1, nsmthg
	    call gausei(1, wk(1), wk(1), ichoic, 1, 
     &                  0.d0, 0.d0, 0.d0, ifbs)
c	    call bpxin(1, wk(1), ichoic)
   50    continue
      else if (igs .eq. 2) then
	    do 45 ii = 1, jl
	       bwk(ii) = wk(ii)
	       wk(ii)  = 0.d0
 45	    continue
            call gssor(1, jl, bwk, wk(1), ichoic, 1, ifbs, 
     &                 0.d0, 0.d0, 0.d0, wk1, nsmthg) 
      end if

      do 70 l = 2, j
         jl  = jd(l+1) - jd(l)
	 call prolng(l, wk(jd(l-1)), y)

	 if (igs .eq. 1) then
	    do 55 i = 1, nsmthg
	       call gausei(l, wk(jd(l)), wk(jd(l)), ichoic, 1,
     &                     0.d0, 0.d0, 0.d0, ifbs)
c	       call bpxin(l, wk(jd(l)), ichoic)
   55       continue
	 else if (igs .eq. 2) then
	    do 51 ii = 1, jl
	       bwk(ii) = wk(jd(l)-1+ii)
	       wk(jd(l)-1+ii) = 0.d0
 51	    continue
            call gssor(l, jl, bwk, wk(jd(l)), ichoic, 1, ifbs, 
     &                 0.d0, 0.d0, 0.d0, wk1, nsmthg) 
	 end if 

         do 60 i = 1, jl
            wk(jd(l)-1+i) = y(i) + wk(jd(l)-1+i)
   60    continue
	 call bpxin(l, wk(jd(l)), ichoic)
   70 continue
	  
      do 120 i = 1, n
         y(i) = wk(jd(j)-1+i) 
  120 continue

      return
      end
c=======================================================================

      subroutine bpxin(l, q, ichoic)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine adjusts the l_th level vector q near the boundary so that
c     a BPX preconditioner can be applied to the appropriate part of the domain.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	

      implicit double precision (a-h, o-z)
      dimension  q((2**l + 1)**2)

      ll = 2**l + 1
      ld = (ll-1)**2

      do 20 j = 1, ll
         kt = (j-1)*ll
         do 10 i = 1, ll
            k = kt + i
	    if (ichoic .eq. 1) then
c*             q(k) = 0 outside the domain and on/near the boundary.
	       if (j.eq.1 .or. j.eq.ll .or. i.eq.1 .or. i.eq.ll) then
	          q(k) = 0.d0
	       else
 	          npos = i*i + j*j
	          if (npos .gt. ld)  q(k) = 0.d0
	       end if
	    else if (ichoic .eq. 3) then
c*             q(k) = 0 on the boundary of the unit square.
	       if (j.eq.1 .or. j.eq.ll) then
	          q(k) = 0.d0
	       else if (i.eq.1 .or. i.eq.ll) then
 	          q(k) = 0.d0
               end if
	    else if (ichoic .eq. 4) then
c*             q(k) = 0 on/near the boundary of the unit square.
	       if (j.eq.1.or.j.eq.2.or.j.eq.(ll-1) .or. j.eq.ll) then
	          q(k) = 0.d0
	       else if(i.eq.1.or.i.eq.2.or.i.eq.(ll-1).or.i.eq.ll) then
 	          q(k) = 0.d0
               end if
            end if
   10    continue
   20 continue

      return
      end
c=======================================================================

      subroutine restrn(l, p, q)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                     l   t               l   t
c     This subroutine computes q = ( I   )  * p, where ( I   )  is the l_th
c                                     l-1                 l-1
c     level restriction operator of the nodal basis preconditioner(BPX).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
      implicit double precision (a-h, o-z)
      dimension p((2**l+1)**2), q((2**(l-1)+1)**2), z(1025, 1025)

      lh = 2**l+1
      ll = 2**(l-1)+1

      do 5 j = 1, ll*ll
	 q(j) = 0.d0
    5 continue

c     Matrix representation(z) of the vector p.

      do 20 j = 1,lh
         kt = (j-1)*lh
	 do 10 i = 1,lh
	    k = kt + i
            z(i,j) = p(k)
   10    continue
   20 continue

      do 50 js = 2, ll-1
         jb = 2*js - 1
	 kt = (js-1)*ll
	 do 40 is = 2, ll-1
            ib = 2*is - 1
	    k  = kt + is
	    q(k) = z(ib,jb) + (z(ib,jb-1) + z(ib,jb+1) + z(ib-1,jb)      
     &           + z(ib+1,jb))/2 + (z(ib-1,jb-1) + z(ib-1,jb+1)  
     &           + z(ib+1,jb-1) + z(ib+1,jb+1))/4 
   40    continue
   50 continue

      return
      end
c=======================================================================

      subroutine prolng(l, p, q)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                     l                   l  
c     This subroutine computes q = ( I   )  * p, where ( I   ) is the l_th
c                                     l-1                 l-1
c     level prolongation operator of the nodal basis preconditioner(BPX).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	
      implicit double precision (a-h, o-z)
      dimension q((2**l+1)**2), p((2**(l-1)+1)**2), z(513, 513)

      ld = 2**(l-1)
      l1 = ld + 1
      l2 = 2**l + 1

      do 5 i = 1, l2*l2
	 q(i) = 0.d0
    5 continue

c     Matrix representation(z) of the vector p.

      do 20 j = 1,l1
         kt = (j-1)*l1
	 do 10 i = 1,l1
	    k = kt + i
            z(i,j) = p(k)
   10    continue
   20 continue

      do 40 j = 2, ld
c        jb = 2*j - 1, kt = (jb-1)*l2
         kt = (2*j-2)*l2
	 do 30 i = 1, ld
	    k      = kt + 2*i
	    q(k-1) = z(i,j)
	    q(k)   = (z(i,j) + z(i+1,j))/2      
   30    continue
c	 q(kt+l2)  = z(l1,j)
   40 continue

      do 60 j = 1, ld
c        jb = 2*j, kt = (jb-1)*l2
         kt = (2*j-1)*l2
	 do 50 i = 1, ld
	    k      = kt + 2*i
	    q(k-1) = (z(i,j) + z(i,j+1))/2
	    q(k)   = (z(i,j) + z(i+1,j) + z(i,j+1) + z(i+1,j+1))/4      
   50    continue
c	 q(kt+l2)  = (z(l1,j) + z(l1,j+1))/2
   60 continue

      return
      end
c=======================================================================

      subroutine gausei(l, x, y, ichoic, kchoic, c, d, f, ifbs)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                          -1          -1  
c     This subroutine computes y = (D - U )  D (D - L )  *x, in the structured 
c                                    l   l    l  l   l 
c     mesh, i.e., the action of the symmetric Gauss-Seidel smoother, where
c     A  = D - L - U  is defined as follows:
c      l    l   l   l
c     kchoic = 1:  A = the stiffness matrix,
c     kchoic = 2:  A = A + c*A_x + d*A_y,
c     kchoic = 3:  A = A + c*A_x + d*A_y + f*M, where M is the mass matrix.
c     ifbs   = 1:  Forward G-S
c     ifbs   = 2:  Backward G-S
c     ifbs   = 3:  Symmetric G-S
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension x((2**l+1)**2), y((2**l+1)**2), xm(1025, 1025)

      ld = 2**l + 1
      m1 = ld - 1
      m2 = m1*m1

      if (ifbs .eq. 0) then
    	 print*, ' G-S iteration is not performed.'
         go to 120
      end if

c     Matrix representation(xm) of some components (corresponding to the 
c     inner nodes) of the vector x.

      do 20 j = 2,ld-1
         kt = (j-1)*ld
	 do 10 i = 2,ld-1
	    k = kt + i
            xm(i,j) = x(k)
 10      continue
 20   continue

      do 30 j = 1,ld
         xm(1,j)  = 0.d0
         xm(ld,j) = 0.d0
         xm(j,1)  = 0.d0
         xm(j,ld) = 0.d0
 30   continue

      do 35 j = 1, ld*ld
         y(j) = 0.d0
 35   continue  

      e1d3   = 1.d0/3
      chd12  = c/(12*m1)
      dhd12  = d/(12*m1)
      fhhd36 = f/(36*m2)
      amm    = - e1d3 - chd12 - dhd12 + fhhd36
      a0m    = - e1d3 - 4*dhd12 + 4*fhhd36
      apm    = - e1d3 + chd12 - dhd12 + fhhd36
      am0    = - e1d3 - 4*chd12 + 4*fhhd36
      a00    = 8*e1d3 + 16*fhhd36
      ap0    = - e1d3 + 4*chd12 + 4*fhhd36
      amp    = - e1d3 - chd12 + dhd12 + fhhd36
      a0p    = - e1d3 + 4*dhd12 + 4*fhhd36
      app    = - e1d3 + chd12 + dhd12 + fhhd36

      if (ifbs .eq. 2) go to 75

c               << Computation of x = (D_l - L_l)^-1 * x >>
      do 50 j = 2, m1                                      !For j = 2,3,...,ld-1
	 npost = j*j
	 do 40 i = 2, m1
 	    npos = npost + i*i
	    if (ichoic .eq. 1 .and. npos .gt. m2) then
               xm(i,j) = 0.d0
	    else if (ichoic.eq.4 .and. 
     &           (j.eq.2 .or. j.eq.m1 .or. i.eq.2 .or. i.eq.m1)) then
               xm(i,j) = 0.d0
	    else
               xm(i,j) = (xm(i,j) - amm*xm(i-1,j-1) - a0m*xm(i,j-1)  
     &                 - apm*xm(i+1,j-1) - am0*xm(i-1,j))/a00
	    end if
 40      continue
 50   continue

      if (ifbs .eq. 1) go to 95

c              << Computation of x = D_l * x >>
      do 70 j = 2, m1 
         do 60 i = 2, m1
            xm(i,j) = xm(i,j)*a00
 60      continue
 70   continue

c               << Computation of x = (D_l - U_l)^-1 * x >>
 75   do 90 j = m1, 2, -1                                  !For j = ld-1,...,3,2
	 npost = j*j
	 do 80 i = m1, 2, -1
 	    npos = npost + i*i
	    if (ichoic .eq. 1 .and. npos .gt. m2) then
               xm(i,j) = 0.d0
	    else if (ichoic.eq.4 .and. 
     &           (j.eq.2 .or. j.eq.m1 .or. i.eq.2 .or. i.eq.m1)) then
               xm(i,j) = 0.d0
	    else
               xm(i,j) = (xm(i,j) - ap0*xm(i+1,j) - amp*xm(i-1,j+1)  
     &                 - a0p*xm(i,j+1) - app*xm(i+1,j+1))/a00
	    end if 
   80    continue
   90 continue

c     Back to the vector y from the matrix xm.
 95   do 110 j = 2,ld-1
         kt = (j-1)*ld
         do 100 i = 2,ld-1
            k = kt + i
            y(k) = xm(i,j)
 100     continue
 110  continue

 120  return
      end
c=======================================================================

      subroutine gsupwd1(l,n,x,y,icdomn,ipde,ifbs,c,d,f) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                          -1          -1  
c     This subroutine computes y = (D - U )  D (D - L )  *x if ifbs=3, in the 
c                                    l   l    l  l   l 
c     structured mesh, i.e., the action of the symmetric Gauss-Seidel smoother,
c     where A  = D - L - U is formed by five point FDS and upwind scheme.
c            l    l   l   l
c     ipde = 1:  A = the stiffness matrix,
c     ipde = 2:  A = A + c*A_x + d*A_y,
c     ipde = 3:  A = A + c*A_x + d*A_y + f*M, where M is the mass matrix.
c     ifbs = 1:  Forward G-S
c     ifbs = 2:  Backward G-S
c     ifbs = 3:  Symmetric G-S
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension y(n), x(n), xm(1025, 1025)

      ld = 2**l + 1            !n = ld*ld
      m1 = ld - 1
      m2 = m1*m1               !h**2 = 1/m2

      if (ifbs .eq. 0) then
    	 print*, ' G-S iteration is not performed.'
         go to 120
      end if

c     Matrix representation(xm) of some components (corresponding to the 
c     inner nodes) of the vector x.

      do 20 j = 2,ld-1
         kt = (j-1)*ld
	 do 10 i = 2,ld-1
	    k = kt + i
            xm(i,j) = x(k)
 10      continue
 20   continue

      do 30 j = 1,ld
         xm(1,j)  = 0.d0
         xm(ld,j) = 0.d0
         xm(j,1)  = 0.d0
         xm(j,ld) = 0.d0
 30   continue

      do 35 j = 1, n
         y(j) = 0.d0
 35   continue  

      ch  = c/m1
      dh  = d/m1
      fhh = f/m2

      if (c .gt. 0.d0) then
	 chm = - ch
	 ch0 =   ch
	 chp = 0.d0
      else
	 chm = 0.d0
	 ch0 = - ch
	 chp =   ch
      end if

      if (d .gt. 0.d0) then
	 dhm = - dh
	 dh0 =   dh
	 dhp = 0.d0
      else
	 dhm = 0.d0
	 dh0 = - dh
	 dhp =   dh
      end if

      a0m    = - 1.d0 + dhm
      am0    = - 1.d0 + chm
      a00    =   4.d0 + ch0 + dh0 + fhh
      ap0    = - 1.d0 + chp
      a0p    = - 1.d0 + dhp

      if (ifbs .eq. 2) go to 75

c               << Computation of x = (D_l - L_l)^-1 * x >>
      do 50 j = 2, m1                                      !For j = 2,3,...,ld-1
	 npost = j*j
	 do 40 i = 2, m1
 	    npos = npost + i*i
	    if (icdomn .eq. 1 .and. npos .gt. m2) then
               xm(i,j) = 0.d0
	    else if (icdomn.eq.4 .and. 
     &           (j.eq.2 .or. j.eq.m1 .or. i.eq.2 .or. i.eq.m1)) then
               xm(i,j) = 0.d0
	    else
               xm(i,j) = (xm(i,j) - a0m*xm(i,j-1) - am0*xm(i-1,j))/a00
	    end if
 40      continue
 50   continue

      if (ifbs .eq. 1) go to 95

c              << Computation of x = D_l * x >>
      do 70 j = 2, m1 
         do 60 i = 2, m1
            xm(i,j) = xm(i,j)*a00
 60      continue
 70   continue

c               << Computation of x = (D_l - U_l)^-1 * x >>
 75   do 90 j = m1, 2, -1                                  !For j = ld-1,...,3,2
	 npost = j*j
	 do 80 i = m1, 2, -1
 	    npos = npost + i*i
	    if (icdomn.eq.4 .and. 
     &           (j.eq.2 .or. j.eq.m1 .or. i.eq.2 .or. i.eq.m1)) then
               xm(i,j) = 0.d0
	    else
               xm(i,j) = (xm(i,j) - ap0*xm(i+1,j) - a0p*xm(i,j+1))/a00
	    end if 
 80      continue
 90   continue

c     Back to the vector y from the matrix xm.
 95   do 110 j = 2,ld-1
         kt = (j-1)*ld
         do 100 i = 2,ld-1
            k = kt + i
            y(k) = xm(i,j)
 100     continue
 110  continue

 120  return
      end
c=======================================================================

      subroutine saxpy(a, x, y, m, z)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes Z = a*X + Y
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision(a-h, o-z)
      dimension x(m), y(m), z(m)

      do 10 i = 1,m
	 z(i) = a*x(i) + y(i)
   10 continue

      return
      end 
c=======================================================================

      double precision function dot(x, y, m)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This function computes inner product of x and y.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision(a-h, o-z)
      dimension x(m), y(m)

      dot = 0.d0
      do 10 i = 1,m
	 dot = dot + x(i)*y(i)
 10   continue

      return
      end
c=======================================================================

      subroutine matvc2(j, x, w, ichoic, kchoic, c, d, f)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes matrix_vector multiplication(W=B*x) in the
c     structured grid, where B is defined as follows:     
c         kchoic = 1 : B = A;
c         kchoic = 2 : B = A + cA_x + dA_y;
c         kchoic = 3 : B = A + cA_x + dA_y + fM;
c         kchoic = 4 : B = cA_x + dA_y;
c         kchoic = 5 : B = cA_x + dA_y + fM;
c         kchoic = 6 : B = fM;
c         kchoic = 7 : B = A + fM;
c     where A is the stiffness matrix, A_x & A_y are lower order matrices, and
c     M is the mass matix.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension x((2**j+1)**2), w((2**j+1)**2), xm(1025, 1025)

      m  = 2**j + 1
      m1 = m - 1
      m2 = m1*m1
c     h = 1.d0/m1

      do 10 i = 1, m*m
	 w(i) = 0.d0
   10 continue

c*    Matrix representation xmat of the vector x.

      do 30 l = 2, m1
	 kt = (l-1)*m 
	 do 20 i = 2, m1
	    k = kt + i
	    xm(i,l) = x(k)
   20    continue
   30 continue

      do 40 i = 1, m
	 xm(1,i) = 0.d0
	 xm(m,i) = 0.d0
         xm(i,1) = 0.d0
         xm(i,m) = 0.d0
   40 continue

      do 60 l = 2, m1
	 kt = (l-1)*m
	 do 50 i = 2, m1
	    k = kt + i
            if (kchoic.ne.4 .and. kchoic.ne.5 .and. kchoic.ne.6) then
	       w(k) = (8*xm(i,l)  - xm(i-1,l) - xm(i+1,l)
     &              - xm(i-1,l-1) - xm(i,l-1) - xm(i+1,l-1)
     &              - xm(i-1,l+1) - xm(i,l+1) - xm(i+1,l+1))/3
            end if
	    if (kchoic.ne.1 .and. kchoic.ne.6 .and. kchoic.ne.7) then
	       w(k) = w(k) + c*((-4*xm(i-1,l) + 4*xm(i+1,l)-xm(i-1,l-1)
     &              + xm(i+1,l-1) - xm(i-1,l+1) + xm(i+1,l+1))/(12*m1))
     &              + d*((-xm(i-1,l-1) - 4*xm(i,l-1) - xm(i+1,l-1) 
     &              + xm(i-1,l+1) + 4*xm(i,l+1) + xm(i+1,l+1))/(12*m1))
	    end if
	    if (kchoic.ne.1 .and. kchoic.ne.2 .and. kchoic.ne.4) then
	       w(k) = w(k) + f*((16*xm(i,l) + 4*xm(i-1,l) + 4*xm(i+1,l)
     &              + xm(i-1,l-1) + 4*xm(i,l-1) + xm(i+1,l-1)
     &              + xm(i-1,l+1) + 4*xm(i,l+1) + xm(i+1,l+1))/(36*m2))
	    end if
   50	 continue
   60 continue

      call bpxin(j, w, ichoic)

      return
      end
c=======================================================================

      subroutine matvc3(j, x, w, ichoic, kchoic, c, d, f)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes matrix_vector multiplication(W=B*x) in the
c     structured grid, where B is defined as follows:     
c         kchoic = 1 : B = A;
c         kchoic = 2 : B = A + cA_x + dA_y;
c         kchoic = 3 : B = A + cA_x + dA_y + fM;
c         kchoic = 4 : B = cA_x + dA_y;
c         kchoic = 5 : B = cA_x + dA_y + fM;
c         kchoic = 6 : B = fM;
c         kchoic = 7 : B = A + fM;
c     where A is the stiffness matrix, A_x & A_y are lower order matrices, and
c     M is the mass matix.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension x((2**j+1)**2), w((2**j+1)**2), xm(1025, 1025)

      m  = 2**j + 1
      m1 = m - 1
      m2 = m1*m1
c     h = 1.d0/m1

      do 10 i = 1, m*m
	 w(i) = 0.d0
   10 continue

c*    Matrix representation xmat of the vector x.

      do 30 l = 2, m1
	 kt = (l-1)*m 
	 do 20 i = 2, m1
	    k = kt + i
	    xm(i,l) = x(k)
cc	    xm(i,l) = 0.d0
   20    continue
   30 continue

cc      xm(3,3) = 1.d0

      do 40 i = 1, m
	 xm(1,i) = 0.d0
	 xm(m,i) = 0.d0
         xm(i,1) = 0.d0
         xm(i,m) = 0.d0
   40 continue

      ch  = c/m1
      dh  = d/m1
      fhh = f/m2

      if (c .gt. 0.d0) then
	 chm = - ch
	 ch0 =   ch
	 chp = 0.d0
      else
	 chm = 0.d0
	 ch0 = - ch
	 chp =   ch
      end if

      if (d .gt. 0.d0) then
	 dhm = - dh
	 dh0 =   dh
	 dhp = 0.d0
      else
	 dhm = 0.d0
	 dh0 = - dh
	 dhp =   dh
      end if

      do 60 l = 2, m1
	 kt = (l-1)*m
	 do 50 i = 2, m1
	    k = kt + i
            if (kchoic.ne.4 .and. kchoic.ne.5 .and. kchoic.ne.6) then
	       w(k) = 4*xm(i,l)  - xm(i-1,l) - xm(i+1,l)
     &              - xm(i,l-1) - xm(i,l+1)
            end if
	    if (kchoic.ne.1 .and. kchoic.ne.6 .and. kchoic.ne.7) then
	       w(k) = w(k) + chm*xm(i-1,l) + ch0*xm(i,l) + chp*xm(i+1,l)
     &              + dhm*xm(i,l-1) + dh0*xm(i,l) + dhp*xm(i,l+1)
	    end if
	    if (kchoic.ne.1 .and. kchoic.ne.2 .and. kchoic.ne.4) then
	       w(k) = w(k) + fhh*xm(i,l)
	    end if
   50	 continue
   60 continue

      call bpxin(j, w, ichoic)

cc      print*, (i, '*', 1.d-3*w(i), i=1,25)
cc      stop

      return
      end
c=======================================================================

      subroutine grflab(ls, n, x, icdomn)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine displays the various 3D views of the computed solution 
c     using a MATLAB tool, 'plot3'. In addition, if the given(original) mesh is
c     uniform on the unit square, then better graphic tools such as 'mesh' and 
c     'surf' in MATLAB can be used to display the computed solution beautifully.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension x(n)

c      print*,' Solution at each node 1: ', (i,'*',x(i),i=1,n)

      k  = 2**ls
      k1 = k - 1
      k2 = k + 1
      kk = k2*k2
      h  = 1.d0/k

      if (icdomn .eq. 3 .or. icdomn .eq. 4) then
         open(unit = 14, file = 'symsol.m', status = 'unknown')
         write(14,*) '% Level =', ls
         write(14,*) 'clear'
         do 30 jj = 1, k2
            kt = (jj-1)*k2
            do 20 ii = 1, k2
               kl = kt + ii
               write(14,*) 'zz(', jj, ',', ii, ') = ', x(kl), ';'
 20         continue
 30      continue
         do 40 jj = 1, k2
            write(14,*) 'xx(', jj, ') = ', (jj-1)*h, ';'
            write(14,*) 'yy(', jj, ') = ', (jj-1)*h, ';'
 40      continue
         write(14,*) 'sslgrf'
      end if
           
      close (14)

      return
      end
c=======================================================================

      subroutine exctxb(j, n, exsoln, b, icdomn, c, d, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     To supply an exact FE solution and a correspoding load vector in order
c     to test the algorithms. (This is a temporary process.)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension exsoln(n), b(n)

      kd = 2**j
      kr = kd - 1
      ke = kd + 1
      h  = 1.d0/kd

c     Find an exact FE solution.
      do 20 jj = 1, ke
         kt = (jj-1)*ke
         y = (jj-1)*h
         do 10 ii = 1, ke
            k = kt + ii
            x = (ii-1)*h
            exsoln(k) = extsol(x, y, icdomn)
 10      continue
 20   continue

      call bpxin(j, exsoln, icdomn)

c     Call MATVC2 to compute b = (A_0 + cA_0x + dA_0y + fM_0)*exsoln
      call matvc2(j, exsoln, b, icdomn, 3, c, d, f)

      return
      end
c=======================================================================

      subroutine error0(j, n, kk, xv, niter, cpu, icdomn, c, d, f,  
     &                  exsoln, ui, errloc, icerr, itbl, dtbl)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     In the structured grid of unit square, it computes (discrete) maximum
c     norm, L^2 norm, H^1 norm, and H^1 seminorm error estimates 
c      icerr = 1: between computed solution(u_h) and the exact finite element
c                 solution(u_e) of the discretized problem; or
c      icerr = 2: between computed solution(u_h) and interpolation(u_I) of
c                 continuous exact solution at grid points, if exact solution 
c                 is known. (User must describe the exact solution, u(x,y).
c                 (See function EXTSOL.)
c      icerr = 3: same as 2, but involved with stopping criteria, tol
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension xv(n), ui(n), errloc(n), exsoln(n)
      dimension itbl(15,2,20), dtbl(15,6,20)
      logical   done1, done2

      tol = 1.d-10
      kd = 2**j
      kr = kd - 1
      ke = kd + 1
      h     = 1.d0/kd
      errmx = 0.d0
      max   = 0
      xmax  = 0.d0
      ymax  = 0.d0
      ersum = 0.d0

      do 5 i = 1, n
	 ui(i)     = 0.d0
	 errloc(i) = 0.d0
 5    continue

      if (icerr .eq. 2 .or. icerr .eq. 3) then
c        Computing maximum error of u_I - u_h, where u is the exact solution,
c        u_I its interpolation, and u_h the computed solution.
         do 20 jj = 1, ke
            kt = (jj-1)*ke
            y = (jj-1)*h
            do 10 ii = 1, ke
               k = kt + ii
               x = (ii-1)*h
c               exsoln(k) = extsol(x, y, icdomn)
               errloc(k) = exsoln(k) - xv(k)
               errabs = dabs(errloc(k))
               if (errabs .gt. errmx) then
                  errmx = errabs
                  max = k
                  xmax = x
                  ymax = y
               end if
 10         continue
 20      continue
      end if

c      print*,' Exact Solution at each node: ', (i,'*',exsoln(i),i=1,n)

c     Computing L^2 error estimate of u_I - u_h.
      do 30 i = 1, n
         ersum = ersum + errloc(i)*errloc(i)
 30   continue
      erl22t = h*h*ersum
      erl2t  = dsqrt(erl22t)

c     Computing L^2 norm error estimate of u_I-u_h, sqrt(errloc^t*M*errloc).
c     Call MATVC2 to compute ui = M_0*errloc
      call matvc2(j, errloc, ui, icdomn, 6, 0.d0, 0.d0, 1.d0)
      erl22 = dot(errloc, ui, n)
      erl2  = dsqrt(erl22)

c     Computing H^1 seminorm error estimate of u_I-u_h, sqrt(errloc^t*A*errloc).
c     Call MATVC2 to compute ui = A_0*errloc
      call matvc2(j, errloc, ui, icdomn, 1, 0.d0, 0.d0, 0.d0)
      erh1s2 = dot(errloc, ui, n)
      erh1s  = dsqrt(erh1s2)

c     Computing H^1 error estimate of u_I - u_h.
      erh1 = dsqrt(erl22 + erh1s2)

cc      print 40, xmax, ymax, errmx
cc 40   format (/' Maximum error ( ||u_I - u_h||_infinity ) at (',f8.5,
cc     &        ',',f8.5,'): ',d12.5)   
cc      print 50, xmax, ymax, exsoln(max)
cc 50   format (/' Exact solution at (',f8.5,',',f8.5,'): ',d12.5) 	
cc      print 60, xmax, ymax, xv(max)
cc 60   format (/' Computed solution at (',f8.5,',',f8.5,'): ',d12.5) 

      ipb = 1

      if (ipb .eq. 1) then
         itbl(kk, 1, 1) = kk
         itbl(kk, 2, 1) = niter
         dtbl(kk, 1, 1) = errmx
         dtbl(kk, 2, 1) = erl2t
         dtbl(kk, 3, 1) = erl2
         dtbl(kk, 4, 1) = erh1s
         dtbl(kk, 5, 1) = erh1
         dtbl(kk, 6, 1) = cpu

      print* 
      print*,' level  niter   maximum    L^2(t)      L^2      H^1 semi
     &     H^1        CPU'
      print*,' --------------------------------------------------------
     &--------------------'
      do 70 i = 1, kk
         print 190, i, itbl(i,2,1),dtbl(i,1,1),dtbl(i,2,1), 
     &   dtbl(i,3,1),dtbl(i,4,1),dtbl(i,5,1),dtbl(i,6,1)
 70   continue
      print*, ' ----------------------------------------
     &------------------------------------'
c      print*
      go to 200

      end if

      if (ipb .eq. 1) niter = 1

      itbl(kk, 1, niter) = kk
      itbl(kk, 2, niter) = niter
      dtbl(kk, 1, niter) = errmx
      dtbl(kk, 2, niter) = erl2t
      dtbl(kk, 3, niter) = erl2
      dtbl(kk, 4, niter) = erh1s
      dtbl(kk, 5, niter) = erh1
      dtbl(kk, 6, niter) = cpu

80    if (ipb .ne. 1 .and. niter .gt. 10) go to 110

      print* 
      print*,' level  niter   maximum    L^2(t)      L^2      H^1 semi
     &     H^1        CPU'
      print*,' --------------------------------------------------------
     &--------------------'
      do 100 i = 1, kk
         print 190, i, itbl(i,2,niter),dtbl(i,1,niter),dtbl(i,2,niter), 
     &   dtbl(i,3,niter),dtbl(i,4,niter),dtbl(i,5,niter),dtbl(i,6,niter)
 100  continue
      print*, ' ----------------------------------------
     &------------------------------------'
c      print*
      go to 200

 110  print* 
      print*,' level  niter   maximum    L^2(t)      L^2      H^1 semi
     &     H^1        CPU'
      print*,' --------------------------------------------------------
     &--------------------'
      do 120 i = 1, kk
         contrl = itbl(i,2,niter) - 10.5d0
         print 180, i, contrl,dtbl(i,1,niter),dtbl(i,2,niter), 
     &   dtbl(i,3,niter),dtbl(i,4,niter),dtbl(i,5,niter),dtbl(i,6,niter)
 120  continue
      print*, ' ----------------------------------------
     &------------------------------------'
c      print*

 180  format(i5, f7.1, e13.4, e11.4, e11.4, e11.4, e11.4, f9.3)
 190  format(i5, i7, e13.4, e11.4, e11.4, e11.4, e11.4, f9.3)

 200  return
      end
c=======================================================================

      double precision function extsol(x, y, ichoic)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Exact solution of the PDE if it is known.  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
 
      pi = dacos(-1.d0)

      if (ichoic .eq. 1) then
c*       Example 1
c	 extsol = (1.d0-x*x-y*y)*x*y

c*       Example 2
c	 extsol = dsin(5*x*y*pi) - dexp(y)

c*       Example 3
c	 extsol = - dsin(x)*(dcos(y) - 1.d0 + 2*y/pi)

c*       Example 4
	 extsol = 5*(x - x**3 - x*y*y)*dsin(3*pi*y)

      else if (ichoic .eq. 3 .or. ichoic .eq. 4) then
         iexmp = 1
         if (iexmp .eq. 2) go to 50

c*       Example 1
	 extsol = 100*(x-x*x)*(y-y*y)

         go to 100

c*       Example 2
c	 extsol = dsin(5*x*y*pi) - dexp(y)

c*       Example 3
c	 extsol = - dsin(pi*x)*(dcos(pi*y) - 1.d0 + 2*y)

50       continue

c*       Example 4
	 extsol =  10*dsin(4*pi*x)*(dexp(y) - dexp(y*y))

c*       Example 5
c	 extsol = 5*(x - x**3 - x*y*y)*dsin(3*pi*y)

      end if

 100  return
      end
c=======================================================================

      double precision function psi(xa, hx, ya, hy, x, y, jchoic)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     psi(j) - 4 element basis funtions on a square with vertices
c             (xa,ya), (xa+hx,ya), (xa,ya+hy), (xa+hx,ya+hy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
 
      if (jchoic .eq. 1) then
	 psi = (x-xa)*(y-ya)/(hx*hy)
      else if (jchoic .eq. 2) then
	 psi = ((x-xa)/hx)*(1.d0 - (y-ya)/hy)
      else if (jchoic .eq. 3) then
	 psi = (1.d0 - (x-xa)/hx)*((y-ya)/hy)
      else if (jchoic .eq. 4) then
	 psi = (1.d0 - (x-xa)/hx)*(1.d0 - (y-ya)/hy)
      end if

      return
      end
c=======================================================================

      double precision function phi(x, y, jchoic)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Element basis functions on the reference triangle with vertices
c     (0., 0.), (1., 0.), (0., 1.).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
 
      if (jchoic .eq. 1) then
	 phi = 1.d0 - x - y
      else if (jchoic .eq. 2) then
	 phi = x
      else 
	 phi = y
      end if

      return
      end
c=======================================================================

      double precision function hxy(x,y)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This function should be defined before running the program, 
c     where hxy is the Dirichlet boundary condition:
c          { -Laplacian(U) = fxy  in Omega                   }
c          {             U = hxy  on the Dirichlet boundary  }
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)

      pi = dacos(-1.d0)

c*    Example 1
      hxy = 0.d0

c*    Example 2
c      if (x .lt. 1.d-15) then
c	 hxy = - dexp(y)
c      else if (y .lt. 1.d-15) then
c	 hxy = -1.d0
c      else
c         hxy = dsin(5*x*y*pi) - dexp(y)
c      end if 

      return
      end
c=======================================================================

      double precision function fxy(x, y, ichoic, c, d, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This function should be defined before running the program,
c     where fxy is given in the Poisson equation:
c        { -Laplacian(U) + c*U_x + d*U_y + f*U = fxy       }
c        {             U = hxy  on the Dirichlet boundary  }
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)

      pi = dacos(-1.d0)

      if (ichoic .eq. 1) then
c*       Example 1
c	 fxy = 12*x*y + c*(1 - 3*x*x - y*y)*y + d*(1 - 3*y*y - x*x)*x
c     &                + f*(1 - x*x - y*y)*x*y

c*       Example 2
c	 fxy = (25*pi*pi*dsin(5*x*y*pi))*(x*x+y*y) + dexp(y)

c*       Example 4      
	 ex = x - x**3 - x*y*y 
	 ed = 1 - 3*x*X - y*y 
         sn = dsin(3*pi*y)
         cs = dcos(3*pi*y)
	 fxy = (40*x + 45*pi*pi*ex + c*5*ed - 10*d*x*y + 5*f*ex)*sn
     &         + (60*pi*x*y + 15*d*pi*ex)*cs

      else if (ichoic .eq. 3 .or. ichoic .eq. 4) then
         iexmp = 1
         if (iexmp .eq. 2) go to 50

c*       Example 1
	 fxy = 2*(x - x*x + y - y*y) 
     &       + c*(1 - 2*x)*(y - y*y) + d*(x - x*x)*(1 - 2*y)
     &       + f*(x - x*x)*(y - y*y)
	 fxy = 100*fxy

         go to 100

c*       Example 2
c	 fxy = (25*pi*pi*dsin(5*x*y*pi))*(x*x+y*y) + dexp(y)

c*       Example 3
c	 fxy = - pi*pi*dsin(pi*x)*(2*dcos(pi*y) - 1.d0 + 2*y)

  50     continue

c*       Example 4
	 p2 = pi*pi
	 sx = dsin(4*pi*x)
	 cx = dcos(4*pi*x)
	 ey = dexp(y)
	 eys = dexp(y*y)
	 ey1 = ey - eys
	 ey2 = ey - 2*y*eys
	 ey3 = ey - (2 + 4*y*y)*eys
	 fxy = (16*p2*ey1-ey3)*sx + c*4*pi*ey1*cx + d*ey2*sx + f*ey1*sx
	 fxy = 10*fxy

c        Example 5
c	 ex = x - x**3 - x*y*y 
c	 ed = 1 - 3*x*X - y*y 
c        sn = dsin(3*pi*y)
c        cs = dcos(3*pi*y)
c	 fxy = (40*x + 45*pi*pi*ex + c*5*ed - 10*d*x*y + 5*f*ex)*sn
c     &         + (60*pi*x*y + 15*d*pi*ex)*cs
c
      end if

 100  return
      end
c=======================================================================
c     END of PROGRAM                                                   c
c=======================================================================
