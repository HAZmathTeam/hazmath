
      program nsauxd

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    <Considering three domains and two preconditioners for each domain>
c
c     NSAUXD uses linear finite elements on quasi-uniform triangle mesh to solve
c     the nonsymmetric/indefinite elliptic Dirichlet boundary value problem by 
c     introducing an auxiliary(structured) space:
c         {  -Laplacian(U) + cU_x + dU_y + fU = fxy   in the domain      }
c         {                                 U = hxy   on the boundary    },
c     where c, d, and f are constants.
c
c     First, the user should describe the functions f(x,y), h(x,y) in the
c     function subprograms 'fxy', 'hxy'.
c
c     Preconditioned conjugate gradient method(PCGM) is applied to solve the 
c     linear system A_0*x = b of corresponding SPD operator, where only nonzero 
c     elements in the lower triangular part of the stiffness matrix A_0 are 
c     stored in the one dimensional array a. As a preconditioner in the PCG 
c     algorithm, a combination of the multilevel nodal basis preconditioner(BPX)
c     B on the structured grid and the matrix T representing the interpolation 
c     from structured mesh onto unstructured mesh is used, namely, R + TBT^t, 
c     where R is a smoother for the unstructured space.
c
c     During execution of the program the user specifies the level(integer
c     value) lvlaux ( 0 < lvlaux < 11 ) of the  so that h = 1/(2^lvlaux), where
c     h is the length of each side of squares (parallel to x,y-axes) in the 
c     uniform mesh of the unit square.
c
c     leng   = number of nodes in the finest structured grid
c              (if lvlaux = 8,  set leng =    66050)
c              (if lvlaux = 9,  set leng =   263170)
c              (if lvlaux = 10, set leng =  1050630)
c              (if lvlaux = 11, set leng =  4198410)
c     lim    = maximum number of nodes & elements in the unstructured mesh
c              (if lvluns = 8, set lim =  66050)
c              (if lvluns = 9, set lim = 263170)
c     lvlaux = level of the auxiliary(structured) space
c     lvluns = level of the unstructured space
c     nodeel = global node numbers for all elements in the unstructured mesh
c     conode = x,y coordinates of all nodes in the unstructured mesh
c     nodedb = node numbers on the Dirichlet boundary in the unstructured mesh
c     noelmt = number of elements in the unstructured mesh
c     nonode = number of nodes in the unstructured mesh
c     nodb   = number of Dirichlet boundary nodes in the unstructured mesh
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     icdomn : choice of domains
c              1 = the quarter circle.
c              3 = the unit square with zeros only on the boundary.
c              4 = the unit square with zeros near/on the boundary.
c     icsolv : choice of solvers of the auxiliary(structured) space
c              1 = considering the BPX preconditioner on the structured mesh; 
c		             B_h = G + T*B*T^t
c              2 = considering an exact solver on the structured mesh, i.e.,
c                  replacing B by A_0^(-1) in the icslov = 1;
c		             B_h = G + T*A_0^(-1)*T^t
c     icalgm : choice of algorithms
c              2 = Algorithm 2 of the paper by H. Kim and J. Xu
c              3 = Algorithm 3 of the paper by H. Kim and J. Xu
c     where
c        A_h : stiffness matrix of the unstructured space
c              (linear elements on triangles)
c        A_0 : stiffness matrix of the uniform mesh in the structured space
c              (bilinear elements on squares)
c        G   : Gauss-Seidel smoother
c        B   : BPX preconditioner
c        B_h : BPX type preconditioner
c        T   : matrix representation of the interpolation from structured mesh
c              to unstructured mesh
c        I   : identity transformation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (leng = 263170, lim = 263170)

      implicit double precision (a-h, o-z)
      dimension nodeel(3,2*lim), conode(2,lim), jcn(15*lim), bwkt(leng), 
     &          nodedb(2500), jre(lim), a(15*lim), b(leng),soludb(2500),
     &          x(leng), ax(15*lim), ay(15*lim), jdiag(lim), am(15*lim),
     &          wk(leng), bwk(leng), xwk(leng), wk1(2*leng),wk2(2*leng),
     &          exsoln(lim), itbl(10,2), dtbl(10,6), bwk1(2*leng), 
     &          jd(11)
      real      time(2), etime, dtime
      logical   done1, done2

      open(unit = 15, file = 'exfesl', status = 'unknown')

 10   continue

      done1 = .false.
      done2 = .false.

      t11 = etime(time)
      t12 = dtime(time)

c     Initial set up (various choices).
c     MG with coarsest level = 1:  c = -24,  d = 26 (-20, 22)
c     MG with coarsest level = 2:  c = -37,  d = 39 (-32, 34)
c     MG with coarsest level = 3:  c = -51,  d = 53 
c     MG with coarsest level = 4:  c = -95,  d = 100 
c     MG with coarsest level = 5:  c = -100, d = 200
c     MGNS with coarsest level = 1:  c = -12, d = 14 (-11, 13)
c     MGNS with coarsest level = 2:  c = -20, d = 22 (-19, 21)
c     MGNS with coarsest level = 3:  c = -30, d = 32
c     MGNS with coarsest level = 4:  c = 
c     MGNS with coarsest level = 5:  c = 

      c      = -5.d0
      d      = -8.d0
      f      =   1.d0

      iexct  = 1               !Choice of the coarsest level, 1, 2.
      icalgm = 2               !Choice of algorithms, 2, 3.
      icsolv = 1               !Choice of solvers of the auxiliary space, 1, 2.
      icdomn = 1               !Choice of domains, 1, 3, 4.
      nstep  = 1               ! 1, 2, 3, 4

      ipde   = 3               !MG parameters - start
      nit    = 150             !        .
      ivbs   = 1               !        .
      jc     = 1               !        .
      nsmthg = 1               !        .
      ifbs1  = 1               !        .
      ifbs2  = 1               !        .
      isolv  = 0               !MG parameters - end


      print 20
 20   format (/' Enter the level ( 2 < lvluns < 10) of the unstructured' 
     &        /' grids so that h =(approx.) 1/(2^lvluns):')
      read*, lvluns

      if (lvluns.le.2 .or. lvluns.ge.10) then
         print*, ' Invalid choice. Reenter!!'
	 go to 10
      end if

      lvlaux = lvluns
      k  = 2**lvlaux 
      k1 = k - 1
      k2 = k + 1
      n  = k2*k2

      jd(1) = 1
      do 31 i = 1, 10
         jd(i+1) = jd(i) + (2**i+1)**2
 31   continue

      print*, ' Please wait a while if the level is large.'
      print*

      if (icdomn .eq. 1) then
	 call meshus(lvluns,nodeel,conode,nodedb,noelmt,nonode,nodb,lim)
      else if (icdomn .eq. 3 .or. icdomn .eq. 4) then
c         c = 0.d0
c         d = 0.d0
         call meshst(lvluns,nodeel,conode,nodedb,noelmt,nonode,nodb,lim)
      end if

c     Call stiffm to compute the stiffness matrix A_h, the lower order matrices
c     A_x and A_y, and the mass matrix M_h of the problem without imposing any 
c     boundary condition (i.e., A = A_h + cA_x + dA_y + fM_h).
      call stiffm(nodeel, conode, noelmt, nonode, a, jcn, jre, jdiag,
     &            ax, ay, am, c, d, f)

c     Call loadbu to compute the load vector b of the linear system A*u = b in 
c     the unstructured mesh.
      call loadbu(nodeel, conode, noelmt, nonode, b, icdomn, c, d, f)

c     Call stifdb to modify matrices A_h, A_x, A_y, M_h and the vector b so that
c     they satisfy the Dirichlet boundary condition.
      if (nodb .ne. 0) then
         call stifdb(a, ax, ay, am, jcn, jre, jdiag, b, conode, nodedb,
     &               soludb, nonode, nodb, c, d, f)
      end if

c     To supply an exact FE solution and a correspoding load vector in order
c     to test the algorithms. (This is a temporary process.)
      call exctxb(nonode, a, ax, ay, am, jcn, jre, conode, nodedb, nodb, 
     &            soludb, exsoln, b, icdomn, c, d, f)

      t21 = etime(time)
      t22 = dtime(time)

c*****Step 1. Solve nonsymmetric problem at auxiliary(structured) space.********

      icb = 1                                       !Choice of load vector b.
      if (icb .eq. 1) then
c        To compute the load vector bwk = (T^t)*b by transpose of T.
         call actntt(conode, nonode, lvlaux, b, bwk, icdomn)
      else if (icb .eq. 2) then
c        To compute the load vector bwk of the linear system
c        (A_0 + cA_0x + dA_0y + fM_0)*u = bwk in the structured mesh directly.
         call loadbs(lvlaux, bwk, icdomn, c, d, f)
      end if 

      icaux = 1              !Choice of solvers in the auxiliary space.

ccc   do 36 nit = 1, 50
ccc      tt11 = etime(time)
ccc      tt21 = dtime(time)
         if (icaux .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by multigrid method in the structured mesh.
            call mgvbs(lvlaux, n, icdomn, ipde, nit, c, d, f, bwk, xwk,
     &           isolv, ivbs,jc, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
         else if (icaux .eq. 2) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by two-grid method in the structured mesh.
            call tgrid1(lvlaux, n, icdomn, iexct, icalgm, nit, c, d, f, 
     &                                    bwk, xwk, 0, jd, wk, wk1, wk2)
         else if (icaux .eq. 3) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by two-grid method in the structured mesh.
            call tgrid3(lvlaux, n, icdomn, iexct, icalgm, nit, c, d, f, 
     &                                    bwk, xwk, 0, jd, wk, wk1, wk2)
         else if (icaux .eq. 4) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Gauss-Seidel method in the structured mesh.
            call gssor(lvlaux, n, bwk, xwk, icdomn,3,3,c,d,f,wk,nit)
         else if (icaux .eq. 5) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by full multigrid method in the structured mesh.
            call mgvf(lvlaux, n, icdomn, iexct, 3, nit, c, d, f, 
     &                          bwk, xwk, jd, wk, wk1, wk2, bwk1)
         else if (icaux .eq. 6) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by full twogrid method in the structured mesh.
            call tgf(lvlaux, n, icdomn, iexct, icalgm, nit, c, d, f,
     &                               bwk, xwk, jd, wk, wk1, wk2, bwk1)
         else if (icaux .eq. 7) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = bwk
c           by Bi-CGSTAB method in the structured mesh.
            call bcgst(lvlaux, n, bwk, xwk, 0.1d0, nit,
     &                  icdomn, 3, c, d, f) 
         end if

ccc      tt12 = etime(time)
ccc      tt22 = dtime(time)
ccc      call error0(lvlaux, n, xwk, nit, tt22, icdomn, c, d, f, exsoln,
ccc  &                  wk, bwkt, 3, itbl, dtbl, done1, done2)
ccc      if (done1 .and. done2)  go to 37
ccc36   continue

c     Call actnt to compute x = T*xwk
      call actnt(conode, nonode, lvlaux, xwk, x, icdomn)

      icount = 0
 40   icount = icount + 1

c*****Step 2. Solve the symmetric part at original(unstructured) space.*********

c     To compute wk = (A_h + cA_x + dA_y + fM)*x on the unstructured grid.
      call matvec(a, ax, ay, am, jcn, jre, nonode, x, wk, 3, c, d,f)

c     Compute bwkt = b - wk.
      do 50 i = 1, nonode
         bwkt(i) = b(i) - wk(i)
 50   continue

c     Call PCGM to solve the SPD linear system A_h*x = b (b = bwk) on the 
c     unstructured grid.
      call pcgm(lvluns, lvlaux, nonode, a, ax, ay, am, jcn, jre, jdiag, 
     &          bwkt, xwk, conode, niter, icdomn, icsolv, c, d, f)

c     Update the solutions.
      do 60 i = 1, nonode
         x(i) = x(i) + xwk(i)
 60   continue

      if (icalgm .eq. 2 .and. icount .eq. nstep) go to 140

c*****Step 3. Solve the nonsymmetric problem at auxiliary(structured) space.****

c     Compute xwk = (A_h + cA_x + dA_y + fM)*x on unstructured grid.
      call matvec(a, ax, ay, am, jcn, jre, nonode, x, xwk, 3, 
     &            c, d, f)

c     Call actntt to compute wk = (T^t)*xwk.
      call actntt(conode, nonode, lvlaux, xwk, wk, icdomn)

c     Compute bwkt = bwk - wk.
      do 90 i = 1, n
         bwkt(i) = bwk(i) - wk(i)
 90   continue

      if (icaux .eq. 1) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*xwk = bwk
c        by multigrid method in the structured mesh.
         call mgvbs(lvlaux, n, icdomn, ipde, nit, c, d, f, bwkt, xwk,
     &           isolv, ivbs,jc, nsmthg, ifbs1, ifbs2, jd, wk, wk1, wk2)
      else if (icaux .eq. 2) then
c        To solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*
c        xwk = bwk by two-grid method in the structured mesh.
         call tgrid1(lvlaux, n, icdomn, iexct, icalgm, nit, c, d, f, 
     &               bwkt, xwk, 0, jd, wk, wk1, wk2)
      end if

c     Call actnt to compute wk = T*xwk.
      call actnt(conode, nonode, lvlaux, xwk, wk, icdomn)

      do 130 i = 1, nonode
         x(i) = x(i) + wk(i)
 130  continue

      if (icalgm .eq. 2) then
         go to 40
      else if (icalgm .eq. 3) then
	 if (icount .eq. nstep) go to 140
         go to 40
      end if

c*****End of Algorithm.*********************************************************

 140  continue

c     Exact solution at the Dirichlet boundary.
      do 160 i = 1, nodb
	 x(nodedb(i)) = soludb(i)
 160  continue

 165  t41 = etime(time)
      t42 = dtime(time)

      print*, 'CPU time = ', t11, t21, t41
      print*, 'CPU time = ', t12, t22, t42
  
c     Presentation of result (computed solution).

 170  Print 180
 180  format(/' 1 Print out the exact finite element solution on the'
     &       /'   structured mesh of the unit square.' 
     &       /'   (Use icaux = 1 or 3 with a small stopping criterion.)'
     &       /' 2 Display the solution vector X with its coordinates.'
     &       /'   (various 3D views using MATLAB tool)'
     &       /' 3 Compute max.norm, L^2 norm, and H^1 (semi)norm errors'
     &	     /'   between computed solution(u_h) and the exact finite'
     &       /'   element solution(u_e) of the discretized problem.'
     &       /' 0 Quit.'
     &      //' Enter choice :')

      read*, iquery
      if (iquery .eq. 1) then
         write(15,185)  lvlaux, n, n/3
 185     format(i3, i8, i8)
         do 195 ii = 1, n/3
            write(15,190) exsoln(3*ii-2), exsoln(3*ii-1), exsoln(3*ii)
 190        format(e23.15, e23.15, e23.15)
 195     continue
         do 197 ii = n/3*3 + 1, n
            write(15,190) exsoln(ii)
 197     continue
	 close (15)
      else if (iquery .eq. 2) then
         if (icalgm .eq. 0)   nonode = n
c        print*,' Solution at each node : ', (i,'*',x(i),i=1, nonode)
c        To display various 3D views of the computed solution using MATLAB tool.
         call grflab(lvluns, lvlaux, conode, nonode, x,icdomn,icalgm)
c         call grflab(lvluns,lvlaux,conode,nonode, exsoln,icdomn,icalgm)
      else if (iquery .eq. 3) then
	 cpu = t22 + t42
c        Compute maximum, L^2, H^1 seminorm, and H^1 norm errors.
         call errors(lvluns, a, ax, ay, am, jcn, jre, nonode, x,
     &        exsoln, niter, cpu, icdomn, c, d, f, itbl, dtbl,wk,xwk)
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
c            close(unit = 15)
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

      subroutine meshus(level, nodeel, conode, nodedb, noelmt,
     &                  nonode, nodb, lim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine reads all information about the triangulation of
c     unstructured mesh from existing data files.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension nodeel(3,2*lim), conode(2,lim), nodedb(2500)

      if (level .eq. 3) then 
         open(unit = 10, file = 'cord3', status = 'old')
         open(unit = 11, file = 'db3', status = 'old')
      else if (level .eq. 4) then 
         open(unit = 10, file = 'cord4', status = 'old')
         open(unit = 11, file = 'db4', status = 'old')
      else if (level .eq. 5) then 
         open(unit = 10, file = 'cord5', status = 'old')
         open(unit = 11, file = 'db5', status = 'old')
      else if (level .eq. 6) then 
         open(unit = 10, file = 'cord6', status = 'old')
         open(unit = 11, file = 'db6', status = 'old')
      else if (level .eq. 7) then 
         open(unit = 10, file = 'cord7', status = 'old')
         open(unit = 11, file = 'db7', status = 'old')
      else if (level .eq. 8) then 
         open(unit = 10, file = 'cord8', status = 'old')
         open(unit = 11, file = 'db8', status = 'old')
      else if (level .eq. 9) then 
         open(unit = 10, file = 'cord9', status = 'old')
         open(unit = 11, file = 'db9', status = 'old')
      end if

      read(10,*) noelmt, nonode

      eps  = .5d-4
c     eps1 = 1.d-5

      do 20 i = 1, nonode
	 read(10,*) conode(1,i), conode(2,i)
c	 if (conode(1,i) .le. eps1)   conode(1,i) = 0.d0
c	 if (conode(2,i) .le. eps1)   conode(2,i) = 0.d0
         xco = conode(1,i)
	 yco = conode(2,i)
	 if (dabs(1.d0 - xco*xco - yco*yco) .lt. eps) then
	    conode(2,i) = dsqrt(1.d0 - xco*xco)
         end if
   20 continue

      do 30 i = 1, noelmt
      	 read(10,*) (nodeel(jj,i), jj = 1,3)
   30 continue
      
      read(11,*) nodb  
      do 50 i = 1, nodb
      	 read(11,*) nodedb(i)
   50 continue

      close (10)
      close (11)

      return
      end
c=======================================================================

      subroutine meshst(level, nodeel, conode, nodedb, noelmt,
     &                  nonode, nodb, lim)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine generates all information about the triangulation of
c     structured(uniform) mesh on the unit square.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension nodeel(3,2*lim), conode(2,lim), nodedb(2500)

      m0     = 2**level - 1
      mh     = m0 + 1
      m      = m0 + 2
      noelmt = 2*mh*mh
      nonode = m*m
      nodb   = 4*mh

      do 30 jj = 1,m
         jt = (jj-1)*m
         do 25 i = 1,m
	    l = jt + i
            conode(1,l) = (i-1)/dble(real(mh))
            conode(2,l) = (jj-1)/dble(real(mh))  
   25    continue  
   30 continue  

      do 40 jj = 1,mh
         jt  = 2*(jj-1)*mh
         jt1 = jt + mh
	 jtv = (jj-1)*m
         do 35 i = 1,mh
	    l  = jt + i
	    l1 = jt1 + i
            lv = jtv + i
            nodeel(1,l)  = lv 
            nodeel(2,l)  = lv + 1 
            nodeel(3,l)  = lv + 1 + m
            nodeel(1,l1) = lv 
            nodeel(2,l1) = lv + 1 + m
            nodeel(3,l1) = lv + m 
   35    continue  
   40 continue  

      do 45 i = 1,m
         nodedb(i)   = i
         nodedb(i+m) = nonode - i + 1
   45 continue  

      do 50 i = 1,m0
         nodedb(2*m+i)    = m*(i+1)
         nodedb(2*m+i+m0) = m*i + 1
   50 continue  

      return
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

      subroutine stiffm(nodeel, conode, l, n, a, jcn, jre, jdiag, ax, 
     &                  ay, am, c, d, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes the stiffness matrix A, the lower order matrices
c     A_x and A_y, and the mass matrix M of linear elements on the given 
c     triangular mesh and stores the matrices in one dimensional arrays a, ax, 
c     ay, and am which are the row-ordered lists that contain only nonzero 
c     entries of the matrices A, A_x, A_y, and M where
c        A is the stiffness matrix;
c        A_x = (a_ij), a_ij = ((phi_j)_x, phi_i);
c        A_y = (a_ij), a_ij = ((phi_j)_y, phi_i);
c        M is the mass matrix.
c     jcn(i)   - indicates that a(i), ax(i), ay(i), or am(i) is a entry in the
c                jcn(i)-th column of the matrices.
c     jdiag(k) - indicates the position of the k_th diagonal entry of A, A_x,
c                A_y, or M in the list a, ax, ay, or am.
c     jre(k)   - indicates the number of nonzero entries in the k_th row of the
c                matrix A, A_x, A_y, or M.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  a(15*n), nodeel(3,l), conode(2,n), jcn(15*n), jre(n),
     &           aij(3,3), co1(2), co2(2), co3(2), ax(15*n), ay(15*n),
     &           axij(3,3), ayij(3,3), amij(3,3), jdiag(n), am(15*n)
      logical done

      m = 15*n

      do 10 i = 1, m
	 a(i)    = 0.d0
	 ax(i)   = 0.d0
	 ay(i)   = 0.d0
	 am(i)   = 0.d0
         jcn(i)  = 0
   10 continue
         
      do 20 i = 1, n
         jre(i)   = 0
         jdiag(i) = 0
   20 continue

      do 120 i = 1, l
         do 30 j = 1, 2
	    co1(j) = conode(j, nodeel(1,i))
	    co2(j) = conode(j, nodeel(2,i))
	    co3(j) = conode(j, nodeel(3,i))
   30    continue

	 call elstif(co1, co2, co3, aij, axij, ayij, amij)

c*       Computation of a, ax, ay, am, jcn, and jre 

         do 110 i1 = 1, 3
	    k1   = nodeel(i1, i)
	    jchk = 15*(k1-1)
	    do 100 i2 = 1, 3
	       k2   = nodeel(i2, i)
	       jpos = jchk + 1
	       done = .false.

	       if (i2 .eq. 1 .and. jre(k1) .eq. 0) then
		  a(jpos)   = aij(i1,i2)
		  ax(jpos)  = axij(i1,i2)
		  ay(jpos)  = ayij(i1,i2)
		  am(jpos)  = amij(i1,i2)
		  jcn(jpos) = k2
		  jre(k1)   = 1
		  done      = .true.
	       end if

   80	       if (.not. done) then
    	          if (jcn(jpos) .eq. k2) then
		     a(jpos)  = a(jpos) + aij(i1,i2)
		     ax(jpos) = ax(jpos) + axij(i1,i2)
		     ay(jpos) = ay(jpos) + ayij(i1,i2)
		     am(jpos) = am(jpos) + amij(i1,i2)
		     done     = .true.
		  else if (jcn(jpos) .gt. k2) then
		     jre(k1) = jre(k1) + 1
		     do 90 i3 = jchk+jre(k1), jpos+1, -1
		        a(i3)   = a(i3-1)
		        ax(i3)  = ax(i3-1)
		        ay(i3)  = ay(i3-1)
		        am(i3)  = am(i3-1)
		        jcn(i3) = jcn(i3-1)
   90		     continue
	             a(jpos)   = aij(i1,i2)
	             ax(jpos)  = axij(i1,i2)
		     ay(jpos)  = ayij(i1,i2)
		     am(jpos)  = amij(i1,i2)
		     jcn(jpos) = k2
		     done      = .true.
	          else
		     jpos = jpos + 1
		  end if

    	          if (jpos .gt. jchk+jre(k1)) then
		     a(jpos)   = aij(i1,i2)
		     ax(jpos)  = axij(i1,i2)
		     ay(jpos)  = ayij(i1,i2)
		     am(jpos)  = amij(i1,i2)
		     jcn(jpos) = k2
		     jre(k1)   = jre(k1) + 1
		     done      = .true.
		  end if
                  go to 80
	       end if
  100	    continue
  110    continue
  120 continue

c      do 140 i = 0, n-1
c	 jchk = i*15
c	 do 130 ii = 1, jre(i+1)
c	    jpos     = jchk + ii
c	    ax(jpos) = c*ax(jpos)
c	    ay(jpos) = d*ay(jpos)
c	    am(jpos) = f*am(jpos)
c  130    continue
c  140 continue

c*    Finding the array jdiag from the list jcn.

      do 160 i = 1, n-1
	 jchk = 15*(i-1)
	 k    = 1
  150	 jpos = jchk + k
	 if (jcn(jpos) .eq. i) then
            jdiag(i) = jpos
	 else
	    k = k + 1
	    go to 150
	 end if
  160 continue
      jdiag(n) = 15*(n-1) + jre(n)

      return
      end
c=======================================================================

      subroutine elstif(co1, co2, co3, aij, axij, ayij, amij) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes the element stiffness matrix aij, the lower order
c     element matrices axij and ayij, and the element mass matrix amij of the 
c     triangle element consisting of three vertices, co1, co2, and co3.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  co1(2), co2(2), co3(2), aij(3,3), axij(3,3), ayij(3,3),
     &           amij(3,3)

      x2m1 = co2(1) - co1(1)
      x3m2 = co3(1) - co2(1)
      x3m1 = co3(1) - co1(1)
      y2m1 = co2(2) - co1(2)
      y3m2 = co3(2) - co2(2)
      y3m1 = co3(2) - co1(2)
      det  = x2m1*y3m1 - x3m1*y2m1
      if (det .gt. 0.d0) then
	 js = 1
      else
	 js  = - 1
	 det = - det
      end if

c*    Computation of aij, axij, ayij, and amij.

      rdet = 0.5d0/det
      aij(1,1) =   rdet*(y3m2*y3m2 + x3m2*x3m2)
      aij(2,2) =   rdet*(y3m1*y3m1 + x3m1*x3m1)
      aij(3,3) =   rdet*(y2m1*y2m1 + x2m1*x2m1)
      aij(2,1) = - rdet*(y3m1*y3m2 + x3m1*x3m2)
      aij(3,1) =   rdet*(y2m1*y3m2 + x2m1*x3m2)
      aij(3,2) = - rdet*(y2m1*y3m1 + x2m1*x3m1)
      aij(1,2) = aij(2,1)
      aij(1,3) = aij(3,1)
      aij(2,3) = aij(3,2)

      axij(1,1) = - js*y3m2/6
      axij(1,2) =   js*y3m1/6
      axij(1,3) = - js*y2m1/6
      ayij(1,1) =   js*x3m2/6
      ayij(1,2) = - js*x3m1/6
      ayij(1,3) =   js*x2m1/6

      do 10 k = 1, 3
         axij(2,k) =   axij(1,k)
         axij(3,k) =   axij(1,k)
         ayij(2,k) =   ayij(1,k)
         ayij(3,k) =   ayij(1,k)
   10 continue

      det24 = det/24
      det12 = 2*det24
      do 30 k = 1, 3
         do 20 kk = 1, 3
            amij(k, kk) = det24
            if (k .eq. kk) amij(k, kk) = det12
 20      continue
 30   continue

      return
      end
c=======================================================================

      subroutine loadbu(nodeel, conode, l, n, b, ichoic, c, d, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes the vector b = (b(j)) in the unstructured
c     mesh, where b(j) = integral_omega(f*phi_j).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  nodeel(3,l), conode(2,n), b(n)

      do 10 i = 1, n
	 b(i) = 0.d0
   10 continue

c*    Integration in Omega

      do 20 i = 1, l
         k1 = nodeel(1,i)
         k2 = nodeel(2,i)
         k3 = nodeel(3,i)
	 ai1 = conode(1,k1)
	 ai2 = conode(2,k1) 
	 x2m1 = conode(1,k2) - ai1
	 x3m1 = conode(1,k3) - ai1
	 y2m1 = conode(2,k2) - ai2
	 y3m1 = conode(2,k3) - ai2
	 det  = dabs(x2m1*y3m1 - x3m1*y2m1)
	 b123 = det*
     &          fxy(ai1+(x2m1+x3m1)/3,ai2+(y2m1+y3m1)/3,ichoic,c,d,f)/6
	 b(k1) = b(k1) + b123
	 b(k2) = b(k2) + b123
	 b(k3) = b(k3) + b123
   20 continue

      return
      end
c=======================================================================

      subroutine stifdb(a, ax, ay, am, jcn, jre, jdiag, b, conode,
     &                  nodedb, soludb, n, nodb, c, d, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine modifies the matrices A, A_x, A_y, M, and the vector b 
c     so that they satisfy the Dirichlet boundary condition.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  a(15*n), b(n), conode(2,n), jcn(15*n),jre(n),ax(15*n),
     &           nodedb(nodb), soludb(nodb), jdiag(n),ay(15*n),am(15*n)

c*    Call sortin to sort the values in the array nodedb into ascending 
c     order using an insertion sort.
      call sortin(nodedb, nodb)

c*    Exact solution at Dirichlet boundary.
      do 10 i = 1, nodb
	 ni        = nodedb(i)
	 soludb(i) = hxy(conode(1,ni), conode(2,ni))
   10 continue

      k = 1
      do 60 i = 1, n
         jchk = 15*(i-1)
	 if (k .gt. nodb) go to 30
	 nk = nodedb(k)

c*       When i is a Dirichlet boundary node;

    	 if (i .eq. nk) then
	    do 20 ii = 1, jre(i)
	       jpos     = jchk + ii
	       a(jpos)  = 0.d0
	       ax(jpos) = 0.d0
	       ay(jpos) = 0.d0
	       am(jpos) = 0.d0
   20       continue
	    a(jdiag(i)) = 1.d0
	    b(i)        = 0.d0
	    k           = k + 1
  	    go to 60
	 end if

c*       When i is not a Dirichlet boundary node;

   30    kk = 1
	 do 50 ii = 1, jre(i)
	    jpos = jchk + ii
   40	    if (kk .gt. nodb)   go to 60
	    kt = nodedb(kk)
	    if (jcn(jpos) .eq. kt) then
	       b(i) = b(i) - soludb(kk)*(a(jpos) + c*ax(jpos) 
     &                              + d*ay(jpos) + f*am(jpos))
	       a(jpos)  = 0.d0
	       ax(jpos) = 0.d0
	       ay(jpos) = 0.d0
	       am(jpos) = 0.d0
	       kk       = kk + 1
 	    else if (jcn(jpos) .gt. kt) then
	       kk = kk + 1
	       go to 40
   	    end if
   50    continue
   60 continue

      return
      end
c=======================================================================

      subroutine sortin(nx, ncount)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     These statements sort the values in the array x into ascending order using
c     an insertion sort. The variable count contains the number of valid data
c     values in x.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  nx(ncount)
      logical done

      do 20 j = 1, ncount - 1
         if (nx(j) .gt. nx(j+1)) then
	    done = .false.
	    k = j
   10       if (.not. done) then
	       nhold   = nx(k)
	       nx(k)   = nx(k+1)
	       nx(k+1) = nhold
	       if (k .eq. 1) then
		  done = .true.
	       else if (nx(k) .ge. nx(k-1)) then
		  done = .true.
	       else
		  k = k - 1
	       end if
	       go to 10
	    end if
	 end if	
   20 continue

      return
      end
c=======================================================================

      subroutine tgrid1(j, n, icdomn, iexct, icalgm, niter,c,d,f,b,x,
     &                                          icit, jd, wk, wk1, wk2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TGRID1 stands for Two-Grid method and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a Two-Grid method.
c          icit = 1: 1, 2, or 4 iterations
c          icit = 0: more than 1 iteration, depending on stopping criteria
c          icit = 2: same as icit = 0, but do not initialize the soln x
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim = 263170)

      implicit double precision (a-h, o-z)
      dimension  x(n), b(n), wk(n), xwk(2*lim), bwk(2*lim), xwkt(2*lim)
      dimension  wk1(2*n), wk2(2*n), jd(11)

      tol   = 1.d-40/(4**j)            !tol - stopping criteria for iterations
c      itmax = 1                      !max. iteration
      itmax = niter
      xnewn = 0.d0
      stopl = 1.d0
      nl    = jd(j) - jd(j-1)             !Dimension of the (j-1)th level.
      iexsol = 1 
      lvlcst = iexct                               !Coarsest level.  
      if (icit .eq. 1) itmax  = 1
c      if (icit .eq. 1 .and. j .eq. lvlcst+1) itmax  = 1

      do 10 i = 1, jd(j)-1
         xwk(i)  = 0.d0
         bwk(i)  = 0.d0
         xwkt(i) = 0.d0
 10   continue

      do 20 i = 1, n
         x(i)           = 0.d0
	 xwkt(jd(j)+i-1)= 0.d0
	 xwk(jd(j)+i-1) = 0.d0
	 bwk(jd(j)+i-1) = b(i)
 20   continue

      if (j .eq. lvlcst) then           !Solving exactly at the coarsest level
         if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c           Compute x = (A_0 + cA_x0 + dA_y0 + fM_0)^(-1)*b.
            call solexc(b, x, icdomn, lvlcst, c, d, f)
         else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c           by Bi-CGSTAB method in the structured mesh.
            nd = jd(lvlcst+1) - jd(lvlcst)
            call bcgst(lvlcst, nd, b, x, 0.1d0, nit,
     &                  icdomn, 3, c, d, f)
         else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c           by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
            nd = jd(lvlcst+1) - jd(lvlcst)
            call gssor(lvlcst, nd, b, x,
     &                 icdomn, 3, 3, c, d, f, wk, nit)
          end if
         go to 555
      end if

      do 150 it = 1, itmax
c********Step 1. Solve the nonsymmetric problem at coarser space.***************

         call restrn(j, bwk(jd(j)), bwk(jd(j-1)))   !Restriction to lower level.
         call bpxin(j-1, bwk(jd(j-1)), icdomn)

         if (j-1 .eq. lvlcst) then
            if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c              Solving exactly at the coarsest level.
               call solexc(bwk(jd(lvlcst)), xwk(jd(lvlcst)), icdomn, 
     &                     lvlcst, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c              by Bi-CGSTAB method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call bcgst(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    0.1d0, nit, icdomn, 3, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c              by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call gssor(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    icdomn, 3, 3, c, d, f, wk, nit)
            end if
            go to 40
         end if

c        Solve (A_0+cA_0x+dA_0y+fM_0)*xwk=bwk by two-grid method at lower level.
         call tgrid2(j-1, nl, icdomn, iexct, icalgm, 4, c, d, f, 
     &               bwk(jd(j-1)), xwk(jd(j-1)), 1, jd, wk, wk1, wk2)

c********Step 2. Solve the symmetric problem at finer space.********************
 40      call prolng(j,xwk(jd(j-1)),xwk(jd(j)))   !Prolongation to upper level.
         call bpxin(j, xwk(jd(j)), icdomn)

         do 50 i = 0, n-1
            xwkt(jd(j)+i) = x(i+1) + xwk(jd(j)+i)
 50      continue

         icwk = 1                               !icit = 1 is better somehow.
         if (icwk .eq. 1) then
c           To compute wk = (A_h + cA_x + dA_y + fM)*xwkt at level j.
            call matvc2(j, xwkt(jd(j)), wk, icdomn, 3, c, d, f)
         else if (icwk .eq. 2) then
c           To compute wk = (cA_x + dA_y + fM)*xwkt at level j.
            call matvc2(j, xwkt(jd(j)), wk, icdomn, 5, c, d, f)
         end if

c        Compute bwk = b - wk = b - (A_h + cA_x + dA_y + fM)*x.
         do 60 i = 1, n
            bwk(jd(j)+i-1) = b(i) - wk(i)
 60      continue

c        To solve the SPD linear system A_h*xwk = bwk.
c         call mgv(j, n, icdomn, 1, 1, nit0, 0.d0, 0.d0, 0.d0, 
c     &            bwk(jd(j)), xwk(jd(j)), jd, wk, wk1, wk2)
         call gs(j, n, bwk(jd(j)), xwk(jd(j)), icdomn, 1, 
     &           0.d0, 0.d0, 0.d0, wk, nit0)
c         call gssor(j, n, bwk(jd(j)), xwk(jd(j)), icdomn, 1, 1,
c     &              0.d0, 0.d0, 0.d0, wk, nit0)

c********Step 3. Correction step.***********************************************

         xdiffn = 0.d0
         do 120 i = 0, n - 1
	    if (icwk .eq. 1) then
               xnew   = xwkt(jd(j)+i) + xwk(jd(j)+i)
	    else if (icwk .eq. 2) then
 	       xnew   = xwk(jd(j)+i)
	    end if
            xd     = xnew - x(i+1)
            xdiffn = xdiffn + xd*xd
            x(i+1) = xnew                    !Correction, i.e., new solution.
 120     continue
c********End of Step 3.*********************************************************

         if (icit .eq. 1 .and. it .ne. itmax) then
	    go to 130
	 else if (icit .eq. 1 .and. it .eq. itmax) then
	    go to 150
	 end if

         xoldn = xnewn
         xnewn = dot(x, x, n)
	 if (xnewn .lt. 1.d-40) then
	    if (xoldn .lt. 1.d-40) then
ccc	       go to 333
	    else
	       stopl = 1.d0
	    end if
	 else
	    stopl  = xdiffn/xnewn
	 end if

c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*x.
 130     call matvc2(j, x, wk, icdomn, 3, c, d, f)

c        Residual: bwk =  b - wk =  b - (A_0 + cA_x0 + dA_y0 + fM_0)*x.
         do 140 i = 1, n
            bwk(jd(j)-1+i) = b(i) - wk(i)
 140     continue
         if (icit .eq. 1)   go to 150

         rnorm = dot(bwk(jd(j)), bwk(jd(j)), n)
c	 print*, ' rnorm in TG1', rnorm

ccc         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 150  continue

c      print*, ' Iteration limit exceeded:', itmax

 333  nit0 = min(it, itmax)
c      print 444, j, nit0
c 444  format(' No. of iterations in TG1 at level ', i2, ':',i3)
c      print*
 
 555  return
      end
c=======================================================================

      subroutine tgrid2(j, n, icdomn, iexct, icalgm, niter,c,d,f,b,x,
     &                                          icit, jd, wk, wk1, wk2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TGRID1 stands for Two-Grid method and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a Two-Grid method.
c          icit = 1: 1, 2, or 4 iterations
c          icit = 0: more than 1 iteration, depending on stopping criteria
c          icit = 2: same as icit = 0, but do not initialize the soln x
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim = 263170)

      implicit double precision (a-h, o-z)
      dimension  x(n), b(n), wk(n), xwk(2*lim), bwk(2*lim), xwkt(2*lim)
      dimension  wk1(2*n), wk2(2*n), jd(11)

      itmax  = niter
      nl     = jd(j) - jd(j-1)             !Dimension of the (j-1)th level.
      lvlcst = iexct                               !Coarsest level.
      iexsol = 1
      if (icit   .eq. 1) itmax  = 1
c      if (icit .eq. 1 .and. j .eq. lvlcst+1) itmax  = 1 

      do 10 i = 1, jd(j)-1
         xwk(i)  = 0.d0
         bwk(i)  = 0.d0
         xwkt(i) = 0.d0
 10   continue

      do 20 i = 1, n
         x(i)           = 0.d0
	 xwkt(jd(j)+i-1)= 0.d0
	 xwk(jd(j)+i-1) = 0.d0
	 bwk(jd(j)+i-1) = b(i)
 20   continue

      do 150 it = 1, itmax
c********Step 1. Solve the nonsymmetric problem at coarser space.***************

         call restrn(j, bwk(jd(j)), bwk(jd(j-1)))   !Restriction to lower level.
         call bpxin(j-1, bwk(jd(j-1)), icdomn)

         if (j-1 .eq. lvlcst) then
            if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c              Solving exactly at the coarsest level.
               call solexc(bwk(jd(lvlcst)), xwk(jd(lvlcst)), icdomn, 
     &                     lvlcst, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c              by Bi-CGSTAB method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call bcgst(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    0.1d0, nit, icdomn, 3, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c              by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call gssor(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    icdomn, 3, 3, c, d, f, wk, nit)
            end if
            go to 40
         end if

c        Solve (A_0+cA_0x+dA_0y+fM_0)*xwk=bwk by two-grid method at lower level.
         call tgrid1(j-1, nl, icdomn, iexct, icalgm, 4, c, d, f, 
     &               bwk(jd(j-1)), xwk(jd(j-1)), 1, jd, wk, wk1, wk2)

c********Step 2. Solve the symmetric problem at finer space.********************
 40      call prolng(j,xwk(jd(j-1)),xwk(jd(j)))   !Prolongation to upper level.
         call bpxin(j, xwk(jd(j)), icdomn)

         do 50 i = 0, n-1
            xwkt(jd(j)+i) = x(i+1) + xwk(jd(j)+i)
 50      continue

         icwk = 1                               !icit = 1 is better somehow.
         if (icwk .eq. 1) then
c           To compute wk = (A_h + cA_x + dA_y + fM)*xwkt at level j.
            call matvc2(j, xwkt(jd(j)), wk, icdomn, 3, c, d, f)
         else if (icwk .eq. 2) then
c           To compute wk = (cA_x + dA_y + fM)*xwkt at level j.
            call matvc2(j, xwkt(jd(j)), wk, icdomn, 5, c, d, f)
         end if

c        Compute bwk = b - wk = b - (A_h + cA_x + dA_y + fM)*x.
         do 60 i = 1, n
            bwk(jd(j)+i-1) = b(i) - wk(i)
 60      continue

c        To solve the SPD linear system A_h*xwk = bwk.
c         call mgv(j, n, icdomn, 1, 1, nit0, 0.d0, 0.d0, 0.d0, 
c     &            bwk(jd(j)), xwk(jd(j)), jd, wk, wk1, wk2)
         call gs(j, n, bwk(jd(j)), xwk(jd(j)), icdomn, 1, 
     &           0.d0, 0.d0, 0.d0, wk, nit0)
c         call gssor(j, n, bwk(jd(j)), xwk(jd(j)), icdomn, 1, 1,
c     &              0.d0, 0.d0, 0.d0, wk, nit0)

c********Step 3. Correction step, i.e., new solution.***************************

         do 120 i = 0, n - 1
	    if (icwk .eq. 1) then
               x(i+1) = xwkt(jd(j)+i) + xwk(jd(j)+i)
	    else if (icwk .eq. 2) then
 	       x(i+1) = xwk(jd(j)+i)
	    end if
 120     continue
c********End of Step 3.*********************************************************

	 if (it .eq. itmax) go to 555

c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*x.
 130     call matvc2(j, x, wk, icdomn, 3, c, d, f)

c        Residual: bwk =  b - wk =  b - (A_0 + cA_x0 + dA_y0 + fM_0)*x.
         do 140 i = 1, n
            bwk(jd(j)-1+i) = b(i) - wk(i)
 140     continue
 150  continue

 555  return
      end
c=======================================================================

      subroutine tgrid3(j, n, icdomn, iexct, icalgm, niter,c,d,f,b,x,
     &                                          icit, jd, wk, wk1, wk2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TGRID1 stands for Two-Grid method and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a Two-Grid method.
c          icit = 1: just 1 iteration
c          icit = 0: more than 1 iteration, depending on stopping criteria
c          icit = 2: same as icit = 0, but do not initialize the soln x
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim = 263170)

      implicit double precision (a-h, o-z)
      dimension  x(n), b(n), wk(n), xwk(2*lim), bwk(2*lim), xwkt(2*lim)
      dimension  wk1(2*n), wk2(2*n), jd(11)

      tol   = 1.d-40/(4**j)            !tol - stopping criteria for iterations
c      itmax = 1                      !max. iteration
      itmax = niter
      iexsol = 1
      xnewn = 0.d0
      stopl = 1.d0
      nl    = jd(j) - jd(j-1)             !Dimension of the (j-1)th level.
 
      lvlcst = iexct                               !Coarsest level.  
      if (icit   .eq. 1) itmax  = 1  
c      if (icdomn .eq. 3 .and. iexct .eq. 1)   lvlcst = 1

      if (j .eq. lvlcst) then           !Solving exactly at the coarsest level
         do 5 i = 1, n
            x(i) = 0.d0
 5       continue

         if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c           Compute x = (A_0 + cA_x0 + dA_y0 + fM_0)^(-1)*b.
            call solexc(b, x, icdomn, lvlcst, c, d, f)
         else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c           by Bi-CGSTAB method in the structured mesh.
            nd = jd(lvlcst+1) - jd(lvlcst)
            call bcgst(lvlcst, nd, b, x, 0.1d0, nit,
     &                  icdomn, 3, c, d, f)
         else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c           by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
            nd = jd(lvlcst+1) - jd(lvlcst)
            call gssor(lvlcst, nd, b, x,
     &                 icdomn, 3, 3, c, d, f, wk, nit)
          end if
         go to 555
      end if

      do 10 i = 1, jd(j+1)-1
         xwk(i)  = 0.d0
         bwk(i)  = 0.d0
         xwkt(i) = 0.d0
 10   continue

      if (icit .eq. 2) then
c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*x.
         call matvc2(j, x, wk, icdomn, 3, c, d, f)
c        Residual: bwk = b - wk =  b - (A_0 + cA_x0 + dA_y0 + fM_0)*x.
         do 20 i = 1, n 
            bwk(jd(j)-1+i) = b(i) - wk(i)
 20      continue
      else
         do 30 i = 1, n
            x(i)           = 0.d0
            bwk(jd(j)+i-1) = b(i)
 30      continue
      end if

      do 150 it = 1, itmax

c********Step 1. Solve the nonsymmetric problem at coarser space.***************
         call restrn(j, bwk(jd(j)), bwk(jd(j-1)))   !Restriction to lower level.
         call bpxin(j-1, bwk(jd(j-1)), icdomn)

         if (j-1 .eq. lvlcst) then
            if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c              Solving exactly at the coarsest level.
               call solexc(bwk(jd(lvlcst)), xwk(jd(lvlcst)), icdomn, 
     &                     lvlcst, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c              by Bi-CGSTAB method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call bcgst(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    0.1d0, nit, icdomn, 3, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c              by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call gssor(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    icdomn, 3, 3, c, d, f, wk, nit)
            end if
            go to 40
         end if

c        Solve (A_0+cA_0x+dA_0y+fM_0)*xwk=bwk by two-grid method at lower level.
         call tgrid4(j-1, nl, icdomn, iexct, icalgm, 1, c, d, f, 
     &               bwk(jd(j-1)), xwk(jd(j-1)), 1, jd, wk, wk1, wk2)

c********Step 2. Solve the symmetric problem at finer space.********************
 40      call prolng(j,xwk(jd(j-1)),xwk(jd(j)))   !Prolongation to upper level.
         call bpxin(j, xwk(jd(j)), icdomn)

         if (icit .ne. 1) then
            do 50 i = 0, n-1
               xwk(jd(j)+i) = x(i+1) + xwk(jd(j)+i)
 50         continue
         end if

c        To compute xwkt = (A_h + cA_x + dA_y + fM)*xwk at level j.
         call matvc2(j, xwk(jd(j)), xwkt(jd(j)), icdomn, 3, c, d, f)

c        Compute bwk = b - xwkt = b - (A_h + cA_x + dA_y + fM)*xwk.
         do 60 i = 0, n-1
            bwk(jd(j)+i) = b(i+1) - xwkt(jd(j)+i)
 60      continue

c        To solve the SPD linear system A_h*xwkt = bwk.
c         call mgv(j, n, icdomn, 1, 1, nit0, 0.d0, 0.d0, 0.d0, 
c     &            bwk(jd(j)), xwkt(jd(j)), jd, wk, wk1, wk2)
         call gs(j, n, bwk(jd(j)), xwkt(jd(j)), icdomn, 1, 
     &           0.d0, 0.d0, 0.d0, wk, nit0)
c         call gssor(j, n, bwk(jd(j)), xwkt(jd(j)), icdomn, 1, 1,
c     &              0.d0, 0.d0, 0.d0, wk, nit0)

         do 70 i = 0, n - 1
            xwkt(jd(j)+i) = xwk(jd(j)+i) + xwkt(jd(j)+i)
 70      continue

c********Step 3. Solve the nonsymmetric problem at coarse space.****************
c        Compute xwk = (A + cA_x + dA_y + fM)*xwkt at level j.
         call matvc2(j, xwkt(jd(j)), xwk(jd(j)), icdomn, 3, c, d, f)

c        Compute bwk = b - xwk = b - (A + cA_x + dA_y + fM)*xwkt.
         do 90 i = 0, n-1
            bwk(jd(j)+i) = b(i+1) - xwk(jd(j)+i)
 90      continue

         call restrn(j, bwk(jd(j)), bwk(jd(j-1)))   !Restriction to lower level.
         call bpxin(j-1, bwk(jd(j-1)), icdomn)

         if (j-1 .eq. lvlcst) then
            if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c              Solving exactly at the coarsest level.
               call solexc(bwk(jd(lvlcst)), xwk(jd(lvlcst)), icdomn, 
     &                     lvlcst, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c              by Bi-CGSTAB method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call bcgst(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    0.1d0, nit, icdomn, 3, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c              by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call gssor(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    icdomn, 3, 3, c, d, f, wk, nit)
            end if
            go to 100
         end if

c        Solve (A_0+cA_0x+dA_0y+fM_0)*xwk=bwk by two-grid method at lower level.
         call tgrid4(j-1, nl, icdomn, iexct, icalgm, 1, c, d, f, 
     &               bwk(jd(j-1)), xwk(jd(j-1)), 1, jd, wk, wk1, wk2)

c********Step 4. Correction step.***********************************************
 100     call prolng(j,xwk(jd(j-1)),xwk(jd(j))) !Prolongation to upper level.
         call bpxin(j, xwk(jd(j)), icdomn)

         xdiffn = 0.d0
         do 130 i = 0, n - 1
            xnew   = xwkt(jd(j)+i) + xwk(jd(j)+i)
            if (icit .eq. 1)   go to 110
            xd     = xnew - x(i+1)
            xdiffn = xdiffn + xd*xd
 110        x(i+1) = xnew       !Correction, i.e., new solution.
 130     continue
c********End of Step 4.*********************************************************

         if (icit .eq. 1)   go to 333
         xoldn = xnewn
         xnewn = dot(x, x, n)
	 if (xnewn .lt. 1.d-40) then
	    if (xoldn .lt. 1.d-40) then
ccc	       go to 333
	    else
	       stopl = 1.d0
	    end if
	 else
	    stopl  = xdiffn/xnewn
	 end if

c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*x.
         call matvc2(j, x, wk, icdomn, 3, c, d, f)

c        Residual: bwk =  b - wk =  b - (A_0 + cA_x0 + dA_y0 + fM_0)*x.
         do 140 i = 1, n
            bwk(jd(j)-1+i) = b(i) - wk(i)
 140     continue

         rnorm = dot(bwk(jd(j)), bwk(jd(j)), n)
c	 print*, ' rnorm in TG1', rnorm

ccc         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 150  continue

c      print*, ' Iteration limit exceeded:', itmax

 333  nit0 = min(it, itmax)
c      print 444, j, nit0
c 444  format(' No. of iterations in TG1 at level ', i2, ':',i3)
c      print*
 
 555  return
      end
c=======================================================================

      subroutine tgrid4(j, n, icdomn, iexct, icalgm, niter, c,d,f,b,x,
     &                                          icit, jd, wk, wk1, wk2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TGRID1 stands for Two-Grid method and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a Two-Grid method.
c          icit = 1: just 1 iteration
c          icit = 0: more than 1 iteration, depending on stopping criteria
c          icit = 2: same as icit = 0, but do not initialize the soln x
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim = 263170)

      implicit double precision (a-h, o-z)
      dimension  x(n), b(n), wk(n), xwk(2*lim), bwk(2*lim), xwkt(2*lim)
      dimension  wk1(2*n), wk2(2*n), jd(11)

      nl     = jd(j) - jd(j-1)             !Dimension of the (j-1)th level.
      lvlcst = iexct                               !Coarsest level.
      iexsol = 1  
c      if (icdomn .eq. 3 .and. iexct .eq. 1)   lvlcst = 1

      do 20 i = 1, jd(j+1)-1
         xwk(i)  = 0.d0
         bwk(i)  = 0.d0
         xwkt(i) = 0.d0
 20   continue

      do 30 i = 1, n
         x(i)           = 0.d0
	 bwk(jd(j)+i-1) = b(i)
 30   continue

c********Step 1. Solve the nonsymmetric problem at coarser space.***************
         call restrn(j, bwk(jd(j)), bwk(jd(j-1)))   !Restriction to lower level.
         call bpxin(j-1, bwk(jd(j-1)), icdomn)

         if (j-1 .eq. lvlcst) then
            if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c              Solving exactly at the coarsest level.
               call solexc(bwk(jd(lvlcst)), xwk(jd(lvlcst)), icdomn, 
     &                     lvlcst, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c              by Bi-CGSTAB method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call bcgst(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    0.1d0, nit, icdomn, 3, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c              by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call gssor(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    icdomn, 3, 3, c, d, f, wk, nit)
            end if
            go to 40
         end if

c        Solve (A_0+cA_0x+dA_0y+fM_0)*xwk=bwk by two-grid method at lower level.
         call tgrid3(j-1, nl, icdomn, iexct, icalgm, 1, c, d, f, 
     &               bwk(jd(j-1)), xwk(jd(j-1)), 1, jd, wk, wk1, wk2)

c********Step 2. Solve the symmetric problem at finer space.********************
 40      call prolng(j,xwk(jd(j-1)),xwk(jd(j)))   !Prolongation to upper level.
         call bpxin(j, xwk(jd(j)), icdomn)

         if (icit .ne. 1) then
            do 50 i = 0, n-1
               xwk(jd(j)+i) = x(i+1) + xwk(jd(j)+i)
 50         continue
         end if

c        To compute xwkt = (A_h + cA_x + dA_y + fM)*xwk at level j.
         call matvc2(j, xwk(jd(j)), xwkt(jd(j)), icdomn, 3, c, d, f)

c        Compute bwk = b - xwkt = b - (A_h + cA_x + dA_y + fM)*xwk.
         do 60 i = 0, n-1
            bwk(jd(j)+i) = b(i+1) - xwkt(jd(j)+i)
 60      continue

c        To solve the SPD linear system A_h*xwkt = bwk.
c         call mgv(j, n, icdomn, 1, 1, nit0, 0.d0, 0.d0, 0.d0, 
c     &            bwk(jd(j)), xwkt(jd(j)), jd, wk, wk1, wk2)
         call gs(j, n, bwk(jd(j)), xwkt(jd(j)), icdomn, 1, 
     &           0.d0, 0.d0, 0.d0, wk, nit0)
c         call gssor(j, n, bwk(jd(j)), xwkt(jd(j)), icdomn, 1, 1,
c     &              0.d0, 0.d0, 0.d0, wk, nit0)

         do 70 i = 0, n - 1
            xwkt(jd(j)+i) = xwk(jd(j)+i) + xwkt(jd(j)+i)
 70      continue

c********Step 3. Solve the nonsymmetric problem at coarse space.****************
c        Compute xwk = (A + cA_x + dA_y + fM)*xwkt at level j.
         call matvc2(j, xwkt(jd(j)), xwk(jd(j)), icdomn, 3, c, d, f)

c        Compute bwk = b - xwk = b - (A + cA_x + dA_y + fM)*xwkt.
         do 90 i = 0, n-1
            bwk(jd(j)+i) = b(i+1) - xwk(jd(j)+i)
 90      continue

         call restrn(j, bwk(jd(j)), bwk(jd(j-1)))   !Restriction to lower level.
         call bpxin(j-1, bwk(jd(j-1)), icdomn)

         if (j-1 .eq. lvlcst) then
            if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c              Solving exactly at the coarsest level.
               call solexc(bwk(jd(lvlcst)), xwk(jd(lvlcst)), icdomn, 
     &                     lvlcst, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c              by Bi-CGSTAB method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call bcgst(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    0.1d0, nit, icdomn, 3, c, d, f)
            else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c              Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c              by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
               nd = jd(lvlcst+1) - jd(lvlcst)
               call gssor(lvlcst, nd, bwk(jd(lvlcst)), xwk(jd(lvlcst)),
     &                    icdomn, 3, 3, c, d, f, wk, nit)
            end if
            go to 100
         end if

c        Solve (A_0+cA_0x+dA_0y+fM_0)*xwk=bwk by two-grid method at lower level.
         call tgrid3(j-1, nl, icdomn, iexct, icalgm, 1, c, d, f, 
     &               bwk(jd(j-1)), xwk(jd(j-1)), 1, jd, wk, wk1, wk2)

c********Step 4. Correction step.***********************************************
 100     call prolng(j,xwk(jd(j-1)),xwk(jd(j))) !Prolongation to upper level.
         call bpxin(j, xwk(jd(j)), icdomn)

         xdiffn = 0.d0
         do 130 i = 0, n - 1
            xnew   = xwkt(jd(j)+i) + xwk(jd(j)+i)
            x(i+1) = xnew       !Correction, i.e., new solution.
 130     continue
c********End of Step 4.*********************************************************
 
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
      dimension  u(n), b(n), wk(n), rhs(2*n), actnB(2*n), jd(11)

      ismth = 2                  !Choice of smoother, 1, 2, 3, 4.
      imatv = 2                  !Choice of matrix, FE or UPWD, 2, 3, 4
      ifbs3 = 1                  !Choice of smoother in exact solver, 1, 2, 3.
c     isolv = 0, 1, 2, 3, 4      !Choice of solver at the coarest level.

c      if (ismth .eq. 3) then
c         tol   = 1.d-12
c      else 
         tol   = 1.d-18/(4**j)      !tol - stopping criteria for iterations
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
	    if (ismth .eq. 1) then 
               call gausei(k, rhs(jd(k)), actnB(jd(k)), icdomn, ipde,
     &                     c, d, f, ifbs1)
c              call bpxin(k, actnB(jd(k)), icdomn)
	       do 80 ii = 2, nsmthg
c                 Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0+fM_0)*actnB.
                  call matvc2(k, actnB(jd(k)), wk, icdomn, ipde,c,d,f)

c                 Residual: wk = rhs - wk = rhs-(A_0+cA_x0+dA_y0+fM_0)*actnB.
                  do 60 i = 1, kd 
                     wk(i) = rhs(jd(k)-1+i) - wk(i)
 60               continue

	          call gausei(k, wk, wk, icdomn, ipde, c, d, f, ifbs1)
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
	    else if (ismth .eq. 4) then
	       if (k .eq. j) then
                  call gssor(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                       ipde, ifbs1, c, d, f, wk, nsmthg)
	       else
                  call gsupwd(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                        ipde, ifbs1, c, d, f, wk, nsmthg)
	       end if 
	    end if 

            if (imatv .eq. 2) then
c              Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
               call matvc2(k, actnB(jd(k)), wk, icdomn, ipde, c, d, f)
            else if (imatv .eq. 3) then
c              Call matvc3 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
               call matvc3(k, actnB(jd(k)), wk, icdomn, ipde, c, d, f)
            else if (imatv .eq. 4) then
	       if (k .eq. j) then
c                 Call matvc2 to compute wk=(A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
                  call matvc2(k, actnB(jd(k)), wk, icdomn, ipde,c,d,f)
               else
c                 Call matvc3 to compute wk=(A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
                  call matvc3(k, actnB(jd(k)), wk, icdomn, ipde,c,d,f)
               end if
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

	    if (ismth .eq. 1) then
               do 140 ii = 1, nsmthg
c                 Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0+fM_0)*actnB.
                  call matvc2(k, actnB(jd(k)), wk, icdomn, ipde, c,d,f)

c                 Residual: wk = rhs - wk =  rhs - (A_0+cA_x0+dA_y0+fM_0)*u.
                  do 120 i = 1, kd 
                     wk(i) = rhs(jd(k)-1+i) - wk(i)
 120              continue

c                 Postsmoothing by the symmetric Gauss-seidel smoother.
                  call gausei(k, wk, wk, icdomn, ipde, c, d, f, ifbs2)
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
	    else if (ismth .eq. 4) then
	       if (k .eq. j) then
                  call gssor(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                       ipde, ifbs2, c, d, f, wk, nsmthg)
	       else
                  call gsupwd(k, kd, rhs(jd(k)), actnB(jd(k)), icdomn,
     &                        ipde, ifbs2, c, d, f, wk, nsmthg)
	       end if 
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

 333  nit = min(kk, itmax)
      print*

      return
      end
c=======================================================================

      subroutine mgv(j, n, icdomn, iexct, ichoic, niter, c, d, f, b, u, 
     &                                               jd, wk, rhs, actnB)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MGV stands for multigrid V-cycle and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a multigrid method using V-cycle
c     algorithm in the structured mesh.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (nsmthg = 1, ifbs1 = 1, ifbs2 = 2)

      implicit double precision (a-h, o-z)
      dimension  u(n), b(n), wk(n), rhs(2*n), actnB(2*n), jd(11)

      tol   = 1.d-50/(4**j)            !tol - stopping criteria for iterations
c      itmax = 50                        !max. iteration
      itmax = niter
      iexsol = 1
      xnewn = 0.d0
      stopl = 1.d0

      lvlcst = iexct
c      if (icdomn .eq. 3 .and. iexct .eq. 1) lvlcst = 1

      do 30 i = 1, jd(j+1) - 1
         rhs(i)    = 0.d0
         actnB(i)  = 0.d0
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

         do 100 k = j, lvlcst+1, -1
            kd = jd(k+1) - jd(k)

c           Presmoothing by the symmetric Gauss-seidel smoother nsmthg times.
            call gausei(k, rhs(jd(k)), actnB(jd(k)), icdomn, ichoic,
     &                  c, d, f, ifbs1)
c           call bpxin(k, actnB(jd(k)), icdomn)
	    do 80 ii = 2, nsmthg
c              Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
               call matvc2(k, actnB(jd(k)), wk, icdomn, ichoic, c, d,f)

c              Residual: wk = rhs - wk = rhs-(A_0+cA_x0+dA_y0+fM_0)*actnB.
               do 60 i = 1, kd 
                  wk(i) = rhs(jd(k)-1+i) - wk(i)
 60            continue

	       call gausei(k, wk, wk, icdomn, ichoic, c, d, f, ifbs1)
c              call bpxin(k, wk, icdomn)

c              Correction step.
               do 70 i = 1, kd 
                  actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 70            continue
 80         continue

c           Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
            call matvc2(k, actnB(jd(k)), wk, icdomn, ichoic, c, d, f)

c           Residual: wk = rhs - wk =  rhs - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
            do 90 i = 1, kd 
               wk(i) = rhs(jd(k)-1+i) - wk(i)
 90         continue

c           Restriction to the lower level.
            call restrn(k, wk, rhs(jd(k-1)))
            call bpxin(k-1, rhs(jd(k-1)), icdomn)
 100     continue

         if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c           Solving exactly at the coarsest level.
            call solexc(rhs(jd(lvlcst)), actnB(jd(lvlcst)), icdomn, 
     &                  lvlcst, c, d, f)
         else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c           by Bi-CGSTAB method in the structured mesh.
            nd = jd(lvlcst+1) - jd(lvlcst)
            call bcgst(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &                 0.1d0, nit, icdomn, ichoic, c, d, f)
         else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c           Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c           by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
            nd = jd(lvlcst+1) - jd(lvlcst)
            call gssor(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &                 icdomn, ichoic, 3, c, d, f, wk, nit)
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

            do 140 ii = 1, nsmthg
c              Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
               call matvc2(k, actnB(jd(k)), wk, icdomn, ichoic, c, d, f)

c              Residual: wk = rhs - wk =  rhs - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
               do 120 i = 1, kd 
                  wk(i) = rhs(jd(k)-1+i) - wk(i)
 120           continue

c              Postsmoothing by the symmetric Gauss-seidel smoother.
               call gausei(k, wk, wk, icdomn, ichoic, c, d, f, ifbs2)
c               call bpxin(k, wk, icdomn)

c              Final correction step.
               do 130 i = 1, kd 
                  actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 130           continue
 140        continue
 150     continue 

c*       New solution:  u = u + B(b - (A_0 + cA_x0 + dA_y0 + fM_0)u).

         do 160 i = 1, n
            u(i) = u(i) + actnB(jd(j)-1+i)
 160     continue

         xoldn = xnewn
         xnewn = dot(u, u, n)
	 if (xnewn .lt. 1.d-50) then
	    if (xoldn .lt. 1.d-50) then
ccc	       go to 333
	    else
	       stopl = 1.d0
	    end if
	 else
	    xdiffn = dot(actnB(jd(j)), actnB(jd(j)) , n)
	    stopl  = xdiffn/xnewn
	 end if

c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         call matvc2(j, u, wk, icdomn, ichoic, c, d, f)

c        Residual: rhs =  b - wk =  b - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         do 170 i = 1, n
            rhs(jd(j)-1+i) = b(i) - wk(i)
 170     continue

         rnorm = dot(rhs(jd(j)), rhs(jd(j)), n)
c	 print*, ' rnorm', rnorm

ccc         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 180  continue

c      print*, ' Iteration limit exceeded:', itmax

 333  nit1 = min(kk, itmax)
c      print 444, j, nit1
c444  format(' No. of iterations in MG V-cycle at level ', i2, ':',i3)
c      print*

      return
      end
c=======================================================================

      subroutine mgv1(j, n, icdomn, iexct, ichoic, niter, c, d, f, b, u, 
     &                                              jd, wk, rhs, actnB)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MGV1 stands for multigrid V-cycle and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a multigrid method using V-cycle
c     algorithm in the structured mesh. It is same as MGV.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim=263170, nsmthg=1, ifbs1=1, ifbs2=2)

      implicit double precision (a-h, o-z)
      dimension  u(n), b(n), wk(n), rhs(2*n), actnB(2*n), jd(11)
      dimension  un(lim), bn(lim)

      tol   = 1.d-50/(4**j)            !tol - stopping criteria for iterations
c      itmax = 50                        !max. iteration
      itmax = niter
      xnewn = 0.d0
      stopl = 1.d0

      lvlcst = iexct
c      if (icdomn .eq. 3 .and. iexct .eq. 1) lvlcst = 1

      do 10 i = 1, n
         u(i)  = 0.d0
         bn(i) = b(i)
 10   continue

c*    Iterate until convergence:  u = u + B(b - (A_0 + cA_x0 + dA_y0 + fM_0)u).

      do 180 kk = 1, itmax
c        Iterator using MG V-cycle, i.e., un = B*bn.
         call premgv(j, n, icdomn, ichoic, c, d, f, bn, un, 
     &               lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, rhs, actnB)

c*       New solution:  u = u + B(b - (A_0 + cA_x0 + dA_y0 + fM_0)u).

         do 160 i = 1, n
            u(i) = u(i) + un(i)
 160     continue

         xoldn = xnewn
         xnewn = dot(u, u, n)
	 if (xnewn .lt. 1.d-50) then
	    if (xoldn .lt. 1.d-50) then
ccc	       go to 333
	    else
	       stopl = 1.d0
	    end if
	 else
	    xdiffn = dot(un, un , n)
	    stopl  = xdiffn/xnewn
	 end if

c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         call matvc2(j, u, wk, icdomn, ichoic, c, d, f)

c        Residual: bn =  b - wk =  b - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         do 170 i = 1, n
            bn(i) = b(i) - wk(i)
 170     continue

         rnorm = dot(bn, bn, n)
c	 print*, ' rnorm', rnorm

ccc         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 180  continue

c      print*, ' Iteration limit exceeded:', itmax

 333  nit1 = min(kk, itmax)
c      print 444, j, nit1
c444  format(' No. of iterations in MG V-cycle at level ', i2, ':',i3)
c      print*

      return
      end
c=======================================================================

      subroutine mgvf(j, n, icdomn, iexct, ichoic, niter, c, d, f, b, u, 
     &                                         jd, wk, rhs, actnB, bwk)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MGVF stands for full multigrid V-cycle and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a full multigrid method using
c     V-cycle algorithm in the structured mesh.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim=263170, nsmthg=1, ifbs1=1, ifbs2=2)

      implicit double precision (a-h, o-z)
      dimension  u(n), b(n), wk(n), rhs(2*n), actnB(2*n), jd(11)
      dimension  un(lim), bn(lim), bwk(2*n)

      iexsol = 1
      lvlcst = iexct
c      if (icdomn .eq. 3 .and. iexct .eq. 1) lvlcst = 1

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

      if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c        Solving exactly at the coarsest level.
         call solexc(bwk(jd(lvlcst)), u, icdomn, lvlcst, c, d, f)
      else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c        by Bi-CGSTAB method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
         call bcgst(lvlcst, nd, bwk(jd(lvlcst)), u,
     &                 0.1d0, nit, icdomn, ichoic, c, d, f)
      else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c        by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
         call gssor(lvlcst, nd, bwk(jd(lvlcst)), u,
     &                 icdomn, ichoic, 3, c, d, f, wk, nit)
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
         call matvc2(k, u, wk, icdomn, ichoic, c, d, f)

c        Residual: bn = bwk - wk =  bwk - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         do 80 ii = 1, kd 
            bn(ii) = bwk(jd(k)-1+ii) - wk(ii)
 80      continue

c*       Iterate until convergence: u = u+B(bwk -(A_0 + cA_x0 + dA_y0 + fM_0)u).

         do 180 kk = 1, itmax
c           Iterator using MG V-cycle, i.e., un = B*bn.
            call premgv(k, kd, icdomn, ichoic, c, d, f, bn, un, 
     &                lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, rhs, actnB)

c*          New solution (final correction step):  u = u + un.
            do 160 i = 1, kd
               u(i) = u(i) + un(i)
 160        continue

            xoldn = xnewn
            xnewn = dot(u, u, kd)
            if (xnewn .lt. 1.d-50) then
               if (xoldn .lt. 1.d-50) then
ccc	          go to 333
               else
                  stopl = 1.d0
               end if
            else
               xdiffn = dot(un, un , kd)
               stopl  = xdiffn/xnewn
            end if

c           Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*u.
            call matvc2(k, u, wk, icdomn, ichoic, c, d, f)

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

      subroutine tgf(j, n, icdomn, iexct, icalgm, niter, c, d, f, b, u, 
     &                                        jd, wk, wk1, wk2, bwk)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TGF stands for full two-grid algorithm  and solves the linear system 
c     (A_0 + cA_x0 + dA_y0 + fM_0)*u = b by a full two-grid method
c     in the structured mesh.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim=263170, nsmthg=1)

      implicit double precision (a-h, o-z)
      dimension  u(n), b(n), wk(n), wk1(2*n), wk2(2*n), jd(11)
      dimension  un(lim), bwk(2*n)

      iexsol = 1
      lvlcst = iexct
c      if (icdomn .eq. 3 .and. iexct .eq. 1) lvlcst = 1

      do 30 i = 1, jd(j) - 1
         bwk(i) = 0.d0
 30   continue

      do 40 i = 1, n
         u(i)           = 0.d0
         un(i)          = 0.d0
         bwk(jd(j)-1+i) = b(i)
 40   continue

      do 60 i = j, lvlcst+1, -1
c        Restriction of b to the lower levels.
         call restrn(i,bwk(jd(i)),bwk(jd(i-1)))
         call bpxin(i-1, bwk(jd(i-1)), icdomn)
 60   continue

      if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c        Solving exactly at the coarsest level.
         call solexc(bwk(jd(lvlcst)), u, icdomn, lvlcst, c, d, f)
      else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c        by Bi-CGSTAB method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
         call bcgst(lvlcst, nd, bwk(jd(lvlcst)), u,
     &                 0.1d0, nit, icdomn, 3, c, d, f)
      else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c        by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
         call gssor(lvlcst, nd, bwk(jd(lvlcst)), u,
     &                 icdomn, 3, 3, c, d, f, wk, nit)
      end if

      do 200 k = lvlcst+1, j                                !Full TG starts.
         kd    = jd(k+1) - jd(k)
         kdl   = jd(k) - jd(k-1)

         do 70 ii = 1, kdl
            un(ii) = u(ii) 
 70      continue

c        Prolongation to the upper level.
         call prolng(k, un, u)
         call bpxin(k, u, icdomn)

c*       Iterate until convergence: u = u+B(bwk -(A_0 + cA_x0 + dA_y0 + fM_0)u).

         call tgrid1(k, kd, icdomn, iexct, icalgm, niter, c, d, f, 
     &               bwk(jd(k)), u, 2, jd, wk, wk1, wk2)

c         print 190, k, niter
c190      format(' No. of iterations of Full TG at level ',i2,':',i3)
c         print*
200   continue

      return
      end
c=======================================================================

      subroutine premgv(j, n, icdomn, ichoic, c, d, f, b, u, 
     &               lvlcst, nsmthg, ifbs1, ifbs2, jd, wk, rhs, actnB)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     PREMGV stands for multigrid V-cycle as a preconditioner, i.e., u=Bb
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension  u(n), b(n), wk(n), rhs(2*n), actnB(2*n), jd(11)

      iexsol = 1

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

c        Presmoothing by the symmetric Gauss-seidel smoother nsmthg times.
         call gausei(k, rhs(jd(k)), actnB(jd(k)), icdomn, ichoic,
     &                  c, d, f, ifbs1)
c        call bpxin(k, actnB(jd(k)), icdomn)
         do 80 ii = 2, nsmthg
c           Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
            call matvc2(k, actnB(jd(k)), wk, icdomn, ichoic, c, d,f)

c           Residual: wk = rhs - wk = rhs-(A_0+cA_x0+dA_y0+fM_0)*actnB.
            do 60 i = 1, kd 
               wk(i) = rhs(jd(k)-1+i) - wk(i)
 60         continue

            call gausei(k, wk, wk, icdomn, ichoic, c, d, f, ifbs1)
c           call bpxin(k, wk, icdomn)

c           Correction step.
            do 70 i = 1, kd 
               actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 70         continue
 80      continue
            
c        Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
         call matvc2(k, actnB(jd(k)), wk, icdomn, ichoic, c, d, f)

c        Residual: wk = rhs - wk =  rhs - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
         do 90 i = 1, kd 
            wk(i) = rhs(jd(k)-1+i) - wk(i)
 90      continue

c        Restriction to the lower level.
         call restrn(k, wk, rhs(jd(k-1)))
         call bpxin(k-1, rhs(jd(k-1)), icdomn)
 100  continue

      if (lvlcst .eq. 1 .or. lvlcst .eq. 2) then
c        Solving exactly at the coarsest level.
         call solexc(rhs(jd(lvlcst)), actnB(jd(lvlcst)), icdomn, 
     &                  lvlcst, c, d, f)
      else if (lvlcst .ge. 3 .and. iexsol .eq. 1) then
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x = b
c        by Bi-CGSTAB method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
         call bcgst(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &                 0.1d0, nit, icdomn, ichoic, c, d, f)
      else if (lvlcst .ge. 3 .and. iexsol .eq. 2) then           
c        Solve the discretized problem (A_0 + cA_0x + dA_0y + fM_0)*x=bwk
c        by Gauss-Seidel or SSOR(SSUR) method in the structured mesh.
         nd = jd(lvlcst+1) - jd(lvlcst)
         call gssor(lvlcst, nd, rhs(jd(lvlcst)), actnB(jd(lvlcst)),
     &                 icdomn, ichoic, 3, c, d, f, wk, nit)
      end if

c     Going upward in the V-cycle.
      do 150 k = lvlcst+1, j 
         kd = (2**k+1)**2

c        Prolongation to the upper level.
         call prolng(k, actnB(jd(k-1)), wk)
         call bpxin(k, wk, icdomn)

c        Correction step.
         do 110 i = 1, kd 
            actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 110     continue

         do 140 ii = 1, nsmthg
c           Call matvc2 to compute wk = (A_0 + cA_x0 + dA_y0 + fM_0)*actnB.
            call matvc2(k, actnB(jd(k)), wk, icdomn, ichoic, c, d, f)

c           Residual: wk = rhs - wk =  rhs - (A_0 + cA_x0 + dA_y0 + fM_0)*u.
            do 120 i = 1, kd 
               wk(i) = rhs(jd(k)-1+i) - wk(i)
 120        continue

c           Postsmoothing by the symmetric Gauss-seidel smoother.
            call gausei(k, wk, wk, icdomn, ichoic, c, d, f, ifbs2)
c           call bpxin(k, wk, icdomn)

c           Final correction step.
            do 130 i = 1, kd 
               actnB(jd(k)-1+i) = actnB(jd(k)-1+i) + wk(i)
 130        continue
 140     continue
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
c      itmax = 2000                     !max. iteration
      itmax = 1
      ld = 2**l + 1
      m1 = ld - 1
      m2 = m1*m1
      xnn   = 0.d0

      do 10 i = 1, n
	 x(i)  = 0.d0
         wk(i) = 0.d0
 10   continue

      do 30 j = 1,ld
	 do 20 i = 1,ld
            xm(i,j) = 0.d0
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

	 if (xnn .lt. 1.d-40) then
	    if (xon .lt. 1.d-40) then
cc	       go to 100
	    else
	       stopl = 1.d0
	    end if
	 else
	    stopl  = xdiffn/xnn
	 end if

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
      dimension b(n), x(n), wk(n), xm(513, 513)

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

 300  format(' Level, no. of iters, relative, residual at GSUPWD: ',
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
      dimension b(n), x(n), wk(n), xm(513, 513)

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
	    print 300, l, kk, stopl, rnorm
	    go to 333
	 end if
 200  continue

 300  format(' Level, no. of iters, relative, residual at GSUPWD: ',
     &       i2, i5, 2e10.3)

      if (itmax .le. 100) go to 555
      print*, ' Iteration limit exceeded:', itmax
 333  nit = min(kk, itmax)
      print*

 555  return
      end
c=======================================================================
 
      subroutine bcgstr(j, n, b, x, tol, nit, icdomn, ipde,c,d,f) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the non-symmetric/indefinite linear system Ax = b 
c     by using Bi-CGSTAB (bi-conjugate gradient stabilized) method.
c     lim = 1089(5), 4225(6), 16641(7), 66049(8), 263169(9), 1050625(10)
c       irst = 0 : start with initial x = 0 
c       irst = 1 : restart with updated x.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim = 4225) 
      implicit double precision (a-h, o-z)
      dimension x(n), b(n), r(lim), p(lim), rs(lim), ap(lim)
      dimension t(lim), ax(lim)

      tol  = 1.d-20/(4**j)       !stopping criterion for iterations
      itmax = 1000               !max iterations
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
         if (k .ge. 200) print*, 'omega', omega

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
	 if (k .ge. 200) print*, 'tau1', tau

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
	 if (k .ge. 200) print*, 'tau2', tau

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

	 if (xnn .lt. 1.d-30) then
	    if (xon .lt. 1.d-30) then
	       go to 333
	    else
	       stopl = 1.d0
	    end if
	 else
	    stopl  = xdn/xnn
	 end if

c        To compute r = r - omega*t.
         call saxpy(-omega, t, r, n, r)

         rnorm = dot(r, r, n)
	 if (k .ge. 200) print*, ' rnorm', rnorm

c         print*, 'x = ', (jj, '*', x(jj), jj=1,n)
         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 200  continue
         
 333  if (irst .eq. 0)  nit1 = min(k, itmax)
      if (irst .eq. 1)  nit = kt + k - 1

      print 444, j, nit
 444  format(' Number of iterations in Bi-CGSTAB at level ',i2,':',i5)
      print*

      return
      end
c=======================================================================

      subroutine bcgst(j, n, b, x, tol, nit, icdomn, ipde, c, d, f) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the non-symmetric/indefinite linear system Ax = b 
c     by using Bi-CGSTAB (bi-conjugate gradient stabilized) method.
c     lim = 1089(5), 4225(6), 16641(7), 66049(8), 263169(9), 1050625(10)
c       ist = 0 : start with initial x = 0 
c       ist = 1 : restart with updated x.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim = 4225) 
      implicit double precision (a-h, o-z)
      dimension x(n), b(n), r(lim), p(lim), rs(lim), ap(lim), t(lim)

      tol  = 1.d-8/(4**j)       !stopping criterion for iterations
      itmax = 5000              !max iterations
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
c	print*
c	if (k .ge. 100) print*, 'omega', omega

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
c	 if (k .ge. 100) print*, 'tau1', tau

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
c	 if (k .ge. 100) print*, 'tau2', tau

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
	    if (xon .lt. 1.d-40) then
	       go to 333
	    else
	       stopl = 1.d0
	    end if
	 else
	    stopl  = xdn/xnn
	 end if

c        To compute r = r - omega*t.
         call saxpy(-omega, t, r, n, r)

         rnorm = dot(r, r, n)
c	 print*, ' rnorm', rnorm

c         print*, 'x = ', (jj, '*', x(jj), jj=1,n)
         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 200  continue

      print*, ' Iteration limit exceeded:', itmax
         
 333  nit1 = min(k, itmax)
c      print 444, j, nit1, rnorm
 444  format(' No. of iters in Bi-CGSTAB at level ',i2,':',i5,e10.3)
c      print*

      return
      end
c=======================================================================
 
      subroutine bcgsto(j, n, b, x, tol, nit, icdomn, ipde, c, d, f) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the non-symmetric/indefinite linear system Ax = b 
c     by using Bi-CGSTAB (bi-conjugate gradient stabilized) method.
c     lim = 1089(5), 4225(6), 16641(7), 66049(8), 263169(9), 1050625(10)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim = 4225) 
      implicit double precision (a-h, o-z)
      dimension x(n), b(n), r(lim), p(lim), rs(lim), ap(lim), t(lim)

      tol  = 1.d-14/(4**j)       !stopping criterion for iterations
      itmax = 1000               !max iterations
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

	 if (dabs(omega) .le. 1.d-40) then
	    print*, 'Bi-CGSTAB breakdown, omega is too small'
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

	 if (dabs(tau) .le. 1.d-40) then
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

	 if (dabs(tau) .le. 1.d-40) then
	    print*, 'Bi-CGSTAB breakdown, tau#2 is too small'
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

	 if (xnn .lt. 1.d-40) then
	    if (xon .lt. 1.d-40) then
	       go to 333
	    else
	       stopl = 1.d0
	    end if
	 else
	    stopl  = xdn/xnn
	 end if

c        To compute r = r - omega*t.
         call saxpy(-omega, t, r, n, r)

         rnorm = dot(r, r, n)
c	 print*, ' rnorm', rnorm

c         print*, 'x = ', (jj, '*', x(jj), jj=1,n)
cc         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 200  continue
         
 333  nit1 = min(k, itmax)
c      print 444, j, nit1
c 444  format(' Number of iterations in Bi-CGSTAB at level ',i2,':',i5)
c      print*

      return
      end
c=======================================================================

      subroutine pcgm(level, j, n, a, ax, ay, am, jcn, jre, jdiag, b, 
     &                xnew, conode, niter, ichoic, jchoic, c, d, f) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the linear system A_h*x = b by using preconditioned
c     conjugate gradient method in the unstructured grid, where A_h is stored in
c     the data structure a.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (lim = 263170)
      implicit double precision (a-h, o-z)
      dimension jcn(n*15), jre(n), a(n*15), b(n), conode(2,n), xnew(n),
     &          r(lim), p(lim), w(lim), Br(lim), xold(lim), jdiag(n),
     &          ax(n*15), ay(n*15), am(n*15)
      equivalence (w, br)

      tol  = 1.d-18/(4**level)        !tol -- stopping criteria for iterations
      itmax = 100                     !max. iteration
      xnewn  = 0.d0
      rdtBrn = 0.d0
      stopl  = 1.d0

      do 10 i = 1, n
	 xnew(i) = 0.d0
	 p(i)    = 0.d0
         w(i)    = 0.d0
	 r(i)    = b(i)
 10   continue

c*    Call function DOT to compute inner product of R and R. 
      rnorm  = dot(r, r, n)

c*    Iterate until convergence:  u = u + B(b - (A_0 + cA_x0 + dA_y0 + fM_0)u).

      do 100 kk = 1, itmax
c*       Call PRECND to compute B*r  qqqq
	 call precnd(r, Br, j, conode, n, a, jcn, jre, jdiag, 
     &               ichoic, jchoic, c, d, f)

	 rdtBro = rdtBrn
	 rdtBrn = dot(Br, r, n)

	 if (kk.eq.1) then
            do 30 i = 1, n
	       p(i) = Br(i)
 30         continue
	 else
	    beta = rdtBrn / rdtBro
c*          Call SAXPY to compute P = B*r + beta*P.
	    call saxpy(beta, p, Br, n, p)
	 end if

c*       Call MATVEC to compute W = a*P (=A_h*P).
	 call matvec(a, ax, ay, am, jcn, jre, n, p, w, 1, c, d, f)

c*       Call function DOT to compute inner product of P and W.
         alpha = rdtBrn/dot(p, w, n)

         do 40 i = 1, n
	    xold(i) = xnew(i)
 40      continue

c*       Call SAXPY to compute Xnew = Xold + alpha*P.
         call saxpy(alpha, p, xold, n, xnew)

         xoldn = xnewn
	 xnewn = dot(xnew, xnew, n)

         do 50 i = 1, n
            xold(i) = alpha * p(i)
 50      continue
         xdiffn = dot(xold, xold, n)
         stopl  = xdiffn/xnewn

c*       Call SAXPY to compute R = R - alpha*W.
         call saxpy(-alpha, w, r, n, r)
         rnorm = dot(r, r, n)

         print 300, level, kk, stopl, rnorm

         if (stopl .lt. tol .and. rnorm .lt. tol)  go to 333
 100  continue

      print*, ' Iteration limit exceeded:', itmax

 300  format(' Level, no. of iters, relative, residual at PCGM:',
     &       i2, i5, 2e10.3)

 333  niter = min(kk, itmax)

      return
      end
c=======================================================================
    
      subroutine precnd(r, p, j, conode, n, a, jcn, jre, jdiag, 
     &                  ichoic, jchoic, c, d, f)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes p = B'*r, where B' is a preconditioner of the
c     form B' = R + TBT^t, where B is a multilevel nodal basis preconditioner
c     (BPX), T is the interpolation matrix from structured mesh to unstructured
c     mesh, and R is a smoother such as Richardson, Jacobi, or symmetric
c     Gauss-Seidel type.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter(length = 263170, lim = 263170)
      implicit double precision (a-h, o-z)
      dimension r(n), p(n), conode(2,n),z(length),zt(length),rs(lim)
      dimension a(n*15), jcn(n*15), jre(n), jdiag(n)

      do 10 i = 1, n
         p(i) = 0.d0
   10 continue

c*    Call actntt to compute z = (T^t)*r
      call actntt(conode, n, j, r, z, ichoic)

      if (jchoic .eq. 1) then
c*       Call bpxg to compute zt = B*z, where B is the BPX preconditioner 
c        with a symmetric Gauss-Seidel smoother.
	 call bpxg(j, z, zt, ichoic)

      else if (jchoic .eq. 2) then
c*       Call pcgbpx to compute zt = (A_h)^-1*z, where A_h is the stiffness
c        matrix of the structured grid inside the domain.
     	 call pcgbpx(j, z, zt, ichoic, c, d, f)
      end if

c*    Call actnt to compute p = T*zt
      call actnt(conode, n, j, zt, p, ichoic)

c*    Call Gause3 to compute the action of symmetric Gauss-Seidel smoother:
c                  rs = (D-U)^(-1)*D*(D-L)^(-1)*r.
      call gause3(a, jcn, jre, jdiag, n, r, rs)

      do 60 i = 1, n
c*        R = I : Richardson smoother
c	  p(i) = p(i) + r(i)
c*        R = D^(-1) : Jacobi smoother
c	  p(i) = p(i) + r(i)/a(jdiag(i))
c*        R = (D-U)^(-1)*D*(D-L)^(-1) : symmetric Gauss-Seidel smoother
	  p(i) = p(i) + rs(i)
   60 continue

      return
      end
c=======================================================================

      subroutine pcgbpx(j, b, xnew, ichoic, c, d, f)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine solves the linear system A*x = b using preconditioned
c     conjugate gradient method, where A is is the stiffness matrix of the
c     structured grid inside the domain and a BPX preconditioner is used.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (length = 263170)
      implicit double precision (a-h, o-z)
      dimension  b((2**j+1)**2), xnew((2**j+1)**2), xold(length),
     &           r(length), p(length), w(length), Br(length)
      equivalence (br, w)

c*    tol -- stopping criteria for iterations
      tol  = 1.d-6/(16**j)
      n    = (2**j+1)**2

      k = 0
      do 10 i = 1, n
	 xnew(i) = 0.d0
	 p(i)    = 0.d0
         w(i)    = 0.d0
	 r(i)    = b(i)
   10 continue

c*    Call function DOT to compute inner product of R and R.
      rnorm  = dot(r, r, n)
      xnewn  = 0.d0
      rdtBrn = 0.d0
      stopl  = 1.d0

   20 continue
c      if (k.eq.4) go to 333

      if (stopl .gt. tol .and. rnorm .gt. tol) then
c*       Call BPXG to compute B*r  
	 call bpxg(j, r, Br, ichoic)
	 rdtBro = rdtBrn
	 rdtBrn = dot(Br, r, n)
         
	 if (k .ge. 200) then
	    print*, 'Over 200 iterations, DIVERGENT!!!'
	    goto 333 
	 end if

	 k = k+1
	 if (k.eq.1) then
            do 30 i = 1, n
	       p(i) = Br(i)
   30       continue
	 else
	    beta = rdtBrn / rdtBro
c*          Call SAXPY to compute P = B*r + beta*P.
	    call saxpy(beta, p, Br, n, p)
	 end if

c*       Call MATVC2 to compute W = (A_0)*P.
	 call matvc2(j, p, w, ichoic, 1, c, d, f)
c*       Call function DOT to compute inner product of P and W.
         alpha = rdtBrn/dot(p, w, n)

         do 40 i = 1, n
	    xold(i) = xnew(i)
   40    continue

c*       Call SAXPY to compute Xnew = Xold + alpha*P.
         call saxpy(alpha, p, xold, n, xnew)

         xoldn = xnewn
	 xnewn = dot(xnew, xnew, n)

	 if (xnewn .lt. 1.d-16) then
	    if (xoldn .lt. 1.d-16) then
	       go to 333
	    else
	       stopl = 1.d0
	    end if
	 else
	    do 50 i = 1, n
 	       xold(i) = alpha * p(i)
   50	    continue
	    xdiffn = dot(xold, xold, n)
	    stopl  = xdiffn/xnewn
	 end if

c*       Call SAXPY to compute R = R - alpha*W.
         call saxpy(-alpha, w, r, n, r)
         rnorm = dot(r, r, n)
	 go to 20
      end if

  333 niter3 = k
      print*, ' Number of iterations in pcgbpx: ', niter3

      return
      end
c=======================================================================

      subroutine bpxg(j, x, y, ichoic)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes y = B*x, where B is a multilevel nodal basis
c     preconditioner (BPX) with a symmetric Gauss-seidel smoother.
c         nsmthg = the number of smoothings, 1, 2, or 3.
c         nwk = 400000, 1402210, depending on j = 9, 10.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	

      parameter( nsmthg = 2, ifbs = 3, nwk = 400000)
      implicit double precision (a-h, o-z)
      dimension x((2**j+1)**2), y((2**j+1)**2), wk(nwk), jd(11)

      jl = (2**j+1)**2

      do i = 1, jl
         y(i) = 0.d0
      end do

      do 10 i = 1, nwk
        wk(i) = 0.d0
   10 continue

      jd(1) = 1
      do 20 i = 1, j-1
	 jd(i+1) = jd(i) + (2**i+1)**2
   20 continue

      do 30 i = 1, jl
         wk(jd(j)-1+i) = x(i)
   30 continue

      do 40 l = j-1, 1, -1
	 call restrn(l+1, wk(jd(l+1)), wk(jd(l)))
c*       Call BPXIN to apply the BPX preconditioner inside the domain only.
	 call bpxin(l, wk(jd(l)), ichoic)
   40 continue

      do 50 i = 1, nsmthg
	 call gausei(1, wk(1), wk(1), ichoic, 1, 0.d0, 0.d0, 0.d0, ifbs)
c	 call bpxin(1, wk(1), ichoic)
   50 continue

      do 70 l = 2, j
         jl  = (2**l+1)**2
	 call prolng(l, wk(jd(l-1)), y)
	 do 55 i = 1, nsmthg
	    call gausei(l, wk(jd(l)), wk(jd(l)), ichoic, 1,
     &                  0.d0, 0.d0, 0.d0, ifbs)
c	    call bpxin(l, wk(jd(l)), ichoic)
   55    continue
         do 60 i = 1, jl
            wk(jd(l)-1+i) = y(i) + wk(jd(l)-1+i)
   60    continue
	 call bpxin(l, wk(jd(l)), ichoic)
   70 continue
	  
      do 120 i = 1, jl
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
      dimension p((2**l+1)**2), q((2**(l-1)+1)**2), z(513, 513)

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
	    if (l .le. 3) q(k) = q(k)/4
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
      dimension x((2**l+1)**2), y((2**l+1)**2), xm(513, 513)

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
  100    continue
  110 continue

  120 return
      end
c=======================================================================

      subroutine gause3(a, jcn, jre, jdiag, n, b, x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                          -1          -1  
c     This subroutine computes x = (D - U )  D (D - L )  * b, i.e., the
c                                    h   h    h  h   h 
c     action of the symmetric Gauss-Seidel smoother, where A  = D - L - U
c                                                           h    h   h   h
c     is the stiffness matrix of the unstructured mesh on the given domain.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	

      implicit double precision (a-h, o-z)
      dimension  a(n*15), jcn(n*15), jre(n), jdiag(n), b(n), x(n)

      do 10 i = 1, n
         x(i)  = b(i)
   10 continue

c*              << Computation of x = (D_l - L_l)^-1 * x >>

c*    For i = 1,2,...,n

      do 30 i = 1, n
	 jchk = (i-1)*15 + 1
	 do 20 j = jchk, jdiag(i) - 1 
            x(i) = x(i) - a(j)*x(jcn(j))
   20    continue
	 x(i) = x(i)/a(jdiag(i))
   30 continue

c*             << Computation of x = D_l * x >>

      do 40 i = 1, n 
         x(i) = x(i)*a(jdiag(i))
   40 continue

c*              << Computation of x = (D_l - U_l)^-1 * x >>

c*    For i = n,...,2,1

      do 80 i = n, 1, -1
	 jchk = (i-1)*15 + jre(i)
	 do 70 j = jdiag(i) + 1, jchk
   	    x(i) = x(i) - a(j)*x(jcn(j))
   70    continue
	 x(i) = x(i)/a(jdiag(i))
   80 continue

      return
      end
c=======================================================================

      subroutine actnt(conode, n, jlevel, w, z, ichoic)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes the action z = Tw, where T is the matrix
c     representation of the interpolation from structured mesh to unstructured
c     mesh.
c     Input:  conode(j,i) - jth coordinates of the node number i
c	      n           - total number of nodes on the unstructured mesh
c	      jlevel      - level of the structured mesh
c	      w           - n_0 vector
c     Output: z           - T*w, where T is the matrix representation of the
c                           interpolation from structured mesh to unstructured
c                           mesh
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension conode(2,n), w((2**jlevel+1)**2), z(n)

      eps = 1.d-5
      m   = 2**jlevel + 1
      m1  = m-1
c     h   = 1.d0/m1

c*    Initialization

      do 10 j = 1,n
	 z(j) = 0.d0
 10   continue

c*    For level = 1

      if (jlevel .eq. 1 )   go to 60
     
c*    For level > 1

c*    w(k) = 0 outside the domain considered(Optional call).
      call bpxin(jlevel, w, ichoic)

c*    Computing z = T*w

      do 50 i = 1, n
	 xco = conode(1,i)
	 yco = conode(2,i)

	 if (xco .lt. eps .or. yco .lt. eps) go to 50

	 if (ichoic .eq. 1) then
	    if (xco*xco + yco*yco .gt. 1.d0 - eps ) go to 50
	 else if (ichoic .eq. 3 .or. ichoic .eq. 4) then
	    if (xco .gt. (1.d0-eps) .or. yco .gt. (1.d0-eps)) go to 50
	 end if

	 ix   = int(xco*m1)
	 iy   = int(yco*m1)
	 k    = iy*m + ix + 1
	 xt   = m1*xco - ix
         yt   = m1*yco - iy
         psi1 = xt          * yt
         psi2 = xt          * (1.d0 - yt)
         psi3 = (1.d0 - xt) * yt
         psi4 = (1.d0 - xt) * (1.d0 - yt)

         z(i) = w(k)*psi4 + w(k+m)*psi3 + w(k+1)*psi2 + w(k+m+1)*psi1
   50 continue

   60 return
      end
c=======================================================================

      subroutine actntt(conode, n, jlevel, w, z, ichoic)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes the action z = (T^t)*w, where T is the matrix
c     representation of the interpolation from structured mesh to unstructured
c     mesh.
c     Input:  conode(j,i) - jth coordinates of the node number i
c	      n           - total number of nodes on the unstructured mesh
c	      jlevel      - level of the structured mesh
c	      w           - n vector
c     Output: z           - (T^t)*w, where T is the matrix representation of the
c                           interpolation from structured mesh to unstructured
c                           mesh
c     h  - mesh size of the structured mesh    
c     n0 - total number of nodes on the structured mesh
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
      implicit double precision (a-h, o-z)
      dimension conode(2,n), w(n), z((2**jlevel+1)**2)

      eps = 1.d-5
      m   = 2**jlevel + 1
      m1  = m - 1
      n0  = m*m
c     h   = 1.d0/(m-1)

c*    Initialization

      do 10 j = 1,n0
	 z(j) = 0.d0
   10 continue

c*    For level = 1
      if (jlevel .eq. 1)  go to 60

c*    For level > 1

c*    Computing z = (T^t)*w

      do 30 i = 1, n
	 xco = conode(1,i)
	 yco = conode(2,i)

	 if (xco .lt. eps .or. yco .lt. eps) go to 30

	 if (ichoic .eq. 1) then
	    if (xco*xco + yco*yco .gt. 1.d0 - eps ) go to 30
	 else if (ichoic .eq. 3 .or. ichoic .eq. 4) then
	    if (xco .gt. (1.d0-eps) .or. yco .gt. (1.d0-eps)) go to 30
	 end if
	 
	 ix   = int(xco*m1)
	 iy   = int(yco*m1)
	 k    = iy*m + ix + 1
	 xt   = m1*xco - ix
         yt   = m1*yco - iy
         psi1 = xt          * yt
         psi2 = xt          * (1.d0 - yt)
         psi3 = (1.d0 - xt) * yt
         psi4 = (1.d0 - xt) * (1.d0 - yt)

         z(k)     = z(k)     + w(i)*psi4 
	 z(k+m)   = z(k+m)   + w(i)*psi3
	 z(k+1)   = z(k+1)   + w(i)*psi2 
	 z(k+m+1) = z(k+m+1) + w(i)*psi1 
   30 continue

c*    z(k) = 0 outside the domain considered.
      call bpxin(jlevel, z, ichoic)

   60 return 
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

      subroutine matvec(a, ax, ay, am, jcn, jre, n, p, w, kchoic, 
     &                  c, d, f)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes matrix_vector multiplication(W=B*p) in the
c     unstructured grid, where B is defined as follows:     
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
      dimension  a(n*15), ax(n*15), ay(n*15), am(n*15), jcn(n*15), 
     &           jre(n), p(n), w(n)

      do 10 i = 1, n
	 w(i) = 0.d0
 10   continue

      if (kchoic.ne.4 .and. kchoic.ne.5 .and. kchoic.ne.6) then
         do 30 i = 1, n
            jchk = (i-1)*15
            do 20 ii = 1, jre(i)
               jpos = jchk + ii
               w(i) = w(i) + a(jpos)*p(jcn(jpos))
 20         continue
 30      continue
      end if

      if (kchoic.ne.1 .and. kchoic.ne.6 .and. kchoic.ne.7) then
         do 50 i = 1, n
            jchk = (i-1)*15
            do 40 ii = 1, jre(i)
               jpos = jchk + ii
               w(i) = w(i) + (c*ax(jpos) + d*ay(jpos))*p(jcn(jpos))
 40         continue
 50      continue
      end if

      if (kchoic.ne.1 .and. kchoic.ne.2 .and. kchoic.ne.4) then
         do 70 i = 1, n
            jchk = (i-1)*15
            do 60 ii = 1, jre(i)
               jpos = jchk + ii
               w(i) = w(i) + f*am(jpos)*p(jcn(jpos))
 60         continue
 70      continue
      end if

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
      dimension x((2**j+1)**2), w((2**j+1)**2), xm(513, 513)

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
      dimension x((2**j+1)**2), w((2**j+1)**2), xm(513, 513)

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

      return
      end
c=======================================================================

      subroutine grflab(lu, ls, conode, n, x, icdomn, icalgm)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine displays the various 3D views of the computed solution 
c     using a MATLAB tool, 'plot3'. In addition, if the given(original) mesh is
c     uniform on the unit square, then better graphic tools such as 'mesh' and 
c     'surf' in MATLAB can be used to display the computed solution beautifully.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension conode(2,n), x(n)

c      print*,' Solution at each node 1: ', (i,'*',x(i),i=1,n)

      k  = 2**ls
      k1 = k - 1
      k2 = k + 1
      kk = k2*k2
      h  = 1.d0/k

      if (icalgm .eq. 0)   go to 15

      open(unit = 13, file = 'nsasol.m', status = 'unknown')
      write(13,*) '% Level =', lu
      write(13,*) 'clear'
      do 10 jj = 1, n
         write(13,*) 'x(', jj, ') = ', conode(1, jj), ';'
         write(13,*) 'y(', jj, ') = ', conode(2, jj), ';'
         write(13,*) 'z(', jj, ') = ', x(jj), ';'
 10   continue
      write(13,*) 'solgrf'

 15   if (icdomn .eq. 3 .or. icdomn .eq. 4) then
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
           
      close (13)
      close (14)

      return
      end
c=======================================================================

      subroutine exctxb(n, a, ax, ay, am, jcn, jre, conode, nodedb,nodb,
     &                  soludb, exsoln, b, icdomn, c, d, f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     To supply an exact FE solution and a correspoding load vector in order
c     to test the algorithms. (This is a temporary process.)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension exsoln(n), b(n),  conode(2,n), a(n*15), nodedb(nodb),
     &          jcn(n*15), jre(n), am(n*15), ax(n*15), ay(n*15),
     &          soludb(nodb)

c     Find an exact FE solution.
      do 20 i = 1,n
	 x = conode(1,i)
	 y = conode(2,i)
         exsoln(i) = extsol(x, y, icdomn)
 20   continue

      do 30 i = 1, nodb
	 exsoln(nodedb(i)) = soludb(i)
 30   continue

c     Call MATVEC to compute b = (A_0 + cA_0x + dA_0y + fM_0)*exsoln
      call matvec(a, ax, ay, am, jcn, jre, n, exsoln, b, 3, c, d, f)

      return
      end
c=======================================================================

      subroutine errors(j, a, ax, ay, am, jcn,jre, n, xv, exsoln,
     &           niter, cpu, icdomn, c, d, f, itbl, dtbl, ui, errloc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This subroutine computes (discrete) maximum, L^2, H^1, and H^1 seminorm 
c     error estimates between computed solution(u_h) and interpolation(u_I) of 
c     the exact solution if exact solution is known. User should describe the 
c     exact solution, u(x,y), before running the program(See function EXTSOL).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
      dimension xv(n), ui(n), errloc(n), a(n*15), exsoln(n), itbl(10,2),
     &          jcn(n*15),jre(n),am(n*15),ax(n*15),ay(n*15), dtbl(10,6)

      h     = 1.d0/(2**j)
      errmx = 0.d0
      max   = 0
      xmax  = 0.d0
      ymax  = 0.d0
      ersum = 0.d0

c     Computing maximum error of u_e - u_h, where u is the exact solution, u_e
c     FE solutioin, and u_h the computed solution.
      do 20 i = 1,n
	 errloc(i) = exsoln(i) - xv(i)
         errabs = dabs(errloc(i))
         if (errabs .gt. errmx) then
	    errmx = errabs
cc	    max = i
cc	    xmax = x
cc	    ymax = y
         end if
 20   continue

c     Computing L^2 error estimate of u_I - u_h.
      do 60 i = 1, n
         ersum = ersum + errloc(i)*errloc(i)
 60   continue
      erl22t = h*h*ersum
      erl2t  = dsqrt(erl22t)

c     Computing L^2 norm error estimate of u_e-u_h, sqrt(errloc^t*M*errloc).
c     Call MATVEC to compute ui = M_h*errloc
      call matvec(a, ax, ay, am, jcn, jre, n, errloc, ui, 6, 
     &            0.d0, 0.d0, 1.d0)
      erl22 = dot(errloc, ui, n)
      erl2  = dsqrt(erl22)

c     Computing H^1 seminorm error estimate of u_e-u_h, sqrt(errloc^t*A*errloc).
c     Call MATVEC to compute ui = A_h*errloc
      call matvec(a, ax, ay, am, jcn, jre, n, errloc, ui, 1, 
     &            0.d0, 0.d0, 0.d0)
      erh1s2 = dot(errloc, ui, n)
      erh1s  = dsqrt(erh1s2)

c     Computing H^1 error estimate of u_I - u_h.
      erh1 = dsqrt(erl22 + erh1s2)

      itbl(j, 1) = j
      itbl(j, 2) = niter
      dtbl(j, 1) = errmx
      dtbl(j, 2) = erl2t
      dtbl(j, 3) = erl2
      dtbl(j, 4) = erh1s
      dtbl(j, 5) = erh1
      dtbl(j, 6) = cpu

      print* 
      print*,' level  niter   maximum    L^2(t)      L^2      H^1 semi
     &     H^1       CPU'
      print*,' --------------------------------------------------------
     &-------------------'
      do 70 i = 3, j
         print 100, i, itbl(i,2), dtbl(i,1), dtbl(i,2), dtbl(i,3),
     &             dtbl(i,4), dtbl(i,5), dtbl(i,6)
70    continue
      print*, ' ----------------------------------------
     &-----------------------------------'
      print*

100   format(i5, i7, e13.4, e11.4, e11.4, e11.4, e11.4, f8.1)

      return
      end
c=======================================================================

      double precision function extsol(x, y, ichoic)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Exact solution of the PDE if it is known.  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h, o-z)
 
      pi = dacos(-1.d0)

      if (ichoic .eq. 1) then
	 iex = 1
	 if (iex .eq. 1) then
c           Example 1
	    extsol = 20*(1.d0-x*x-y*y)*x*y
	 else if (iex .eq. 2) then
c           Example 2
	    extsol = dsin(5*x*y*pi) - dexp(y)
	 else if (iex .eq. 3) then
c	    Example 3
	    extsol = - dsin(x)*(dcos(y) - 1.d0 + 2*y/pi)
	 else if (iex .eq. 4) then
c	    Example 4
	    extsol = 5*(x - x**3 - x*y*y)*dsin(3*pi*y)
	 end if

      else if (ichoic .eq. 3 .or. ichoic .eq. 4) then
c*       Example 1
	 extsol = 100*(x-x*x)*(y-y*y)

c*       Example 2
c	 extsol = dsin(5*x*y*pi) - dexp(y)

c*       Example 3
c	 extsol = - dsin(pi*x)*(dcos(pi*y) - 1.d0 + 2*y)

c*       Example 4
cc	 extsol =  10*dsin(4*pi*x)*(dexp(y) - dexp(y*y))

c*       Example 5
c	 extsol = 5*(x - x**3 - x*y*y)*dsin(3*pi*y)

      end if

      return
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
	 iex = 1
	 if (iex .eq. 1) then
c           Example 1
	    fxy = 12*x*y + c*(1 - 3*x*x - y*y)*y + d*(1 - 3*y*y - x*x)*x
     &                + f*(1 - x*x - y*y)*x*y
	    fxy = 20*fxy
	 else if (iex .eq. 2) then
c           Example 2
	    fxy = (25*pi*pi*dsin(5*x*y*pi))*(x*x+y*y) + dexp(y)
	 else if (iex .eq. 4) then
c           Example 4      
	    ex = x - x**3 - x*y*y 
	    ed = 1 - 3*x*X - y*y 
	    sn = dsin(3*pi*y)
	    cs = dcos(3*pi*y)
	    fxy = (40*x + 45*pi*pi*ex + c*5*ed - 10*d*x*y + 5*f*ex)*sn
     &          + (60*pi*x*y + 15*d*pi*ex)*cs
	 end if

      else if (ichoic .eq. 3 .or. ichoic .eq. 4) then
c*       Example 1
	 fxy = 2*(x - x*x + y - y*y) 
     &       + c*(1 - 2*x)*(y - y*y) + d*(x - x*x)*(1 - 2*y)
     &       + f*(x - x*x)*(y - y*y)
	 fxy = 100*fxy

c*       Example 2
c	 fxy = (25*pi*pi*dsin(5*x*y*pi))*(x*x+y*y) + dexp(y)

c*       Example 3
c	 fxy = - pi*pi*dsin(pi*x)*(2*dcos(pi*y) - 1.d0 + 2*y)

c*       Example 4
cc	 p2 = pi*pi
cc	 sx = dsin(4*pi*x)
cc	 cx = dcos(4*pi*x)
cc	 ey = dexp(y)
cc	 eys = dexp(y*y)
cc	 ey1 = ey - eys
cc	 ey2 = ey - 2*y*eys
cc	 ey3 = ey - (2 + 4*y*y)*eys
cc	 fxy = (16*p2*ey1-ey3)*sx + c*4*pi*ey1*cx + d*ey2*sx + f*ey1*sx
cc	 fxy = 10*fxy

c        Example 5
c	 ex = x - x**3 - x*y*y 
c	 ed = 1 - 3*x*X - y*y 
c        sn = dsin(3*pi*y)
c        cs = dcos(3*pi*y)
c	 fxy = (40*x + 45*pi*pi*ex + c*5*ed - 10*d*x*y + 5*f*ex)*sn
c     &         + (60*pi*x*y + 15*d*pi*ex)*cs
c
      end if

      return
      end
c=======================================================================
c     END of PROGRAM                                                   c
c=======================================================================
