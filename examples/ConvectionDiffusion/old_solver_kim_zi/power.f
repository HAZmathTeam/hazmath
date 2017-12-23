C=====================================================================
      subroutine power(eigval, eigvec, j, n, a, jcn, jrs, au, jcnu, 
     &                 jrsu, na, conode, nodedb, nodb, ichoic, jchoic)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine computes the dominant eigenpair(eigenvalue, eigenvector)
c     of the matrix BA using the power method, where A is the stiffness matrix
c     of the unstructured mesh in the quarter circle and B is a BPX type
c     preconditioner.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      parameter(lim = 270000)
      implicit double precision (a-h, o-z)
      dimension eigvec(n), xo(lim), xn(lim), conode(2,n),a(na),jcn(na)
      dimension jrs(n+1), au(na), jcnu(na), jrsu(n+1), nodedb(nodb)
C
      tol    = 1.d-10
      nmax   = 9999
      k      = 1
      eig0   = 0.d0
      eig1   = 0.d0
      eigval = 3.d0

c*    Initial guess

      do 10 i = 1, n
         eigvec(i) = 1.d0
   10 continue  

      do 15 i = 1, nodb
	 eigvec(nodedb(i)) = 0.d0
   15 continue

      xnnorm = dsqrt(dot(eigvec,eigvec,n))

      do 20 i = 1, n
         eigvec(i) = eigvec(i)/xnnorm
   20 continue

   30 if (k .lt. nmax) then
c*       Call MATVEC to compute xo = A*eigvec
         call matvec(a, jcn, jrs, n, na, eigvec, xo)
c*       Call PRECND to compute xn = B*xo
         call precnd(xo, xn, j, conode, n, a, jcn, jrs, au, jcnu,
     &               jrsu, na, ichoic, jchoic)

	 eig2 = dot(eigvec, xn, n)
         eig3 = eigval

         if (dabs(eig2-2*eig1+eig0) .ge. 1.d-16)  then 
c*          Aitken acceleration to get the eigenvalue fast
            eigval = eig0 - (eig1-eig0)**2/(eig2-2*eig1+eig0)
         endif

         print 40,  k, eig2, eigval
   40    format('k, eigenvalues in power method :',i4,2(1x,f16.10))

         xnnorm = dsqrt(dot(xn,xn,n))
	 if (xnnorm .le. 1.d-16) then
	    print*, ' Eigenvalue 0, select new initial guess!!'
	    go to 60
         endif

         do 50 i = 1, n
            eigvec(i) = xn(i)/xnnorm
   50    continue

c         print*, 'eigvector(',k,') = ', (i, '*', eigvec(i), i=1,n) 

         if (k .ge. 20) then
c	    if (abs((eigval-eig3)/eigval) .lt. tol)  go to 60
	    if (abs((eig2-eig1)/eig2) .lt. tol)  go to 60
         end if

         k      = k + 1
         eig0   = eig1
         eig1   = eig2
         go to 30
      else
	 print*, 'Maximum number of iterations exceeded.'
      end if

   60 return
      end
C=====================================================================
      subroutine inv_pwr(eigval,m,ia,ja,idir,ifree,
     >     r,ksol,krhs,kra,krm,keig,kfree,
     >     n,nx,ny,nnz,nel,lf,lc,tol)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(1),r(1)
C=====================================================================
C...  Computes the smallest eigenpair(eigenvalue, eigenvector) of the
C...  form Au = s*Mu by using the inverse power method, where A is
C...  the stiffness matrix, M the mass matrix, u = eigvec, s = eigval.
C...  To solve A*y = x we use a V-cycle MG method.
C=====================================================================
cc      tol    = 1.d-12
      nmax   = 10*n
      k      = 1
      eig0   = 0.d0
      eig1   = 0.d0
      eigval = 1.d0
C
C...  Initial guess with zero on the Dirichlet boundary.
C
ccc      call ones(r(keig),n)
ccc      call zero_dir(r(keig),m(idir),n)
C
      do 15 i = 1,n
         p1 = m(ia+i-1)
         p2 = m(ia+i) - 1
         if(p2-p1 .lt. 0) stop 10
         if(p2-p1 .eq. 0) then
            r(keig+i-1) = 0.d0
ccc            r(kra+p1-1) = 1.d0
         else
           r(keig+i-1) = 1.d0
         end if
 15   continue
C
C...  To compute H^1 seminorm, i.e., \sqrt(eigvec^t*A*eigvec)
C
      call abyvg(m(ia),m(ja),r(kra),r(keig),n,r(krhs))
      call scpro(r(keig),r(krhs),eig_norm,n)
      eig_norm = dsqrt(eig_norm)
C
C...  Normalization.
C
      eig_norm_inv = 1.d0/eig_norm
      call usmultu(r(keig),n,eig_norm_inv)
C
      do while (k .lt. nmax)
C...     Compute sol = A^(-1)M*eigvec
C
         tol_mg = 1.d-13
         call abyvg(m(ia),m(ja),r(krm),r(keig),n,r(krhs))
         call mg_s(m(1),ia,ja,idir,ifree,r(1),ksol,krhs,kra,kfree,
     >        n,nx,ny,nnz,nel,lf,lc,tol_mg) 
C
	 call scpro(r(keig),r(ksol),eig2,n)
         call scpro(r(keig),r(keig),eig_l2,n)
         eig2 = eig2/eig_l2
         eig3 = eigval
C
         if (dabs(eig2-2*eig1+eig0) .ge. 1.d-16)  then 
C...        Aitken acceleration to get the eigenvalue fast.
            eigval = eig0 - (eig1-eig0)**2/(eig2-2*eig1+eig0)
         endif
C
         if (k .ge. 2) print 40,  k, 1/eig2, 1/eigval
 40      format('k, eigenvalues in inverse power method :', i6,
     >        2(3x,f13.10))
C
C...     To compute H^1 seminorm, i.e., \sqrt(sol^t*A*sol)
C
         call abyvg(m(ia),m(ja),r(kra),r(ksol),n,r(krhs))
         call scpro(r(ksol),r(krhs),eig_norm,n)
         eig_norm = dsqrt(eig_norm)
C
C...     Normalization.
C
	 if (eig_norm .le. 1.d-16) then
	    write(*,*) ' Eigenvector 0, select new initial guess!!'
	    stop 19
         endif
C
         eig_norm_inv = 1.d0/eig_norm
         call usmultv(r(keig),r(ksol),n,eig_norm_inv)
C
ccc         print*, 'eigvector(',k,') = ',(i,'*',r(keig+i-1), i=1,n)
C
         if (k .ge. 5) then
            if (dabs((eigval-eig3)/eigval) .lt. tol)  go to 60
ccc            if (dabs((eig2-eig1)/eig2) .lt. tol)  go to 60
         end if
C
         k      = k + 1
         eig0   = eig1
         eig1   = eig2
      end do
C
      write(*,*) ' Maximum number of iterations exceeded.'
C
  60  return
      end
C=====================================================================
