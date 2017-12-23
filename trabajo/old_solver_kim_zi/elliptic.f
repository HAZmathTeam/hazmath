C=====================================================================
      program FEM_MG
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension m(20 000 000), r(10 000 000)
      character*100 meshf(0:20), mtrxf(20)
      common /quad/ wg(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /quad_gs_ntrl/ z_ntrl(7),zw(5,7)
      common /quad_gs/ w(7),z(7)
      common /menu/ ich(200)
      common /fname/ meshf,mtrxf
C---------------------------------------------------------------------
C...  This program solves the 2nd order elliptic linear PDEs of the
C...  form
C...
C...     - div (alpha(x,y) grad(u)) + gamma(x,y) u = f(x,y).
C...
C...  where ALPHA is a symmetric 2x2 matrix. Dirichlet, Neumann, or
C...  Robin type boundary condition is imposed. The standard Galerkin 
C...  finite element is solved by Multigrid.
C...
C...  Parameter:
C---------------------------------------------------------------------
      call input
      call quad_elt
      call quad_data
      call quad_data1
C 
      mt     = 16
      lstiff = ich(1)
      jdf    = ich(3)
      lread  = ich(5)
      lfmg   = ich(36)
      lcmg   = ich(37)
C
      write(*,10) lstiff,lread,lfmg,lcmg
 10   format(2x,' lstiff = ', i7 / 2x,
     >          ' lread  = ', i7 / 2x,
     >          ' lfmg   = ', i7 / 2x,
     >          ' lcmg   = ', i7)
C
 20   continue
C
      lf = lfmg
C
      go to (100, 200, 200) lstiff
C
 100  continue
      if(lread .eq. 1) then
         open(mt,file=meshf(lf),status='unknown',form='formatted')
         read(mt,*) n1,n2
         read(mt,*) nel,n,ned
      else
         open(mt,file=meshf(lf),status='unknown',form='unformatted')
         read(mt) n1,n2
         read(mt) nel,n,ned
      end if
      write(*,*) ' No. of elements, nodes, & bdry edges:',nel,n,ned
C
C...  Compute the addresses of each array in the arrays m and r.
C
      ie     = 1
      je     = ie + nel + 1
      idir   = je + nel*3
      ia     = idir + n + 5
      ja     = ia + n + 1
      inedg  = ja + 2*(nel + n + 5) + n
      iet    = inedg + ned * 2
      jet    = iet + n + 1
      jat    = jet + nel*3
      lasti  = jat + 2*(nel + n + 5) + n
C     
      kx     = 1
      ky     = kx + n
      ksol   = ky + n
      krhs   = ksol + n
      kra    = krhs + n
      krat   = kra + 2*(nel + n + 5) + n
      lastr  = krat + 2*(nel + n + 5) + n
C
      call rdmesh(lread,mt,m(ie),m(je),ned,m(inedg),
     >     m(idir),r(kx),r(ky),nel,n,jdf)
C     
      call iit(m(ie),m(je),nel,n,m(iet),m(jet))
C     
      nnz = 0
      call smbasg(m(ie),m(je),m(iet),m(jet),m(ia),m(ja),n,nnz,m(idir))
C
      iip = iet
      call mtrx_sym(
     I     nel,n,jdf,m(ie),m(je),r(kx),r(ky),
     O     m(ia),m(ja),r(kra),nnz,r(krhs),
     W     ned,m(inedg),m(idir),m(iip))
C
      go to 400
C
 200  continue
      if (lstiff .eq. 2) then
         open(mt,file=mtrxf(lf),status='unknown',form='formatted')
         read(mt,*) n
      else if (lstiff .eq. 3) then
         open(mt,file=mtrxf(lf),status='unknown',form='unformatted')
         read(mt) n
      end if
C
      idir = 1
      ia   = idir + n
C
      if (lstiff .eq. 2) then
         read(mt,*) (m(ia+i), i=0,n)
      else if (lstiff .eq. 3) then
         read(mt) (m(ia+i), i=0,n)
      end if
C
      nnz   = m(ia+n) - 1
      ja    = ia + n + 1
      lasti = ja + nnz
C
      ksol  = 1
      krhs  = ksol + n
      kra   = krhs + n     
      lastr = kra + nnz
C
      nnz1  = nnz - 1
      n1    = n - 1
C
      if (lstiff .eq. 2) then
         read(mt,*) (m(ja+i), i=0,nnz1)
         read(mt,*) (r(kra+i), i=0,nnz1)
         read(mt,*) (r(krhs+i), i=0,n1)
      else if (lstiff .eq. 3) then
         read(mt) (m(ja+i), i=0,nnz1)
         read(mt) (r(kra+i), i=0,nnz1)
         read(mt) (r(krhs+i), i=0,n1)
      end if
C
      write(*,*) ' No. of nodes & no. of nonzero entries:',n,nnz
      go to 400
C
 400  continue
C     
C...  Solve by MG.
C
      img = ia
      kmg = kra
C
      call mg(m(idir),m(img),r(ksol),r(krhs),r(kmg),n,nnz,lf,lcmg)
C
C...  OUTPUT.
C
      if (lstiff .eq. 1 .and. n .lt. 2000) then
         call output(r(kx),r(ky),r(ksol),n)
         call output_sol(m(ie),m(je),r(kx),r(ky),r(ksol),n,nel)
      end if 
C
      close(mt)
      endfile(800)
      close(800)
      endfile(100)
      close(100)
      endfile(199)
      close(199) 
C
 990  stop
      end
C=====================================================================
