C=====================================================================
      program PTEST
C=====================================================================
      parameter (nlast = 5 000 000, nrlast = 3 000 000)
cc      parameter (nlast = 20 000 000, nrlast = 10 000 000)
cc      parameter (nlast = 60 000 000, nrlast = 30 000 000)
      implicit real*8(a-h,o-z)
      common /top/ m(nlast)              !To use the machine witt
      common /top1/ r(nrlast)            !To use the machine witt
cc      dimension m(nlast), r(nrlast)
      character*100 meshf(0:20)
      common /quad/ wg(7),quadpt(2,7),phi(3,7),phix(3,7),phiy(3,7)
      common /quad_gs_ntrl/ z_ntrl(7),zw(5,7)
      common /quad_gs/ w(7),z(7)
      common /menu/ ich(200)
      common /fname/ meshf
C---------------------------------------------------------------------
C...  This program tests something.
C---------------------------------------------------------------------
      n     = 0
      nnz   = 0
      nd    = 10000
C
      ia    = 1
      ja    = ia + nd
      ir    = ja + 20*nd
      ic    = ir + 20*nd
      ilast = ic + 20*nd
C
      kra   = 1
      kaij  = kra + 20*nd
      klast = kaij + 20*nd
C
      call mat_conv(m(ia),m(ja),r(kra),n,nnz,m(ir),m(ic),r(kaij))
      call outmat1(m(ia),m(ja),r(kra),n,nnz)

      write(*,*), 'ia****** ', (m(k), k = 1,n)
      write(*,*), 'ja****** ', (m(ja+k-1), k = 1,nnz)
      write(*,*), 'aa****** ', (r(kra+k-1), k = 1,nnz)

      stop
      end
C=====================================================================
