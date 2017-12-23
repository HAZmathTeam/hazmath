C=====================================================================
      program Match_MG
C=====================================================================
      implicit none
      real*8 ra(*)
      integer*4 ia(*),mt,n,iget1,iget,malloc
      pointer (mia,ia),(mra,ra)
      external malloc
C---------------------------------------------------------------------
C...  Initializing the memory allocation for RA and IA.
C---------------------------------------------------------------------
      mt = 16
      open(mt,file='matrices/matrix_rhs.unf5',
     >     status='unknown',form='unformatted')
C
      read(mt) n
      rewind mt
C
      iget1 = 26*n
      iget = iget1 * 4
      mia = malloc(iget)
C
      write(*,*) mia,n,iget1,iget
      read(*,*) 
C
      iget = iget1 * 8
      mra = malloc(iget)
C     
      write(*,*) mra,n,iget1,iget
      read(*,*)
C
cc      call solv00(ia,ra,mt)
      call solv00(ia,ra,iget1,iget1)

      call free(mia)
      call free(mra)
C
      stop
      end
C=====================================================================
      subroutine solv00(ia,ra,irend,iaend)
C=====================================================================
      real*8 ra(1)
      integer*4 ia(1)
C---------------------------------------------------------------------
C...  Testing the program.
C---------------------------------------------------------------------
      do k = irend-99, irend
         ra(k) = k
         write(*,*) k,ra(k)
      end do
      read(*,*)
      do k = iaend-99, iaend
         ia(k) = k
         write(*,*) k,ia(k)
      end do
C
      return
      end
C=====================================================================

