C====================================================================
      subroutine wcycle_lvl(lvl,kend,mcont)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension mcont(0:128)
C--------------------------------------------------------------------
C...  To control the levels of W-cycle type multigrid version of
C...  two-grid method.
C--------------------------------------------------------------------
      kend = 2**(lvl-2)
C
      do 30 i = 1, lvl-2
         k =  kend/(2**i)
         do 20 j = 1, kend/2
            kk = k + (j-1)*2*k
            if (kk .lt. kend) then
               mcont(kk) = lvl + 1 - i
            else
               go to 30
            end if
 20      continue
 30   continue
      mcont(kend) = 0
      mcont(0)    = lvl
C
ccc      write(*,*) 'mcont:', (mcont(kk), kk = 1, kend)
C
      return
      end
C====================================================================
