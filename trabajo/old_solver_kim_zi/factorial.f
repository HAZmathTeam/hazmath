C=====================================================================
cc      program FACTORIAL
      subroutine dummy
C=====================================================================
      parameter (nlast = 5 000 000, nrlast = 2 500 000) !level = 8
      parameter (nsize = 200)
      implicit real*8(a-h,o-z)
      dimension nb(nsize)
      real*8 ntotal
C
      ipt = 1
      nb(ipt) = 10
      nl = 9
      ntotal = 1.d0
C
      if (nb(ipt) .gt. nl) then
         nb(ipt+1) = nb(ipt) - 1
         ipt = ipt + 1
         call recursive1(nb,nl,ntotal,nsize,ipt)
         ipt = ipt - 1
         ntotal = ntotal * nb(ipt) 
      else if (nb(ipt) .eq. nl) then
         ntotal = nb(ipt)
      else 
         write(*,*) 'NB is less than NL.'
         ntotal = 0
      end if
C
      write(*,*) 'ipt, nb, nl, ntotal at main:', 
     >     ipt, nb(ipt), nl, ntotal
C
      stop
      end
C=====================================================================
      subroutine recursive1(nb,nl,ntotal,nsize,ipt)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension nb(nsize)
      real*8 ntotal
C
      write(*,*) 'ipt, nb, nl, ntotal at AA11:', 
     >     ipt, nb(ipt), nl, ntotal
C
      if (nb(ipt) .gt. nl) then
         nb(ipt+1) = nb(ipt) - 1
         ipt = ipt + 1
         call recursive2(nb,nl,ntotal,nsize,ipt)
         ipt = ipt - 1
         ntotal = ntotal * nb(ipt) 
      else if (nb(ipt) .eq. nl) then
         ntotal = nb(ipt)
      else 
         write(*,*) 'NB is less than NL.'
      end if
C
      write(*,*) 'ipt, nb, nl, ntotal at BB11:', 
     >     ipt, nb(ipt), nl, ntotal
C
      return
      end   
C=====================================================================
      subroutine recursive2(nb,nl,ntotal,nsize,ipt)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension nb(nsize)
      real*8 ntotal
C
      write(*,*) 'ipt, nb, nl, ntotal at AAAA2222:', 
     >     ipt, nb(ipt), nl, ntotal
C
      if (nb(ipt) .gt. nl) then
         nb(ipt+1) = nb(ipt) - 1
         ipt = ipt + 1
         call recursive1(nb,nl,ntotal,nsize,ipt)
         ipt = ipt - 1
         ntotal = ntotal * nb(ipt) 
      else if (nb(ipt) .eq. nl) then
         ntotal = nb(ipt)
      else 
         write(*,*) 'NB is less than NL.'
      end if
C
      write(*,*) 'ipt, nb, nl, ntotal at BBBB2222:', 
     >     ipt, nb(ipt), nl, ntotal
C
      return
      end   
C=====================================================================
C=====================================================================
      program FACTORIAL_new
C=====================================================================
      parameter (nlast = 5 000 000, nrlast = 2 500 000) !level = 8
      parameter (nsize = 200)
      implicit real*8(a-h,o-z)
      dimension nb(nsize)
      real*8 ntotal
C
      ipt = 1
      nb(ipt) = 10
      nl = 1
      ntotal = 1.d0
C
      if (nb(ipt) .gt. nl) then
         call recursive1_new(nb,nl,ntotal,nsize,ipt)
      else if (nb(ipt) .eq. nl) then
         ntotal = ntotal*nb(ipt)
      else 
         write(*,*) 'NB is less than NL.'
         ntotal = 0
      end if
C
      write(*,*)
      write(*,*) 'ipt, nb, nl, ntotal at main:', 
     >     ipt, nb(ipt), nl, ntotal
C
      stop
      end
C=====================================================================
      subroutine recursive1_new(nb,nl,ntotal,nsize,ipt)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension nb(nsize)
      real*8 ntotal
C
      write(*,*) 'ipt, nb, nl, ntotal at AA11:', 
     >     ipt, nb(ipt), nl, ntotal
C
      nb(ipt+1) = nb(ipt) - 1
      ipt = ipt+1
      if (nb(ipt) .eq. nl) then
         ipt = ipt - 1
         ntotal = ntotal*nb(ipt)*nb(ipt+1)
      else if (nb(ipt) .gt. nl) then
         call recursive2_new(nb,nl,ntotal,nsize,ipt)
         ipt = ipt - 1
         ntotal = ntotal * nb(ipt) 
      end if
C
      write(*,*) 'ipt, nb, nl, ntotal at BB11:', 
     >     ipt, nb(ipt), nl, ntotal
C
      return
      end   
C=====================================================================
      subroutine recursive2_new(nb,nl,ntotal,nsize,ipt)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension nb(nsize)
      real*8 ntotal
C
      write(*,*) 'ipt, nb, nl, ntotal at AA22:', 
     >     ipt, nb(ipt), nl, ntotal
C
      nb(ipt+1) = nb(ipt) - 1
      ipt = ipt+1
      if (nb(ipt) .eq. nl) then
         ipt = ipt - 1
         ntotal = ntotal*nb(ipt)*nb(ipt+1)
      else if (nb(ipt) .gt. nl) then
         call recursive1_new(nb,nl,ntotal,nsize,ipt)
         ipt = ipt - 1
         ntotal = ntotal*nb(ipt) 
      end if
C
      write(*,*) 'ipt, nb, nl, ntotal at BB22:', 
     >     ipt, nb(ipt), nl, ntotal
C
      return
      end   
C=====================================================================
