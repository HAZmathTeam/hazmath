C=====================================================================
      subroutine sgauss(ia,ja,a,b,maxa,au,n)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),a(1),maxa(1),au(1),b(1)
C---------------------------------------------------------------------
C...  Gaussian elimination.
C---------------------------------------------------------------------
cc      write(*,*) ' Gauss elimination started.'
C
      call Upper(ia,ja,a,maxa,au,n,nwk)
      call Lower(ia,ja,a,maxa,au(nwk+1),n,nwk)
      np = n + 1
cc      write(*,*) (b(k),k=1,n)
cc      write(*,*) (maxa(kk),kk=1,np)
cc      read(*,*)
cc      write(*,*) (au(kk),kk=1,maxa(np)-1)
cc      write(*,*) (au(kk+nwk),kk=1,maxa(np)-1)
      call decons(au,au(nwk+1),maxa,n)
      call redbns(au,au(nwk+1),b,maxa,n)
CCCCCCCCCCCCCCCCCC If symmetric:
!      call decomp(au,maxa,ish,n)
!      call redbak(au,b,maxa,n)
C
cc      write(*,*) ' Gauss elimination ended.'
      return
      end
C=====================================================================
      subroutine Upper(xadj,adjncy,a,maxa,au,n,nwk)
C=====================================================================
      integer*4 adjncy(1),xadj(1),maxa(1)
      real*8 a(1),au(1)
C---------------------------------------------------------------------
C...  This subroutine transforms the adjacency graph structure into 
C...  envelope one. The goal is to invert them into AL and AU by gauss 
C...  elimination. 
C---------------------------------------------------------------------
      maxa(1) = 1
      do 20 k = 2, n+1
         nodei = k-1
         icol  = 0
         do 10 jk = xadj(nodei), xadj(nodei+1) - 1
            icol = max0(icol,k-adjncy(jk))
 10      end do
         maxa(k) = maxa(k-1)+icol
 20   end do
C
      nnm = n+1
      nwk = maxa(nnm) - 1
cc      write(*,*) nwk, ' nonzeroes in the envelope'
C
      call nullv(au,nwk)
C
      kkz = 0
      do k = 1, n
         nodei = k
         do 50 jk = xadj(nodei),xadj(nodei+1)-1
            nodej = adjncy(jk)
            idist = iabs(nodej - k)
            if(nodej-k) 50, 30, 40
 30         au(maxa(k)) = a(jk)
            go to 50
 40         au(maxa(nodej)+idist) = a(jk)
            go to 50
 50      end do
      end do
      return
      end
C=====================================================================
      subroutine Lower(xadj,adjncy,a,maxa,al,n,nwk)
C=====================================================================
      integer*4 adjncy(1),xadj(1),maxa(1)
      real*8 a(1),al(1)
C---------------------------------------------------------------------
C...  This subroutine transforms the adjacency graph structure into 
C...  envelope one. The goal is to invert them into AL and AU by gauss 
C...  elimination. 
C---------------------------------------------------------------------
      nnm = n+1
cc      write(*,*) nwk, ' nonzeroes in the envelope'
C
      call nullv(al,nwk)
C
      do 50 k = 1 , n
         nodei = k
         do 40 jk = xadj(nodei),xadj(nodei+1)-1
            nodej = adjncy(jk)
            idist = iabs(nodej - k)
            if(nodej-k) 10, 20, 40
 10         al(maxa(nodei)+idist) = a(jk)
            go to 40
 20         al(maxa(k)) = 1.0d00
 40      end do
 50   end do
      return
      end
C=====================================================================
      subroutine decons (au,al,maxa,nn)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension au(1),al(1),maxa(1)
C---------------------------------------------------------------------
C...  To calculate (l)*(u) factorization of nonsymmetric
C...  stiffness matrix using upper and lower envelopes
C---------------------------------------------------------------------
      if (nn.eq.1) return
      nlow=0
      do 100 n=1,nn
      kn=maxa(n)
      ku=maxa(n+1)-1
      kh=ku-kn
C
      al(kn) = 1.d0
      if (kh) 100,100,110
C
 110  k = n - kh
      ki=maxa(k)
      al(ku) = al(ku)/au(ki)
cc    au(ku) = au(ku) !!!!!!
      kh1 = kh - 1
C
      if (kh1) 120,120,130
 130  k=n-kh1
      ic=0
      klt=ku
      do 140 j=1,kh1
         ic=ic+1
         klt=klt-1
         ki=maxa(k)
         nd=maxa(k+1)-ki-1
C
         if (nd) 141,145,145
 145     kk=min(ic,nd)
         cl=0.0d00
         cu=0.0d00
         do 150 l=1,kk
            kind = ki  + l
            lind = klt + l
            cl=cl + au(kind)*al(lind)
            cu=cu + al(kind)*au(lind)
 150     continue
         al(klt)=(al(klt)-cl)/au(ki)
         au(klt)= au(klt)-cu
 141     k=k+1
 140  continue
C
 120  b=0.d0
      kl = kn + 1
      do 160 kk=kl,ku
         cu=au(kk)
         cl=al(kk)
         b=b+cl*cu
 160  continue
      au(kn)=au(kn) - b
      if((au(kn)) .lt. 1.d-13) then
         klll = n
 111     write(*,1111) klll,au(kn)
      end if
 100  continue
C
 1111 format(//10x,'          **** General error ****   '/
     >     10x,'***Zero/negative diagonal entry in stiffnes matrix***'/
     >     /10x,'      i = ',i10,/10x,'     a(i,i) =',e15.7/)
      return
      end
C=====================================================================
      subroutine redbns(au,al,v,maxa,nn)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension au(1),al(1),v(1),maxa(1)
C---------------------------------------------------------------------
C...  To reduce and back-substitute iteration vectors
C...  for nonsymmetrical matrix  
C---------------------------------------------------------------------
      do 400 n = 1 , nn
         kl = maxa(n) + 1
         ku = maxa(n+1) - 1
         if (ku - kl) 400,410,410
 410     k = n
         c = 0.d0
         do 420 kk = kl , ku
            k = k - 1
            c = c + al(kk) * v(k)
 420     continue
         v(n) = v(n) - c
 400  continue
C
cc      read (mt1) (a(j),j=1,nwk)
C
      k = maxa(nn)
      v(nn) = v(nn) / au(k)
C
      n = nn
      do 500 l = 2 , nn
         kl = maxa(n) + 1
         ku = maxa(n+1) - 1
         if (ku - kl) 530,510,510
 510     k = n
         do 520 kk = kl , ku
            k = k - 1
            v(k) = v(k) - au(kk) * v(n)
 520     continue
 530     nl = n - 1
         ln = maxa(nl)
         v(nl) = v(nl) / au(ln)
         n = n - 1
 500  continue
C
      return
      end
C=====================================================================
      subroutine decomp (a,maxa,nn,ish)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension a(1),maxa(1)
C---------------------------------------------------------------------
C...  To calculate (l)*(d)*(l)(i) factorization of stiffness matrix
C---------------------------------------------------------------------
      if (nn.eq.1) return
      nlow=0
      do 200 n=1,nn
         kn=maxa(n)
         kl=kn+1
         ku=maxa(n+1)-1
         kh=ku-kl
C
         if (kh) 304,240,210
 210     k=n-kh
         ic=0
         klt=ku
         do 260 j=1,kh
            ic=ic+1
            klt=klt-1
            ki=maxa(k)
            nd=maxa(k+1)-ki-1
            if (nd) 260,260,270
 270        kk=min(ic,nd)
            c=0.
            do 280 l=1,kk
               c=c+a(ki+l)*a(klt+l)
 280        continue
            a(klt)=a(klt)-c
            k=k+1
 260     continue
C
 240     k=n
         b=0.
         do 300 kk=kl,ku
            k=k-1
            ki=maxa(k)
            c=a(kk)/a(ki)
            if (dabs(c).lt.1.e07) go to 290
            write (*,2010) n,c
            stop
 290        b=b+c*a(kk)
            a(kk)=c
 300     continue
         a(kn)=a(kn)-b
C
 304     if (a(kn)) 310,310,200
 310     nlow=nlow+1
         if (ish.eq.0) go to 320
         if (a(kn).eq.0) a(kn)=-1.e-16
         go to 200
 320     write(*,2000) n,a(kn)
         stop
 200  continue
C
      write(*,2020) nlow    
      write(*,'(/a/)') ' Decomp ended'
C
 2000 format(//48h stop - stiffness matrix not positive definite  ,//
     >         32h nonpositive pivot for equation ,i4,//
     >         10h pivot =  ,e20.12)
 2010 format (//'   stop - sturm sequense check failed because of
     >  multiplier growth for column number',i4,//'multiplier=',e20.8)
 2020 format(/5x,'there are ',i4,' roots lower then shift')
C
      return
      end
C=====================================================================
      subroutine redbak (a,v,maxa,nn)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension a(1),v(1),maxa(1)
C---------------------------------------------------------------------
C...  To reduce and back-substitute iteration vectors
C---------------------------------------------------------------------
      do 400 n=1,nn
         kl=maxa(n)+1
         ku=maxa(n+1)-1
         if (ku-kl) 400,410,410
 410     k=n
         c=0.
         do 420 kk=kl,ku
            k=k-1
            c=c+a(kk)*v(k)
 420     continue
         v(n)=v(n) - c
 400  continue
C
      do 480 n=1,nn
         k=maxa(n)
         v(n)=v(n)/a(k)
 480  continue
      if (nn.eq.1) return
      n=nn
      do 500 l=2,nn
         kl=maxa(n)+1
         ku=maxa(n+1)-1
         if (ku-kl) 500,510,510
 510     k=n
         do 520 kk=kl,ku
            k=k-1
            v(k)=v(k)-a(kk)*v(n)
 520     continue
         n=n-1
 500  continue
C
      return
      end
C=====================================================================
