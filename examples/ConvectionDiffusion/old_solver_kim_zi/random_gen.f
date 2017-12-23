C=====================================================================
      program   random_gen
C=====================================================================
      parameter (imax = 1 000 000 000, jmax = 12 000 000)
      implicit real*8(a-h,o-z)
      integer*4 iperm(jmax)
C---------------------------------------------------------------------
C...  Generating a permutaion (quasi-random order) of N integers 
C...  {1,2,...,N}.
C...  
C...  Parameter
C...    IPERM - a permutation of N integers {1,2,...,N}.
C---------------------------------------------------------------------
      write(*,*) 'Enter the number of nodes:'
      read(*,*) n
      k  = 0
      n1 = n + 1
C
      do 10 i = 1, n
         iperm(i) = 0
 10   continue
C
      do 40 i = 1, imax
         x  = dble(i)
         x  = dsin(x)
         x  = dabs(x)
         k1 = n1*x
         k1 = mod(k1, n) + 1
         if (iperm(k1) .eq. 0) then
            k         = k + 1 
            iperm(k1) = k
            if (k .eq. n) go to 50
         end if
 40   continue
C 
 50   continue
      write(*,*) (i,'*',iperm(i), i=1,n)
C
cc      return
      end
C=====================================================================
      subroutine random_gen1
C=====================================================================
      parameter (imax = 1 000 000 000, jmax = 12 000 000)
      implicit real*8(a-h,o-z)
      IMPLICIT integer*8 (I-N)
C---------------------------------------------------------------------
C...  Generating random numbers on (0,1).
C...  
C...  Parameter
C...    I
C---------------------------------------------------------------------
      i1=  2
      do i=1,100
ccc        r2 =  rand(i1)
        i1=iabs(irand(i))
cccccccc        i2 = 130*((rng(i1)+1.d00)*0.5d00)
        i2 = 1.d17*((rng(i1)+1.d00)*0.5d00)
        i2 = mod(i2,127)+1
cccccccccc        print *,'i1=',i1,' random is = ',i2
        print *,' random is = ',i2
      enddo
      stop
      end
C=====================================================================
      REAL*8 FUNCTION RNG(ISEED)
C=====================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
C---------------------------------------------------------------------
C     RANDOM NUMBER GENERATOR
C---------------------------------------------------------------------
      ISEED=MOD(ISEED*24+23,1361)
      RNG=ISEED/1359.*.99923+.00038
      RETURN
      END
C=====================================================================


