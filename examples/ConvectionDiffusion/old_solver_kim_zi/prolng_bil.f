C=====================================================================
      subroutine prolng_symb_bil(nxf,nyf,nf,nxc,nyc,nc,ipp,jpp,np)
C=====================================================================
      dimension ipp(1),jpp(1)
C---------------------------------------------------------------------
C...  Symbolic assembly of a prolongation matrix for uniform meshes
C...  on a rectangular domain. Ordering of the grid points is given as 
C...  follows: start from left bottom corner, go upward in y direction,
C...  and do the same way for the following columns.
C...
C...  Input:
C...    NXF,NYF  - dimension of the fine mesh with NXF being the 
C...               number of grid points in x direction and NYF in
C...               y direction. The following relations hold:
C...               NXF = 2*NXC - 1, NYF = 2*NYC - 1.
C...
C...  Working Parameters:
C...    NXC,NYC  - dimension of the coarse mesh with NXC being the 
C...               number of grid points in x direction and NYC in
C...               y direction. The following relations hold:
C...               NXC = (NXF + 1)/2, NYC = (NYF + 1)/2.
C...
C...  Output:
C...    IPP,JPP  - structure of the prolongation matrix of dimension
C...               NFxNC, in RRCO form.
C...
C...  Note:
C...    NF+1     - the dimension of IPP.
C...    NF,NC    - dimension of the prolongation matrix; NF rows and
C...               NC columns. Here NF = NXF*NYF and NC = NXC*NYC.
C...    NP       - dimension of JPP; NP = IPP(NF+1)-1
C---------------------------------------------------------------------
      nxc  = (nxf + 1)/2
      nyc  = (nyf + 1)/2
      nc   = nxc*nyc
      nxc1 = nxc - 1
      nyc1 = nyc - 1
      m1   = nyf + nyc1
      m2   = 2*nyf + 2*nyc1
C
C...  Prolongate first for the odd columns of the rectangular grid.
C
      m  = 1 - m2
      kf = 1 - nyf
      kc = 0
      do 20 i = 1,nxc
         m  = m + m2
         kf = kf + nyf
         do 10 j = 1, nyc1
            ipp(kf)   = m
            ipp(kf+1) = m + 1
            kc        = kc + 1
            jpp(m)    = kc
            jpp(m+1)  = kc
            jpp(m+2)  = kc + 1
            m         = m + 3
            kf        = kf + 2
 10      continue
         kc      = kc + 1
         ipp(kf) = m
         jpp(m)  = kc
         m       = m + 1
         kf      = kf + 1
 20   continue
C
      ipp(nf+1) = m
      np        = m - 1
C
C...  Prolongate for the even columns of the rectangular grid.
C
      m  = 1
      kf = 1
      kc = 0
      do 40 i = 1,nxc1
         m  = m + m1
         kf = kf + nyf
         do 30 j = 1, nyc1
            ipp(kf)   = m
            ipp(kf+1) = m + 2
            kc        = kc + 1
            kcc       = kc + nyc
            jpp(m)    = kc
            jpp(m+1)  = kcc
            jpp(m+2)  = kc 
            jpp(m+3)  = kc + 1 
            jpp(m+4)  = kcc
            jpp(m+5)  = kcc + 1 
            m         = m + 6
            kf        = kf + 2
 30      continue
         kc       = kc + 1
         ipp(kf)  = m
         jpp(m)   = kc
         jpp(m+1) = kcc + 1
         m        = m + 2
         kf       = kf + 1
 40   continue
C
cc      do kkk = 1 , nf
cc         write(*,*) kkk
cc         write(*,*) (jpp(k123),k123=ipp(kkk),ipp(kkk+1)-1)
cc         read(*,*)
cc      end do

      return
      end
C=====================================================================
      subroutine pbyvg_bil(ipp,jpp,b,nc,c,nf)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ipp(1),jpp(1),b(nc),c(nf)
C--------------------------------------------------------------------
C...  PRODUCT--- Prolongation: c = P*b.
C--------------------------------------------------------------------
      do 10 i = 1, nf
         ipa = ipp(i)
         ipb = ipp(i+1) - 1
         ipd = ipb - ipa
         if(ipd .eq. 0) then
            c(i) = b(jpp(ipa)) 
         else if(ipd .eq. 1) then
            c(i) = (b(jpp(ipa)) + b(jpp(ipb)))*0.5d0
         else if(ipd .eq.3) then
            c(i) = (b(jpp(ipa)) + b(jpp(ipa+1)) 
     >           + b(jpp(ipa+2)) + b(jpp(ipb)))*0.25d0
         end if
 10   continue
      return
      end
C=====================================================================
      subroutine pbyvcs_bil(ipp,jpp,b,nc,c,nf)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ipp(1),jpp(1),b(nc),c(nf)
C--------------------------------------------------------------------
C...  PRODUCT--- Vector + Prolongation: c = c + P*b.
C--------------------------------------------------------------------
      do 10 i = 1, nf
         ipa = ipp(i)
         ipb = ipp(i+1) - 1
         ipd = ipb - ipa
         if(ipd .eq. 0) then
            c(i) = c(i) + b(jpp(ipa)) 
         else if(ipd .eq. 1) then
            c(i) = c(i) + (b(jpp(ipa)) + b(jpp(ipb)))*0.5d0
         else if(ipd .eq.3) then
            c(i) = c(i) + (b(jpp(ipa)) + b(jpp(ipa+1)) 
     >           + b(jpp(ipa+2)) + b(jpp(ipb)))*0.25d0
         end if
 10   continue
      return
      end
C=====================================================================
      subroutine pbyvcm_bil(ipp,jpp,b,nc,c,nf)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ipp(1),jpp(1),b(nc),c(nf)
C--------------------------------------------------------------------
C...  PRODUCT--- Vector - Prolongation: c = c - P*b.
C--------------------------------------------------------------------
      do 10 i = 1, nf
         ipa = ipp(i)
         ipb = ipp(i+1) - 1
         ipd = ipb - ipa
         if(ipd .eq. 0) then
            c(i) = c(i) - b(jpp(ipa)) 
         else if(ipd .eq. 1) then
            c(i) = c(i) - (b(jpp(ipa)) + b(jpp(ipb)))*0.5d0
         else if(ipd .eq.3) then
            c(i) = c(i) - (b(jpp(ipa)) + b(jpp(ipa+1)) 
     >           + b(jpp(ipa+2)) + b(jpp(ipb)))*0.25d0
         end if
 10   continue
      return
      end
C=====================================================================
      subroutine pbyvas_bil(a,ipp,jpp,b,nc,c,nf)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ipp(1),jpp(1),b(nc),c(nf),a(nf)
C--------------------------------------------------------------------
C...  PRODUCT--- Vector + Prolongation: c = a + P*b.
C--------------------------------------------------------------------
      do 10 i = 1, nf
         ipa = ipp(i)
         ipb = ipp(i+1) - 1
         ipd = ipb - ipa
         if(ipd .eq. 0) then
            c(i) = a(i) + b(jpp(ipa)) 
         else if(ipd .eq. 1) then
            c(i) = a(i) + (b(jpp(ipa)) + b(jpp(ipb)))*0.5d0
         else if(ipd .eq.3) then
            c(i) = a(i) + (b(jpp(ipa)) + b(jpp(ipa+1)) 
     >           + b(jpp(ipa+2)) + b(jpp(ipb)))*0.25d0
         end if
 10   continue
      return
      end
C=====================================================================
      subroutine pbyvam_bil(a,ipp,jpp,b,nc,c,nf)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ipp(1),jpp(1),b(nc),c(nf),a(nf)
C--------------------------------------------------------------------
C...  PRODUCT--- Vector - Prolongation: c = a - P*b.
C--------------------------------------------------------------------
      do 10 i = 1, nf
         ipa = ipp(i)
         ipb = ipp(i+1) - 1
         ipd = ipb - ipa
         if(ipd .eq. 0) then
            c(i) = a(i) - b(jpp(ipa)) 
         else if(ipd .eq. 1) then
            c(i) = a(i) - (b(jpp(ipa)) + b(jpp(ipb)))*0.5d0
         else if(ipd .eq.3) then
            c(i) = a(i) - (b(jpp(ipa)) + b(jpp(ipa+1)) 
     >           + b(jpp(ipa+2)) + b(jpp(ipb)))*0.25d0
         end if
 10   continue
      return
      end
C=====================================================================
      subroutine rvbyp_bil(ipp,jpp,b,nf,c,nc)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ipp(1),jpp(1),b(nf),c(nc)
C--------------------------------------------------------------------
C...  PRODUCT--- Restriction: c = b^t*P = P^t*b.
C--------------------------------------------------------------------
      do 10 i = 1, nc
         c(i) = 0.0d00
 10   continue
      do 30 i = 1, nf
         ipa = ipp(i)
         ipb = ipp(i+1) - 1
         ipd = ipb - ipa
         z = b(i)
         if(ipd .eq. 0) then
            j = jpp(ipa)
            c(j) = c(j) + z
         else if(ipd .eq. 1) then
            do 20 k = ipa, ipb
               j = jpp(k)
               c(j) = c(j) + 0.5d0*z
 20         continue
         else if(ipd .eq.3) then
            do 25 k = ipa, ipb
               j = jpp(k)
               c(j) = c(j) + 0.25d0*z
 25         continue
         end if
 30   continue
      return
      end
C=====================================================================
      subroutine abyp_bil(ia,ja,ipp,jpp,np,ic,jc,an,cn,x,nq)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),ipp(1),jpp(1),ic(1),jc(1)
      dimension an(1),cn(1),x(nq)
C--------------------------------------------------------------------
C...  Multiplication of a prolongation matrix P by a general sparse 
C...  matrix A : C = A*P.
C--------------------------------------------------------------------
      do i  = 1 , np
         ica = ic(i) 
         icb = ic(i+1) - 1
         if(icb .ge. ica) then
            do j = ica,icb
               x(jc(j)) = 0.0d00
            end do
            iaa = ia(i)
            iab = ia(i+1) - 1
            do jp = iaa,iab
               j = ja(jp)
               a = an(jp)
               ipa = ipp(j)
               ipb = ipp(j+1) - 1
               ipd = ipb - ipa
               if(ipd .eq. 0) then
                  k    = jpp(ipa)
                  x(k) = x(k) + a
               else if(ipd .eq. 1) then
                  a5   = a*0.5d0
                  do k1 = ipa, ipb
                     k = jpp(k1)
                     x(k) = x(k) + a5
                  end do
               else if(ipd .eq.3) then
                  ad4   = a*0.25d0
                  do k1 = ipa, ipb
                     k = jpp(k1)
                     x(k) = x(k) + ad4
                  end do
               end if
            end do
            do j = ica,icb
               cn(j) = x(jc(j))
            end do
         end if
      end do
      return
      end
C=====================================================================
      subroutine formdir(idir_fine,ip,jp,n,idir,nc)
C=====================================================================
      integer idir_fine(1),ip(1),jp(1),idir(1)
C---------------------------------------------------------------------
C...  Generates Dirichlet boundary nodes for the lower level using the
C...  Dirichlet boundary data in the upper level and prolongation
C...  matrix.
C---------------------------------------------------------------------
      call inullv(idir,nc)
C
      do i = 1, n
         if ((ip(i+1)-ip(i).gt.1) .or. (idir_fine(i).le.n)) go to 10
         icoarse  = jp(ip(i))
         idir(icoarse) = nc + 1
 10      continue
      end do
C
      return
      end
C=====================================================================
      subroutine adir_new(ia,ja,an,idir,n)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),idir(1)
C---------------------------------------------------------------------
C...  Modify a general matrix A in order to take care of the Dirichlet
C...  boundary node. All the entries of A which correspond to the
C...  Dirichlet nodes (either column part or row part) will be 
C...  assigned to be 0, except the diagonal entries, in this case, it
C...  is assigned to be 1.
C...  This subroutine does not store the entries of new zeros.
C...  (This gives the difference from the subroutine ADIR_OLD.)
C---------------------------------------------------------------------
      do 30 i = 1, n
         if (idir(i) .gt. n) then
            iaa = ia(i)
            iab = ia(i+1) - 1
            do 10 j = iaa,iab
               k = ja(j)
               if (k .eq. i) then
                  an(j) = 1.d0
               else 
                  ja(j) = 0
               end if
 10         end do
         else
            iaa = ia(i)
            iab = ia(i+1) - 1
            do 20 j = iaa,iab
               k = ja(j)
               if (idir(k) .gt. n) ja(j) = 0
 20         end do
         end if
 30   end do
      new_nnz = 0
      iaa_old = ia(1)
      do i = 1, n
         iab_old = ia(i+1) - 1
         do j = iaa_old,iab_old
            k = ja(j)
            if (k .ne. 0) then
               new_nnz = new_nnz + 1
               an(new_nnz) = an(j)
               ja(new_nnz) = ja(j)
            end if
         end do
         ia(i+1) = new_nnz + 1
         iaa_old = iab_old + 1
      end do
      return
      end
C=====================================================================
      subroutine adir_old(ia,ja,an,idir,n)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(1),ja(1),an(1),idir(1)
C---------------------------------------------------------------------
C...  Modify a general matrix A in order to take care of the Dirichlet
C...  boundary node. All the entries of A which correspond to the
C...  Dirichlet nodes (either column part or row part) will be 
C...  assigned to be 0, except the diagonal entries, in this case, it
C...  is assigned to be 1.
C---------------------------------------------------------------------
      do 30 i = 1, n
         if (idir(i) .gt. n) then
            iaa = ia(i)
            iab = ia(i+1) - 1
            do 10 j = iaa,iab
               k = ja(j)
               if (k .eq. i) then
                  an(j) = 1.d0
               else 
                  an(j) = 0.d0
               end if
 10         end do
         else
            iaa = ia(i)
            iab = ia(i+1) - 1
            do 20 j = iaa,iab
               k = ja(j)
               if (idir(k) .gt. n) an(j) = 0.d0
 20         end do
         end if
 30   end do
      return
      end
C=====================================================================
