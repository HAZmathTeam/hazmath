C=====================================================================
      subroutine ddata(lvl,nsub,n1x,n2x,n1y,n2y, 
     >     n1midx,n2midx,n1midy,n2midy) 
C=====================================================================
      integer n1x(1),n2x(1),n1y(1),n2y(1)
      integer n1midx(1),n2midx(1),n1midy(1),n2midy(1)
C---------------------------------------------------------------------
C...  Form the node information of each subdomain.
C...  Both nonoverlapping and overlapping subdomains are considered.
C...  
C...  Parameters
C...    LVL        - the level of the finest mesh.
C...    NSUB       - the number of subdomains
C...    N1X,N1MIDX - global node numbering of left end point in the 
C...                 x-direction of each subdomain 
C...    N2X,N2MIDX - global node numbering of right end point in the
C...                 x-direction of each subdomain 
C...    N1Y,N1MIDY - global node numbering of bottom end point in the
C...                 y-direction of each subdomain 
C...    N2Y,N2MIDY - global node numbering of top end point in the
C...                 y-direction of each subdomain
C---------------------------------------------------------------------
      nxy  = 2**lvl + 1
      nmid = nxy/2 + 1
      nlb  = 3 * 2**(lvl-3) + 1
      nrt  = 5 * 2**(lvl-3) + 1
C...
      nmidx = nmid
      nmidy = nmidx
C     
      if (nsub .eq. 1) then
         n1midx(1) = 1
         n2midx(1) = nxy
         n1midy(1) = 1
         n2midy(1) = nxy
      else if (nsub .eq. 2) then
         n1midx(1) = 1
         n2midx(1) = nmidx
         n1midy(1) = 1
         n2midy(1) = nxy
         n1midx(2) = nmidy
         n2midx(2) = nxy
         n1midy(2) = 1
         n2midy(2) = nxy
      else if (nsub .eq. 4) then
         n1midx(1) = 1
         n2midx(1) = nmidx
         n1midy(1) = 1
         n2midy(1) = nmidy
         n1midx(2) = 1
         n2midx(2) = nmidx
         n1midy(2) = nmidy
         n2midy(2) = nxy
         n1midx(3) = nmidx
         n2midx(3) = nxy
         n1midy(3) = 1
         n2midy(3) = nmidy
         n1midx(4) = nmidx
         n2midx(4) = nxy
         n1midy(4) = nmidy
         n2midy(4) = nxy
      end if
C
      if (nsub .eq. 1) then
         n1x(1) = 1
         n2x(1) = nxy
         n1y(1) = 1
         n2y(1) = nxy
      else if (nsub .eq. 2) then
         n1x(1) = 1
         n2x(1) = nrt
         n1y(1) = 1
         n2y(1) = nxy
         n1x(2) = nlb
         n2x(2) = nxy
         n1y(2) = 1
         n2y(2) = nxy
      else if (nsub .eq. 4) then
         n1x(1) = 1
         n2x(1) = nrt
         n1y(1) = 1
         n2y(1) = nrt
         n1x(2) = 1
         n2x(2) = nrt
         n1y(2) = nlb
         n2y(2) = nxy
         n1x(3) = nlb
         n2x(3) = nxy
         n1y(3) = 1
         n2y(3) = nrt
         n1x(4) = nlb
         n2x(4) = nxy
         n1y(4) = nlb
         n2y(4) = nxy
      end if
C
      return
      end
C=====================================================================
      subroutine formsub(lvl,n,n1,n2,nsub,ns,n1s,n2s,m,mfree,
     >     isub_all,jsub_all,jsub_mid,idirs,ielsub,nel)
C=====================================================================
      parameter (nsubmax = 10)
      integer n1x(nsubmax),n2x(nsubmax),n1y(nsubmax),n2y(nsubmax)
      integer n1midx(nsubmax),n2midx(nsubmax),
     >        n1midy(nsubmax),n2midy(nsubmax)
      integer idirs(1),isub_all,jsub_all,jsub_mid
      integer m(1),ns(1),n1s(1),n2s(1),ielsub(1)
C---------------------------------------------------------------------
C...  Generate mesh information about each subdomain.
C...  
C...  Parameters
C...    N,N1,N2 - Mesh information on the global domain
C...    NSUB    - the number of subdomains
C...    N1X     - global node numbering of left end point in the 
C...              x-direction of each subdomain 
C...    N2X     - global node numbering of right end point in the
C...              x-direction of each subdomain 
C...    N1Y     - global node numbering of bottom end point in the
C...              y-direction of each subdomain 
C...    N2Y     - global node numbering of top end point in the
C...              y-direction of each subdomain 
C...    NS(k)   - the dimension of subdomain k
C...    IDIRS   - pointers to arrays, which identify the Dirichlet nodes 
C...              of each subdomain; it contains ns(k)+1 for a Dirichlet
C...              node of k-th subdomain, 0 otherwise.
C...    ISUB_ALL,JSUB_ALL - the global reference numbering correspoding
C...              to the local numbering of each overlapping subdomain
C...    ISUB_ALL,JSUB_MID - identifies nodes which are in the
C...              nonoverlapping region for each subdomain.
C...    IELSUB  - To each global element (triangle) assign the subdomain
C...              number to which the element belongs
C---------------------------------------------------------------------
C
C...  ERROR STOP
C
      call ddata(lvl,nsub,n1x,n2x,n1y,n2y,n1midx,n2midx,n1midy,n2midy) 
C
      nodesub = 0
      do k = 1, nsub
         n1s(k)  = (n2x(k)-n1x(k)+1)
         n2s(k)  = (n2y(k)-n1y(k)+1)
         ns(k)   = n1s(k)*n2s(k)
         nodesub = nodesub + ns(k)
      end do
C     
      isub_all = 1
      jsub_all = isub_all + nsub + 1
      idirs(1) = jsub_all + nodesub
      do k = 1, nsub - 1
         idirs(k+1) = idirs(k) + ns(k)
      end do
C      
      jsub_mid = idirs(nsub) + ns(nsub) 
      imfree = jsub_mid + nodesub
C
      write(*,*) n1,n2,n
      if(n1*n2 - n .ne.  0) stop 10
C
      call inullv(m(jsub_mid),nodesub)
C
      m(isub_all) = 1
      ipoint = 1
      do k = 1, nsub
         k1x = n1x(k)
         k2x = n2x(k)
         k1y = n1y(k)
         k2y = n2y(k)
         do jx = k1x, k2x
            do jy = k1y, k2y
C...           Get the  node number 
               nodenum = nomxy(jx,jy,n1,n2)
               m(jsub_all+ipoint-1) = nodenum
               if  (jx .ge. n1midx(k) .and. jx .le. n2midx(k) .and.
     >              jy .ge. n1midy(k) .and. jy .le. n2midy(k)) then
                  m(jsub_mid+ipoint-1) = m(jsub_mid+ipoint-1) + 1
               end if
C
               if  (jx.eq.k1x .or. jx.eq.k2x .or. 
     >              jy.eq.k1y .or. jy.eq.k2y) then
                  m(idirs(1)+ipoint-1) = ns(k) + 1
               else
                  m(idirs(1)+ipoint-1) = 0
               end if
               ipoint = ipoint + 1
            end do
         end do
         m(isub_all+k) = ipoint
      end do
C
      mfree = mfree - 1
C
      isub_all = isub_all + mfree
      jsub_all = jsub_all + mfree
      jsub_mid = jsub_mid + mfree
C
      do k = 1 , nsub
         idirs(k)  = idirs(k) + mfree
      end do
      mfree = imfree + mfree
C
      do k = 1 , nsub
         k1x = n1midx(k)
         k2x = n2midx(k)
         k1y = n1midy(k)
         k2y = n2midy(k)
         do i = k1x, k2x-1
            do j = k1y, k2y-1
               nelnum = nomxy(i,j,n1-1,n2-1)
               iel1 = 2*nelnum-1 
               iel2 = 2*nelnum
               ielsub(iel1) = k 
               ielsub(iel2) = k 
            end do
         end do
      end do
cc      write(*,*) (ielsub(kk),kk=1,nel)
      return
      end
C====================================================================
      subroutine get_blocks(n,ia,ja,a,ias,jas,kas,nsub,ns,
     >     isub_all,jsub_all,iwork,idirs,m,r,mfree,kfree)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),a(1),isub_all(1),jsub_all(1)
      dimension iwork(1),ias(1),jas(1),kas(1),m(1),r(1)
      dimension idirs(1),ns(1)
C---------------------------------------------------------------------
C...  Generate the stiffness matrix of each subdomain.
C...
C...  Parameters:
C...    N, NS(k)    - the dimension of global system and k-th subdomain
C...    IA, JA, A   - the stiffness matrix for the global mesh
C...    IAS,JAS,KAS - the pointer of each subdomain's stiffness matrix
C...    NSUB        - the number of subdomains 
C...    IWORK       - the temporary working space and later it restores
C...                  the information of the Dirichlet nodes.
C...    IDIRS       - array which identifies Dirichlet nodes of 
C...                  each subdomain; it contains ns(k)+1 for a 
C...                  Dirichlet node of k-th subdomain, 0 otherwise.
C...    ISUB_ALL,JSUB_ALL - the global reference numbering correspoding
C...                        to the local numbering of each subdomain
C---------------------------------------------------------------------
      mfree = mfree - 1
      kfree = kfree - 1
C
      isubt = 1
      jsubt = isubt + n + 1
C
      call iit(isub_all,jsub_all,nsub,n,m(isubt),m(jsubt))
      call iit(m(isubt),m(jsubt),n,nsub,isub_all,jsub_all)
C
      call inullv(iwork,n)
C
      ias(1) = 1
      jas(1) = ias(1) + isub_all(2) - isub_all(1) + 1
      kas(1) = 1
C
      kdir_point = 1
C
      do i = 1 , nsub
         ipoint = 1
         do ksub = isub_all(i),isub_all(i+1)-1
            iwork(jsub_all(ksub)) = ipoint
            ipoint = ipoint + 1
         end do
C
         nsp = ns(i) + 1
         ipoint = 1
         m(ias(i)) = 1
C
         do ksub = isub_all(i),isub_all(i+1)-1
            k    = jsub_all(ksub)
            kblk = iwork(k)
            kdir = idirs(kdir_point+kblk-1)
            if (kdir .eq. nsp) then
               m(jas(i) + ipoint - 1) = kblk
               r(kas(i) + ipoint - 1) = 1.d0
               ipoint = ipoint + 1
               go to 10
            end if
C
            do jk = ia(k),ia(k+1) - 1
               l = ja(jk)
               lblk = iwork(l)
               if (lblk .ne. 0 .and. 
     >              idirs(kdir_point+lblk-1) .eq. 0) then
                  m(jas(i) + ipoint - 1) = lblk
                  r(kas(i) + ipoint - 1) = a(jk)
                  ipoint = ipoint + 1
               end if
            end do
 10         continue
            m(ias(i) + kblk) = ipoint
         end do         
C
         nnz = ipoint - 1
C
         do ksub = isub_all(i),isub_all(i+1)-1
            iwork(jsub_all(ksub)) = 0
         end do
C
         if(i .ne. nsub) then
            ias(i+1) = jas(i) + nnz
            jas(i+1) = ias(i+1) + ns(i+1) + 1
            kas(i+1) = kas(i) + nnz
         else
            imfree = jas(i) + nnz
            ikfree = kas(i) + nnz
         end if
         kdir_point = kdir_point + ns(i)
      end do
C...
      do i = 1 , nsub
         ias(i) = ias(i) + mfree
         jas(i) = jas(i) + mfree
         kas(i) = kas(i) + kfree
      end do
C
      mfree = mfree + imfree
      kfree = kfree + ikfree
C
C...  Restoring the Dirichlet nodes.
C
      do k = 1 , n
         if(ia(k+1)-ia(k) .le. 1) then
            iwork(k) = n + 1
         else
            iwork(k) = 0
         end if
      end do
      return
      end
C=====================================================================
      subroutine formlast(node_overlap,ififth_inv,ififth,n,nsfifth)
C=====================================================================
      integer*4 node_overlap(1),ififth(1),ififth_inv(1),n,nsfifth
C---------------------------------------------------------------------
C...  Get mesh information about the last domain which is overlapped by 
C...  at least two subdomains.
C...
C...  Parameters:
C...    NODE_OVERLAP  - controls the number of times overlapped in each node.
C...    IFIFTH        - global numbering of nodes in the overlapped domain.
C...    IFIFTH_INV    - the inverse of IFIFTH.
C---------------------------------------------------------------------
      call inullv(ififth_inv,n)
C
      nsfifth = 0
      do k = 1 , n
         if(node_overlap(k) .gt. 1) then
            nsfifth = nsfifth + 1
            ififth(nsfifth) = k
            ififth_inv(k) = nsfifth
         end if
      end do
      return
      end
C=====================================================================
      subroutine get_fifth_block(ia,ja,a,nsfifth,
     >        ififth_inv,ififth,ia5,ja5,a5,nnz5)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ififth(1),ia5(1),ja5(1),a5(1),ififth_inv(1)
      dimension ia(1),ja(1),a(1)
C---------------------------------------------------------------------
C...  Generate the stiffness matrix for the all overlapping region.
C...
C...  Parameters:
C...    IFIFTH     - global numbering of nodes in the overlapped domain.
C...    IFIFTH_INV - the inverse of IFIFTH.
C---------------------------------------------------------------------
      ipoint = 1
      ia5(1) = 1
C
      do i = 1 , nsfifth
         k    = ififth(i)
         kblk = ififth_inv(k)
         if (kblk .gt. 0) then
            do jk = ia(k),ia(k+1) - 1
               l = ja(jk)
               lblk = ififth_inv(l)
               if (lblk .gt. 0) then
                  ja5(ipoint) = lblk
                  a5(ipoint) = a(jk)
                  ipoint = ipoint + 1
               end if
            end do
            ia5(kblk+1) = ipoint
         end if
      end do         
C
      nnz5 = ipoint - 1
C
      return
      end
C=====================================================================
      subroutine overlap_all(nsub,isub_all,jsub_all,ns,
     >     idir,node_overlap)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension isub_all(1),jsub_all(1),idir(1),node_overlap(1),ns(1)
C---------------------------------------------------------------------
C...  Compute the number of times overlapped in each node. Here all
C...  the nodes overlapped are considered.
C...  
C...  Parameters:
C...    NODE_OVERLAP  - controls the number of times overlapped
C...                    in each node.
C---------------------------------------------------------------------
      do i = 1, nsub
         do k = isub_all(i), isub_all(i+1)-1
            if(idir(k) .le. ns(i)) then
               node_overlap(jsub_all(k)) = node_overlap(jsub_all(k))+1
            end if
         end do
      end do
      return
      end
C=====================================================================
      subroutine overlap_mid(nsub,isub_all,jsub_all,ns,
     >     jsub_mid,node_overlap)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension isub_all(1),jsub_all(1),jsub_mid(1),node_overlap(1)
      dimension ns(1)
C---------------------------------------------------------------------
C...  Compute the number of times overlapped in each node. Here only
C...  midpoints of overlapped region are considered.
C...  
C...  Parameters:
C...    NODE_OVERLAP  - controls the number of times overlapped
C...                    in each node.
C---------------------------------------------------------------------
      do i = 1, nsub
         do k = isub_all(i), isub_all(i+1)-1
            if(jsub_mid(k) .gt. 0) then
               node_overlap(jsub_all(k)) = node_overlap(jsub_all(k))+1
            end if
         end do
      end do
      return
      end
C=====================================================================
      function nomxy(i,j,nx,ny)
C=====================================================================
C---------------------------------------------------------------------
C...  NOMXY gives the global number of the node (i,j)
C---------------------------------------------------------------------
      NOMXY = (i-1)*ny + j
C
      if(i .le. 0 .or. j .le. 0) then
         nomxy = 0
      end if
      return
      end
C====================================================================
