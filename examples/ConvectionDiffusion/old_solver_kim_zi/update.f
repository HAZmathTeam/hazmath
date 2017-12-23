C=====================================================================
      subroutine update(u,err,n,node_overlap)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension u(1),node_overlap(1),err(1)
C---------------------------------------------------------------------
C...  Update the global solution U by adding the solutions ERR  of the
C...  residual equations of each subdomain. For the overlapping regions
C...  the numerical values are averaged.
C...
C...  Parameter:
C...    NODE_OVERLAP - contains the number of times overlapped for each node.
C---------------------------------------------------------------------
      do k = 1, n
         if(node_overlap(k) .lt. 1) go to 10
         u(k) = u(k) + err(k)/dble(node_overlap(k))
 10      continue
      end do
      return
      end
C=====================================================================
      subroutine update1(u,err,n,node_overlap,u_exact)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension u(1),node_overlap(1),err(1),u_exact(1)
C---------------------------------------------------------------------
C...  Update the global solution U by adding the solutions ERR  of the
C...  residual equations of each subdomain. For the overlapping regions
C...  the numerical values are averaged.
C...
C...  Parameter:
C...    NODE_OVERLAP - contains the number of times overlapped for each node.
C---------------------------------------------------------------------
      do k = 1, n
         if(node_overlap(k) .eq. 1) then
            u(k) = u(k) + err(k)
         else
            u(k) = u_exact(k)
            if(node_overlap(k) .gt. 1) stop 66
         end if
      end do
      return
      end
C=====================================================================
      subroutine scatter(jsub_all,u,us,ns,idir)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension jsub_all(1),us(1),u(1),idir(1)
C---------------------------------------------------------------------
C...  Distribute the residual to each subdomain.
C---------------------------------------------------------------------
      do k = 1, ns
         if (idir(k) .le. ns) then
            us(k) = u(jsub_all(k))
         else
            us(k) = 0.0d00
         end if
      end do
      return
      end
C=====================================================================
      subroutine scatter1(jsub_all,u,us,ns)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension jsub_all(1),us(1),u(1)
C---------------------------------------------------------------------
C...  Distribute the residual to the subdomain without considering
C...  the boudary conditions.
C---------------------------------------------------------------------
      do k = 1, ns
         us(k) = u(jsub_all(k))
      end do
      return
      end
C=====================================================================
      subroutine gather_all(jsub_all,us,u,ns,idir)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension jsub_all(1),idir(1),us(1),u(1)
C---------------------------------------------------------------------
C...  Collect the coarse grid solutions US to combine them into the 
C...  global working array U. On the overlapped regions the averaged 
C...  values will be computed everywhere in the subroutine UPDATE.
C...  
C...  Parameters:
C...    US, U  - the coarse and fine grid solutions
C---------------------------------------------------------------------
      do k = 1, ns
         if(idir(k) .le. ns) then
            u(jsub_all(k)) = u(jsub_all(k)) + us(k)
         end if
      end do
      return
      end
C=====================================================================
      subroutine gather_mid(jsub_all,us,u,ns,jsub_mid)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension jsub_all(1),us(1),u(1),jsub_mid(1)
C---------------------------------------------------------------------
C...  Collect the coarse grid solutions US to combine them into the 
C...  global working array U. On the overlapped regions the averaged 
C...  values will be computed only on the common boundaries of 
C...  subdomains in the subroutine UPDATE.
C...  
C...  Parameters:
C...    US, U  - the coarse and fine grid solutions
C---------------------------------------------------------------------
      do k = 1,ns
         if(jsub_mid(k) .gt. 0) then
            u(jsub_all(k)) = u(jsub_all(k)) + us(k)
         end if
      end do
      return
      end
C=====================================================================
      subroutine gather_all_old(jsub_all,us,u,ns,idir,node_overlap)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension jsub_all(1),idir(1),us(1),u(1),node_overlap(1)
C---------------------------------------------------------------------
C...  Collect the coarse grid solutions US to combine them into the 
C...  global working array U. On the overlapped regions the averaged 
C...  values will be computed everywhere in the subroutine UPDATE.
C...  
C...  Parameters:
C...    US, U        - the coarse and fine grid solutions
C...    NODE_OVERLAP - controls the number of times added in each node.
C---------------------------------------------------------------------
      do k = 1,ns
         if(idir(k) .le. ns) then
            u(jsub_all(k)) = u(jsub_all(k)) + us(k)
            node_overlap(jsub_all(k)) = node_overlap(jsub_all(k)) + 1
         end if
      end do
      return
      end
C=====================================================================
      subroutine gather_mid_old(jsub_all,us,u,ns,jsub_mid,node_overlap)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension jsub_all(1),us(1),u(1),jsub_mid(1),node_overlap(1)
C---------------------------------------------------------------------
C...  Collect the coarse grid solutions US to combine them into the 
C...  global working array U. On the overlapped regions the averaged 
C...  values will be computed only on the common boundaries of 
C...  subdomains in the subroutine UPDATE.
C...  
C...  Parameters:
C...    US, U        - the coarse and fine grid solutions
C...    NODE_OVERLAP - controls the number of times added in each node.
C---------------------------------------------------------------------
      do k = 1,ns
         if(jsub_mid(k) .gt. 0) then
            u(jsub_all(k)) = u(jsub_all(k)) + us(k)
            node_overlap(jsub_all(k)) = node_overlap(jsub_all(k))  + 1
         end if
      end do
      return
      end
C=====================================================================
      subroutine dir_two(v,ip,jp,nf,nc,idir)
C=====================================================================
      integer ip(1),jp(1),idir(1)
      real*8 v(1)
C---------------------------------------------------------------------
C...  Generates Dirichlet boundary nodes for the lower level using the
C...  Dirichlet boundary data in the upper level and prolongation
C...  matrix. The values of V at the Dirichlet nodes are assigned to 
C...  be zero. 
C---------------------------------------------------------------------
      do i = 1 , nf
         if ((ip(i+1)-ip(i) .gt. 1) .or. 
     >        (idir(i) .le. nf)) go to 10
         icoarse  = jp(ip(i))
         v(icoarse) = 0.0d00
         idir(i) = -icoarse
 10      continue
         ip(i) = 0
      end do
      do i = 1 , nf
         if(idir(i) .lt. 0) then
            ip(iabs(idir(i))) = nc + 1
         end if
      end do
      call icopyv(ip,idir,nc)
      return
      end
C=====================================================================
