C=====================================================================
      subroutine depth_first_search_block(n,xadj,adjncy,
     >     p_vstack,lowlink,nblk,subtree_blk)
C=====================================================================
C=====================================================================
cc      subroutine depth_first_search_block(n,xadj,adjncy,p_vstack)
C=====================================================================
      integer*4 xadj(1),adjncy(1),p_vstack(1),subtree_blk(1)
      integer*4 n,np1,nb,nblk,cnt,lfp1,sp,vp,v,wp,w,qvp1,i
      integer*4 lowlink(1),edgeptr(900),lstack_number(900)
cc      integer*4 q(1)
C---------------------------------------------------------------------
C...  This subroutine performs the Tarjan's algorithm for finding the
C...  strong components of a diggraph. 
C...  This is the fred gustavson's implementation of the algorithm. 
C...  This is an implementation without q(). q is commented out.
C---------------------------------------------------------------------
C...  Initialization.
C
      nblk = 0
      nb = 0
C
      cnt  = nb+1
      lfp1 = cnt
      np1 = n+1
      vp = np1
      sp = np1
C
      do i = 1 , n
         edgeptr(i) = xadj(i)
         lstack_number(i) = 0
         lowlink(i) = 0
      end do
C
      lstack_number(n+1) = 0
C
C...  Get out when the output renumbering is full;
C...  Otherwise begin the search at vertex LFP1 for another 
C...  tree in the forest.
C
 10   if(cnt .eq. np1) then
C
C...     Reverse the ordering and exit.
C
         do i = 1,n/2
            ipp = p_vstack(i)
            p_vstack(i) = p_vstack(n-i+1)
            p_vstack(n-i+1) = ipp
         end do
C
C...     Now the starting addresses for the blocks are in lowlink.
C
         lowlink(1) = 1
         do i = 1, nblk
            lowlink(i+1) = lowlink(i)+
     >                subtree_blk(nblk-i+2)-subtree_blk(nblk-i+1)
         enddo
         return
      end if
C
      do i = lfp1,n
         if(lstack_number(i) .eq. 0) go to 30
      end do
      write(*,*) ' There is an error :::: '
      stop 1
C
 30   continue
      v = i
      lfp1 = v + 1
      go to 50
C     
C...  Recursive call of sico or whatever:::
C     
 40   continue
      vp = vp - 1
      subtree_blk(vp) = v
      v=w
C
C...  Add tree branch to the current tree::::
C
 50   continue
      nb = nb + 1
      lstack_number(v) = nb
      lowlink(v) = lstack_number(v)
      sp = sp - 1
      p_vstack(sp) = v
      qvp1 = v+1
cc    qvp1 = q(v)+1
C
C...  Examine all the edges leading out of v
C
 60   continue
      wp = edgeptr(v)
      w = adjncy(wp)
      edgeptr(v) = wp+1
      if(lstack_number(w) .ge. lstack_number(v)) go to 70
      if(lstack_number(w) .eq. 0) go to 40
      lowlink(v) = min0(lowlink(v),lstack_number(w))
 70   continue
      if(edgeptr(v) .lt. xadj(qvp1)) go to 60
C
C...  Check if v is a strong component root:::
C...  Gather the component if it is:::
C  
      if(lowlink(v) .lt. lstack_number(v)) go to 90
C
      nblk = nblk + 1
      subtree_blk(nblk) = cnt
C
 80   continue
      w = p_vstack(sp)
      lstack_number(w) = np1
      sp = sp + 1
      p_vstack(cnt) = w
      cnt = cnt + 1
      if(v .ne. w) go to 80
C
C...  Go to the begining. The present tree is finished...
C
      if(sp .eq. np1) go to 10
 90   continue
      w = v
      v = subtree_blk(vp)
      vp = vp + 1
      qvp1 = v + 1
cc    qvp1 = q(v) + 1
      lowlink(v) = min0(lowlink(v),lowlink(w))
      go to 70
C
      end
C====================================================================


