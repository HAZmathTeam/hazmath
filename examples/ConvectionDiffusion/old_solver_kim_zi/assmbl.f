C====================================================================
      subroutine smbass(ie,je,iet,jet,ia,ja,n,nnz)
C====================================================================
      dimension ia(1),ie(1),je(1),iet(1),jet(1),ja(1)
C--------------------------------------------------------------------
C...  Symbolic assembly of a symmetric nodal assembly matrix.
C...
C...  Input:
C...    IE, JE   - mesh connectivity matrix in RRCU.
C...    IET, JET - transpose of IE and JE in RRCO.
C...    N        - number of nodes in the mesh.
C...    IA       - array of dimension N+1 which contains the value  
C...               N+1 in the positions corresponding to Dirichlet 
C...               nodes and 0 elsewhere.
C...
C...  Output:
C...    IA, JA   - structure of the nodal assembly matrix of 
C...               dimension NxN, in RRUU form.
C...
C...  Note:
C...    N+1 is the dimension of IA.
C...    NNZ is the dimension of JA.
C--------------------------------------------------------------------
      jp = 1
      nm = n - 1
      np = n + 1
      do 40 i = 1, nm
         jpi = jp
         if (ia(i) .eq. np) go to 30
         ieta = iet(i)
         ietb = iet(i+1) - 1
         do 20 ip = ieta, ietb
            j = jet(ip)
            iea = ie(j)
            ieb = ie(j+1) - 1
            do 10 kp = iea,ieb
               k = je(kp)
               if (k .le. i) go to 10
               if(ia(k) .ge. i) go to 10
               ja(jp) = k
               jp = jp + 1
               ia(k) = i
 10         continue
 20      continue
 30      ia(i) = jpi
 40   continue
      ia(n) = jp
      ia(np) = jp
      nnz = ia(np)-1
      return
      end
C====================================================================
      subroutine smbasg(ie,je,iet,jet,ia,ja,n,nnz,iwk)
C====================================================================
      dimension ia(1),ie(1),je(1),iet(1),jet(1),ja(1),iwk(1)
C--------------------------------------------------------------------
C...  Symbolic assembly of a general nodal assembly matrix.
C...
C...  Input:
C...    IE, JE   - mesh connectivity matrix in RRCU.
C...    IET, JET - transpose of IE and JE in RRCO.
C...    N        - number of nodes in the mesh.
C...    IWK      - array of dimension N+1 which contains the value  
C...               N+1 in the positions corresponding to Dirichlet 
C...               nodes and 0 elsewhere.
C...
C...  Output:
C...    IA, JA   - structure of the nodal assembly matrix of 
C...               dimension NxN, in RRCU form.
C...
C...  Note:
C...    N+1 is the dimension of IA.
C...    NNZ is the dimension of JA.
C--------------------------------------------------------------------
      np = n + 1
      nnz = 0
      jp = 1
      do 50 i = 1, n
         jpi = jp
         if(iwk(i) .eq. np) then
            ja(jp) = i
            jp = jp + 1
            go to 40
         end if
         ieta = iet(i)
         ietb = iet(i+1) - 1
         do 30 ip = ieta,ietb
            j = jet(ip)
            iea = ie(j)
            ieb = ie(j+1) - 1
            do 20 kp = iea,ieb
               k = je(kp)
               if (iwk(k) .ge. i) go to 20
               ja(jp) = k
               jp = jp + 1
               iwk(k) = i
 20         continue
 30      continue
 40      ia(i) = jpi
 50   continue
      ia(np) = jp
      nnz = ia(np)-1
      return
      end
C====================================================================
      subroutine assmbs(ia,ja,idir,ae,be,jep,nn,an,ad,b,ip,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),jep(1),ae(1),be(1),an(1),b(1),ad(1)
      dimension ip(1),idir(1)
C--------------------------------------------------------------------
C...  Numerical assembly of an element matrix AE and vector BE into 
C...  the nodal assembly matrix A and right-hand vector B: Symmetric 
C...  case.
C...
C...  Input:  
C...    IA, JA - structure of matrix A in RRUU. The order of A is N.
C...    IDIR   - array which identifies Dirichlet nodes; it contains
C...             n+1 for a Dirichlet node, 0 otherwise.
C...    AE, BE - element matrix and vector, respectively, to be
C...             assembled in AN, AD, and B.
C...    JEP    - reference numbers of the nodes associated with the
C...             rows and columns of AE, i.e., global reference numbers 
C...             of nodes of a triangle.
C...    NN     - order of AE and BE.
C...
C...  Output: 
C...    AN, AD - numerical values of nonzeros of A in RRUU.
C...    B      - right-hand vector of system of linear equations.
C...
C...  Working space:
C...    IP     - of demension N, initialized to 0. IP is used and 
C...             then reset to 0. IP is the expanded integer array of
C...             pointers.
C...
C...  Note:
C...    The prescribed values of the Dirichlet unknowns should be 
C...    stored in the corresponding positions of B. In other words if
C...    i is a Dirichlet node with prescribed value C_i, then, before
C...    using this algorithm, set IDIR(i) = n+1, B(i) = C_i.
C--------------------------------------------------------------------
      do 40 l = 1, nn
         i = jep(l)
         if (idir(i) .gt. n) go to 40
         k = l - nn
         ad(i) = ad(i) + ae(k+l*nn)
         b(i) = b(i) + be(l)
         kk = 0
C
         do 20 ll = 1, nn
            k = k + nn
            if (ll .eq. l) go to 20
            j = jep(ll)
            if (idir(j) .gt. n) go to 10
            if (j .lt. i) go to 20
            ip(j) = k
            kk = 1
            go to 20
 10         b(i) = b(i) - ae(k)*b(j)
 20      continue
C
         if (kk .eq. 0) go to 40
         iaa = ia(i)
         iab = ia(i+1) - 1
         do 30 j = iaa, iab
            k = ip(ja(j))
            if (k .eq. 0) go to 30
            an(j) = an(j) + ae(k)
            ip(ja(j)) = 0
 30      continue
 40   continue
C     
      return
      end
C====================================================================
      subroutine assmbg(ia,ja,idir,ae,be,jep,nn,an,b,ip,n)
C====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),jep(1),ae(1),be(1),an(1),b(1)
      dimension ip(1),idir(1)
C--------------------------------------------------------------------
C...  Numerical assembly of an element matrix AE and vector BE into 
C...  the nodal assembly matrix A and right-hand vector B: General 
C...  case.
C...
C...  Input:  
C...    IA, JA - structure of matrix A in RRCU. The order of A is N.
C...    IDIR   - Array which identifies Dirichlet nodes; it contains
C...             n+1 for a Dirichlet node, 0 otherwise.
C...    AE, BE - element matrix and vector, respectively, to be
C...             assembled in AN, and B.
C...    JEP    - reference numbers of the nodes associated with the
C...             rows and columns of AE, i.e., global reference numbers 
C...             of nodes of a triangle.
C...    NN     - order of AE and BE.
C...
C...  Output: 
C...    AN     - numerical values of nonzeros of A in RRCU.
C...    B      - right-hand vector of system of linear equations.
C...
C...  Working space:
C...    IP     - of demension N, initialized to 0. IP is used and 
C...             then reset to 0. IP is the expanded integer array of
C...             pointers.
C...
C...  Note:
C...    The prescribed values of the Dirichlet unknowns should be 
C...    stored in the corresponding positions of B. In other words if
C...    i is a Dirichlet node with prescribed value C_i, then, before
C...    using this algorithm, set IDIR(i) = n+1, B(i) = C_i.
C--------------------------------------------------------------------
      do 40 l = 1, nn
         i = jep(l)
         if (idir(i) .gt. n) go to 40
         k = l - nn
         b(i) = b(i) + be(l)
C
         do 20 ll = 1, nn
            k = k + nn
            j = jep(ll)
            if(idir(j) .gt. n) go to 10
            ip(j) = k
            go to 20
 10         b(i) = b(i) - ae(k)*b(j)
 20      continue
C
         iaa = ia(i)
         iab = ia(i+1) - 1
         kkk = 0
         do 30 j = iaa, iab
            k = ip(ja(j))
            if(k .eq. 0) go to 30
            an(j) = an(j) + ae(k)
            ip(ja(j)) = 0
            kkk = kkk + 1
            if(kkk .eq. nn) go to 40
 30      continue
 40   continue
C     
      return
      end
C====================================================================


