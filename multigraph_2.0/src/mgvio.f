c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine vutl(ncolor,red,green,blue,isock,fname)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            integer(kind=iknd) :: vioctr,viodtr,vioutl
            real(kind=rknd), dimension(*) :: red,green,blue
            character(len=4) :: type
            character(len=8) :: keywrd
            character(len=80) :: host,host0,fname,fname0
            common /gph0/id
cy
        if(ncolor<=0) then
            call viostr(id,'putl',4)
            call vioint(id,ncolor,1)
            ii=vioutl(id,'c')
            ii=viodtr(id)
            call viostp()
            return
        else
c
            if(rknd==rsngl) then
                keywrd='bhsingle'
            else if(rknd==rdble) then
                keywrd='bhdouble'
            else
                keywrd='bhdouble'
            endif
c
            call viosta()
            if(isock<0) then
                call fstr(fname0,lenf,fname,0_iknd)
                type='FILE'
                host='localhost'
            else
                call sint(fname0,lenf,isock)
                type='INET'
                host=fname
            endif
            call fstr(host0,lenh,host,0_iknd)
            id=vioctr(type,'XDR',host0,lenh,fname0,lenf,'w')
            ii=vioutl(id,'o')
            call viostr(id,keywrd,8)
            call viostr(id,'putl',4)
            call vioint(id,ncolor,1)
            if(rknd==rsngl) then
                call vioflt(id,red,ncolor)
                call vioflt(id,green,ncolor)
                call vioflt(id,blue,ncolor)
            else if(rknd==rdble) then
                call viodbl(id,red,ncolor)
                call viodbl(id,green,ncolor)
                call viodbl(id,blue,ncolor)
            else
                call vioqud(id,red,ncolor)
                call vioqud(id,green,ncolor)
                call vioqud(id,blue,ncolor)
            endif
        endif
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine vframe(iframe)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            common /gph0/id
cy
        call viostr(id,'list',4)
        call vioint(id,iframe,1)
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine vline(x,y,z,n,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(n) :: x,y,z
            common /gph0/id
cy
        call viostr(id,'line',4)
        call vioint(id,icolor,1)
        call vioint(id,n,1)
        if(rknd==rsngl) then
            call vioflt(id,x,n)
            call vioflt(id,y,n)
            call vioflt(id,z,n)
        else if(rknd==rdble) then
            call viodbl(id,x,n)
            call viodbl(id,y,n)
            call viodbl(id,z,n)
        else
            call vioqud(id,x,n)
            call vioqud(id,y,n)
            call vioqud(id,z,n)
        endif
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine vfill(x,y,z,n,icolor)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(n) :: x,y,z
            common /gph0/id
cy
        call viostr(id,'fill',4)
        call vioint(id,icolor,1)
        call vioint(id,n,1)
        if(rknd==rsngl) then
            call vioflt(id,x,n)
            call vioflt(id,y,n)
            call vioflt(id,z,n)
        else if(rknd==rdble) then
            call viodbl(id,x,n)
            call viodbl(id,y,n)
            call viodbl(id,z,n)
        else
            call vioqud(id,x,n)
            call vioqud(id,y,n)
            call vioqud(id,z,n)
        endif
        return
        end
c-----------------------------------------------------------------------
c
c            piecewise lagrange triangle multi grid package
c
c                    edition 11.0 - - - june, 2012
c
c-----------------------------------------------------------------------
        subroutine vioqud(id,x,n)
cx
            use mthdef
            implicit real(kind=rknd) (a-h,o-z)
            implicit integer(kind=iknd) (i-n)
            real(kind=rknd), dimension(*) :: x
            real(kind=rdble), dimension(n) :: t
cy
        do i=1,n
            t(i)=x(i)
        enddo
        call viodbl(id,t,n)
        return
        end
