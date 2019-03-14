 program main
    implicit none
    integer,parameter::np=16
    integer::m,n,nrot
    real(kind=8),dimension(np,np)::a=0,b=0,v
    real(kind=8),dimension(np)::d

!需要根据矩阵的维数修改输出格式
100 FORMAT(16(f10.5))
200 FORMAT(16(f20.16,/))
    !定义Hilbert矩阵
    do m=1,np
    a(m,m)=(2.d0*real(m)-1.d0)
    end do

    do m=1,np-1
    b(m,m+1)=0.5d0
    b(m+1,m)=0.5d0
    end do

    write(*,*)np,"维A矩阵:"
    do m=1,np
        write(*,100) a(m,:)
    end do


    call jacobi(A,np,np,d,v,nrot)
    call eigsrt(d,v,np,np)


    write(*,*)"本征值"
    write(*,200) d


    write(*,*)"本征向量"
    do m=1,np
        write(*,100) v(m,:)
    end do

end

    SUBROUTINE jacobi(a,n,np,d,v,nrot)

! Computes all eigenvalues and eigenvectors of a real symmetric matrix a, 
! which is of size n by n, stored in a physical np by np array. 
! On output, elements of a above the diagonal are destroyed. 
! d returns the eigenvalues of a in its first n elements. 
! v is a matrix with the same logical and physical dimensions as a, 
! whose columns contain, on output, the normalized eigenvectors of a. 
! nrot returns the number of Jacobi rotations that were required.
! Please notice that the eigenvalues are not ordered on output. 
! If the sorting is desired, the addintioal routine "eigsrt" 
! can be invoked to reorder the output of jacobi.

      INTEGER::n,np,nrot,NMAX
      REAL(kind=8)::a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL(kind=8)::c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+&
     &g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END







     SUBROUTINE eigsrt(d,v,n,np)
!
! Given the eigenvalues d and eigenvectors v as output from jacobi or tqli,
! this routine sorts the eigenvalues into descending order,
! and rearranges the columns of v correspondingly.
!
      INTEGER n,np
      REAL(kind=8)::d(np),v(np,np)
      INTEGER i,j,k
      REAL(kind=8) p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END
