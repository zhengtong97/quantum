Module KroneckerProduct_Mod
Implicit None
Integer , parameter , private :: DP = Selected_Real_Kind( 9 )
contains

Subroutine KroneckerProduct( A , B , H )
Real( Kind = DP ) , Intent( IN )  :: A(:,:) , B(:,:)
Real( Kind = DP ) , Intent( OUT ) :: H(:,:)
Integer :: i , j , m , n , p , q
m = size( A , dim = 1 )
n = size( A , dim = 2 )
p = size( B , dim = 1 )
q = size( B , dim = 2 )
Do i = 1 , m
Do j = 1 , n
H( p*(i-1)+1 : p*i , q*(j-1)+1 : q*j ) = B * A(i,j)
End Do
End Do
End Subroutine KroneckerProduct

End Module KroneckerProduct_Mod



    PROGRAM quantum_usingJCB
use KroneckerProduct_Mod
    implicit none
    integer,parameter::np=10,axnp=20
    integer::i,j,k,info,time_int,L
    real(kind=8)::g1,g2,emp,del,time_real,EPS
    real(kind=8),dimension(np,np)::a=0.d0,adrag=0.d0,II=0.d0,adraga=0.d0
    real(kind=8),dimension(2,2)::sigma_x,sigma_z,sigma_xbar
    real(kind=8),dimension(axnp,axnp)::ax,adragx,xhat,H,zhat
    real(kind=8),dimension(axnp,axnp)::V
    real(kind=8),dimension(2)::up,down
    real(kind=8),dimension(axnp)::lam
    real(kind=8),dimension(axnp)::wave,c,hec
    real(kind=8),dimension(axnp,1)::hec1
    complex(kind=8)::epsx,epsz,change,com
    com=(0,1)
    
!测试赋值
 g1=0.0
 g2=0.0
 emp=1.0
 v=0.0
 del=1.0
    
!对波函数初始态赋值
    wave(1)=0.5**0.5
    wave(2)=0.5**0.5
    wave(3:axnp)=0

!定义泡利矩阵
sigma_x(1,1)=1.d0
sigma_x(2,2)=1.d0
sigma_z(1,1)=1.d0
sigma_z(2,2)=-1.d0
sigma_xbar(1,2)=1.d0
sigma_xbar(2,1)=1.d0

up=(/1.d0,0.d0/)
down=(/0.d0,1.d0/)


!统一定义输出格式，可以直接在此处修改

100 FORMAT(10(f5.2))
200 FORMAT(20(f5.2))

!定义a，adrag ,I


write(*,*)"Matrix I"

do j=1,np
II(j,j)=1.d0
end do

do j=1,np
write(*,100) II(j,:)
end do

do j=1,np-1
a(j,j+1)=real(j**0.5)
end do

write(*,*)"Matrix a"

do j=1,np
write(*,100) a(j,:)
end do


write(*,*)"Matrix adrag"
do j=1,np-1
adrag(j+1,j)=real(j**0.5)
end do

do j=1,np
write(*,100) adrag(j,:)
end do

call KroneckerProduct(a,sigma_x,ax)

write(*,*)"Matrix ax"

do j=1,axnp
write(*,200) ax(j,:)
end do


call KroneckerProduct(adrag,sigma_x,adragx)

write(*,*)"Matrix adragx"

do j=1,axnp
write(*,200) adragx(j,:)
end do

write(*,*)"Matrix adragx*ax"
write(*,200)matmul(adragx,ax)

call KroneckerProduct(II,sigma_xbar,xhat)

write(*,*)"Matrix xhat"
do j=1,axnp
write(*,200) xhat(j,:)
end do

call KroneckerProduct(II,sigma_z,zhat)

write(*,*)"Matrix zhat"
do j=1,axnp
write(*,200) zhat(j,:)
end do

H=del*zhat+emp*xhat+v*matmul(adragx,ax)+g1*Matmul((ax+adragx),xhat)+g2*Matmul(Matmul((ax+adragx),(ax+adragx)),xhat)
write(*,*)"Matrix H"
do j=1,axnp
write(*,200) H(j,:)
end do

EPS=1.e-30 !雅可比方法求矩阵特征值，特征向量，先设置精度，减小误差


call CJCBI(H,axnp,EPS,V,L)
!!调用函数后H对角元素为特征值，V为特征向量组成的矩阵，EPS为人工定义的精度
write(*,*)"对角线上为特征值："
write(*,200)H

write(*,*)"提取特征值,组成向量"
do i=1,axnp
    lam(i)=H(i,i)
    write(*,*)lam(i)
end do




write(*,*)"Matrix composed of special hectors"
write(*,200) V
write(*,*)"special value"
write(*,*)lam

do j=1,axnp
    do i=1,axnp
        hec(i)=V(i,j)
    end do
    c(j)=dot_product(wave,hec)
end do

write(*,*)"wave"
write(*,*)wave
write(*,*)"C(i)"
write(*,200)c

300 FORMAT(A4,10X,A8,20X,A8)
write(*,300)"t","epsx","epsz"

do time_int=0,40,1
time_real=real(time_int)/10.d0
epsx=0
epsz=0
do j=1,axnp
    do k=1,axnp
        hec1(k,1)=V(k,j)
    end do
    do i=1,axnp
        do k=1,axnp
            hec(k)=V(k,i)
        end do
        change=com*(lam(j)-lam(i))*time_real
        epsx=epsx+c(i)*c(j)*exp(change)*sum(matmul(hec,matmul(xhat,hec1)))!把x，z写到同一个
        epsz=epsz+c(i)*c(j)*exp(change)*sum(matmul(hec,matmul(zhat,hec1)))!循环中，节约时间
    end do
end do

400 FORMAT(F5.2,5X,("("F8.3,",",F8.3,2X,")"),5X,("(",F10.5,",",F10.5,2X,")"))!使输出更加美观
write(*,400)time_real,epsx,epsz

end do

END




SUBROUTINE CJCBI(A,N,EPS,V,L)

REAL(kind=8)::A(N,N),V(N,N)
REAL(kind=8)::EPS
DOUBLE PRECISION FM,CN,SN,OMEGA,X,Y
INTEGER P,Q

L=1
DO 20 I=1,N
V(I,I)=1.0
DO 10 J=1,N
IF (I.NE.J) V(I,J)=0.0
10      CONTINUE
20    CONTINUE
25    FM=0.0
DO 30 I=2,N
DO 30 J=1,I-1
IF (ABS(A(I,J)).GT.FM) THEN
FM=ABS(A(I,J))
P=I
Q=J
END IF
30    CONTINUE
IF (FM.LT.EPS) THEN
L=1
RETURN
END IF
IF (L.GT.1000) THEN
L=0
RETURN
END IF
L=L+1
X=-A(P,Q)
Y=(A(Q,Q)-A(P,P))/2.0
OMEGA=X/SQRT(X*X+Y*Y)
IF (Y.LT.0.0) OMEGA=-OMEGA
SN=1.0+SQRT(1.0-OMEGA*OMEGA)
SN=OMEGA/SQRT(2.0*SN)
CN=SQRT(1.0-SN*SN)
FM=A(P,P)
A(P,P)=FM*CN*CN+A(Q,Q)*SN*SN+A(P,Q)*OMEGA
A(Q,Q)=FM*SN*SN+A(Q,Q)*CN*CN-A(P,Q)*OMEGA
A(P,Q)=0.0
A(Q,P)=0.0
DO 60 J=1,N
IF ((J.NE.P).AND.(J.NE.Q)) THEN
FM=A(P,J)
A(P,J)=FM*CN+A(Q,J)*SN
A(Q,J)=-FM*SN+A(Q,J)*CN
END IF
60    CONTINUE
DO 70 I=1,N
IF ((I.NE.P).AND.(I.NE.Q)) THEN
FM=A(I,P)
A(I,P)=FM*CN+A(I,Q)*SN
A(I,Q)=-FM*SN+A(I,Q)*CN
END IF
70    CONTINUE
DO 80 I=1,N
FM=V(I,P)
V(I,P)=FM*CN+V(I,Q)*SN
V(I,Q)=-FM*SN+V(I,Q)*CN
80    CONTINUE
GOTO 25
END




