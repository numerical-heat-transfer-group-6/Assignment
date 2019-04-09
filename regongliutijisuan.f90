module publicVariables
integer,parameter :: N=10 ! 划分网格个数
real :: Fw,Fe,Dw,De,F,D !F是对流强度，D是扩散强度,下角标是w,e界面
real fai0,faiL ! fai0是左边界初始值，faiL是右边界初始值
real :: fai(N) ! 迭代过程中的节点值储存数组
real :: faiq(N) ! 迭代节点初始值设置
real :: r ! 迭代收敛残差
end module


program main

use publicVariables
implicit none

integer :: option ! 控制迭代格式
integer :: i !i是迭代字母
real :: gamma, rou, deltx, u, L ! gamma扩散系数；rou密度；deltx相邻网格距离；u速度；L长度

write(*,*) "Please input the velocity u: "
read(*,*) u
write(*,*) "Plase choose differential scheme: 1-centerdiff, 2-in-wind,3-Mixed: "
read(*,*) option


gamma=0.1
r=1
rou=1.0
L=1.0
deltx=L/N 


Fw=rou*u
Fe=rou*u

if (Fw==Fe) then
  F=Fw
end if


Dw=gamma/deltx
De=gamma/deltx
  
if (Dw==De) then
  D=Dw
end if


fai0=1.0  !左边界条件

faiL=0.0  !右边界条件

do i=1,N
  faiq(i)= 1.0  ! 高斯-赛德尔迭代初值
  fai(i) = 0.0  ! 高斯-赛德尔迭代初值
end do

if (option == 1) then 
   call centerdiff
end if

if (option == 2) then
   call in_wind
end if 

if (option == 3) then
   call  Mixed
end if

! 数据输出
write(*,fmt="(10(1XF10.5))") (fai(i),i=1,N,1)
write(*,fmt="(10(1XF10.5))") (faiq(i),i=1,N,1)
write(*,*) r 
pause

end program

! CENTER_DIFF SUBROUTINE  调用中心差分格式子程序
subroutine centerdiff  ! 中心差分格式

use publicVariables

real :: r_no, r_de
integer :: i

  do while(r>1e-6)  ! 高斯-赛德尔迭代主体以及收敛判别标准

     fai(1)=((2*D+F)*fai0+(D-F/2)*faiq(2))/(3*D+F/2)
     
     do  i=2,N-1

     fai(i)=(faiq(i+1)*(D-F/2)+fai(i-1)*(D+F/2))/(2*D)
       
     end do

     fai(N)=(fai(N-1)*(D+F/2)+faiL*(2*D-F))/(3*D-F/2)

     r_no = 0.0
     r_de = 0.0
     r = 0.0
     do i = 1,N,1
        r_no = (fai(i)-faiq(i))**2 + r_no
        r_de = fai(i)**2 + r_de
     end do
     r = sqrt(r_no)/sqrt(r_de)

   do i=1,N,1
      faiq(i) = fai(i)
   end do
  end do

end subroutine
    
! IN_WIND SUBROUTINE  调用迎风格式子程序
subroutine in_wind    ! 迎风格式

use publicVariables

real :: r_no, r_de
integer :: i
	
	do while(r>1e-6)

	 fai(1)=((2*D+F)*fai0+D*faiq(2))/(3*D+F)
     
     do  i=2,N-1

     fai(i)=(D*faiq(i+1)+fai(i-1)*(D+F))/(2*D+F)
       
     end do

     fai(N)=(fai(N-1)*(D+F)+faiL*(2*D))/(3*D+F)

     r_no = 0.0
     r_de = 0.0
     r = 0.0
     do i = 1,N,1
        r_no = (fai(i)-faiq(i))**2 + r_no
        r_de = fai(i)**2 + r_de
     end do
     r = sqrt(r_no)/sqrt(r_de)

   do i=1,N,1
      faiq(i) = fai(i)
   end do
  end do

end subroutine    

! Mixed_SUBROUTINE  调用混合格式子程序

subroutine Mixed   ! 混合格式

use publicVariables

real :: r_no, r_de
integer :: i
	
	do while(r>1e-6)

	 fai(1)=((2*D+F)*fai0+D*faiq(2))/(3*D+F)
     
     do  i=2,N-1

     fai(i)=(D*faiq(i+1)+fai(i-1)*(D+F))/(2*D+F)
       
     end do

     fai(N)=(fai(N-1)*(D+F)+faiL*(2*D))/(3*D+F)

     r_no = 0.0
     r_de = 0.0
     r = 0.0
     do i = 1,N,1
        r_no = (fai(i)-faiq(i))**2 + r_no
        r_de = fai(i)**2 + r_de
     end do
     r = sqrt(r_no)/sqrt(r_de)

   do i=1,N,1
      faiq(i) = fai(i)
   end do
  end do

end subroutine    