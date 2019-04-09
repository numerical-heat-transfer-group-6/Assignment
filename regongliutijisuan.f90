module publicVariables
integer,parameter :: N=10 ! �����������
real :: Fw,Fe,Dw,De,F,D !F�Ƕ���ǿ�ȣ�D����ɢǿ��,�½Ǳ���w,e����
real fai0,faiL ! fai0����߽��ʼֵ��faiL���ұ߽��ʼֵ
real :: fai(N) ! ���������еĽڵ�ֵ��������
real :: faiq(N) ! �����ڵ��ʼֵ����
real :: r ! ���������в�
end module


program main

use publicVariables
implicit none

integer :: option ! ���Ƶ�����ʽ
integer :: i !i�ǵ�����ĸ
real :: gamma, rou, deltx, u, L ! gamma��ɢϵ����rou�ܶȣ�deltx����������룻u�ٶȣ�L����

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


fai0=1.0  !��߽�����

faiL=0.0  !�ұ߽�����

do i=1,N
  faiq(i)= 1.0  ! ��˹-���¶�������ֵ
  fai(i) = 0.0  ! ��˹-���¶�������ֵ
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

! �������
write(*,fmt="(10(1XF10.5))") (fai(i),i=1,N,1)
write(*,fmt="(10(1XF10.5))") (faiq(i),i=1,N,1)
write(*,*) r 
pause

end program

! CENTER_DIFF SUBROUTINE  �������Ĳ�ָ�ʽ�ӳ���
subroutine centerdiff  ! ���Ĳ�ָ�ʽ

use publicVariables

real :: r_no, r_de
integer :: i

  do while(r>1e-6)  ! ��˹-���¶����������Լ������б��׼

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
    
! IN_WIND SUBROUTINE  ����ӭ���ʽ�ӳ���
subroutine in_wind    ! ӭ���ʽ

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

! Mixed_SUBROUTINE  ���û�ϸ�ʽ�ӳ���

subroutine Mixed   ! ��ϸ�ʽ

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