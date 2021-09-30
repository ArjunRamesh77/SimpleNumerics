!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Subroutine for Runge-Kutta 4th order
!
!****************************************************************************

! RK4 Subroutine
subroutine arjRK4(func,num,numrange,init,rangex,ystore)
implicit none

!Variable declaration
real,dimension(num):: init
real,dimension(numrange):: rangex
real,allocatable,dimension(:):: k1,k2,k3,k4,y
real,dimension(num,numrange):: ystore
real:: h,x
integer:: num,numrange,i

!Function interface
interface 
function func(x,y)
real,intent(in):: x
real,dimension(2),intent(in):: y
real,dimension(2):: func
end function func
end interface

!Allocate 
allocate(k1(num))
allocate(k2(num))
allocate(k3(num))
allocate(k4(num))
allocate(y(num))

!Calculate step size
h = (rangex(numrange) - rangex(1))/(numrange)

!Initial Values
ystore(1:num,1) = init
x = rangex(1)
y = init

!Main Loop - RK4
do i=2,size(rangex)

k1 = func(x,y)

k2 = func(x + 0.5*h,y + (0.5*h)*k1)

k3 = func(x + 0.5*h,y + (0.5*h)*k2)

k4 = func(x + h,y + h*k3)

y = y + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
x = x + h

ystore(1:num,i) = y
end do

!Deallocate memory 
deallocate(k1)
deallocate(k2)
deallocate(k3)
deallocate(k4)
deallocate(y)

end subroutine arjRK4