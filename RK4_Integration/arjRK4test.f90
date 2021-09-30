!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Program to test arjRK4.dll
!
!****************************************************************************

program arjunRK4
implicit none

!Variable declaration
real,dimension(2):: init
real,dimension(200):: rangex
real,dimension(2,200):: ystore
real:: h
integer:: num,i,numrange

!Function declaration
interface 
function func(x,y)
real,intent(in):: x
real,dimension(2),intent(in):: y
real,dimension(2):: func
end function func
end interface

!Parameters
num = 2
numrange = 200

!Step size and range
h = (10.0 - 0.0)/numrange

rangex(1) = 0.0
do i=2,numrange
  rangex(i) = rangex(i-1) + h
end do

!Initial Value
init = (/0.2,0.2/)

call arjRK4(func,num,numrange,init,rangex,ystore)

print *,ystore(1,1:numrange)

end program arjunRK4

!**********************************************************
function func(x,y)

!Variable declaration
real,intent(in):: x
real,dimension(100),intent(in):: y
real,dimension(100):: func

!Function definition
func(1) = y(1)*(1 - y(1))
func(2) = y(2)*(1 - y(2))

end function func
!**********************************************************

