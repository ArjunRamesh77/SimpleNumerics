!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Program to test arjreg.dll
!
!****************************************************************************

program arjunreg
implicit none

! Variable declaration
real,dimension(9):: y
real,dimension(9):: x
real,dimension(4):: c
integer:: order,num

! Data
x = (/0.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/)
y = (/2.0,6.0,9.0,12.0,15.0,22.0,30.0,43.0,50.0/)

num = 9
order = 3

call arjreg(x,y,num,order,c)

print *,c

end program arjunreg
