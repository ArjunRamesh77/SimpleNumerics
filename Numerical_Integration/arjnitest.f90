!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		01/01/2015
!  DESCRIPTION: Program to test arjni.dll
!
!****************************************************************************

program arjunsimpson
implicit none

!Variable declaration
real:: h,int
real,dimension(6):: y
integer:: num

!Parameters
num = 6
int = 1

!Data
h = (1.0/6.0)
y = (/1.0,0.8571,0.75,0.6666,0.6,0.5454/)

! Call integration routine
call arjsimp(y,h,num,int)
print *,int

end program arjunsimpson

