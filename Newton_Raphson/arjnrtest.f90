!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Test Program for arjnr.dll
!
!****************************************************************************

program arjunnewton
implicit none

!Variable declaration
real,dimension(2):: xo
real:: delx,tol
integer:: num

!Interface to vector function
interface
function fvector(x)
real,dimension(2),intent(in)::x
real,dimension(2)::fvector
end function fvector
end interface

!Set Parameters
num = 2
delx = 0.00001
tol = 0.000001

!Initial Guess
xo(1) = -1
xo(2) = -1

! Call Solver
call arjnr(fvector,num,xo,delx,tol)

print *,xo

end program arjunnewton

!******************************************************************
! VECTOR FUNCTION
function fvector(x)
implicit none

!Function variable declaration
real,dimension(2),intent(in)::x
real,dimension(2)::fvector

!Function definition
fvector(1) = x(1)**2 + 2*x(2)**2 - 12
fvector(2) = x(1)*x(1) - 30

end function fvector
!******************************************************************

