!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Subroutine to perform Broydens method
!
!****************************************************************************

program arjunnewton
implicit none

!Variable declaration
real,dimension(4):: xo
real:: delx,tol
integer:: num

!Interface to vector function
interface
function fvector(x)
real,dimension(4),intent(in)::x
real,dimension(4)::fvector
end function fvector
end interface

!Set Parameters
num = 4
delx = 0.01
tol = 0.0001

!Initial Guess
xo = (/0.0,0.0,0.0,0.0/)

! Call Solver
call arjnewton(fvector,num,xo,delx,tol)

print *,xo

return 
end program arjunnewton

!******************************************************************
! VECTOR FUNCTION
function fvector(x)
implicit none

!Function variable declaration
real,dimension(4),intent(in)::x
real,dimension(4)::fvector

! DEFINE EUATIONS
fvector(1) = x(1) - 3*x(2)**2 + 10*x(3) + 20
fvector(2) = x(1)*x(2) + x(3)*3 + x(4)**3 - 10*sin(x(1))
fvector(3) = -x(1) + 123*x(2) + 3*x(3)**2 - 0.001*x(4)*exp(x(3)) - 400
fvector(4) = x(1)*sin(x(2)) - x(2) + x(3) + x(4)**3

end function fvector
!******************************************************************

!$$$$$$ ! VECTOR FUNCTION
!$$$$$$ function fvector(x)
!$$$$$$ implicit none
!$$$$$$ 
!$$$$$$ !Function variable declaration
!$$$$$$ real,dimension(7),intent(in)::x
!$$$$$$ real,dimension(7)::fvector
!$$$$$$ real:: l1,l2,l3,l4,l5,l6,l7,pi
!$$$$$$ 
!$$$$$$ l1 = 100.0; 
!$$$$$$ l2 = 100.0; 
!$$$$$$ l3 = 200.0; 
!$$$$$$ l4 = 75.0; 
!$$$$$$ l5 = 100.0; 
!$$$$$$ l6 = 75.0; 
!$$$$$$ l7 = 50.0;
!$$$$$$ pi = 3.141; 
!$$$$$$ 
!$$$$$$ ! DEFINE EUATIONS
!$$$$$$ fvector(1) = l3*x(3)**2 - l4*x(4)**2 + 0.00000001;
!$$$$$$ fvector(2) = l2*x(2)**2 + l4*x(4)**2 + l5*x(5)**2 - l6*x(6)**2 + 0.00000001; 
!$$$$$$ fvector(3) = l1*x(1)**2 + l6*x(6)**2 + l7*x(7)**2 - ((5.2e5)*(pi**2)*((0.2)**5))/(8.0*0.02*998.0);
!$$$$$$ fvector(4) = x(1) - x(2) - x(6) + 0.00000001;
!$$$$$$ fvector(5) = x(2) - x(4) - x(3) + 0.00000001;
!$$$$$$ fvector(6) = x(4) + x(3) - x(5) + 0.00000001;
!$$$$$$ fvector(7) = x(6) + x(5) - x(7) + 0.00000001; 
!$$$$$$ 
!$$$$$$ end function fvector