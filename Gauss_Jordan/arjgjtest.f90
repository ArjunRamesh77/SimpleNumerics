!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Test program for arjgj.dll
!
!****************************************************************************

program arjungausselim
implicit none

!Varible declarations
real,dimension(3,4):: A
real,dimension(3):: b,sol 
integer:: num,i,j

!Parameters
num = 3

!Set Matrix
read *,A
read *,b

!Call Routine 
call arjgj(A,b,num,sol)

! Print Gauss Eliminated Matrix
do i=1,3
  write(*,11)  (A(i,j),j=1,4)
end do
11 format(10f10.2)

print *,'  '
print *,sol

end program arjungausselim
