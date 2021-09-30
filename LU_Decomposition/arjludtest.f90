!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Program to test arjlud.dll
!
!****************************************************************************

program arjunlud
implicit none

!Varible declarations
real,dimension(3,3)::A 
real,dimension(3,4):: L,U
real,dimension(3):: b,x,y 
integer:: num,i,j

!Parameters
num = 3

!Set Matrix
read *,A
read *,b

!Call Routine 
call arjlud(A,L,U,num)

! Print L Matrix
do i=1,3
  write(*,11)  (L(i,j),j=1,3)
end do
11 format(10f10.2)

print *,'  '

! Print U Matrix
do i=1,3
  write(*,12)  (U(i,j),j=1,3)
end do
12 format(10f10.2)

! Ly = b
call arjge(L,b,num,y)

! Ux = y
call arjge(U,y,num,x)

! Print solution
print *,x

end program arjunlud