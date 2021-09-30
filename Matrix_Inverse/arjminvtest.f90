!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Program to test arjminv.dll
!
!****************************************************************************

program matrixinv
implicit none

!Variable declaration
real,dimension(3,3):: A,invsol,sol
integer:: num,i,j

num = 3

read *,A

call arjminv(A,num,invsol)

! Print inverse
do i=1,num
write(*,10)  (invsol(i,j),j=1,num)
end do
10 format(10f10.3)   

end program matrixinv