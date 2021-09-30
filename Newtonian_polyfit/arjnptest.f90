!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Program to test arjnp.dll
!
!****************************************************************************

program arjnewpoly
implicit none

!Variable declaration
real,dimension(5):: X,Y
real,dimension(10):: val,ans
integer:: num,valsize,i

!Data
X = (/1.0,2.0,3.0,4.0,5.0/)
Y = (/2.0,6.0,10.0,15.0,20.0/)

!Parameters
num = 5
valsize = 10

!Create val array
val(1) = 1.0
do i=2,10
  val(i) = val(i-1) + (5.0 - 1.0)/9.0
end do


call arjnpoly(X,Y,num,valsize,val,ans)

print *,ans

end program arjnewpoly