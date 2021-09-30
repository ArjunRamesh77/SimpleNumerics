!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Program to perform Newtonian Polynomial interpolation
!
!****************************************************************************


! Subroutine to determine Newton Polynomials
subroutine arjnpoly(X,Y,num,valsize,val,ans)
implicit none

!Variable declaration
real,dimension(num):: X,Y
real:: F,P
real,allocatable,dimension(:):: a,coef
real,dimension(valsize):: val,ans
integer:: num,valsize,i,j,nl

!Allocate memory
allocate(a(num))
allocate(coef(num))

a = Y

!Create coeff array
do i=1,num
  coef(i) = 0.0
end do

!Determine Coefficients
coef(1) = a(1)
do i=1,size(X)-1
   do j=1,size(X) - i
      a(j) = (a(j) - a(j+1))/(X(j) - X(j+i))
      if(j == 1) coef(i+1) = a(j)
   end do
end do


!Main Loop
do nl = 1,size(val)
    
    F = coef(1)

do i=1,size(X) - 1
   P = coef(i+1)
   do j=1,i
     P = P*(val(nl) - X(j))
   end do
   F = F + P
end do

ans(nl) = F
end do

!dallocate memory
deallocate(a)
deallocate(coef)

end subroutine arjnpoly
