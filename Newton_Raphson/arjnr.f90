!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Subroutine to solve system of non-linear equations using
!				Newton-Raphson Technique
!
!****************************************************************************

! NEWTON RAPHSON SOLVER
subroutine arjnr(fvector,num,xo,delx,tol)
implicit none

!Variable declarations
integer:: num,i,j,k,l,inv,maxiter
real,dimension(num):: xo
real,allocatable,dimension(:):: x,sol,b,ff,fb,fprev,faft,D
real,allocatable,dimension(:,:):: Jac1,A,invsol
real:: delx,Lm,smallnum,tol

!Interface to vector function
interface
function fvector(x)
real,dimension(1000),intent(in)::x
real,dimension(1000)::fvector
end function fvector
end interface

!Allocate memmory 
allocate(Jac1(num,num))
allocate(A(num,num + 1))
allocate(x(num))
allocate(sol(num))
allocate(b(num))
allocate(ff(num))
allocate(fb(num))
allocate(invsol(num,num))
allocate(fprev(num))
allocate(faft(num))
allocate(D(num))


!Main Loop
maxiter = 1000
do l = 1,maxiter

!Determine Jacobian
do i = 1,num
    do j = 1,num
    x = xo
    x(j) = x(j) + delx       
    ff = fvector(x)

    fb = fvector(xo)
    Jac1(i,j) = (ff(i) - fb(i))/delx
    
    end do
end do

!$$$$$$ ! Partial Pivoting
!$$$$$$ do i=1,num
!$$$$$$     if(Jac1(i,i) == 0) then
!$$$$$$       D = Jac1(i,1:num)
!$$$$$$         if(i==1) then 
!$$$$$$           Jac1(i,1:num) = Jac1(i+1,1:num)
!$$$$$$           Jac1(i+1,1:num) = D 
!$$$$$$         end if
!$$$$$$           
!$$$$$$         if(i==num) then
!$$$$$$           Jac1(i,1:num) = Jac1(i-1,1:num)
!$$$$$$           Jac1(i+1,1:num) = D
!$$$$$$         end if
!$$$$$$         
!$$$$$$         if((i>1).AND.(i<num)) then
!$$$$$$         Jac1(i,1:num) = Jac1(i+1,1:num)
!$$$$$$         Jac1(i+1,1:num) = D
!$$$$$$         end if
!$$$$$$     end if
!$$$$$$ end do

!Determine first inverse
do inv = 1,num 

!Identity vector generator
do i = 1,num
  if(i == inv) then
    b(i) = 1.0
   else
    b(i) = 0.0
   end if
end do

!Create A Matrix
do i = 1,num
  do j = 1,num
    A(i,j) = Jac1(i,j)
  end do
end do

!Create Augmented Matrix
do i = 1,num
  A(i,num+1) = b(i)
end do

!Gauss Elimination Solution
do i = 1,num
    do j = (i+1),num
      
        Lm = A(j,i)/A(i,i) 
        
        if (.NOT.(A(i,i) == 0)) then
        do k=1,(num+1)
           A(j,k) = A(j,k) - Lm*A(i,k)
        end do   
        end if
    end do
end do

!Back Substitution
do i=1,num
  b(i) = A(i,num+1)
end do

do i = num,1,-1
    sol(i) = b(i)
    
    do j = num,i,-1
        if(j /= i) sol(i) = sol(i) - sol(j)*A(i,j)
    end do
    
    sol(i) = sol(i)/A(i,i)
end do

!Add inverse solution to final solution matrix
do i = 1,num
invsol(i,inv) = sol(i)
end do
end do

! Newton Update Function
fprev = fvector(xo)
xo = xo - matmul(invsol,fvector(xo))
faft = fvector(xo)

! Stopping criterion
if(sqrt(DOT_PRODUCT(fprev - faft,fprev - faft)) <= tol) exit

end do

!Deallocate memory
deallocate(Jac1)
deallocate(A)
deallocate(x)
deallocate(sol)
deallocate(b)
deallocate(ff)
deallocate(fb)
deallocate(invsol)
deallocate(fprev)
deallocate(faft)
deallocate(D)

end subroutine arjnr