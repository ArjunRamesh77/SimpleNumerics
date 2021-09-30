!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Subroutine for Linear regression
!
!****************************************************************************

subroutine arjreg(x,y,num,order,c)
implicit none

!Variable declaration
real,dimension(num):: x
real,dimension(num):: y
real:: Lm
real,allocatable,dimension(:,:):: A
real,allocatable,dimension(:,:):: AtA,invsol
integer:: order,num,i,j
real,dimension(order + 1):: c

!Allocate memory
allocate(A(num,order + 1))
allocate(AtA(order+1,order+1))
allocate(invsol(order+1,order+1))

!Create A matrix
do i=1,num
  do j = 0,order

    if(j==0) then
      A(i,j+1) = 1.0
      else 
    A(i,j+1) = x(i)**(j)
    end if
  end do
end do

!Determine AtA
AtA = MATMUL(TRANSPOSE(A),A)

!Determine inverse
call arjminv(AtA,order + 1,invsol)

!Determine coefficients
c = MATMUL(MATMUL(invsol,TRANSPOSE(A)),y)

!Deallocate memory
deallocate(A)
deallocate(AtA)
deallocate(invsol)

end subroutine arjreg


!MATRIX INVERSION
subroutine arjminv(Aa,num,invsol)
implicit none

!Variable declaration
integer:: num,inv,i,j,k
real:: Lm
real,allocatable,dimension(:):: b,sol,D
real,dimension(num,num):: Aa,invsol
real,allocatable,dimension(:,:):: A

!Allocate memory
allocate(b(num))
allocate(sol(num))
allocate(A(num,num+1))
allocate(D(num))


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

! Partial Pivoting
do i=1,num
    if(Aa(i,i) == 0) then
      D = Aa(i,1:num)
        if(i==1) then 
          Aa(i,1:num) = Aa(i+1,1:num)
          Aa(i+1,1:num) = D 
        end if
          
        if(i==num) then
          Aa(i,1:num) = Aa(i-1,1:num)
          Aa(i+1,1:num) = D
        end if
        
        if((i>1).AND.(i<num)) then
        Aa(i,1:num) = Aa(i+1,1:num)
        Aa(i+1,1:num) = D
        end if
    end if
end do

!Create A Matrix
do i = 1,num
  do j = 1,num
    A(i,j) = Aa(i,j)
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

!deallocate memory
deallocate(b)
deallocate(sol)
deallocate(A)
deallocate(D)

end subroutine arjminv