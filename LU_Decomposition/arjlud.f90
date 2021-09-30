!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Subroutine to perform arjlud.dll
!
!****************************************************************************

! LU-Decomposition routine
subroutine arjlud(A,L,U,num)
implicit none

!Variable declaration
integer:: num,i,j,k,n
real,dimension(num,num):: A,L,U
real,allocatable,dimension(:):: Lm,D

!Allocate memory
allocate(Lm(num))
allocate(D(num))

!Partial Pivoting
do i=1,num
    if(A(i,i) == 0) then
      D = A(i,1:num)
        if(i==1) then 
          A(i,1:num) = A(i+1,1:num)
          A(i+1,1:num) = D 
        end if
          
        if(i==num) then
          A(i,1:num) = A(i-1,1:num)
          A(i+1,1:num) = D
        end if
        
        if((i>1).AND.(i<num)) then
        A(i,1:num) = A(i+1,1:num)
        A(i+1,1:num) = D
        end if
    end if
end do

! Upper Triangular Matrix
n = 1
do i=1,num
    do j=i+1,num
        Lm(n) = A(j,i)/A(i,i)
        do k = 1,num
           A(j,k) = A(j,k) - Lm(n)*A(i,k)
        end do
            n = n+1;
    end do
end do

!Extract Lower Triangular Matrix
do i=1,num
  do j=1,num
    L(i,j) = 0.0
  end do
end do

! Force diagonals elements to 1
do i=1,num
    do j=1,num
        if(i==j) L(i,j) = 1.0
    end do
end do

! Force Lower Triangular elements to Lambda
n = 1
do i=1,num
    do j=i+1,num
        L(j,i) = Lm(n)
        n=n+1;
    end do
end do

! Extract Lower Triangular Matrix
do i=1,num
    do j=i+1,num
        A(j,i) = 0.0
    end do
end do

U = A
!Aeallocate memory
deallocate(Lm)
deallocate(D)

end subroutine arjlud

!Gauss Elimintation
subroutine arjge(A,b,num,sol)
implicit none

!Variable declarations
real,dimension(num,num+1):: A
real,dimension(num):: b,sol,D
real:: Lm
integer:: num,i,j,k

!Partial Pivoting
do i=1,num
    if(A(i,i) == 0) then
      D = A(i,1:num)
        if(i==1) then 
          A(i,1:num) = A(i+1,1:num)
          A(i+1,1:num) = D 
        end if
          
        if(i==num) then
          A(i,1:num) = A(i-1,1:num)
          A(i+1,1:num) = D
        end if
        
        if((i>1).AND.(i<num)) then
        A(i,1:num) = A(i+1,1:num)
        A(i+1,1:num) = D
        end if
    end if
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

end subroutine arjge