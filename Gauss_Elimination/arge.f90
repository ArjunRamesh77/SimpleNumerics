!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Subroutine to peroform Gauss Elimination
!
!****************************************************************************

! GAUSS ELIMINATION ROUTINE
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