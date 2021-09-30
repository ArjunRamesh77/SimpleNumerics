!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		12/22/2014
!  DESCRIPTION: Subroutine to peroform Gauss-Jordan Elimination
!
!****************************************************************************

! GAUSS-JORDAN ELIMINATION ROUTINE
subroutine arjgj(A,b,num,sol)
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

!Gauss Elimination 1
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

!Gauss Elimination 2
do i=num,1,-1
    do j= i-1,1,-1
        
		Lm = A(j,i)/A(i,i)
        
		if (.NOT.(A(i,i) == 0)) then
        do k = 1,(num+1)
           A(j,k) = A(j,k) - Lm*A(i,k)
        end do   
        end if
    end do
end do

!Back Substitution
do i=num,1,-1
	sol(i) = A(i,num + 1)/A(i,i)
end do


end subroutine arjgj