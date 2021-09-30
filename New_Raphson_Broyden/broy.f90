! NEWTON RAPHSON SOLVER
subroutine arjnewton(fvector,num,xo,delx,tol)
implicit none

!Variable declarations
integer:: num,i,j,k,l,inv,maxiter
real,dimension(num):: xo
real,allocatable,dimension(:):: x,sol,b,ff,fb,fprev,faft,D,xprev,dfvector,delxo,p1
real,allocatable,dimension(:,:):: Jac1,A,invsol,op
real:: delx,Lm,tol

!Interface to vector function
interface
function fvector(x)
real,dimension(100),intent(in)::x
real,dimension(100)::fvector
end function fvector
end interface

!Allocate memmory 
allocate(Jac1(num,num))
allocate(x(num))
allocate(ff(num))
allocate(fb(num))
allocate(invsol(num,num))
allocate(fprev(num))
allocate(faft(num))
allocate(xprev(num))
allocate(dfvector(num))
allocate(delxo(num))
allocate(p1(num))
allocate(op(num,num))

!Determine First Jacobian
do i = 1,num
    do j = 1,num
    x = xo
    x(j) = x(j) + delx       
    ff = fvector(x)

    fb = fvector(xo)
    Jac1(i,j) = (ff(i) - fb(i))/delx
    
    end do
end do


!Determine first inverse
call arjminv(Jac1,num,invsol)

maxiter = 100
xprev = xo - delx

! Main Loop
do l = 1,maxiter

delxo = xo - xprev
dfvector = fvector(xo) - fvector(xprev)

! Sherman-Morrison Inverse Update
!Intermittent Product
p1 = ((delxo - matmul(invsol,dfvector))/dot_product(matmul(delxo,invsol),dfvector))

do i=1,size(p1)
  do j=1,size(delxo)
    op(i,j) = p1(i)*delxo(j)
end do
end do

invsol = invsol + matmul(op,invsol)

! Newton Update Function
xprev = xo
fprev = fvector(xo)
xo = xo - matmul(invsol,fvector(xo))
faft = fvector(xo)

! Stopping criterion
if(sqrt(DOT_PRODUCT(faft - fprev,faft - fprev)) <= tol) exit

end do

!Deallocate memory
deallocate(Jac1)
deallocate(x)
deallocate(ff)
deallocate(fb)
deallocate(invsol)
deallocate(fprev)
deallocate(faft)

end subroutine arjnewton

! Determine Inverse
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

