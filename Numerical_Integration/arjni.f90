!****************************************************************************
!
!  AUTHOR: 		ARJUN RAMESH
!  DATE: 		01/01/2015
!  DESCRIPTION: Subroutine for Numerical integration
!
!****************************************************************************

subroutine arjsimp(y,h,num,int)
implicit none

!Variable declaration
integer:: num,numpsc,max38,max13,i,j,set,numspc
real,dimension(num):: x,y
real:: h,int38,int13,int,inttr

!Determine parameters
numspc = num - 1
max38 = numspc - mod(numspc,3) 

!Main loop - (3/8)th rule
int38 = 0.0
do i = 1,max38 + 1

if((i == 1).or.(i == (max38 + 1))) then
  int38 = int38 + y(i)
else
if(mod(i-1,3) == 0) then
  int38 = int38 + 2*y(i)
else
  int38 = int38 + 3*y(i)
end if
end if
end do

! Add (3/8)th rule solution  
int = ((3.0*h)/8.0)*int38

! Simpsons (1/3)rd rule solution
int13 = 0.0
if(mod(numspc,3) == 2) then
set = 1
do j = max38 + 1,num

if((j == (max38 + 1)).or.(j == num)) then
  int13 = int13 + y(j)
else
if(set==1) then
  int13 = int13 + 4*y(j)
  set = 0
else
  int13 = int13 + 2*y(j)
  set = 1
end if
end if
end do
end if
 
int13 = (h/3.0)*int13
int = int + int13

! Trapezoidal rule
inttr = 0.0
if(mod(numspc,3) == 1) inttr = (h/2.0)*(y(num) + y(num - 1))
int = int + inttr

end subroutine arjsimp