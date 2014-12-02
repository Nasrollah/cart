!===================================================================!
!> \brief
!! This subroutine solves a block tridiagonal system of equations
!! using the Thomas algorithm
!!
!! Perhaps time using matmul vs loops
!!
!! Versions:\par
!!    - Leffell 09/25/2014
!!
!! Uses:\par
!!    
!!
!! Source code:\par
!!   \include blockThomas.f90
!!
!====================================================================!
subroutine blockThomas(a,b,c,f,nq,n)
!
implicit none
!
integer, intent(in)       :: nq                                  !< dimension of the block
integer, intent(in)       :: n                                   !< number of blocks
real*8,  intent(inout)    :: a(nq,nq,n),b(nq,nq,n),c(nq,nq,n)    !< left, central and right block arrays
real*8,  intent(inout)    :: f(nq,n)                             !< rhs array on entry and solution on exit
!
! local variables
!
integer :: i,j,k,l
real*8  :: b11,b21,b22,b31,b32,b33,b41,b42,b43,b44,b51,b52,b53,b54,b55
real*8  :: u12,u13,u14,u15,u23,u24,u25,u34,u35,u44,u45
real*8  :: d1,d2,d3,d4,d5
real*8  :: c1,c2,c3,c4,c5
real*8  :: t(5,5)
i = 1
   
! Decompose b(i) into L and U
b11 = 1.0d0/b(1,1,i)
u12 = b(1,2,i)*b11
u13 = b(1,3,i)*b11
u14 = b(1,4,i)*b11
u15 = b(1,5,i)*b11

b21 = b(2,1,i)
b22 = 1.0d0/(b(2,2,i)-b21*u12)
u23 = (b(2,3,i)-b21*u13)*b22
u24 = (b(2,4,i)-b21*u14)*b22
u25 = (b(2,5,i)-b21*u15)*b22

b31 = b(3,1,i)
b32 = b(3,2,i)-b31*u12
b33 = 1.0d0/(b(3,3,i)-b31*u13-b32*u23)
u34 = (b(3,4,i)-b31*u14-b32*u24)*b33
u35 = (b(3,5,i)-b31*u15-b32*u25)*b33

b41 = b(4,1,i)
b42 = b(4,2,i)-b41*u12
b43 = b(4,3,i)-b41*u13-b42*u23
b44 = 1.0d0/(b(4,4,i)-b41*u14-b42*u24-b43*u34)
u45 = (b(4,5,i)-b41*u15-b42*u25-b43*u35)*b44

b51 = b(5,1,i)
b52 = b(5,2,i)-b51*u12
b53 = b(5,3,i)-b51*u13-b52*u23
b54 = b(5,4,i)-b51*u14-b52*u24-b53*u34
b55 = 1.0d0/(b(5,5,i)-b51*u15-b52*u25-b53*u35-b54*u45)

d1 = f(1,i)*b11
d2 = (f(2,i)-b21*d1)*b22
d3 = (f(3,i)-b31*d1-b32*d2)*b33
d4 = (f(4,i)-b41*d1-b42*d2-b43*d3)*b44
d5 = (f(5,i)-b51*d1-b52*d2-b53*d3-b54*d4)*b55

f(5,i) = d5
f(4,i) = d4-u45*d5
f(3,i) = d3-u34*f(4,i)-u35*d5
f(2,i) = d2-u23*f(3,i)-u24*f(4,i)-u25*d5
f(1,i) = d1-u12*f(2,i)-u13*f(3,i)-u14*f(4,i)-u15*d5


! c(i) = Inv(b(i))*c(i)
do j = 1,nq
   c1 =  c(1,j,i)*b11
   c2 = (c(2,j,i)-b21*c1)*b22
   c3 = (c(3,j,i)-b31*c1-b32*c2)*b33
   c4 = (c(4,j,i)-b41*c1-b42*c2-b43*c3)*b44
   c5 = (c(5,j,i)-b51*c1-b52*c2-b53*c3-b54*c4)*b55
   
   c(5,j,i) = c5;
   c(4,j,i) = c4-u45*c5;
   c(3,j,i) = c3-u34*c(4,j,i)-u35*c5
   c(2,j,i) = c2-u23*c(3,j,i)-u24*c(4,j,i)-u25*c5
   c(1,j,i) = c1-u12*c(2,j,i)-u13*c(3,j,i)-u14*c(4,j,i)-u15*c5
enddo

do i = 2,n-1

   b(:,:,i) = b(:,:,i) - matmul(a(:,:,i),c(:,:,i-1))
   f(:,i)   = f(:,i)   - matmul(a(:,:,i),f(:,i-1))

   
   ! Decompose b(i) into L and U
   b11 = 1.0d0/b(1,1,i)
   u12 = b(1,2,i)*b11
   u13 = b(1,3,i)*b11
   u14 = b(1,4,i)*b11
   u15 = b(1,5,i)*b11
   
   b21 = b(2,1,i)
   b22 = 1.0d0/(b(2,2,i)-b21*u12)
   u23 = (b(2,3,i)-b21*u13)*b22
   u24 = (b(2,4,i)-b21*u14)*b22
   u25 = (b(2,5,i)-b21*u15)*b22

   b31 = b(3,1,i)
   b32 = b(3,2,i)-b31*u12
   b33 = 1.0d0/(b(3,3,i)-b31*u13-b32*u23)
   u34 = (b(3,4,i)-b31*u14-b32*u24)*b33
   u35 = (b(3,5,i)-b31*u15-b32*u25)*b33

   b41 = b(4,1,i)
   b42 = b(4,2,i)-b41*u12
   b43 = b(4,3,i)-b41*u13-b42*u23
   b44 = 1.0d0/(b(4,4,i)-b41*u14-b42*u24-b43*u34)
   u45 = (b(4,5,i)-b41*u15-b42*u25-b43*u35)*b44

   b51 = b(5,1,i)
   b52 = b(5,2,i)-b51*u12
   b53 = b(5,3,i)-b51*u13-b52*u23
   b54 = b(5,4,i)-b51*u14-b52*u24-b53*u34
   b55 = 1.0d0/(b(5,5,i)-b51*u15-b52*u25-b53*u35-b54*u45)

   d1 = f(1,i)*b11
   d2 = (f(2,i)-b21*d1)*b22
   d3 = (f(3,i)-b31*d1-b32*d2)*b33
   d4 = (f(4,i)-b41*d1-b42*d2-b43*d3)*b44
   d5 = (f(5,i)-b51*d1-b52*d2-b53*d3-b54*d4)*b55

   f(5,i) = d5
   f(4,i) = d4-u45*d5
   f(3,i) = d3-u34*f(4,i)-u35*d5
   f(2,i) = d2-u23*f(3,i)-u24*f(4,i)-u25*d5
   f(1,i) = d1-u12*f(2,i)-u13*f(3,i)-u14*f(4,i)-u15*d5

   
   ! c(i) = Inv(b(i))*c(i)
   do j = 1,nq
      c1 =  c(1,j,i)*b11
      c2 = (c(2,j,i)-b21*c1)*b22
      c3 = (c(3,j,i)-b31*c1-b32*c2)*b33
      c4 = (c(4,j,i)-b41*c1-b42*c2-b43*c3)*b44
      c5 = (c(5,j,i)-b51*c1-b52*c2-b53*c3-b54*c4)*b55
      
      c(5,j,i) = c5;
      c(4,j,i) = c4-u45*c5;
      c(3,j,i) = c3-u34*c(4,j,i)-u35*c5
      c(2,j,i) = c2-u23*c(3,j,i)-u24*c(4,j,i)-u25*c5
      c(1,j,i) = c1-u12*c(2,j,i)-u13*c(3,j,i)-u14*c(4,j,i)-u15*c5
   enddo
enddo

i = n
b(:,:,i) = b(:,:,i) - matmul(a(:,:,i),c(:,:,i-1))
f(:,i)   = f(:,i)   - matmul(a(:,:,i),f(:,i-1))

   
! Decompose b(i) into L and U
b11 = 1.0d0/b(1,1,i)
u12 = b(1,2,i)*b11
u13 = b(1,3,i)*b11
u14 = b(1,4,i)*b11
u15 = b(1,5,i)*b11

b21 = b(2,1,i)
b22 = 1.0d0/(b(2,2,i)-b21*u12)
u23 = (b(2,3,i)-b21*u13)*b22
u24 = (b(2,4,i)-b21*u14)*b22
u25 = (b(2,5,i)-b21*u15)*b22

b31 = b(3,1,i)
b32 = b(3,2,i)-b31*u12
b33 = 1.0d0/(b(3,3,i)-b31*u13-b32*u23)
u34 = (b(3,4,i)-b31*u14-b32*u24)*b33
u35 = (b(3,5,i)-b31*u15-b32*u25)*b33

b41 = b(4,1,i)
b42 = b(4,2,i)-b41*u12
b43 = b(4,3,i)-b41*u13-b42*u23
b44 = 1.0d0/(b(4,4,i)-b41*u14-b42*u24-b43*u34)
u45 = (b(4,5,i)-b41*u15-b42*u25-b43*u35)*b44

b51 = b(5,1,i)
b52 = b(5,2,i)-b51*u12
b53 = b(5,3,i)-b51*u13-b52*u23
b54 = b(5,4,i)-b51*u14-b52*u24-b53*u34
b55 = 1.0d0/(b(5,5,i)-b51*u15-b52*u25-b53*u35-b54*u45)

d1 = f(1,i)*b11
d2 = (f(2,i)-b21*d1)*b22
d3 = (f(3,i)-b31*d1-b32*d2)*b33
d4 = (f(4,i)-b41*d1-b42*d2-b43*d3)*b44
d5 = (f(5,i)-b51*d1-b52*d2-b53*d3-b54*d4)*b55

f(5,i) = d5
f(4,i) = d4-u45*d5
f(3,i) = d3-u34*f(4,i)-u35*d5
f(2,i) = d2-u23*f(3,i)-u24*f(4,i)-u25*d5
f(1,i) = d1-u12*f(2,i)-u13*f(3,i)-u14*f(4,i)-u15*d5

! Backward Sweep
do i = n-1,1,-1
   ! f(i) = f(i) - c(i+1)*f(i+1)
   f(:,i) = f(:,i) - matmul(c(:,:,i),f(:,i+1))
enddo
!
return
end subroutine blockThomas
