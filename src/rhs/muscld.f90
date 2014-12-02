!==========================================================================
!>
!!  Koren's Differentiable limiter for 3rd order accuracy
!!  at cell face points
!!
!!  Modified for 2nd order at boundaries
!!
!!  J. Sitaraman 10/12/08
!!
!!
subroutine muscld(f,ql,qr,is,imax,th,qt,eps,ibmin,ibmax,mdim,nq)
!
implicit none
!      
integer, intent(in) :: nq            !< number of q-variables
integer, intent(in) :: mdim          !< dimension of variables
integer, intent(in) :: is,imax       !< start-end dimensions for reconstruction
integer, intent(in) :: ibmin,ibmax   !< boundary processing, =-2 implies do 2nd order extrapolation, = 0 implies 1st order
!      
real*8, intent(in) ::  th,qt,eps     !< parameters for limiter, th=1/3, qt=1/4 = 3rd order
real*8, intent(in) :: f(nq,mdim)     !< cell centered variables
real*8, intent(out) :: ql(nq,mdim),qr(nq,mdim) !< left and right state variables 
!
!  local variables
!
real*8, allocatable ::  f2(:,:)
integer i,n,im
real*8 thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt
!
!****************
!
im=imax-1
!
allocate(f2(nq,mdim))
!
! first order if qt=0
!
if(qt.eq.0.0)then
   do n=1,nq
      do i=is,imax
         ql(n,i)=f(n,i)
         qr(n,i)=f(n,i)
      enddo
   enddo
else
!
! do 3rd order otherwise
!
   thm    = 1.0 - th
   thp    = 1.0 + th
   
   do n=1,nq
      do i=is,im
         f2(n,i) = f(n,i+1) - f(n,i)
      enddo
      
      do I=IS+1,IM
         f2i    = f2(n,i)
         f2i1   = f2(n,i-1)
         a1     = 3.0*(f2i*f2i1+eps)
         a2     = 2.0*(f2i-f2i1)**2 + a1
         f3i    =  a1/a2 
         f3qt   = qt*f3i
         ql(n,i)= f(n,i)+f3qt*( thm*f2i1 + thp*f2i )
         qr(n,i)= f(n,i)-f3qt*( thp*f2i1 + thm*f2i )
      enddo
      !
      !..boundaries at present only 2nd order
      !
      ql(n,is)=f(n,is)
      qr(n,is)=f(n,is)
      ql(n,imax)=f(n,imax)
      qr(n,imax)=f(n,imax)
      !
      if(ibmin.eq.2) then
         f2i =f(n,is+1)-f(n,is)
         f2i1=f(n,is+2)-f(n,is+1)
         ql(n,is)=f(n,is)+0.5*(f2i*4./3.-f2i1*1./3.)
      endif

      if(ibmax.eq.2) then
         f2i  =f(n,imax)-f(n,imax-1)
         f2i1 =f(n,imax-1)-f(n,imax-2)
         qr(n,imax)=f(n,imax)-0.5*(f2i*4./3.-f2i1*1./3.)
      endif
      
   enddo
!
!...  Check for -ve pressure / density
!
   if(ql(1,is).le.0..or.ql(nq,is).le.0.) then
      ql(:,is)=f(:,is)
   endif
   
   if(qr(1,imax).le.0..or.qr(nq,imax).le.0.) then
      qr(:,imax)=f(:,imax)
   endif
endif
!
deallocate(f2)
!
return
end subroutine muscld

