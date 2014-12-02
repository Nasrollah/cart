!
! 5th order weno scheme
!
subroutine fweno(f,ql,qr,is,ie,th,qt,eps,ibmin,ibmax,mdim,nvar)
!      
implicit none
integer, intent(in) :: is,ie,ibmin,ibmax,mdim,nvar
      
real*8, intent(in) :: th,qt,eps
real*8, intent(in) ::  f(nvar,mdim)
real*8, intent(out) :: ql(nvar,mdim),qr(nvar,mdim)
!
!..   local variables
!
integer i,n,im
real*8, allocatable :: f0(:),f1(:),f2(:),f2m(:),slope(:)
real*8 ammd,ammd4,weno5
real*8 f0bmin,f1bmin,f2bmin,f0bmax,f1bmax,f2bmax
real*8 at,s1,ati1,t1,at1
!
allocate(f0(mdim),f1(mdim),f2(mdim),f2m(mdim),slope(mdim))
!
!..this is just 1st order upwind
!
if(qt.eq.0.0)then
   do n=1,nvar
      do i=is,ie
         ql(n,i)=f(n,i)
         qr(n,i)=f(n,i)
      enddo
   enddo
   return
endif

!
! 5th order weno scheme follows
!
im=ie-1
!
do n=1,nvar
!
!..let's load up a few difference arrays
!
   do i=is,ie
      f0(i) = f(n,i)
   enddo
   !
   !..1st difference at i+1/2
   !
   do i=is,ie-1
      f1(i) = f0(i+1) - f0(i)
   enddo
   !
   !..2nd difference at i
   !
   do i=is+1,ie-1
      f2(i)  = f1(i) - f1(i-1)
   enddo
   !
   !..extrapolate at the boundaries
   !
   f2(is) = 2.*f2(is+1)-f2(is+2)
   f2(ie) = 2.*f2(ie-1)-f2(ie-2)
   !
   !..modify at boundaries, if needed
   !
   if(ibmin.eq.-2) then
      f0bmin = f0(is) !fmin(n)
      f1bmin = f0(is) - f0bmin
      f2(is) = f1(is) - f1bmin
      f2bmin = 2.*f2(is)-f2(is+1)
   endif
   if(ibmax.eq.-2) then
      f0bmax = f0(ie) !fmax(n)
      f1bmax = f0bmax - f0(ie)
      f2(ie) = f1bmax - f1(ie-1)
      f2bmax = 2.*f2(ie)-f2(ie-1)
   endif
   !
   !..limit 2nd difference to i+1/2
   !
   do i = is,ie-1
      f2m(i) = ammd(f2(i),f2(i+1))
   enddo
   !
   !..now combine everything to get ql and qr
   !
   !..first-order at boundary?
   !
   at    = f1(is) - 0.5*f2m(is)
   slope(is) = ammd(at,2.*f1(is))
   ql(n,is) = f0(is) + 0.5*slope(is)
   !
   if(ibmin.eq.-2) then
      i = is
      s1    = ammd(f1(i),f1bmin)
      at    = f1(i)  - 0.5*f2m(i)
      ati1  = f1bmin + 0.5*ammd(f2(i),f2bmin)
      t1    = ammd(at,ati1)
      slope(i)  =  sign(1d0,t1)* &
               dmin1(0.5*abs(at+ati1),dmax1(2.*abs(s1),abs(t1)))
      ql(n,i) = f0(i) + 0.5*slope(i)
   endif
   !
   at1   = f1(ie-1) + 0.5*f2m(ie-1)
   slope(ie) = ammd(at1,2.*f1(ie-1))
   qr(n,ie) = f0(ie) - 0.5*slope(ie)
   !
   if(ibmax.eq.-2) then
      i = ie
      s1    = ammd(f1bmax,f1(i-1))
      at    = f1bmax  - 0.5*ammd(f2(i-1),f2bmax)
      ati1  = f1(i-1) + 0.5*f2m(i-1)
      t1    = ammd(at,ati1)
      slope(i)  =  sign(1d0,t1)* &
              dmin1(0.5*abs(at+ati1),dmax1(2.*abs(s1),abs(t1)))
      qr(n,i) = f0(i) - 0.5*slope(i)
   endif
   !
   !..sonic-a near the boundary?
   !
   do i=is+1,ie-1,max(ie-is-2,1)
      !
      !limited curvatures to calculate new slopes
      !
      at    = f1(i)   - 0.5*f2m(i)
      ati1  = f1(i-1) + 0.5*f2m(i-1)
      !
      !..limit slopes at i
      !
      s1    = ammd(f1(i),f1(i-1))
      t1    = ammd(at,ati1)
      !
      !..now find appropriate slope
      !
      slope(i)  =  sign(1d0,t1)* &
           dmin1(0.5*abs(at+ati1),dmax1(2.*abs(s1),abs(t1)))
      !
      !..use slope to calculate ql and qr
      !
      ql(n,i) = f0(i) + 0.5*slope(i)
      qr(n,i) = f0(i) - 0.5*slope(i)
   enddo
   !
   !..suresh at interior
   !
   do i=is+2,ie-2
      ql(n,i) = weno5(f0(i-2),f0(i-1),f0(i),f0(i+1),f0(i+2))
      qr(n,i) = weno5(f0(i+2),f0(i+1),f0(i),f0(i-1),f0(i-2))
   enddo
enddo
deallocate(f0,f1,f2,f2m,slope)
return
end subroutine fweno
!
! weno5 function
!
function weno5(a,b,c,d,e)
!
implicit none
!
real*8 a,b,c,d,e
real*8 b1,b2,epsw,djm1,ejm1,dj,ej
real*8 djp1,ejp1,dis0,dis1,dis2,q30,q31,q32
real*8 d01,d02,a1ba0,a2ba0,w0,w1,w2,weno5
!     
b1 = 13./12. 
b2 = 1./6.
epsw = 1.e-6
djm1 = a-2.*b+c
ejm1 = a-4.*b+3.*c
dj   = b-2.*c+d
ej   = b-d
djp1 = c-2.*d+e
ejp1 = 3.*c-4.*d+e
dis0 = b1*djm1*djm1+0.25*ejm1*ejm1+epsw
dis1 = b1*dj*dj+0.25*ej*ej+epsw
dis2 = b1*djp1*djp1+0.25*ejp1*ejp1+epsw
q30 = 2.*a-7.*b+11.*c
q31 = -b+5.*c+2.*d
q32 = 2.*c+5.*d-e
d01 = dis0/dis1
d02 = dis0/dis2
a1ba0 = 6.*d01
a2ba0 = 3.*d02
w0 = 1./(1.+a1ba0+a2ba0)
w1 = a1ba0*w0
w2 = 1.-w0-w1
weno5 = b2*(w0*q30+w1*q31+w2*q32)

return
end function weno5
!
! min mod limiter
!
function ammd(fl,fr)
  real*8 fl,fr,ammd
  
  ammd=0.5*(sign(1d0,fl)+sign(1d0,fr))*dmin1(abs(fl),abs(fr))
  
  return
end function ammd
