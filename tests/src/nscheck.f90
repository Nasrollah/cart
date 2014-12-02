!==============================================================================
program nscheck
!===============================================================================
! This is the main program for verifying accuracy of the
! the inviscid and viscous RHS
!
use NS_mms
implicit none
!include 'mpif.h'
!
real*8, allocatable :: q(:,:,:,:),x(:,:,:,:),qwork(:,:,:,:)
real*8, allocatable :: spec(:,:,:)
real*8, allocatable :: rhs(:,:,:,:),rhsv(:,:,:,:)
real*8, allocatable :: rhsexact(:,:,:,:),rhsvexact(:,:,:,:)
!
integer :: j,k,l,n
integer :: jmax,kmax,lmax
integer :: nq,nvar
integer :: fluxOrder,dissOrder,viscOrder
integer :: p
integer :: ej(4)
!
real*8  :: xlim,ylim,zlim,dx,dy,dz,vol,fsmach,rey,pr,prtr,gamma
real*8  :: t1,t2,dissCoef,e1,e2,t3,t4,dxold
real*8  :: tm(3),error(3),errold(3),slope(3)
!
character*16 :: scheme,reconstruction,istor
!
namelist /inputs/ fluxOrder,dissOrder,dissCoef,scheme,&
        reconstruction,viscOrder,istor
!
! begin
! 
read(5,inputs)
write(6,inputs)
!
gamma=1.4d0
fsmach=0.5d0
rey=1000
pr=0.72d0
prtr=0.9d0
call set_mms_params(gamma,pr,prtr,rey)
!
jmax=40
kmax=40
lmax=40
nq=6
nvar=5
!
do p=1,10
   !
   if (istor=='column') then
      allocate(q(jmax,kmax,lmax,nq),rhs(jmax,kmax,lmax,nq),x(jmax,kmax,lmax,3))
      allocate(qwork(jmax,kmax,lmax,9))
      allocate(rhsexact(jmax,kmax,lmax,nq))
      allocate(spec(jmax,kmax,lmax))
      allocate(rhsv(jmax,kmax,lmax,nq))
      allocate(rhsvexact(jmax,kmax,lmax,nq))
   else
      allocate(q(nq,jmax,kmax,lmax),rhs(nq,jmax,kmax,lmax),x(3,jmax,kmax,lmax))
      allocate(qwork(9,jmax,kmax,lmax))
      allocate(rhsexact(nq,jmax,kmax,lmax))
      allocate(spec(jmax,kmax,lmax))
      allocate(rhsv(nq,jmax,kmax,lmax))
      allocate(rhsvexact(nq,jmax,kmax,lmax))
   endif
   !
   xlim=1d0
   ylim=1d0
   zlim=1d0
   !
   dx=xlim/(jmax)
   dy=ylim/(kmax)
   dz=zlim/(lmax)
   tm=(/0.,0.,0./)
   vol=dx*dy*dz
   !
   ! define the points in the grid
   !
   do l=1,lmax
      do k=1,kmax
         do j=1,jmax
            if (istor=='column') then
               x(j,k,l,1)=(j-1)*dx!-xlim*0.5d0
               x(j,k,l,2)=(k-1)*dy!-ylim*0.5d0
               x(j,k,l,3)=(l-1)*dz!-zlim*0.5d0
               call getqreal(q(j,k,l,1:5),x(j,k,l,1),x(j,k,l,2),x(j,k,l,3))
               q(j,k,l,6)=0d0
               rhs(j,k,l,:)=0d0
               rhsexact(j,k,l,:)=0d0
               rhsv(j,k,l,:)=0d0
               rhsvexact(j,k,l,:)=0d0
            else
               x(1,j,k,l)=(j-1)*dx!-xlim*0.5d0
               x(2,j,k,l)=(k-1)*dy!-ylim*0.5d0
               x(3,j,k,l)=(l-1)*dz!-zlim*0.5d0
               call getqreal(q(1:5,j,k,l),x(1,j,k,l),x(2,j,k,l),x(3,j,k,l))
               q(6,j,k,l)=0d0
               rhs(:,j,k,l)=0d0
               rhsexact(:,j,k,l)=0d0
               rhsv(:,j,k,l)=0d0
               rhsvexact(:,j,k,l)=0d0
            endif
         enddo
      enddo
   enddo
   !
   if (scheme=='central') then
      call cpu_time(t1)
      call inviscidRHSunified(nq,nvar,gamma,q,rhs,spec,tm,dx,dy,dz,jmax,kmax,lmax,&
           fluxOrder,dissOrder,dissCoef,istor)
      call cpu_time(t2)
   elseif (scheme=='upwind') then
      call cpu_time(t1)
      call inviscidRHSupwind(nq,nvar,gamma,q,rhs,spec,tm,dx,dy,dz,jmax,kmax,lmax,&
        reconstruction,istor)
      call cpu_time(t2)
   endif
   !
   !
   call inviscidDivergence(x,rhsexact,nq,jmax*kmax*lmax,istor)
   !
   call cpu_time(t3)
   call viscousRHS(rey,pr,prtr,nq,nvar,gamma,q,qwork,rhsv,dx,dy,dz,jmax,kmax,lmax,&
        min(4,viscOrder),istor)
   !call viscousRHSfast(rey,pr,prtr,nq,nvar,gamma,q,rhsv,dx,dy,dz,jmax,kmax,lmax)
   !
   call cpu_time(t4)
   !
   call viscousDivergence(x,rhsvexact,nq,jmax*kmax*lmax,istor)
   !   
   error=0d0
   !
   ! check the error between discrete and exact divergence of
   ! inviscid and viscous fluxes
   !
   do l=4,lmax-3
      do k=4,kmax-3
         do j=4,jmax-3
            do n=1,5
               if (istor=='column') then
                  e1=(rhs(j,k,l,n)/vol-rhsexact(j,k,l,n))
                  e2=(rhsv(j,k,l,n)/vol-rhsvexact(j,k,l,n))
               else
                  e1=(rhs(n,j,k,l)/vol-rhsexact(n,j,k,l))
                  e2=(rhsv(n,j,k,l)/vol-rhsvexact(n,j,k,l))
               endif
               error(1)=max(error(1),abs(e1)) ! inf-norm of inviscid flux divergence
               if (abs(e2) > error(2)) then
                  ej=(/j,k,l,n/)
               endif
               error(2)=max(error(2),abs(e2)) ! inf-norm of viscous flux divergence
               error(3)=max(error(3),abs(e1+e2)) ! overall inf-norm
            enddo
         enddo
      enddo
   enddo
   !   
   if (p > 1) then
      slope=(log(errold)-log(error))/(log(dxold)-log(dx))
      write(6,"(3(1x,F8.4),6(1x,E10.3))") slope,error,(t2-t1)/(nvar*jmax*kmax*lmax),&
        (t4-t3)/(nvar*jmax*kmax*lmax),(t4-t3+t2-t1)/(nvar*jmax*kmax*lmax)
   else
      write(6,*)
      write(6,*) 'Accuracy check for Navier-Stokes discretization,&
                   storage=',istor
      write(6,*) 'Timing is reported as time per node per field'
      write(6,*)
      write(6,"(100(A1))") ('=',j=1,100)
      write(6,"(9(A8,2x))") 'invisc','visc','overall','invisc','visc',' overall',&
           'invisc','visc',' overall'
      write(6,"(9(A8,2x))") 'order','order','order','error','error','error',&
           'time','time','time'
      write(6,"(100(A1))") ('=',j=1,100)
   endif
   !
   dxold=dx
   errold=error
   !
   deallocate(q,qwork,rhs,x,rhsexact,spec,rhsv,rhsvexact)
   !
   jmax=jmax*1.1
   kmax=kmax*1.1
   lmax=lmax*1.1
   !
enddo
!
end program
