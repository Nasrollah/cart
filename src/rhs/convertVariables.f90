subroutine convertVariables(gamma,nq,jmax,kmax,lmax,q,imode,istor)
!
implicit none
!
real*8, intent(in) :: gamma                      !< ratio of specific heats
integer, intent(in) :: jmax                      !< coordinate 1 dimension
integer, intent(in) :: kmax                      !< coordinate 2 dimension
integer, intent(in) :: lmax                      !< coordinate 3 dimension
integer, intent(in) :: nq                        !< number of q-variables
real*8, intent(inout) :: q(nq*jmax*kmax*lmax)    !< q-data (conservative or primitive)
integer, intent(in) :: imode                     !< imode (0= cons to prim, 1= prim to cons)
character*(*), intent(in) :: istor                !< storage type 'row' or 'column'
                                                  !< row=> q(:,j,k,l), column=>q(j,k,l,:
!
integer :: j,k,l,n,iq,iloc
integer :: qstride,qmult
real*8  :: gm1,ggm1,ggm1i
real*8  :: qprim(5)
real*8  :: rinv,usq
!
gm1=gamma-1.0d0
ggm1=gamma*gm1
ggm1i=1d0/ggm1
!
! use storage string to set the variable
! location order
!
if (istor=='row') then
   qstride=1
   qmult=nq
else
   qstride=jmax*kmax*lmax
   qmult=1
endif
!
! compute velocities and temperature
! from conservative variables
! and replace q(2:5) with those
!
if (imode==0) then
   iq=0
   do l=1,lmax
      do k=1,kmax
         do j=1,jmax
            !
            ! get conservative variables here
            !
            iloc=iq*qmult+1
            
            do n=1,5
               qprim(n)=q(iloc)
               iloc=iloc+qstride
            end do
            !
            usq=0.
            !
            ! use qwork(1:4) as temporary storage for 
            ! primitives
            !
            !
            rinv=1./qprim(1)
            qprim(2)=qprim(2)*rinv !< u-velocity
            qprim(3)=qprim(3)*rinv !< v-velocity
            qprim(4)=qprim(4)*rinv !< w-velocity
            usq=qprim(2)**2+qprim(3)**2+qprim(4)**2
            qprim(5)=ggm1*(qprim(5)*rinv-0.5*usq) !< temperature
            !
            iloc=iq*qmult+1
            iloc=iloc+qstride
            !
            do n=2,5
               q(iloc)=qprim(n)
               iloc=iloc+qstride
            end do
            !
            iq=iq+1
            !
         enddo
      enddo
   enddo
else
   ggm1=1./ggm1
   iq=0
   do l=1,lmax
      do k=1,kmax
         do j=1,jmax
            !
            ! get conservative variables here
            !
            iloc=iq*qmult+1
            
            do n=1,5
               qprim(n)=q(iloc)
               iloc=iloc+qstride
            end do
            !
            usq=0.
            !
            ! use qwork(1:4) as temporary storage for 
            ! primitives
            !
            usq=qprim(2)**2+qprim(3)**2+qprim(4)**2
            qprim(2)=qprim(2)*qprim(1) !< u-momentum
            qprim(3)=qprim(3)*qprim(1) !< v-momentum
            qprim(4)=qprim(4)*qprim(1) !< w-momentum
            qprim(5)=qprim(5)*qprim(1)*ggm1i+0.5*usq*qprim(1) !< total energy
            !
            iloc=iq*qmult+1
            iloc=iloc+qstride
            !
            do n=2,5
               q(iloc)=qprim(n)
               iloc=iloc+qstride
            end do
            !
            iq=iq+1
            !
         enddo
      enddo
   enddo
endif
!
return
end subroutine convertVariables

