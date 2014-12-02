!==============================================================================
program linearizationTest
  !===============================================================================
  ! This program verifies the linearization of the Central + Diss
  ! scheme invoked by inviscidRHSunified.f90
  !
  implicit none
  !
  real*8, allocatable :: q(:,:,:,:),x(:,:,:,:),spaceMetric(:,:,:,:)
  real*8, allocatable :: timeMetric(:,:,:,:)
  real*8, allocatable :: spec(:,:,:)
  real*8, allocatable :: tscale(:,:,:)
  real*8, allocatable :: qdq(:,:,:,:),dq(:,:,:,:),ddq(:,:,:,:)
  real*8, allocatable :: rhs1(:,:,:,:)
  real*8, allocatable :: rhs(:,:,:,:),rhsv(:,:,:,:),qstar(:,:,:,:) 
  real*8, allocatable :: rhs_q(:,:,:,:),rhs_qdq(:,:,:,:),res(:,:)
  real*8, allocatable :: res_r(:),res_l(:)
  !
  integer :: j,k,l,i,n,nq,nvar,n_dq
  integer :: jmax,kmax,lmax,nsteps
  integer :: fluxOrder,dissOrder,bctyp
  integer, parameter :: out_unit=20
  !
  real*8  :: tf,t,h,cfl,a1,a2,a3,a4,h1
  real*8  :: dx,dy,dz,second,fourth,sixth,vol
  real*8  :: fsmach,rey,totime,viscous
  real*8  :: alpha=0.
  real*8  :: beta=0.
  real*8  :: gamma=1.4
  real*8  :: amp,freq
  real*8  :: pi=acos(-1.)
  real*8  :: norm,two_norm,dissCoef,s_dq,slope(3)
  real*8  :: max_res,max_res_r,max_res_l,res_rl
  !
  character*80 :: integer_string,xfile,qfile,fmt
  integer :: ii,jj,kk,ll,iloc,myid,numprocs,ierr,max_j,max_k,max_l

  !
  ! begin
  ! 
  n_dq = 30
  !
  nq=5
  fsmach=0.5
  jmax=40
  kmax=40
  lmax=40
  !
  allocate(q(5,jmax,kmax,lmax),rhs(5,jmax,kmax,lmax),x(3,jmax,kmax,lmax))
  allocate(ddq(5,jmax,kmax,lmax),qdq(5,jmax,kmax,lmax))
  allocate(dq(5,jmax,kmax,lmax))
  allocate(rhs1(5,jmax,kmax,lmax),rhs_q(5,jmax,kmax,lmax),rhs_qdq(5,jmax,kmax,lmax))
  allocate(spaceMetric(10,jmax,kmax,lmax),timeMetric(3,jmax,kmax,lmax))
  allocate(spec(jmax,kmax,lmax))
  allocate(tscale(jmax,kmax,lmax))
  allocate(res(n_dq,4))
  allocate(res_r(jmax*kmax*lmax),res_l(jmax*kmax*lmax))
  !
  dx=10./(jmax)
  dy=10./(kmax)
  dz=10./(lmax)
  nvar=5
  vol=dx*dy*dz
  !
  ! define the points in the grid
  !
  do l=1,lmax
     do k=1,kmax
	do j=1,jmax
           x(1,j,k,l)=(j-1)*dx-5.0
           x(2,j,k,l)=(k-1)*dy-5.0
           x(3,j,k,l)=(l-1)*dz-5.0
	enddo
     enddo
  enddo
  !
  ! initialize data
  !
  call init_data(q,fsmach,alpha,beta,gamma,jmax,kmax,lmax)
  call random_seed()
  call random_number(ddq)
!  q = q + 5e-4*ddq
  
  call init_metrics(spaceMetric,timeMetric,dx,dy,dz,jmax,kmax,lmax)
  !
  t=0
  !
  ! amplitude and frequency of oscillation
  !
  amp=0.01
  freq=1
  !
  t=(2*pi/freq)
  h = dx*cfl/(1.0+fsmach)
  tscale=h
  bctyp = 0

  ! Numerical Scheme Inputs (to be inputs)
  dissOrder = 2
  dissCoef  = 0.5d0!0.5d0
  write(*,*) "rho = ", q(1,20,20,20)
  do i=1,3
     
     fluxOrder = 2*i

     ! Evaluate R(q)
     call bc_per1(q,  x,fsmach,gamma,amp,freq,t,jmax,kmax,lmax,bctyp,fluxOrder,dissOrder)
     call inviscidRHSunified2(nq,nvar,gamma,q,  rhs_q,  spec,timeMetric,dx,dy,dz,jmax,kmax,lmax,fluxOrder,dissOrder,dissCoef,'row')

    ! write(8,*) "rhs_q(5,20,20,20): ", rhs_q(5,20,20,20)

     ! Cycle over range of dq values
     s_dq = 1.0d0

     do j=1,n_dq
        
        s_dq = 0.5*s_dq
        
        ! Preturb only single q value
        dq = 0.0d0
        dq(1,20,20,20) = s_dq
        !dq  = s_dq
 
        ! Update qdq = q + dq 
        qdq = q + dq
        
        !call bc_per1(qdq,x,fsmach,gamma,amp,freq,t,jmax,kmax,lmax,bctyp,fluxOrder,dissOrder)
        
        ! Evaluate R(q+dq)
        call inviscidRHSunified2(nq,nvar,gamma,qdq,rhs_qdq,spec,timeMetric,&
	   dx,dy,dz,jmax,kmax,lmax,fluxOrder,dissOrder,dissCoef,'row')
        
        !write(8,*) "rhs_qdq(5,20,20,20): ", rhs_qdq(5,20,20,20)

        ! Evaluate (dR/dq)*delta_q               
        call fluxJacobianMatVec(nq,nvar,gamma,q,rhs,spec,tscale,timeMetric,dx,dy,dz,&
	  jmax,kmax,lmax,fluxOrder,dissOrder,dissCoef,'row',dq)
        
        !write(8,*) "Adq(5,20,20,20): ", rhs(5,20,20,20)

        ! Store R(q+dq)-R(q) and dRdQ*dQ for coarsest and second order                                   
        if (i.eq.1.and.j.eq.1) then
           iloc = 0
           max_res   = 0.0d0
           max_res_r = 0.0d0
           max_res_l = 0.0d0
           do ll = 1,lmax
              do kk = 1,kmax
                 do jj = 1,jmax
                    iloc = iloc + 1
                    if (jj.gt.2.and.kk.gt.2.and.ll.gt.2.and.jj.lt.(jmax-1).and.kk.lt.(kmax-1).and.ll.lt.(lmax-1)) then
                       res_r(iloc) = rhs_qdq(1,jj,kk,ll)-rhs_q(1,jj,kk,ll)
                       if (abs(res_r(iloc)).gt.max_res_r) then
                          max_res_r = abs(res_r(iloc))
                       endif
                       
                       res_l(iloc) = rhs(1,jj,kk,ll)
                       if (abs(res_l(iloc)).gt.max_res_l) then
                          max_res_l = abs(res_l(iloc))
                       endif
                       
                       res_rl = res_r(iloc)-res_l(iloc)
                       if (abs(res_rl).gt.max_res) then
                          max_res = abs(res_rl)
                          max_j = jj
                          max_k = kk
                          max_l = ll
                       endif
                       if (abs(res_rl).gt.0.0d0) then
                          write(*,*) "Nonzero residual at: "
                          write(*,*) "j = ",jj
                          write(*,*) "k = ",kk
                          write(*,*) "l = ",ll
                          write(*,*) "with value: ", res_rl
                          write(*,*) "rhs: ", res_r(iloc)
                          write(*,*) "lhs: ", res_l(iloc)
                       endif
                    endif
                 enddo
              enddo
           enddo
           
           write(*,*) "Checking Second Order Linearization"
           write(*,*) "DeltaQ = ", s_dq
           write(*,*) "Maximum Residual in R(q+dq)-R(q) is: ", max_res_r
           write(*,*) "Maximum Residual in dRdQ*dQ is: ", max_res_l
           write(*,*) "Maximum Total Residual is: ", max_res
           write(*,*) "Location of Maximum residual: (j,k,l) = (",max_j,",",max_k,",",max_l,")"

           open(unit=out_unit,file="res_array.txt",action="write",status="replace")
           do jj = 1,iloc
              write(out_unit,*) jj,res_r(jj),res_l(jj)
           enddo
           close(out_unit)
        endif
        
        ! Evaluate residual                                                            
        rhs = rhs_qdq-rhs_q-rhs

        ! Evaluate norm of residual                                                                                                                                       
        call eval_two_norm(two_norm,rhs,jmax,kmax,lmax)

        ! Output norm of residual
        if (i.eq.1) then
           res(j,1) = s_dq
        endif
        res(j,i+1) = two_norm
     enddo

     ! Determine Convergence Rate
     slope(i) = (log(res(20,i+1))-log(res(10,i+1)))/(log(res(20,1))-log(res(10,1)))

     write(*,*) "Convergence rate for flux order ", fluxOrder, ": ", slope(i)
  enddo

  ! Write Output
  fmt = "(E18.8,E18.8,E18.8,E18.8)"
  open(unit=out_unit,file="res_conv.txt",action="write",status="replace")
  do i = 1,n_dq
     write(out_unit,fmt) res(i,1),res(i,2),res(i,3),res(i,4)
  enddo
  close(out_unit)

  deallocate(q,qdq,rhs,rhs_q,rhs_qdq,res,res_r,res_l)
end program linearizationTest

subroutine eval_two_norm(norm,rhs,jmax,kmax,lmax)
  ! 
  implicit none
  !                                                                                                                                                                      
  integer, intent(in) :: jmax,kmax,lmax               
  real*8, intent(in) :: rhs(5,jmax,kmax,lmax)                    
  real*8, intent(inout) :: norm                                  
  real*8 :: tnorm                                                
  integer :: ierr                                                
  integer :: n,j,k,l                                               
  !                                                                                                                                                                               
  tnorm=0.
  do l=4,lmax-3
     do k=4,kmax-3
        do j=4,jmax-3
           do n = 1,5
              tnorm=tnorm+rhs(n,j,k,l)**2
           enddo
        enddo
     enddo
  enddo
  !                                                                                                                                                                    
  norm=sqrt( tnorm/(jmax*kmax*lmax*n) )
                              
  return
end subroutine eval_two_norm
