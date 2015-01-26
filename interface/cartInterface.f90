module cartInterface
  !
  use spatialCommunication, only : initSpatialComm,partitionGrid,&
	                           initDataBuffers,initiateDataSend,&
                                   initiateDataRecv,finalizeCommunication
  include 'mpif.h'
  real*8, allocatable :: q(:,:,:,:),dq(:,:,:,:),x(:,:,:,:),spec(:,:,:)
  real*8, allocatable :: qn(:,:,:,:),qnm1(:,:,:,:)
  real*8, allocatable :: rhs(:,:,:,:),qstar(:,:,:,:) ,qwork(:,:,:,:)
  real*8, allocatable :: spaceMetric(:,:,:,:),timeMetric(:,:,:,:)
  real*8, allocatable :: tscale(:,:,:)
  !
  integer :: j,k,l,i,n,bctyp
  integer :: jmax,kmax,lmax,nsteps
  integer :: irhs,ilhs
  !
  real*8  :: tf,t,h,cfl,a1,a2,a3,a4,h1,dti
  real*8  :: dx,dy,dz,second,fourth,sixth,vol
  real*8  :: fsmach,rey,totime,pr,prtr
  real*8  :: alpha=0.d0
  real*8  :: beta=0.d0
  real*8  :: gamma=1.4d0
  real*8  :: amp,freq
  real*8  :: norm,norm_two,norm_inf,norm_rhs,norm_dq
  integer :: numprocs,myid,ierr
  real*8 :: tm(3)
  integer :: nq,nvar
  real*8  :: t1,t2,dissCoef
  integer :: fluxOrder,dissOrder
  real*8  :: dt=1e6
  integer :: ivisc
  integer :: nf
  integer :: viscorder
  integer :: nsubiter
  integer :: nsave
  real*8  :: tscal
  real*8  :: xx0(3),xx1(3)
  character*16 :: timeIntegrator,icase
  character*6 :: istor
  real*8  :: t_total,s_total,tcomp_solver,tcomp_rhs,tcomp_lhs,tcomp_ts,tcomm_ts
  integer :: ninstances,nsweep
  real,parameter :: pi=acos(-1.)
  integer :: nplanes
  integer :: iplanes_send(6),iplanes_recv(6)
  integer :: iperiodic(3)
contains
  subroutine cart_set_defaults
    !
    ! default inputs
    !
    nsteps=20
    fsmach=0.5
    fluxorder=2
    dissorder=4
    disscoef=0.5
    CFL=1.0
    dt=1e6
    nq=5
    nvar=5
    rey=1e6
    ivisc=0
    viscorder=4
    timeIntegrator='ts'
    nsubiter=5
    jmax=40
    kmax=40
    lmax=40
    nsave=100
    istor='row'
    icase='taylor-green'
    ilhs=1
    irhs=0
    nsweep=1
  end subroutine cart_set_defaults
    
  subroutine cart_param_input
    namelist /inputs/ nsteps,fluxOrder,dissOrder,dissCoef,CFL,dt,ivisc,viscorder,&
         jmax,kmax,lmax,fsmach,rey,pr,prtr,nq,nsubiter,timeIntegrator,nsave,istor,icase,&
	 ilhs,irhs
    !
    call cart_set_defaults
    ! 
    open(unit=1,file='cart.input',form='formatted',err=1000)
    read(1,inputs)
    close(1)
    if (myid==0) write(6,inputs)
    !write(6,inputs)
    return
1000 continue
    write(6,*) '#################################################'
    write(6,*) 'File cart.input not found will use default values'
    write(6,*) 'File cart.input not found will use default values'
    write(6,*) '#################################################'
    return
  end subroutine cart_param_input
  !
  !!
  !> initialize mpi
  !!
  subroutine cart_mpi_init
    implicit none
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world,numprocs,ierr)
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    ninstances=numprocs
  end subroutine cart_mpi_init
  !!
  !> initialize data
  !!
  subroutine cart_init_data
    !
    implicit none
    !
    ! begin
    ! 
    t_total = 0.0
    call cpu_time(s_total)
    tcomp_solver = 0.0
    tcomp_rhs    = 0.0
    tcomp_lhs    = 0.0
    tcomp_ts     = 0.0
    tcomm_ts     = 0.0
    !
    !
    if (timeIntegrator.ne.'ts') then
       if (icase=='vortex') then
          dx=20./(jmax)
          dy=15./(kmax)
          dz=10./(lmax)
          xx0=5d0
       else
          !pi=acos(-1.)
          dx=2*pi/(jmax)
          dy=2*pi/(kmax)
          dz=2*pi/(lmax)
          xx0=pi
       endif
    else
       dx=10./(jmax)
       dy=10./(kmax)
       dz=10./(lmax)
       xx0=5d0
    endif
    !
    nvar=5
    nf=max(max(fluxorder,dissorder),viscorder)/2
    !
    call mpi_barrier(mpi_comm_world,ierr)
    call initSpatialComm(mpi_comm_world)
    call partitionGrid(jmax,kmax,lmax,dx,dy,dz,xx1,nf)
    call initdatabuffers(istor,jmax,kmax,lmax,nq,nf)
    !write(1000+myid,"(3(I10,1x),3(F10.4,1x))") jmax,kmax,lmax,xx1
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if (istor=='row') then
       allocate(q(nq,jmax,kmax,lmax),dq(nq,jmax,kmax,lmax),rhs(nq,jmax,kmax,lmax),x(3,jmax,kmax,lmax))
       allocate(qn(nq,jmax,kmax,lmax),qnm1(nq,jmax,kmax,lmax))
       allocate(spec(jmax,kmax,lmax))
       allocate(qstar(nq,jmax,kmax,lmax))
       allocate(qwork(9,jmax,kmax,lmax))
       allocate(spaceMetric(10,jmax,kmax,lmax))
       allocate(timeMetric(3,jmax,kmax,lmax))
       allocate(tscale(jmax,kmax,lmax))
    else
       allocate(q(jmax,kmax,lmax,nq),dq(jmax,kmax,lmax,nq),rhs(jmax,kmax,lmax,nq),x(jmax,kmax,lmax,3))
       allocate(qn(jmax,kmax,lmax,nq),qnm1(jmax,kmax,lmax,nq))
       allocate(spec(jmax,kmax,lmax))
       allocate(qstar(jmax,kmax,lmax,nq))
       allocate(qwork(jmax,kmax,lmax,9))
       allocate(spaceMetric(jmax,kmax,lmax,10))
       allocate(timeMetric(jmax,kmax,lmax,3))
       allocate(tscale(jmax,kmax,lmax))   
    endif
    !
    vol=dx*dy*dz
    !
    ! define the points in the grid
    !
    if (istor=='row') then
       do l=1,lmax
          do k=1,kmax
             do j=1,jmax         
                x(1,j,k,l)=(j-1)*dx-xx0(1)+xx1(1)
                x(2,j,k,l)=(k-1)*dy-xx0(2)+xx1(2)
                x(3,j,k,l)=(l-1)*dz-xx0(3)+xx1(3)
             enddo
          enddo
       enddo
    else if (istor=='column') then
       do l=1,lmax
          do k=1,kmax
             do j=1,jmax
                x(j,k,l,1)=(j-1)*dx-xx0(1)+xx1(1)
                x(j,k,l,2)=(k-1)*dy-xx0(2)+xx1(2)
                x(j,k,l,3)=(l-1)*dz-xx0(3)+xx1(3)
             enddo
          enddo
       enddo
    endif
    !
    !
    ! initialize data
    !
    freq=1.
    call init_data(q,fsmach,alpha,beta,gamma,jmax,kmax,lmax,nq,istor)
    if (timeIntegrator .ne. 'ts') then
       if (icase=='vortex') then
          t=0d0
          call vortex_exact_cart(fsmach,gamma,q,x,&
               nf+1,jmax-nf+1,1,kmax,nf+1,lmax-nf+1,jmax,kmax,lmax,nq,istor)
          call bc_case(q,nq,jmax,kmax,lmax,nf,icase,istor)
       elseif (icase=='taylor-green') then
          t=0d0
          call taylor_green(fsmach,gamma,q,x,&
            nf+1,jmax-nf,nf+1,kmax-nf,nf+1,lmax-nf,jmax,kmax,lmax,nq,istor)
          !call bc_case(q,nq,jmax,kmax,lmax,nf,icase,istor)
       endif
       
       if (dt > 1e5) then
          h=CFL*min(dz,min(dx,dy))/(1d0+fsmach)
       else
          h=dt
       endif
    else 
       t=(2*pi/freq)*(1.0*myid)/numprocs
       if (myid == 0) then
          if (bctyp == 2) then
             write(6,*) "Boundary Condition Type 2: Nonharmonic forcing"
          else
             write(6,*) "Boundary Condition Type 0: Default"
          endif
       endif
       amp=0.01
       call bc_per1(q,x,fsmach,gamma,amp,freq,t,jmax,kmax,lmax,bctyp,fluxOrder,dissOrder)       
       h =cfl*min(dz,min(dx,dy))/(1d0+fsmach)
       if (myid==0) then
          write(6,*) "Timestep = ", h
       endif
    endif
    !
    call init_metrics(spaceMetric,timeMetric,dx,dy,dz,jmax,kmax,lmax,istor)
    !
    iperiodic=(/1,1,1/)
    call updateAllFringes(q,iperiodic,nf,jmax,kmax,lmax,nq)
    !
    !call storep3di(x,q,myid,fsmach,alpha,rey,totime,jmax,kmax,lmax,nq,nf,istor)
    !
    qn=q
    qnm1=q
    !
    tscale=h
    tscal=h
    vol=dx*dy*dz
    a1=1./4/vol
    a2=8./15/vol
    a3=5./12/vol
    a4=3./4/vol
    tm=(/0.,0.,0./)
    dti=0.5*vol/h
    n=1
    !
  end subroutine cart_init_data
  !!
  !> Read in input parameters
  !!
  subroutine cart_init_param(cfl1,bctyp1,jmax1,kmax1,lmax1,irhs1,ilhs1)
    implicit none
    !
    real*8,  intent(in) :: cfl1
    integer, intent(in) :: bctyp1,jmax1,kmax1,lmax1,irhs1,ilhs1
    !
    call cart_set_defaults
    !
    cfl   = cfl1
    bctyp = bctyp1
    jmax  = jmax1
    kmax  = kmax1
    lmax  = lmax1
    irhs  = irhs1
    ilhs  = ilhs1
  end subroutine cart_init_param
  !!
  !> Evaluate RHS (Inviscid)
  !!
  subroutine cart_rhs_inviscid
    implicit none
    !
    real*8 :: t_start,t_end
    !
    call cpu_time(t_start)
    if (irhs.eq.0) then
       
       call inviscidRHSunified(nq,nvar,gamma,q,rhs,spec,tm,dx,dy,dz,jmax,kmax,lmax,&
              fluxOrder,dissOrder,dissCoef,istor)
    elseif (irhs.eq.1) then
       call inviscidRHSupwind(nq,nvar,gamma,q,rhs,spec,timeMetric,dx,dy,dz,jmax,kmax,lmax,&
            'muscl','row')
    endif
    call cpu_time(t_end)
    !
    tcomp_rhs    = tcomp_rhs    + t_end - t_start
  end subroutine cart_rhs_inviscid
  !!
  !> Evaluate RHS (viscous)
  !!
  subroutine cart_rhs_viscous
    implicit none
    !
    real*8 :: t_start,t_end
    call cpu_time(t_start)
    call viscousRHS(rey,pr,prtr,nq,nvar,gamma,qstar,qwork,rhs,dx,dy,dz,jmax,kmax,lmax,&
         min(4,viscOrder),istor)
    call cpu_time(t_end)
    !
    tcomp_rhs    = tcomp_rhs    + t_end - t_start
  end subroutine cart_rhs_viscous
  !!
  !! bdf_source term for dual time stepping
  !!
  subroutine cart_bdf_source
    implicit none
    rhs=rhs-(3d0*q-4d0*qn+qnm1)*dti !> bdf source term
  end subroutine cart_bdf_source
  !!
  !! update time
  !!
  subroutine cart_update_time
    !call bc_case(q,nq,jmax,kmax,lmax,nf,icase,istor)
    n=n+1
    t=t+h
    qnm1=qn
    qn=q      
    tscal=(2d0/3d0)*h
  end subroutine cart_update_time
  !!
  !> integrate in time using RK3
  !!
  subroutine cart_step
    !
    n=n+1
    t=t+h
    !
    call inviscidRHSunified(nq,nvar,gamma,q,rhs,spec,tm,dx,dy,dz,jmax,kmax,lmax,&
         fluxOrder,dissOrder,dissCoef,istor)
    if (ivisc==1) &
         call viscousRHS(rey,pr,prtr,nq,nvar,gamma,q,qwork,rhs,dx,dy,dz,jmax,kmax,lmax,&
         min(4,viscOrder),istor)
    call compute_norm(norm,rhs,jmax,kmax,lmax,nq,istor)
    !
    qstar=q+a2*h*(rhs)
    q=q+a1*h*(rhs)
    !
    call bc_case(q,nq,jmax,kmax,lmax,nf,icase,istor)
    call bc_case(qstar,nq,jmax,kmax,lmax,nf,icase,istor)
    !
    !  RK stage2
    !
    call inviscidRHSunified(nq,nvar,gamma,qstar,rhs,spec,tm,dx,dy,dz,jmax,kmax,lmax,&
         fluxOrder,dissOrder,dissCoef,istor)
    if (ivisc==1) &
         call viscousRHS(rey,pr,prtr,nq,nvar,gamma,qstar,qwork,rhs,dx,dy,dz,jmax,kmax,lmax,&
         min(4,viscOrder),istor)
    !
    qstar=q+a3*h*(rhs)
    call bc_case(qstar,nq,jmax,kmax,lmax,nf,icase,istor)
    !
    !  RK stage3
    !
    call inviscidRHSunified(nq,nvar,gamma,qstar,rhs,spec,tm,dx,dy,dz,jmax,kmax,lmax,&
         fluxOrder,dissOrder,dissCoef,istor)
    if (ivisc==1) &
         call viscousRHS(rey,pr,prtr,nq,nvar,gamma,qstar,qwork,rhs,dx,dy,dz,jmax,kmax,lmax,&
         min(4,viscOrder),istor)
    !
    q=q+a4*h*(rhs)
    call bc_case(q,nq,jmax,kmax,lmax,nf,icase,istor)
    write(6,*) n,t,norm
    !
  end subroutine cart_step
  !!
  !> LHS 
  !! solve the linearized problem 
  !! A\delta q = -R(q)
  !!
  subroutine cart_lhs(it)
    implicit none
    integer, intent(in) :: it
    !
    real*8 :: t_start,t_end,tsum_ts,tsum
    !> Perform Implicit Solve (LU-SGS or ADI or DADI)
    call cpu_time(t_start)
    !
    call updateAllFringes(rhs,iperiodic,nf,jmax,kmax,lmax,nq)
    !
    if (ilhs.eq.0) then
       if (istor=='row') then
          call lusgs_hyb(gamma,rey,q,rhs,tscal,timeMetric,dx,dy,dz,nq,nvar,jmax,kmax,lmax)
       else
          call lusgs_hyb_col(gamma,rey,q,rhs,tscal,timeMetric,dx,dy,dz,nq,nvar,jmax,kmax,lmax)
       endif
    elseif (ilhs.eq.1) then
       call adi(nq,nvar,gamma,q,rhs,spec,tscal,timeMetric,dx,dy,dz,jmax,kmax,lmax,&
            fluxOrder,dissOrder,dissCoef,istor)
    elseif (ilhs.eq.2) then
       call ddadi(nq,nvar,gamma,q,rhs,spec,tscal,timeMetric,dx,dy,dz,jmax,kmax,lmax,&
            fluxOrder,dissOrder,dissCoef,istor)
    elseif (ilhs.eq.3) then
       call gslr(nq,nvar,gamma,q,rhs,spec,tscal,timeMetric,dx,dy,dz,jmax,kmax,lmax,&
            fluxOrder,dissOrder,dissCoef,istor,nsweep,myid)
    elseif (ilhs.eq.4) then
       call gslr_storage(nq,nvar,gamma,q,rhs,spec,tscal,timeMetric,dx,dy,dz,jmax,kmax,lmax,&
            fluxOrder,dissOrder,dissCoef,istor,nsweep,myid)
    endif
    call cpu_time(t_end)
    !
    !> Update Solution  
    !
    dq = rhs/vol
    call updateAllFringes(dq,iperiodic,nf,jmax,kmax,lmax,nq)
    q=q+dq
    !
    !
    if (timeIntegrator.eq.'ts') then
       call compute_norm_two(norm,dq,jmax,kmax,lmax,ninstances)
       if (myid==0) write(6,*) numprocs,it,norm
    else
       call compute_norm(norm,rhs,jmax,kmax,lmax,nq,istor)
       if (myid==0) write(6,*) n,it,norm
    endif
    !
    tcomp_lhs    = tcomp_lhs    + t_end - t_start
    !
  end subroutine cart_lhs
  !!
  !> compute norm on dq and rhs                                                                      
  !!                                                                                                 
  subroutine ts_compute_norm(rhs_res_two,dq_res_two,rhs_res_inf,dq_res_inf)
    !use tsglobalvar
    implicit none
    real*8, intent(inout) :: rhs_res_two,dq_res_two,rhs_res_inf,dq_res_inf

    !> Call Spatial Residual                                                                         
    if (irhs.eq.0) then
       call inviscidRHSunified(nq,nvar,gamma,q,rhs,spec,timeMetric,dx,dy,dz,jmax,kmax,lmax,&
            fluxOrder,dissOrder,dissCoef,istor)
    elseif (irhs.eq.1) then
       call inviscidRHSupwind(nq,nvar,gamma,q,rhs,spec,timeMetric,dx,dy,dz,jmax,kmax,lmax,&
            'muscl','row')
    endif
    
    !> Add Temporal                                                                                  
    call ts_source_term(q,rhs,vol,myid,ninstances,jmax,kmax,lmax)
    
    !> Evaluate norm of rhs (2 and Inf)                                                              
    call compute_norm_two(rhs_res_two,rhs,jmax,kmax,lmax,ninstances)
    call compute_norm_inf(rhs_res_inf,rhs,jmax,kmax,lmax,ninstances)
    
    !> Evaluate norm of dq (2 and Inf)                                                                                                                                                                                                    
    call compute_norm_two(dq_res_two,dq,jmax,kmax,lmax,ninstances)
    call compute_norm_inf(dq_res_inf,dq,jmax,kmax,lmax,ninstances)


  end subroutine ts_compute_norm
  !!
  !> output data in plot3d format
  !!
  subroutine cart_output
    !
    implicit none
    !
    ! write output plot3d files
    !
   if (timeIntegrator=='ts') then
       call storep3d(x,q,myid*1000+n,fsmach,alpha,rey,totime,jmax,kmax,lmax,nq,nf,istor)
    else
       call storep3di(x,q,myid,fsmach,alpha,rey,totime,jmax,kmax,lmax,nq,nf,istor)
    endif
    !
  end subroutine cart_output
  !!
  !> clean up all the arrays
  !> and exit gracefully
  !!
  subroutine cart_cleanup
    implicit none
    !
    deallocate(q,dq,rhs,x,spaceMetric,spec,tscale)
    call mpi_finalize(ierr)
    !
  end subroutine cart_cleanup
end module cartInterface
