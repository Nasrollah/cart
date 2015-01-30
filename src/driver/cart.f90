program cart
  use cartInterface
  implicit none
  integer :: istep,it
  real*8 :: t_start,t_end,t_iter
  !
  call cart_param_input
  call cart_mpi_init
  !
  call cart_init_data
  !
  write(*,*) "about to enter loop"
  t_iter = 0.0d0
  do istep=1,nsteps
     if (myid==0) write(6,*) '# timestep :',istep
     do it=1,nsubiter
        call cart_rhs_inviscid
        call cart_rhs_viscous
        call cart_bdf_source
        call cpu_time(t_start)
        call cart_lhs(it)
        call cpu_time(t_end)
        t_iter = t_iter + t_end - t_start
     enddo
     call cart_update_time
  enddo
  call cpu_time(t_end)
  t_iter = t_iter/real(nsubiter*nsteps)
  if (myid==0) write(*,*) "Average iteration time: ", t_iter
  !
  call cart_output
  call cart_cleanup
  !
end program cart
