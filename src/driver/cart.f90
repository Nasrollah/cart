program cart
  use cartInterface
  implicit none
  integer :: istep,it
  !
  call cart_param_input
  call cart_mpi_init
  !
  call cart_init_data
  !
  do istep=1,nsteps
     write(6,*) '# timestep :',istep
     do it=1,nsubiter
        call cart_rhs_inviscid
        !call cart_rhs_viscous
        call cart_bdf_source
        call cart_lhs(it)
     enddo
     call cart_update_time
  enddo
  !
  call cart_output
  call cart_cleanup
  !
end program cart
