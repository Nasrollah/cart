subroutine get_qcons(q,qcons,nq,nvar,jmax,kmax,lmax,j,k,l,qmult,qskip,dim1,dim2,jstride,iqp)
!
implicit none
!
integer,                              intent(in)    :: nq,nvar,jmax,kmax,lmax,j,k,l
integer,                              intent(in)    :: qmult,qskip,dim1,dim2,jstride
real*8, dimension(nq*jmax*kmax*lmax), intent(in)    :: q
real*8, dimension(nq),                intent(inout) :: qcons
integer,                              intent(out)   :: iqp
!
integer :: n,iq,iloc
!
call get_iloc(j,k,l,qmult,dim1,dim2,jstride,iloc,iqp)
!
do n=1,nvar
   qcons(n) =q(iloc)
   iloc=iloc+qskip
enddo

return
end subroutine get_qcons

subroutine get_iloc(j,k,l,qmult,dim1,dim2,jstride,iloc,iqp)
!                                                                                                                                                                                                            
implicit none
!                                                                                                                                                                                                       
integer, intent(in)    :: j,k,l
integer, intent(in)    :: qmult,dim1,dim2,jstride
integer, intent(out)   :: iloc,iqp
!                                                                                                                               
integer :: iq
!                                                                                                                                                                                                            
iq   = (l-1)*dim2+(k-1)*dim1+(j-1)*jstride
iqp  = iq + 1
iloc = iq*qmult+1
!
return
end subroutine get_iloc
