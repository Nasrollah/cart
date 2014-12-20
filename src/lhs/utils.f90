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

subroutine store_fj(q,mdim,nq,nvar,k,l,dim0,dim1,dim2,qmult,qskip,idir,faceSpeed,gm1,efac,&
     dfdq,pressure,spec,sigma,flux_order,diss_order,jstride,jmax,kmax,lmax)
!
implicit none
!
real*8, dimension(nq*jmax*kmax*lmax), intent(in) :: q
integer, dimension(3), intent(in) :: dim0
integer,                           intent(in) :: nq,nvar,mdim,k,l
integer, intent(in) :: dim1,dim2,qmult,qskip,idir
real*8,                            intent(in) :: faceSpeed,gm1,efac
real*8, dimension(nq,nq,mdim), intent(inout) :: dfdq
real*8, dimension(jmax*kmax*lmax), intent(in) :: pressure,spec
real*8, dimension(mdim),           intent(inout) :: sigma
integer, intent(in) :: flux_order,diss_order,jstride,jmax,kmax,lmax
!
integer :: j,iq,iqp,iloc,n
real*8 :: vel,rvel,qcons(nq)
!
iq=(l-1)*dim2+(k-1)*dim1
!
do j=1,dim0(idir)
   iloc=iq*qmult+1
   iqp=iq+1
   do n=1,nvar
      qcons(n)=q(iloc)
      iloc=iloc+qskip
   enddo
   !                                                                                                
   vel = qcons(idir+1)/qcons(1)
   rvel= vel-faceSpeed
   
   call fluxJacobian(nq,gm1,qcons,vel,rvel,pressure(iqp),&
        flux_order,diss_order,efac,idir,dfdq(:,:,j))
   !                                                                                                
   sigma(j)  = spec(iqp)+abs(vel)
   iq        = iq+jstride
   !                        
enddo
!
return
end subroutine store_fj
