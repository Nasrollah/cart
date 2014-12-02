!=========================================================================
!> 
!! Roe's approximate Riemann Solver
!! 
!!
!! FLUX = 0.5*((FL+FR)-|A|*(QR-QL))
!
!! This implementation is from
! 
!! Vatsa, V. N., Thomas, J., L. and Wedan, B.,W.,"Navier-Stokes Computations 
!! of Prolate Spheroids at Angle of Attack," AIAA 87-2627
!!
!! Note:
!!
!! Includes grid speed terms
!! The left and right states are primitive variables
!!
!!  J. Sitaraman 06/09/09
!!
subroutine roeflxcart(gamma, f,ql,qr,idir,faceArea,faceSpeed,spec,is,ie,mdim )
!
implicit none
!
! subroutine arguments
!  
real*8, intent(in) :: gamma       !< ratio of specific heats
integer, intent(in) :: is,ie,mdim !< start-end limits and arrayd dimension
integer, intent(in) :: idir       !< coordinate direction 
real*8, intent(in)  :: faceArea   !< face area
real*8, intent(in)  :: faceSpeed  !< time metric
real*8, intent(in)  :: ql(5,mdim),qr(5,mdim) !< left and right state variables
real*8, intent(out) :: f(5,mdim)   !< reconstructed fluxes at cell faces
real*8, intent(inout) :: spec(mdim)  !< spectral radius scaled by cell face area at cell faces
!
! local variables
!
real*8 gm1
real*8 eps,rlft,ulft,vlft,wlft,plft
real*8 rlfti,rulft,rvlft,rwlft,uvwl,elft,hlft,clft
real*8 rrht,urht,vrht,wrht,prht
real*8 rrhti,rurht,rvrht,rwrht,uvwr,erht,hrht,crht
real*8 rat,rati,rav,uav,vav,wav,hav,uvw,cav
real*8 aq1,aq2,aq3,aq4,aq5,rr,r4
real*8 uu,c2,c2i,auu,aupc,aumc,uulft,uurht,upclft,upcrht
real*8 umclft,umcrht,dauu,dauus,daupc,daumc,daumcs,rcav,aquu
real*8 daupcs,c2ih,ruuav,b1,b2,b3,b4,b5,b6,b7,b8,aj
real*8 plar,eplft,eprht,fssub
integer i,i1
!
! first executable statement
!
eps = 1.0d-6
gm1=gamma-1.0
gm1=0.4
if (idir==1) then
   do  i = is,ie
      i1    = i + 1
      rlft = ql(1,i)
      ulft = ql(2,i)
      vlft = ql(3,i)
      wlft = ql(4,i)
      plft = ql(5,i)
      rlfti = 1.0/rlft
      rulft = rlft*ulft
      rvlft = rlft*vlft
      rwlft = rlft*wlft
      uvwl = 0.5*( ulft*ulft + vlft*vlft + wlft*wlft )
      elft = plft/gm1 + rlft*uvwl
      hlft = ( elft + plft )*rlfti
      clft = sqrt( gm1*( hlft - uvwl ) )

      rrht = qr(1,i1)
      urht = qr(2,i1)
      vrht = qr(3,i1)
      wrht = qr(4,i1)
      prht = qr(5,i1)
      rrhti = 1.0/rrht
      rurht = rrht*urht
      rvrht = rrht*vrht
      rwrht = rrht*wrht
      uvwr = 0.5*( urht*urht + vrht*vrht + wrht*wrht )
      erht = prht/gm1 + rrht*uvwr
      hrht = ( erht + prht )*rrhti
      crht = sqrt( gm1*( hrht - uvwr ) )

      rat  = sqrt( rrht*rlfti )
      rati = 1.0/( rat + 1. )
      rav  =   rat*rlft
      uav  = ( rat*urht + ulft )*rati
      vav  = ( rat*vrht + vlft )*rati
      wav  = ( rat*wrht + wlft )*rati
      hav  = ( rat*hrht + hlft )*rati
      uvw  = 0.5*( uav*uav + vav*vav + wav*wav )
      cav  = sqrt( gm1*( hav - uvw ) )

      aq1  = rrht - rlft
      aq2  = urht - ulft
      aq3  = vrht - vlft
      aq4  = wrht - wlft
      aq5  = prht - plft

      rr=faceArea
      r4=faceSpeed/rr

      uu  = uav !r1*uav + r2*vav + r3*wav + r4
      c2  = cav*cav
      c2i = 1.0/c2

      auu   = abs( uu    )
      aupc  = abs( uu+cav )
      aumc  = abs( uu-cav )

      uulft = ulft + r4 !r1*ulft + r2*vlft + r3*wlft + r4
      uurht = urht + r4 !r1*urht + r2*vrht + r3*wrht + r4
      upclft= uulft + clft
      upcrht= uurht + crht
      umclft= uulft - clft
      umcrht= uurht - crht

      dauu = 4.*(uurht-uulft)+eps
      dauus = dmax1(dauu,0.0)
      if (auu.le.0.5*dauus) auu = auu**2/dauu+0.25*dauu
      daupc = 4.*(upcrht-upclft)+eps
      daupcs = dmax1(daupc,0.0)
      if (aupc.le.0.5*daupcs) aupc = aupc**2/daupc+0.25*daupc
      daumc = 4.*(umcrht-umclft)+eps
      daumcs = dmax1(daumc,0.0)
      if (aumc.le.0.5*daumcs) aumc = aumc**2/daumc+0.25*daumc

      spec(i)=max(auu,aupc)
      spec(i)=max(spec(i),aumc)
      spec(i)=spec(i)*rr

      rcav = rav*cav
      aquu = uurht - uulft
      c2ih = 0.5*c2i
      ruuav= auu*rav
      b1   = auu*( aq1 - c2i*aq5 )
      b2   = c2ih*aupc*( aq5 + rcav*aquu )
      b3   = c2ih*aumc*( aq5 - rcav*aquu )
      b4   = b1 + b2 + b3
      b5   = cav*( b2 - b3 )
      b6   = ruuav*( aq2 - aquu )
      b7   = ruuav*( aq3 ) !- r2*aquu )
      b8   = ruuav*( aq4 ) !- r3*aquu )

      aq1 = b4
      aq2 = uav*b4  + b5 + b6
      aq3 = vav*b4  + b7 ! + r2*b5
      aq4 = wav*b4  + b8 ! + r3*b5
      aq5 = hav*b4 + ( uu-r4 )*b5 + uav*b6 + vav*b7 + wav*b8 &
           - c2*b1/gm1
      aj    = 0.5*rr
      plar  = plft + prht
      eplft = elft + plft
      eprht = erht + prht
      fssub = rr*r4
      fssub = 0.0

      f(1,i) = aj*(  rlft*uulft +  rrht*uurht           - aq1 )
      f(2,i) = aj*( rulft*uulft + rurht*uurht + plar - aq2 )
      f(3,i) = aj*( rvlft*uulft + rvrht*uurht  - aq3 )  !+ r2*plar
      f(4,i) = aj*( rwlft*uulft + rwrht*uurht  - aq4 )  !+ r3*plar
      f(5,i) = aj*( eplft*uulft + eprht*uurht - r4*plar - aq5 )
   enddo
elseif (idir==2) then
      do  i = is,ie
      i1    = i + 1
      rlft = ql(1,i)
      ulft = ql(2,i)
      vlft = ql(3,i)
      wlft = ql(4,i)
      plft = ql(5,i)
      rlfti = 1.0/rlft
      rulft = rlft*ulft
      rvlft = rlft*vlft
      rwlft = rlft*wlft
      uvwl = 0.5*( ulft*ulft + vlft*vlft + wlft*wlft )
      elft = plft/gm1 + rlft*uvwl
      hlft = ( elft + plft )*rlfti
      clft = sqrt( gm1*( hlft - uvwl ) )

      rrht = qr(1,i1)
      urht = qr(2,i1)
      vrht = qr(3,i1)
      wrht = qr(4,i1)
      prht = qr(5,i1)
      rrhti = 1.0/rrht
      rurht = rrht*urht
      rvrht = rrht*vrht
      rwrht = rrht*wrht
      uvwr = 0.5*( urht*urht + vrht*vrht + wrht*wrht )
      erht = prht/gm1 + rrht*uvwr
      hrht = ( erht + prht )*rrhti
      crht = sqrt( gm1*( hrht - uvwr ) )

      rat  = sqrt( rrht*rlfti )
      rati = 1.0/( rat + 1. )
      rav  =   rat*rlft
      uav  = ( rat*urht + ulft )*rati
      vav  = ( rat*vrht + vlft )*rati
      wav  = ( rat*wrht + wlft )*rati
      hav  = ( rat*hrht + hlft )*rati
      uvw  = 0.5*( uav*uav + vav*vav + wav*wav )
      cav  = sqrt( gm1*( hav - uvw ) )

      aq1  = rrht - rlft
      aq2  = urht - ulft
      aq3  = vrht - vlft
      aq4  = wrht - wlft
      aq5  = prht - plft

      rr=faceArea
      r4=faceSpeed/rr

      uu  = vav+r4 !r1*uav + r2*vav + r3*wav + r4
      c2  = cav*cav
      c2i = 1.0/c2

      auu   = abs( uu    )
      aupc  = abs( uu+cav )
      aumc  = abs( uu-cav )

      uulft = vlft+r4 !r1*ulft + r2*vlft + r3*wlft + r4
      uurht = vrht+r4 !r1*urht + r2*vrht + r3*wrht + r4
      upclft= uulft + clft
      upcrht= uurht + crht
      umclft= uulft - clft
      umcrht= uurht - crht

      dauu = 4.*(uurht-uulft)+eps
      dauus = dmax1(dauu,0.0)
      if (auu.le.0.5*dauus) auu = auu**2/dauu+0.25*dauu
      daupc = 4.*(upcrht-upclft)+eps
      daupcs = dmax1(daupc,0.0)
      if (aupc.le.0.5*daupcs) aupc = aupc**2/daupc+0.25*daupc
      daumc = 4.*(umcrht-umclft)+eps
      daumcs = dmax1(daumc,0.0)
      if (aumc.le.0.5*daumcs) aumc = aumc**2/daumc+0.25*daumc

      spec(i)=max(auu,aupc)
      spec(i)=max(spec(i),aumc)
      spec(i)=spec(i)*rr

      rcav = rav*cav
      aquu = uurht - uulft
      c2ih = 0.5*c2i
      ruuav= auu*rav
      b1   = auu*( aq1 - c2i*aq5 )
      b2   = c2ih*aupc*( aq5 + rcav*aquu )
      b3   = c2ih*aumc*( aq5 - rcav*aquu )
      b4   = b1 + b2 + b3
      b5   = cav*( b2 - b3 )
      b6   = ruuav*( aq2 ) !- r1*aquu )
      b7   = ruuav*( aq3 - aquu )
      b8   = ruuav*( aq4 ) !- r3*aquu )

      aq1 = b4
      aq2 = uav*b4 +  b6     ! + r1*b5
      aq3 = vav*b4  + b7 + b5      
      aq4 = wav*b4  + b8 ! + r3*b5
      aq5 = hav*b4 + ( uu-r4 )*b5 + uav*b6 + vav*b7 + wav*b8 &
           - c2*b1/gm1
      aj    = 0.5*rr
      plar  = plft + prht
      eplft = elft + plft
      eprht = erht + prht
      fssub = rr*r4
      fssub = 0.0
      f(1,i) = aj*(  rlft*uulft +  rrht*uurht           - aq1 )
      f(2,i) = aj*( rulft*uulft + rurht*uurht  - aq2 )     ! +r1*plar
      f(3,i) = aj*( rvlft*uulft + rvrht*uurht  +plar- aq3 )  
      f(4,i) = aj*( rwlft*uulft + rwrht*uurht  - aq4 )  !+ r3*plar
      f(5,i) = aj*( eplft*uulft + eprht*uurht - r4*plar - aq5 )
   enddo
else
   do  i = is,ie
      i1    = i + 1
      rlft = ql(1,i)
      ulft = ql(2,i)
      vlft = ql(3,i)
      wlft = ql(4,i)
      plft = ql(5,i)
      rlfti = 1.0/rlft
      rulft = rlft*ulft
      rvlft = rlft*vlft
      rwlft = rlft*wlft
      uvwl = 0.5*( ulft*ulft + vlft*vlft + wlft*wlft )
      elft = plft/gm1 + rlft*uvwl
      hlft = ( elft + plft )*rlfti
      clft = sqrt( gm1*( hlft - uvwl ) )

      rrht = qr(1,i1)
      urht = qr(2,i1)
      vrht = qr(3,i1)
      wrht = qr(4,i1)
      prht = qr(5,i1)
      rrhti = 1.0/rrht
      rurht = rrht*urht
      rvrht = rrht*vrht
      rwrht = rrht*wrht
      uvwr = 0.5*( urht*urht + vrht*vrht + wrht*wrht )
      erht = prht/gm1 + rrht*uvwr
      hrht = ( erht + prht )*rrhti
      crht = sqrt( gm1*( hrht - uvwr ) )

      rat  = sqrt( rrht*rlfti )
      rati = 1.0/( rat + 1. )
      rav  =   rat*rlft
      uav  = ( rat*urht + ulft )*rati
      vav  = ( rat*vrht + vlft )*rati
      wav  = ( rat*wrht + wlft )*rati
      hav  = ( rat*hrht + hlft )*rati
      uvw  = 0.5*( uav*uav + vav*vav + wav*wav )
      cav  = sqrt( gm1*( hav - uvw ) )

      aq1  = rrht - rlft
      aq2  = urht - ulft
      aq3  = vrht - vlft
      aq4  = wrht - wlft
      aq5  = prht - plft

      rr=faceArea
      r4=faceSpeed/rr

      uu  = wav+r4 !r1*uav + r2*vav + r3*wav + r4
      c2  = cav*cav
      c2i = 1.0/c2

      auu   = abs( uu    )
      aupc  = abs( uu+cav )
      aumc  = abs( uu-cav )

      uulft = wlft+r4 !r1*ulft + r2*vlft + r3*wlft + r4
      uurht = wrht+r4 !r1*urht + r2*vrht + r3*wrht + r4
      upclft= uulft + clft
      upcrht= uurht + crht
      umclft= uulft - clft
      umcrht= uurht - crht

      dauu = 4.*(uurht-uulft)+eps
      dauus = dmax1(dauu,0.0)
      if (auu.le.0.5*dauus) auu = auu**2/dauu+0.25*dauu
      daupc = 4.*(upcrht-upclft)+eps
      daupcs = dmax1(daupc,0.0)
      if (aupc.le.0.5*daupcs) aupc = aupc**2/daupc+0.25*daupc
      daumc = 4.*(umcrht-umclft)+eps
      daumcs = dmax1(daumc,0.0)
      if (aumc.le.0.5*daumcs) aumc = aumc**2/daumc+0.25*daumc

      spec(i)=max(auu,aupc)
      spec(i)=max(spec(i),aumc)
      spec(i)=spec(i)*rr

      rcav = rav*cav
      aquu = uurht - uulft
      c2ih = 0.5*c2i
      ruuav= auu*rav
      b1   = auu*( aq1 - c2i*aq5 )
      b2   = c2ih*aupc*( aq5 + rcav*aquu )
      b3   = c2ih*aumc*( aq5 - rcav*aquu )
      b4   = b1 + b2 + b3
      b5   = cav*( b2 - b3 )
      b6   = ruuav*( aq2 ) !- r1*aquu )
      b7   = ruuav*( aq3 ) !- r2*aquu
      b8   = ruuav*( aq4 - aquu) 

      aq1 = b4
      aq2 = uav*b4 +  b6     ! + r1*b5
      aq3 = vav*b4  + b7     ! + r2*b5
      aq4 = wav*b4  + b8 + b5
      aq5 = hav*b4 + ( uu-r4 )*b5 + uav*b6 + vav*b7 + wav*b8 &
           - c2*b1/gm1
      aj    = 0.5*rr
      plar  = plft + prht
      eplft = elft + plft
      eprht = erht + prht
      fssub = rr*r4
      fssub = 0.0
      f(1,i) = aj*(  rlft*uulft +  rrht*uurht           - aq1 )
      f(2,i) = aj*( rulft*uulft + rurht*uurht  - aq2 )     ! +r1*plar
      f(3,i) = aj*( rvlft*uulft + rvrht*uurht  - aq3 )     ! +r2*plar  
      f(4,i) = aj*( rwlft*uulft + rwrht*uurht  +plar- aq4 ) 
      f(5,i) = aj*( eplft*uulft + eprht*uurht - r4*plar - aq5 )
   enddo
endif
!
return
end subroutine roeflxcart

