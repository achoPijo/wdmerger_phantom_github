!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine to obtain temperature density angular velocity and keplerian velocity 
!  from a given dumpfile
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel Blanco
!
!  $Id: $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: 
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'Density,Temperature and other radial profiles and relevant data'
 logical :: firstcall = .true.

 public  :: do_analysis

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use io,           only: fatal
 use dim,          only: maxp,maxvxyzu
 use centreofmass, only: reset_centreofmass,get_centreofmass
 use part,         only: igas,iamtype,iphase,maxphase,rhoh
 use prompting,    only: prompt
! use setbinary,    only: L1_point
 use units,        only: umass,udist,utime,unit_density,unit_pressure
 use vectorutils,  only: cross_product3D
! use eos_helmholtz,only: xmass,speciesmax
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer, parameter  :: nrpoints = 1000000
 integer             :: i,j,ierr,ncountx,ncounty,ncountz,npoints
 real                :: rmax,rTmax,Tmax,dr,mtot,r,rin,rout
 real                :: rtab(nrpoints),Ttabx(nrpoints),Ttaby(nrpoints),Ttabz(nrpoints),Ttab(nrpoints),rhotab(nrpoints),rhotabx(nrpoints),rhotaby(nrpoints),rhotabz(nrpoints)
 real                :: omegatab(nrpoints),keplertab(nrpoints),rsample
 real                :: vec(3),xcom(3),vcom(3),v
 character(len=200)  :: fileout

 !
 !-- Initialization
 !
 rtab(:)          = 0.
 Ttabx(:)          = 0.
 rhotabx(:)        = 0.
 Ttabz(:)          = 0.
 rhotabz(:)        = 0.
 omegatab(:)      = 0.
 keplertab(:)     = 0.
 rsample = 0.0005
 call prompt('rsample in code units?',rsample)
 npoints = 500
 call prompt('number of points', npoints)
 mtot         = npart*particlemass
 !--------------
 
 call reset_centreofmass(npart,xyzh(:,:),vxyzu(:,:))
 !call get_centreofmass(xcom,vcom,40000,xyzh(1:40000,:),vxyzu(1:40000,:))
 !do i=1,npart
 !   xyzh(1:3,i)  = xyzh(1:3,i)  - xcom(1:3)
 !   vxyzu(1:3,i) = vxyzu(1:3,i) - vcom(1:3)
 !enddo
 !-- Loop over all particles to obtain Radius of star (POSSIBLE PROBLEM EJECTED PARTICLES)
 rmax  = 0.
 do i=1,npart
    r = sqrt(xyzh(1,i)**2+xyzh(2,i)**2+xyzh(3,i)**2)
    rmax = max(rmax,r)
 enddo

 rmax = log10(rmax)
 rin  = log10(1e-7)
 dr   = (rmax-rin)/npoints
 rout = rin + dr

 do i=1,npoints
    ncountx = 0
    ncounty = 0
    ncountz = 0
    rtab(i) = 10**rin
    keplertab(i) = sqrt(mtot/rtab(i))
    do j=1,npart
       r = sqrt(xyzh(1,j)**2+xyzh(2,j)**2)
       v = sqrt(vxyzu(1,j)**2 + vxyzu(2,j)**2)

       !if (log10(r) > rin .and. log10(r) <= rout .and. abs(xyzh(3,j)) < rsample) then !(xyzh(3,j)**2 + xyzh(2,j)**2)
       !
       !   rhotabx(i)   = rhotabx(i) + rhoh(xyzh(4,j),particlemass)
       !   Ttabx(i)     = Ttabx(i) + vxyzu(5,i)
       !   call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),vec(:))
       !   omegatab(i)  = omegatab(i) + v!vec(3)
       !   ncountx = ncountx + 1
       !endif
       if (log10(xyzh(1,j)) > rin .and. log10(xyzh(1,j)) < rout .and. (xyzh(2,j)**2 + xyzh(3,j)**2) < rsample) then

          rhotabx(i)   = rhotabx(i) + rhoh(xyzh(4,j),particlemass)
          Ttabx(i)     = Ttabx(i) + vxyzu(5,i)
          ncountx = ncountx + 1
       endif       
       if (log10(xyzh(2,j)) > rin .and. log10(xyzh(2,j)) < rout .and. (xyzh(1,j)**2 + xyzh(3,j)**2) < rsample) then

          rhotaby(i)   = rhotabz(i) + rhoh(xyzh(4,j),particlemass)
          Ttabz(i)     = Ttabz(i) + vxyzu(5,i)
          ncountz = ncountz + 1
       endif
       if (log10(xyzh(3,j)) > rin .and. log10(xyzh(3,j)) < rout .and. (xyzh(1,j)**2 + xyzh(2,j)**2) < rsample) then

          rhotabz(i)   = rhotabz(i) + rhoh(xyzh(4,j),particlemass)
          Ttabz(i)     = Ttabz(i) + vxyzu(5,i)
          ncountz = ncountz + 1
       endif
    enddo
    
    if (ncountx /= 0) then
       rhotabx(i)   = rhotabx(i)/ncountx
       Ttabx(i)     = Ttabx(i)/ncountx
       omegatab(i)  = omegatab(i)/ncountx
    endif
    if (ncountz /= 0) then
       rhotabz(i)   = rhotabz(i)/ncountz
       Ttabz(i)     = Ttabz(i)/ncountz
    endif

    rin = rout
    rout = rout + dr

 enddo

 !
 !-- Write informtaion to file
 !

 fileout = trim(dumpfile)//'_radialprofiles.dat'
 open(iunit,file=fileout,status='replace')
 write(iunit,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'r',  &
        2,'Temperature x [K]',   &
        3,'Density x [g/cm3]', &
        4,'Temperature z [K]',   &
        5,'Density z [g/cm3]',  &
        6,'omegatab',   &
        7,'keplertab'

 do j=1,npoints
      write(iunit,'(7(1pe18.10,1x))') rtab(j),Ttabx(j),rhotabx(j)*unit_density,Ttabz(j),rhotabz(j)*unit_density,omegatab(j),keplertab(j)
 enddo
 close(iunit)


 !
 !-- INCLUDE A FILE OUTPUT WITH TMAX,
 ! 


end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
