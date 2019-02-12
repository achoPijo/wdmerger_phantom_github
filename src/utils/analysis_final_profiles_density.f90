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

 integer, parameter  :: nrpoints = 2000
 integer             :: i,j,ierr,ncountx,ncountz
 real                :: rmax,rTmax,Tmax,dr,mtot,r
 real                :: rtab(nrpoints),Ttabx(nrpoints),Ttabz(nrpoints),Ttab(nrpoints),rhotab(nrpoints),rhotabx(nrpoints),rhotabz(nrpoints)
 real                :: omegatab(nrpoints),keplertab(nrpoints),ncountxtab(nrpoints),ncountztab(nrpoints),macumtab(nrpoints),rsample,maxdens,macum
 real                :: halfgravacctab(nrpoints),centrifugalacctab(nrpoints)
 real                :: vec(3),xcom(3),vcom(3),xdens(3)
 character(len=200)  :: fileout

 !
 !-- Initialization
 !
 rtab(:)          = 0.
 Ttabx(:)         = 0.
 rhotabx(:)       = 0.
 Ttabz(:)         = 0.
 rhotabz(:)       = 0.
 omegatab(:)      = 0.
 keplertab(:)     = 0.
 ncountxtab(:)    = 0.
 ncountztab(:)    = 0.
 macumtab(:)      = 0.
 macum            = 0.
 rsample = 0.001
 !nrpoints = 2000
 call prompt('rsample in code units?', rsample)
 !call prompt('number of points', nrpoints)
 mtot         = npart*particlemass
 !--------------
 
 !call reset_centreofmass(npart,xyzh(:,:),vxyzu(:,:))
 call get_centreofmass(xcom,vcom,npart,xyzh(:,:),vxyzu(:,:))

 !Find point of maximum density
 maxdens = 0.0
 do i=1,npart
    if (rhoh(xyzh(4,i),particlemass) > maxdens ) then
       xdens(:) = xyzh(1:3,i)
       maxdens  = rhoh(xyzh(4,i),particlemass)
    endif
 enddo
 print *, npart
 print *, xdens

 do i=1,npart
    xyzh(1:3,i)  = xyzh(1:3,i)  - xdens(1:3)!xcom(1:3)
    vxyzu(1:3,i) = vxyzu(1:3,i) - vcom(1:3)
 enddo
 !-- Loop over all particles to obtain Radius of star (POSSIBLE PROBLEM EJECTED PARTICLES)
 rmax  = 0.
 do i=1,npart
    r = sqrt(xyzh(1,i)**2+xyzh(2,i)**2+xyzh(3,i)**2)
    rmax = max(rmax,r)
 enddo

 rmax = 0.1
 dr = rmax/nrpoints
 rtab(1)=dr
 do i=2,nrpoints
    rtab(i)  = rtab(i-r)+dr
    macum    = 0.
    do j = 1,npart
       r = sqrt(xyzh(1,j)**2+xyzh(2,j)**2+xyzh(3,j)**2)
       if (r < rtab(i)) then
          macum = macum + particlemass
       endif
    enddo
    keplertab(i) = sqrt(macum/rtab(i))/rtab(i) !!!!KEPLEEEER SI WITH ACCUMULATED MASS!!!!!!!
    macumtab(i)  = macum
    halfgravacctab(i) = (1/2)*macumtab(i)/(rtab(i)**2) !helf gravitational acceleration
 enddo

 do i=1,nrpoints
    ncountx = 0
    ncountz = 0

    do j=1,npart
       r = sqrt(xyzh(1,j)**2+xyzh(2,j)**2+xyzh(3,j)**2)

       !if (xyzh(1,j) < rtab(i) + dr/2 .and. xyzh(1,j) > rtab(i) - dr/2 .and. sqrt(xyzh(3,j)**2 + xyzh(2,j)**2) < rsample) then
       if (sqrt(xyzh(1,j)**2 + xyzh(2,j)**2) < rtab(i) + dr/2 .and. sqrt(xyzh(1,j)**2 + xyzh(2,j)**2) > rtab(i) - dr/2 .and. abs(xyzh(3,j)) < rsample) then

          rhotabx(i)   = rhotabx(i) + rhoh(xyzh(4,j),particlemass)
          Ttabx(i)     = Ttabx(i) + vxyzu(5,j)
          call cross_product3D(xyzh(1:3,j),vxyzu(1:3,j),vec(:))
          omegatab(i)  = omegatab(i) + vec(3)/(r**2)
          ncountx = ncountx + 1
       endif
       if (xyzh(3,j) < rtab(i) + dr/2 .and. xyzh(3,j) > rtab(i) - dr/2 .and. sqrt(xyzh(1,j)**2 + xyzh(2,j)**2) < rsample) then

          rhotabz(i)   = rhotabz(i) + rhoh(xyzh(4,j),particlemass)
          Ttabz(i)     = Ttabz(i) + vxyzu(5,j)
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

    ncountxtab(i) = ncountx
    ncountztab(i) = ncountz
    centrifugalacctab(i)  = (omegatab(i)**2)/rtab(i)

 enddo

 !
 !-- Write informtaion to file
 !

 fileout = trim(dumpfile)//'_radialprofiles.dat'
 open(iunit,file=fileout,status='replace')
 write(iunit,"('#',12(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'r',  &
        2,'Temperature x [K]',   &
        3,'Density x [g/cm3]', &
        4,'Temperature z [K]',   &
        5,'Density z [g/cm3]',  &
        6,'omegatab',   &
        7,'keplertab',  &
        8,'macumtab',   &
        9,'ncountx',    &
       10,'ncountz',    &
       11,'halfgravacctab',    &
       12,'ncountz',    &

 do j=1,nrpoints
      write(iunit,'(12(1pe18.10,1x))') rtab(j),Ttabx(j),rhotabx(j)*unit_density,Ttabz(j),rhotabz(j)*unit_density,omegatab(j)/utime,keplertab(j)/utime,macumtab(j),ncountxtab(j),ncountztab(j),halfgravacctab(j),centrifugalacctab(j)
 enddo
 close(iunit)


 !
 !-- INCLUDE A FILE OUTPUT WITH TMAX omegamax locations
 ! 


end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
