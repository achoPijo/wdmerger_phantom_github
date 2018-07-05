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
 use setbinary,    only: L1_point
 use units,        only: umass,udist,utime,unit_density,unit_pressure
 use vectorutils,  only: cross_product3D
 use eos_helmholtz,only: xmass,speciesmax
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer                      :: i,j,ierr,ncount
 integer                      :: nrpoints = 10000
 real                         :: rmax,rTmax,Tmax,dr,mtot,r
 real                         :: rtab(nrpoints),Ttab(nrpoints),rhotab(nrpoints)
 real                         :: omegatab(nrpoints),keplertab(nrpoints)
 real                         :: compAverage(speciesmax-1,nrpoints)
 real                         :: vec(3)
 character(len=200)           :: fileout

 !
 !-- Initialization
 !
 rtab(:)          = 0.
 Ttab(:)          = 0.
 rhotab(:)        = 0.
 omegatab(:)      = 0.
 keplertab(:)     = 0.
 compAverage(:,:) = 0.

 mtot         = npart*particlemass
 !--------------

 call reset_centreofmass(npart,xyzh(:,:),vxyzu(:,:))
 
 !-- Loop over all particles to obtain Radius of star (POSSIBLE PROBLEM EJECTED PARTICLES)
 rmax  = 0.
 rTmax = 0.
 Tmax  = 0.
 do i=1,npart
    r = sqrt(xyzh(1,i)**2+xyzh(2,i)**2+xyzh(3,i)**2)
    rmax = max(rmax,r)
    Tmax = max(Tmax,vxyzu(5,i))
    if (Tmax == vxyzu(5,i)) then 
       rTmax = r
    endif
 enddo

 dr = rmax/nrpoints
 rtab(1)=dr
 do i=2,nrpoints
    rtab(i)      = rtab(i-r)+dr
    keplertab(i) = sqrt(mtot/rtab(i))
 enddo
 T = 0
 do i=1,nrpoints
    ncount = 0

    do j=1,npart
       r = sqrt(xyzh(1,j)**2+xyzh(2,j)**2+xyzh(3,j)**2)

       if (r < rtab(i) + dr/2 .and. r > rtab(i) - dr/2) then

          rhotab(i)    = rhotab(i) + rhoh(xyzh(4,j),particlemass)
          Ttab(i)      = Ttab(i) + vxyzu(5,i)
          call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),vec(:))
          omegatab(i)  = omegatab(i) + vec(3)
          compAverage(:,i) = compAverage(:,i) + xmass(1:15,j)

          ncount = ncount + 1
       endif
    enddo
    
    if (ncount =/ 0) then
       rhotab(i)   = rhotab(i)/ncount
       Ttab(i)     = Ttab(i)/ncount
       omegatab(i) = omegatab(i)/ncount
       compAverage(:,i) =compAverage(:,i)/ncount
    endif
 enddo

 !
 !-- Write informtaion to file
 !

 fileout = trim(dumpfile)//'_radialprofiles.dat'
 open(iunit,file=fileout,status='replace')
 write(iunit,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'r',&
        2,'Temperature',   &
        3,'Density', &
        4,'Omega',   &
        5,'Kepler V',&
        6,'Composition 15 columns'
 do j=1,nrpoints
      write(iunit,'(20(1pe18.10,1x))') rtab(j),Ttab(j),rhotab(j),omegatab(j),keplertab(j),compAverage(:,j)
 enddo
 close(iunit)


 !
 !-- INCLUDE A FILE OUTPUT WITH TMAX,
 ! 


end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
