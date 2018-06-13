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
!  Analysis routine check whether ddensity profile corresponds to that of 
!  a polytrope
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
 character(len=20), parameter, public :: analysistype = 'Burning'
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
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer                      :: i,iloc
 real                         :: rTmax,Tmax,r
 logical                      :: iexist
 character(len=120)           :: fileprefix,fileout





 !
 !-- Initialization

 !--------------
 !--------------
print *, vxyzu(1,49999)
print *, vxyzu(2,49999)
print *, vxyzu(3,49999)
print *, vxyzu(4,49999)
print *, vxyzu(5,49999)

print *, "------------"

print *, vxyzu(1,50000)
print *, vxyzu(2,50000)
print *, vxyzu(3,50000)
print *, vxyzu(4,50000)
print *, vxyzu(5,50000)


!print *, xyzh(1,49999)
!print *, xyzh(2,49999)
!print *, xyzh(3,49999)
!print *, xyzh(4,49999)
!
!print *, "------------"
!
!print *, xyzh(1,50000)
!print *, xyzh(2,50000)
!print *, xyzh(3,50000)
!print *, xyzh(4,50000)



end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
