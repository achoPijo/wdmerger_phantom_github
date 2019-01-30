!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION: Resets centre of mass
!              
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel
!
!  $Id:  $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, part
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'resetcenterofmass'

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass
 use centreofmass, only: reset_centreofmass
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in)    :: particlemass,time
 real                            :: xpos(3),vpos(3)
 logical                         :: iexist
 character(len=200)              :: fileout

 !
 ! Call reset CoM of mass subroutine
 !

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

end subroutine do_analysis

end module
