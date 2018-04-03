!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!   Take two relaxed stars from phanom_moddump_binary_system.f90
!   and place them in orbit at a given distance and with given 
!   initial conditions (irrotational - co-rotating)
!   Author: Jose Miguel Blanco (supervisor: Pablo Lorén-Aguilar)
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel Blanco
!
!  $Id: -
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, extern_gwinspiral, externalforces, io,
!    options, part, physcon, prompting, timestep, units
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

 logical, parameter :: use_defaults = .false.  ! if .true. will automatically use default values
                                               ! if .false., will ask user to prompt new values
 !--The default values
 real,    private   :: separation   = 50.
 logical, private   :: useirrinit   = .false.

 public             :: modify_dump
 private

contains
!-----------------------------------------------------------------------

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use eos,            only: relflag
 use io,             only: iprint,fatal
 use prompting,      only: prompt
 use options,        only: iexternalforce,nfulldump,damp
 use part,           only: igas
 use units,          only: unit_velocity,utime
 use physcon,        only: c,pi
 use timestep,       only: tmax, dtmax
 use centreofmass,   only: get_centreofmass,reset_centreofmass                 
 use externalforces, only: iext_gwinspiral
 use extern_gwinspiral, only: Nstar
 integer, intent(inout)    :: npart
 integer, intent(inout)    :: npartoftype(:)
 real,    intent(inout)    :: massoftype(:)
 real,    intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                   :: i,ierr
 real                      :: com(3),com_star1(3),com_star2(3),vcom(3)
 real                      :: rad1,rad2,mstar1,mstar2,mtotal, omega, omega2    !CHNGCODE added omega
 real                      :: xcm1,xcm2,ycm1,ycm2,r1,r2
 real                      :: c_code,tmax0 
 character(len=120)        :: setupfile

 !--Reset centre of mass location 
 call reset_centreofmass(npart,xyzh(:,:),vxyzu(:,:))

 !
 !--Restart velocity                     
 vxyzu(1:3,:) = 0.0                             ! reset velocity
 !
 !--Set new runtime parameters
 tmax           =    2.                 !2 time units
 dtmax          =    1./utime           !1 second
 damp           =    0.
 nfulldump      =    1
 iexternalforce =    0
 relflag        =   .false.
 !
 return
end subroutine modify_dump
!-----------------------------------------------------------------------
end module moddump
