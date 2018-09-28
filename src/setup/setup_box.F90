!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  This module sets up a box with a given lattice structure, pressure  
!  and density
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel Blanco
!
!  $Id: fed1eb61c5cdf5aa93e91dcf6c500fc4baefa0aa $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: physcon, units, kernel, parts
!+
!--------------------------------------------------------------------------
module setup
 use part,              only: maxvxyzu
 use units,             only: utime
 implicit none

 real(kind=8)       :: udist,umass
 character(len=20)  :: dist_unit,mass_unit
 logical            :: use_prompt, iexist, binary
 integer            :: np

 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  Setup routine for White Dwarfs
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use boundary,      only: set_boundary
 use centreofmass,  only: reset_centreofmass 
 use eos,           only: ieos, equationofstate, init_eos,finish_eos
 use dim,           only: maxp
 use io,            only: master
 use kernel,        only: hfact_default
 use options,       only: iexternalforce,nfulldump,damp,alphau
 use part,          only: igas,rhoh
 use physcon,       only: solarm,solarr,pi,planckh,mass_electron_cgs,mass_proton_cgs
 use prompting,     only: prompt
 use timestep,      only: tmax, dtmax
 use units,         only: set_units, select_unit, unit_pressure, unit_density
 use io,            only: warning

 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real,              parameter     :: mue=2.0d0
 real                             :: inMatrix(128000,10)
 real                             :: xboundmin,xboundmax,yboundmin,yboundmax,zboundmin,zboundmax,xbound,ybound,zbound
 character(len=120)               :: setupfile,inname
 logical                          :: write_setup
 integer                          :: i, ierr
 real                             :: rhozero,presszero,mtotal


 !
 ! Set default units
 !
 call set_units(dist='cm',mass='g',time='s')
 
 npartoftype(:) = 0
 massoftype(:)  = 0.d0
 !
 !--Initializations
 !
 npart = 128000
 npartoftype(igas) = npart
 rhozero = 1.
 presszero = 900.

 !
 ! Set equation of state
 !
 ieos = 2
 !
 ! set gamma
 !
 gamma = 5./3.
 polyk = 0.

 time  = 0.
 hfact = hfact_default

 ! set_boundary
 xboundmin=0.
 xboundmax=1.
 xbound   =xboundmax-xboundmin
 yboundmin=0.
 yboundmax=1.
 ybound   =yboundmax-yboundmin
 zboundmin=0.
 zboundmax=1.
 zbound   =zboundmax-zboundmin
 call set_boundary(xboundmin,xboundmax,yboundmin,yboundmax,zboundmin,zboundmax)

 !call set_unifids(npart,hfact,mstar,xyzh)
 mtotal = rhozero*xbound*ybound*zbound
 massoftype(igas) = mtotal/npart

 xyzh(4,:)    = hfact*((massoftype(igas)/rhozero)**(1./3.))
 vxyzu(1:3,:) = 0.
 vxyzu(4,i)   = rhozero*(gamma-1)/presszero

 ! RETRIEVE lattice

 open(10,file='Grid3DBCC', form='formatted')
 do i=1,npart
    read(10,*) inMatrix(i,:)
 enddo
 do i=1,npart
    xyzh(1,i)=inMatrix(i,1)
    xyzh(2,i)=inMatrix(i,2)
    xyzh(3,i)=inMatrix(i,3)
 enddo
 close(10)

 !
 !--Set new runtime parameters
 tmax           =    2.                 !2 time units
 dtmax          =    1./utime           !1 second
 nfulldump      =    1
 iexternalforce =    0
 alphau         =    0

end subroutine setpart

end module setup

