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
 real,    private   :: separation           = 0.05
 logical, private   :: use_irr_init         = .false.
 logical, private   :: use_corotating_frame = .true.
 logical, private   :: binary               = .true.
 logical, private   :: use_relocation       = .false.

 public             :: modify_dump
 private

contains
!-----------------------------------------------------------------------

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use eos,                     only: relflag
 use io,                      only: iprint,fatal
 use prompting,               only: prompt
 use options,                 only: iexternalforce,nfulldump,damp
 use part,                    only: igas
 use units,                   only: unit_velocity,utime
 use physcon,                 only: c,pi
 use timestep,                only: tmax, dtmax
 use centreofmass,            only: get_centreofmass,reset_centreofmass                 
 use externalforces,          only: iext_corotate
 use extern_corotate,         only: omega_corotate,dynfac
 use extern_gwinspiral,       only: Nstar
 use nuc_reactions,           only: nuc_burn
 integer, intent(inout)    :: npart
 integer, intent(inout)    :: npartoftype(:)
 real,    intent(inout)    :: massoftype(:)
 real,    intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                   :: i,ierr
 real                      :: xcom1(3),xcom2(3),vcom1(3),vcom2(3)
 real                      :: rad1,rad2,mstar1,mstar2,mtotal, omega, omega2    !CHNGCODE added omega
 real                      :: xcm1,xcm2,ycm1,ycm2,r1,r2
 character(len=120)        :: setupfile

 !--Reset centre of mass location 
 call reset_centreofmass(npart,xyzh(:,:),vxyzu(:,:))

 return
end subroutine modify_dump
!-----------------------------------------------------------------------

subroutine read_Nstar_from_setup(filename,ierr)
 use extern_gwinspiral, only: Nstar
 use infile_utils,      only: open_db_from_file,inopts,close_db,read_inopt
 use io,                only: error
 use part,              only: maxvxyzu
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)
 real                          :: dummyr
 integer                       :: dummyint
 logical                       :: binary
 !
 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(*, '(1x,2a)') 'Setup_wd: Reading setup Nstar options from ',trim(filename)
 !
 nerr = 0
 call read_inopt(dummyr,'dist_unit',db,ierr)
 nerr= nerr + ierr
 call read_inopt(dummyr,'mass_unit',db,ierr)
 nerr= nerr + ierr
 call read_inopt(dummyr,'utime',db,ierr)
 nerr= nerr + ierr
 if (maxvxyzu == 5) then
     call read_inopt(dummyr,'Tin',db,ierr)
     nerr= nerr + ierr
 endif
 call read_inopt(binary,'binary',db,ierr)
 nerr= nerr + ierr
 call read_inopt(dummyint,'ntotal',db,ierr)
 nerr= nerr + ierr
 call read_inopt(dummyr,'mstar1',db,ierr)
 nerr= nerr + ierr
 call read_inopt(dummyr,'mstar2',db,ierr)
 nerr= nerr + ierr
 call read_inopt(Nstar(1),'Nstar1',db,ierr)
 nerr= nerr + ierr
 call read_inopt(Nstar(2),'Nstar2',db,ierr)
 nerr= nerr + ierr

if (nerr > 0) then
    call error ('setup_wd','there were some errors in reading the setup file')
    print "(1x,a,i2,a)",'Setup_wd: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_Nstar_from_setup

!------------------------------------
! Compute estimate of separation for 
! start of Roche Lobe overflow
! Eggleton (1983) ApJ 268, 368-369
!------------------------------------
real function separation_from_Rochelobe(m1,m2,R1,R2)
 real, intent(in) :: m1,m2,R1,R2
 real :: q,q13,q23,RL

 if (m1 <= m2) then
    q = m1/m2
    RL = R1
 else
    q = m2/m1
    RL = R2
 endif
 q13 = q**(1./3.)
 q23 = q13*q13 
 separation_from_Rochelobe = RL / (0.49*q23/(0.6*q23 + log(1. + q13)))
end function separation_from_Rochelobe

end module moddump
