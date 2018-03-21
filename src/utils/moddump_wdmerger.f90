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
 use units,          only: unit_velocity
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

 ! ADD HERE  A PROMPT TO ASK FOR THE SETUP FILE TO CHECK THE VALUES OF NSTAR1 AND NSTAR2
 call prompt('Enter file name for the setup: ', setupfile)
 call read_Nstar_from_setup(setupfile,ierr)

 !
 !--Check particle numbers
 if (Nstar(1) <= 0 .or. Nstar(2) <= 0) call fatal('moddump','Require particle numbers in both stars')
 !
 !--Request parameters (unless hardcoded to use defaults)
 ! now determine the parameters of their new orbit
 if (.not.use_defaults) then
    !call prompt('Enter desired separation:',separation,0.)
    call prompt('Use irrotational initial conditions?',useirrinit,.false.)              ! CHANGE to implement our desired initial conditions
 endif
 !
 !--Reset centre of mass location 
 call reset_centreofmass(Nstar(1),xyzh(:,1:Nstar(1)),vxyzu(:,1:Nstar(1)))
 call reset_centreofmass(Nstar(2),xyzh(:,Nstar(1)+1:npart),vxyzu(:,Nstar(1)+1:npart))

 mstar1 = Nstar(1) * massoftype(igas)
 mstar2 = Nstar(2) * massoftype(igas)
 mtotal = npart  * massoftype(igas)
 r1 = maxval(xyzh(1:3,1:Nstar(1)))
 r2 = maxval(xyzh(1:3,Nstar(1)+1:npart))

 separation = separation_from_Rochelobe(mstar1,mstar2,r1,r2)


 !
 !--Calcuate the new orbital parameters
 rad1   =  separation * mstar2/mtotal           ! distance of star 1 from the CoM
 rad2   =  separation * mstar1/mtotal           ! distance of star 2 from the CoM
 omega  =  sqrt(mtotal/separation**3)           ! orbital rotational speed
 if (useirrinit) then
    omega2 = -1.0 * omega                       ! spin velocity of the stars, assuming sinchronaized spin speeds
 else
    omega2 = -0.0 * omega
 endif
 
 !
 !--Place stars on new orbits
 !  for simplicity, assume stars are on the x-axis                             !Take here into account both options co-rotating and irrotational
 vxyzu(1:3,:) = 0.0                             ! reset velocity
 xcm1 = - rad1
 ycm1 = 0
 xcm2 = rad2
 ycm2 = 0
 do i=1,Nstar(1)
    xyzh(1,i)  =  xyzh(1,i) + xcm1
    vxyzu(1,i) =  -xyzh(2,i)*(omega+omega2)+omega2*ycm1
    vxyzu(2,i) =  xyzh(1,i)*(omega+omega2)-omega2*xcm1
    vxyzu(3,i) =  0.0d0
 enddo
 do i=Nstar(1)+1,npart
    xyzh(1,i)  = xyzh(1,i) + xcm2
    vxyzu(1,i) =  -xyzh(2,i)*(omega+omega2)+omega2*ycm2
    vxyzu(2,i) =  xyzh(1,i)*(omega+omega2)-omega2*xcm2
    vxyzu(3,i) =  0.0d0 
 enddo
 !
 !--Set new runtime parameters
 tmax           =   10.*2.*pi/omega           !Ten orbits
 dtmax          =   2.*pi/omega/50.           !50 timesteps per orbit
 damp           =    0.
 nfulldump      =    1
 iexternalforce =    0
 relflag        =   .false.
 !
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
