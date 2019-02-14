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
 real,    private   :: separation           = 0.06
 logical, private   :: use_irr_init         = .false.
 logical, private   :: use_corotating_frame = .true.
 logical, private   :: binary               = .true.
 logical, private   :: use_relocation       = .true.

 public             :: modify_dump
 private

contains
!-----------------------------------------------------------------------

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use eos,                     only: relflag, equationofstate,init_eos
 use io,                      only: iprint,fatal
 use prompting,               only: prompt
 use options,                 only: iexternalforce,nfulldump,damp,alphau
 use part,                    only: igas,rhoh
 use units,                   only: unit_velocity,utime
 use physcon,                 only: c,pi
 use timestep,                only: time,tmax, dtmax
 use centreofmass,            only: get_centreofmass,reset_centreofmass                 
 use externalforces,          only: iext_corotate
 use extern_corotate,         only: omega_corotate,dynfac
 use extern_gwinspiral,       only: Nstar
 use nuc_reactions,           only: nuc_burn
 use eos_helmholtz,           only: xmass,speciesmax
 integer, intent(inout)    :: npart
 integer, intent(inout)    :: npartoftype(:)
 real,    intent(inout)    :: massoftype(:)
 real,    intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                   :: i,ierr
 real                      :: xcom1(3),xcom2(3),vcom1(3),vcom2(3)
 real                      :: rad1,rad2,mstar1,mstar2,mtotal, omega, omega2    !CHNGCODE added omega
 real                      :: xcm1,xcm2,ycm1,ycm2,r1,r2,densi,Tnew,dummyponrhoi,dummyspsoundi,cvi
 character(len=120)        :: setupfile

 call prompt('Is this a binary setup?', binary)

 if (binary) then

    setupfile='wd.setup'
    !-- PROMPT TO ASK FOR THE SETUP FILE TO CHECK THE VALUES OF NSTAR1 AND NSTAR2
    call prompt('Enter file name for the setup: ', setupfile)

    call read_Nstar_from_setup(setupfile,ierr)
   
    !--Check particle numbers
    if (Nstar(1) <= 0 .or. Nstar(2) <= 0) call fatal('moddump','Require particle numbers in both stars')
    !
    !--Request parameters (unless hardcoded to use defaults)
    ! 
    if (.not.use_defaults) then
       !call prompt('Enter desired separation:',separation,0.)
       call prompt('Use irrotational initial conditions?',use_irr_init)
       !-- Relocation?
       call prompt('Relocate centers of mass?',use_relocation)
       !-- Desired separation between centers of mass?
       if (use_relocation) call prompt('Enter desired separtation',separation)
       !-- Use a corrotating frame?
       call prompt('Use a corotating frame?',use_corotating_frame)
       !--
       if (use_corotating_frame) call prompt ('Enter desired dynfac',dynfac)

       call prompt('Alpha u',alphau)
    endif

    !-- Compute masses
    mstar1 = Nstar(1) * massoftype(igas)
    mstar2 = Nstar(2) * massoftype(igas)
    mtotal = npart  * massoftype(igas)
    
    !-- Set temperature of stars to a new values

    !call init_eos(15,ierr)

    !do i=1,npart
    !   Tnew = 100000
    !   densi = rhoh(xyzh(4,i),massoftype(igas))
    !   vxyzu(5,i) = Tnew
    !   call equationofstate(15,dummyponrhoi,dummyspsoundi,densi,xyzh(1,i),xyzh(2,i),xyzh(3,i), &
    !                        tempi=Tnew,xmassi=xmass(:,i),cvi=cvi)
    !   vxyzu(4,i) = Tnew*cvi

    !enddo

   
   
    if (use_relocation) then
       !--Reset centre of mass location 
       call reset_centreofmass(Nstar(1),xyzh(:,1:Nstar(1)),vxyzu(:,1:Nstar(1)))
       call reset_centreofmass(Nstar(2),xyzh(:,Nstar(1)+1:npart),vxyzu(:,Nstar(1)+1:npart))

       !-- Obtain new separation from Roche Lobe
       !r1 = maxval(xyzh(1:3,1:Nstar(1)))
       !r2 = maxval(xyzh(1:3,Nstar(1)+1:npart))
       !separation = separation_from_Rochelobe(mstar1,mstar2,r1,r2)

       !--Calcuate the new orbital parameters
       rad1   =  separation * mstar2/mtotal           ! distance of star 1 from the CoM
       rad2   =  separation * mstar1/mtotal           ! distance of star 2 from the CoM
       omega  =  sqrt(mtotal/separation**3)           ! orbital rotational speed
       if (use_irr_init) then
          omega2 = -1.0 * omega                       ! spin velocity of the stars, assuming sinchronaized spin speeds
       else
          omega2 = -0.0 * omega
       endif

       xcm1 = - rad1
       ycm1 = 0
       xcm2 = rad2
       ycm2 = 0

    else
       call get_centreofmass(xcom1,vcom1,Nstar(1),xyzh(:,1:Nstar(1)),vxyzu(:,1:Nstar(1)))
       call get_centreofmass(xcom2,vcom2,Nstar(2),xyzh(:,Nstar(1)+1:npart),vxyzu(:,Nstar(1)+1:npart))
       separation = sqrt((xcom1(1)-xcom2(1))**2+(xcom1(2)-xcom2(2))**2+(xcom1(3)-xcom2(3))**2)
       print *, separation
       print *, mtotal
       omega  =  sqrt(mtotal/separation**3) 
       print *, omega
       xcm1=xcom1(1)
       xcm2=xcom2(1)
       ycm1=xcom1(2)
       ycm2=xcom2(2)
    endif

    !-- Determine individual star spin
    if (use_irr_init) then
       omega2 = -1.0 * omega  !-- spin velocity of the stars, assuming sinchronized spin speeds
    else
       omega2 = -0.0 * omega
    endif

    
    if (use_corotating_frame) then

       vxyzu(1:3,:) = 0.0                             ! reset velocity

       do i=1,Nstar(1)
          if (use_relocation) xyzh(1,i)  =  xyzh(1,i) + xcm1
          vxyzu(1,i) =  omega2*(-xyzh(2,i)+ycm1)
          vxyzu(2,i) =  omega2*(xyzh(1,i)-xcm1)
          vxyzu(3,i) =  0.0d0
       enddo
       do i=Nstar(1)+1,npart
          if (use_relocation) xyzh(1,i)  =  xyzh(1,i) + xcm2
          vxyzu(1,i) =  omega2*(-xyzh(2,i)+ycm2)
          vxyzu(2,i) =  omega2*(xyzh(1,i)-xcm2)
          vxyzu(3,i) =  0.0d0 
       enddo  
      
       omega_corotate = omega
       iexternalforce = iext_corotate
       !
       !--Set new runtime parameters
       tmax           =   50.*2.*pi/omega                !50 orbits large enough for RLO to occur
       dtmax          =   2.*pi/(sqrt(mtotal/(separation/3)**3))/50.   !50 timesteps per aprox final orbit

    else
                           
       vxyzu(1:3,:) = 0.0                             ! reset velocity

       do i=1,Nstar(1)
          if (use_relocation) xyzh(1,i)  =  xyzh(1,i) + xcm1
          vxyzu(1,i) =  -xyzh(2,i)*(omega+omega2)+omega2*ycm1
          vxyzu(2,i) =  xyzh(1,i)*(omega+omega2)-omega2*xcm1
          vxyzu(3,i) =  0.0d0
       enddo
       do i=Nstar(1)+1,npart
          if (use_relocation) xyzh(1,i)  =  xyzh(1,i) + xcm2
          vxyzu(1,i) =  -xyzh(2,i)*(omega+omega2)+omega2*ycm2
          vxyzu(2,i) =  xyzh(1,i)*(omega+omega2)-omega2*xcm2
          vxyzu(3,i) =  0.0d0 
       enddo
       !
       !--Set new runtime parameters
       tmax           =  0.800 !50.*2.*pi/omega   !50 orbits large enough for RLO to occur
       dtmax          =  0.001 !2.*pi/omega/50.   !50 timesteps per  orbit
       
       iexternalforce = 0
       damp           = 0.
    endif
      

    !
   
 else
    !
    !--Restart velocity                     
    vxyzu(1:3,:) = 0.0
    !--Reset centre of mass location 
    call reset_centreofmass(npart,xyzh(:,:),vxyzu(:,:))
   
                                ! reset velocity
    !
    !--Set new runtime parameters
    tmax           =    2.0                       !2 time units
    dtmax          =    0.1   !1./utime           !1 second
    iexternalforce =    0
    damp           =    0.
    !
 endif
 
 !damp           =    0.
 if (use_defaults) then
    alphau      =    0.000
 endif
 nfulldump      =    1
 if (use_corotating_frame) then
    relflag        =    1  
 else
    relflag        =    4
 endif
 nuc_burn       =    0

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
