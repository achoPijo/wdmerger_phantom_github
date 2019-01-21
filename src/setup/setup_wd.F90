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
!  This module sets up a white dwarf star with a given mass and  
!  number of particles
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

 integer            :: Nstar(2) = 0  ! give default value in case dump header not read
 real(kind=8)       :: udist,umass
 character(len=20)  :: dist_unit,mass_unit
 logical            :: iexist
 logical            :: use_prompt = .true.
 logical            :: binary     = .true.
 integer            :: np         = 100000
 real               :: mstar      = 0.6d0
 real               :: mstar2     = 0.6d0
 real               :: Tin        = 1.d7
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  Setup routine for White Dwarfs
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use centreofmass,  only: reset_centreofmass 
 use eos,           only: ieos, equationofstate, init_eos,finish_eos,relflag
 use eos_helmholtz, only: tmaxhelmeos,tminhelmeos,xmass,speciesmax
 use dim,           only: maxp
 use io,            only: master
 use kernel,        only: hfact_default
 use options,       only: iexternalforce,nfulldump,damp,alphau
 use part,          only: igas,rhoh
 use physcon,       only: solarm,solarr,pi,planckh,mass_electron_cgs,mass_proton_cgs
 use prompting,     only: prompt
 use timestep,      only: tmax, dtmax
 use units,         only: set_units, select_unit, unit_pressure, unit_density
 use white_dwarf,   only: set_wd
 use io,            only: warning
 use nuc_reactions, only: nuc_burn

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
 real                             :: K,K_cgs,cvi,densi,dummyponrhoi,dummyspsoundi,tff,R1,R2
 character(len=120)               :: setupfile,inname
 logical                          :: write_setup
 integer                          :: i, ierr

 !
 !--Initializations
 !
 R1=0.
 R2=0.
 !
 ! Set default units
 !
 call set_units(dist=solarr,mass=solarm,G=1.d0)
 
 !
 ! Set equation of state
 !
 if (maxvxyzu == 5) then
    ieos = 15

 else
    if (maxvxyzu == 4) then
      ieos = 2
    else
      ieos = 1
    endif
 endif
 !
 ! set gamma
 !
 gamma = 5./3.
 !
 ! Polytropic constant K calculation for a white dwarf !REVISE cgs units!!
 !

 K_cgs = ((3.)**(2./3.)*(planckh**2))/(20.*(pi**(2./3.))*             &
               mass_electron_cgs*((mass_proton_cgs*mue)**(gamma)))
 K     = K_cgs*(unit_density**gamma)/(unit_pressure)       
 time  = 0.
 polyk = K
 hfact = hfact_default
 npart = 0
 npartoftype(:) = 0
 massoftype(:)  = 0.d0
 

 write_setup  = .false.
 !
 ! determine if an .in file exists !REVISE, dont think this is needed
 !
 inname=trim(fileprefix)//'.in'
 inquire(file=inname,exist=iexist)
 !
 ! determine if a .setup file exists
 !
 setupfile = trim(fileprefix)//'.setup'
 call read_setupfile(setupfile,ierr)
 if ( (ierr /= 0 .or. .not.iexist) .and. id==master) then
    ! setup file does not exist or is incomplete
    !
    ! Using prompts, determine the parameters the users wishes:
    !         mass unit, distance unit, binary setup?, mass of star(s), number of particles 
    ! MODIFIED TO ENTER AUTOMATICALLY SOME FIXED INPUT
    !ierr = 1
    !do while (ierr /= 0)
    !   call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
    !   call select_unit(mass_unit,umass,ierr)
    !   if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
    !enddo
    mass_unit = "solarm"
    umass = solarm
    !ierr = 1
    !do while (ierr /= 0)
    !   call prompt('Enter distance unit (e.g. solarr,au,pc,kpc,0.1pc)',dist_unit)
    !   call select_unit(dist_unit,udist,ierr)
    !   if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    !enddo
    dist_unit = "solarr"
    udist = solarr

    if (use_prompt) then 
       call prompt('Initial Temperature',Tin,tminhelmeos,tmaxhelmeos)
       call prompt('Set up a binary system?',binary)
       call prompt('Enter the total number of particles',np,0,maxp)
       if (binary) then
          call prompt('Enter the mass of star 1(code units)', mstar,0.0d0)
          call prompt('Enter the mass of star 2(code units)', mstar2,0.0d0)
       else
          call prompt('Enter the mass of the star (code units)', mstar,0.d0)
       endif
    endif
   


    
    write_setup = .true.
 endif
 !
 ! set units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)

 
 npart = np
 npartoftype(igas) = npart

 if (binary) then
    !
    ! Set the number of particles for each star
    !
    Nstar(1) = int(npart * mstar/(mstar+mstar2))
    Nstar(2) = npart - Nstar(1)
    !
    ! Generate each star
    !
    call set_wd(Nstar(1),hfact,mstar,xyzh(:,1:Nstar(1)))
    call set_wd(Nstar(2),hfact,mstar2,xyzh(:,Nstar(1)+1:npart))    
    !
    ! Separate the stars so they are isolated and can be relaxed without 
    ! interacting with each other
    !
    R1=maxval(xyzh(1,1:Nstar(1)))
    R2=maxval(xyzh(1,Nstar(1)+1:npart))
    do i=1,Nstar(1)
       xyzh(1,i)=xyzh(1,i)-100.
    enddo
    do i=Nstar(1)+1,npart
       xyzh(1,i)=xyzh(1,i)+100.
    enddo
    massoftype(igas) = (mstar+mstar2)/npart
 else
    Nstar(1) = npart
    call set_wd(npart,hfact,mstar,xyzh)
    massoftype(igas) = mstar/npart
    R1=maxval(xyzh(1,1:Nstar(1)))
 endif

 ! REVISE need to build an if clause depending on whether ISOTHERMAL is 
 ! declared or not. Our default for now is ISOTHERMAL = YES
 vxyzu(1:3,:) = 0.
 ! REVISE maybe for binary want to put both stars in orbit so if relaxation is too slow they wont fall together
 
 !
 ! call reset_centreofmass(npart,xyzh,vxyzu) ! We dont want to reset center of mass as this way we know where to look for the stars (+-1000 in x)
 !

 !
 ! add energies
 ! 
 
 call init_eos(ieos,ierr)

 ! set the mass weightings of each species
 ! currently depends on the star's mass
 ! TODO: update this to be set by user at runtime
 if (binary) then 
    call star_comp(xmass(:,1:Nstar(1)),Nstar(1),mstar)
    call star_comp(xmass(:,Nstar(1)+1:npart),Nstar(2),mstar2)
 else
    call star_comp(xmass(:,1:Nstar(1)),Nstar(1),mstar)
    
 endif

 do i=1,npart
    if (sum(xmass(1:speciesmax-1,i)) > 1.0+tiny(xmass) .or. sum(xmass(1:speciesmax-1,i)) < 1.0-tiny(xmass)) then
      call warning('eos_helmholtz', 'mass fractions total != 1')
      ierr = 1
      return
    endif
 enddo




 do i=1,npart
    densi = rhoh(xyzh(4,i),massoftype(igas))
    if (maxvxyzu==4) then
       if (gamma < 1.00001) then
          vxyzu(4,i) = polyk
       else
          !
          !  Note: Running the polytrope with u stored is not quite
          !  the same as using P = K rho^gamma because we really
          !  should use the actual rho, not rhozero.
          !           
          vxyzu(4,i) = polyk*densi**(gamma-1.)/(gamma-1.)
       endif
    endif
    if (maxvxyzu == 5) then
       vxyzu(5,i) = Tin
       call equationofstate(ieos,dummyponrhoi,dummyspsoundi,densi,xyzh(1,i),xyzh(2,i),xyzh(3,i), &
                            tempi=Tin,xmassi=xmass(:,i),cvi=cvi)
       vxyzu(4,i) = Tin*cvi
    endif
 enddo

 !
 ! Write setup file as exact mass will depend on the lane-emden solution iteration
 !
 call write_setupfile(setupfile)

 !
 !--Compute freefall time
 !
 !--If binary choose the frefall time from the less massive star
 !
 if (binary) then
    tff=max((pi/2.)*(R1**(3./2.))/sqrt(2.*mstar),(pi/2.)*(R2**(3./2.))/sqrt(2.*mstar2))
 else
    tff=(pi/2.)*(R1**(3./2.))/sqrt(2.*mstar)
 endif
 !
 !--Set new runtime parameters
 tmax           =    0.1                !2 time units
 dtmax          =    0.01   !1./utime   !1 second
 damp           =    1/tff              !
 nfulldump      =    1
 iexternalforce =    0
 alphau         =    0
 relflag        =    1
 nuc_burn       =    0

end subroutine setpart
!-----------------------------------------------------------------------
!+
!  Write setup parameters to input setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 use dim,          only: tagline
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20

 write(*,"(a)") ' Writing '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# '//trim(tagline)
 write(iunit,"(a)") '# input file for Phantom whitedwarf setup'

 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
 call write_inopt(utime,'time_unit','time unit (in seconds)',iunit)

 if (maxvxyzu == 5) then
    write(iunit,"(/,a)") '# Initial Temperature'
    call write_inopt(Tin,'Temperature','[K]',iunit)
 endif

 write(iunit,"(/,a)") '# binary setup?'
 call write_inopt(binary,'binary','Produce a binary star system',iunit)

 write(iunit,"(/,a)") '# number of particles and mass of star(s)'
 if (binary) then
    call write_inopt(Nstar(1)+Nstar(2),'ntotal','Total number of particles',iunit)
    call write_inopt(mstar,'mstar1','Mass of star 1 [code units]',iunit)
    call write_inopt(mstar2,'mstar2','mass of star 2 [code units]',iunit)
    call write_inopt(Nstar(1),'Nstar1','number of particles of star 1',iunit)
    call write_inopt(Nstar(2),'Nstar2','number of particles of star 2',iunit)
 else
    call write_inopt(Nstar(1),'ntotal','Total number of particles of the white dwarf',iunit)
    call write_inopt(mstar,'mstar','Mass of star 1 [code units]',iunit)
 endif

 close(iunit)

end subroutine write_setupfile
!-----------------------------------------------------------------------
!+
!  Read setup parameters from input setup file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use dim,          only: maxp
 use eos_helmholtz,only: tmaxhelmeos,tminhelmeos
 use infile_utils, only: open_db_from_file,inopts,close_db,read_inopt
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)
 real                          :: dummy
 !
 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(*, '(1x,2a)') 'Setup_wd: Reading setup options from ',trim(filename)
 !
 nerr = 0
 call read_inopt(dist_unit,'dist_unit',db,ierr)
 nerr= nerr + ierr
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 nerr= nerr + ierr
 call read_inopt(dummy,'utime',db,ierr)

 if (maxvxyzu == 5) then
    call read_inopt(Tin,'Temperature',db,ierr)
 endif

 call read_inopt(binary,'binary',db,ierr)
 
 if (binary) then
    call read_inopt(np,'ntotal',db,ierr)
    nerr= nerr + ierr
    call read_inopt(mstar,'mstar1',db,ierr)
    nerr= nerr + ierr
    call read_inopt(mstar2,'mstar2',db,ierr)
    nerr= nerr + ierr
 else
    call read_inopt(np,'ntotal',db,ierr)
    nerr= nerr + ierr
    call read_inopt(mstar,'mstar',db,ierr)
    nerr= nerr + ierr
 endif

 !
 ! parse units
 !
 call select_unit(mass_unit,umass,ierr)
 if (ierr /= 0) then
    call error('setup_wd','mass unit not recognised')
 endif
 call select_unit(dist_unit,udist,ierr)
 if (ierr /= 0) then
    call error('setup_wd','length unit not recognised')
 endif
 if (mstar < 0. .or. mstar2 < 0.) then
    call error('setup_wd','mass in setup file is smaller than 0.')
    ierr = 1 
    nerr = nerr + ierr
 endif
 if (np < 0 .or. np > maxp) then
    call error('setup_wd','number of particles in setup file out of range: 0 to maxp')
    ierr = 1
    nerr = nerr + ierr
 endif
 if (Tin < tminhelmeos .or. Tin > tmaxhelmeos) then
    call error('setup_wd','Temperature in input file outside valid range. Tmin: 1e3 K, Tmax:5e9 K')
    ierr = 1
    nerr = nerr + ierr
 endif

 if (nerr > 0) then
    print "(1x,a,i2,a)",'Setup_wd: ',nerr,' error(s) during read of setup file.'
    print *,"Please, insert next your setup options "
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile
!-----------------------------------------------------------------------
!+
!  Return star composition based on star mass
!+
!-----------------------------------------------------------------------
subroutine star_comp(composition,nstar,massstar)
 real, intent(inout) :: composition(:,:)
 integer, intent(in) :: nstar
 real,    intent(in) :: massstar
 integer             :: i

    if (massstar <= 0.45) then
       do i=1,nstar
          composition(:,i) = 0.0
          composition(2,i) = 1.0
       enddo
    else if (massstar > 0.45 .and. massstar <= 1.1 ) then
       do i=1,nstar
          composition(:,i) = 0.0
          composition(3,i) = 0.4
          composition(4,i) = 0.6
       enddo
    else
       do i=1,nstar
          composition(:,i) = 0.0
          composition(4,i) = 0.8
          composition(5,i) = 0.2
       enddo
    endif
 

end subroutine star_comp




end module setup

