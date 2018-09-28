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

 !
 ! Set default units
 !
 call set_units(dist=cm,mass=g,time=s)
 
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
 write(*, '(1x,2a)') 'Setup_box: Reading setup options from ',trim(filename)
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
    call error('setup_box','mass unit not recognised')
 endif
 call select_unit(dist_unit,udist,ierr)
 if (ierr /= 0) then
    call error('setup_box','length unit not recognised')
 endif
 if (np < 0 .or. np > maxp) then
    call error('setup_box','number of particles in setup file out of range: 0 to maxp')
    ierr = 1
    nerr = nerr + ierr
 endif

 if (nerr > 0) then
    print "(1x,a,i2,a)",'Setup_box: ',nerr,' error(s) during read of setup file.'
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

