!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION: Determines the centre of mass of the system and writes it
!               to file
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel Blanco
!
!  $Id: 710fd2514e754b4af37500156be0a1839dc2c485 $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, part
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'centerofmass'
 logical, private :: firstcall = .true.

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,                 only: nptmass,xyzmh_ptmass,vxyz_ptmass
 use centreofmass,         only: get_centreofmass
 use extern_gwinspiral,    only: Nstar
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 real                         :: xpos1(3),vpos1(3),xpos2(3),vpos2(3)
 logical                      :: iexist
 character(len=200)           :: fileout,setupfile
 integer                      :: ierr
 
 !-- PROMPT TO ASK FOR THE SETUP FILE TO CHECK THE VALUES OF NSTAR1 AND NSTAR2
 call prompt('Enter file name for the setup: ', setupfile)
 !setupfile='wd.setup'
 call read_Nstar_from_setup(trim(setupfile),ierr)

 !--Check particle numbers
 if (Nstar(1) <= 0 .or. Nstar(2) <= 0) call fatal('moddump','Require particle numbers in both stars')

 !
 ! Open file (appendif exists)
 !
 fileout = trim(dumpfile(1:INDEX(dumpfile,'_')-1))//'distances.dat'
 inquire(file=fileout,exist=iexist)
 if ( .not.iexist .or. firstcall ) then
    firstcall = .false.
    open(iunit,file=fileout,status='replace')
    write(iunit,"('#',2(1x,'[',i2.2,1x,a11,']',2x))") &
          1,'time', &
          2,'distance'

 else
    open(iunit,file=fileout,position='append')
 endif
 !
 ! Call centre of mass subroutine
 !
 call get_centreofmass(xpos1,vpos1,Nstar(1),xyzh(1:Nstar(1)),vxyzu(1:Nstar(1)),nptmass,xyzmh_ptmass,vxyz_ptmass)
 call get_centreofmass(xpos2,vpos2,Nstar(2),xyzh(Nstar(1)+1:npart),vxyzu(Nstar(1)+1:npart),nptmass,xyzmh_ptmass,vxyz_ptmass)
 !
 ! compute distance between centers of mass of both stars
 !
 distance = sqrt((xpos1(1)-xpos2(1))**2+(xpos1(2)-xpos2(2))**2+(xpos1(3)-xpos2(3))**2)
 !
 ! Write results to file
 !
 write(iunit,'(2(es18.10,1x))') time,distance

 close(iunit)

end subroutine do_analysis

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

end module
