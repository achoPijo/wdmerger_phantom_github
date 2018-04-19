!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: corot_binary_relaxation
!
!  DESCRIPTION:
!  Module containing utils for the relaxation of the binary system
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel Blanco
!
!  $Id: $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io,centreofmass,options,part,setubinary
!+
!--------------------------------------------------------------------------
module corot_binary_relaxation
 implicit none

 public  :: reduce_separation,compute_omega

 private :: read_Nstars_from_setup

contains
!--------------------------------------------------------------------------
subroutine reduce_separation(xyzh,vxyzu,npart,dt)
 use io,              only: fatal
 use centreofmass,    only: reset_centreofmass,get_centreofmass
 use options,         only: damp
 use part,            only: igas,massoftype
 use setbinary,       only: L1_point
 integer,          intent(in)    :: npart
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in)    :: dt

 integer                      :: i,ierr,ncount,Nstar1,Nstar2
 real                         :: mstar1,mstar2,separation
 real                         :: newseparation,separationdiff,dynfac

 real                         :: xcom1(3),xcom2(3),vcom1(3),vcom2(3),dr(3)
 logical                      :: iexist
 character(len=120)           :: setupfile

 !
 !-- Initialization
 !
 !--------------

! -- Read setup file for Nstar1 and Nstar2
 setupfile = 'wd.setup'
 inquire(file=trim(setupfile),exist=iexist)
 if (iexist) then
    call read_Nstars_from_setup(setupfile,ierr,Nstar1,Nstar2)
    if (ierr==1) call fatal('relax_binary_corotate','problem obtaining data from setup file')
 else
     call fatal('relax_binary_corotate','setup file does not exist')
 endif

 mstar1 = Nstar1 * massoftype(igas)
 mstar2 = Nstar2 * massoftype(igas)

 !-- Find the current separation between the centers of mass of both stars 
 call get_centreofmass(xcom1,vcom1,Nstar1,xyzh(:,1:Nstar1),vxyzu(:,1:Nstar1))

 call get_centreofmass(xcom2,vcom2,Nstar2,xyzh(:,Nstar1+1:npart),vxyzu(:,Nstar1+1:npart))

 separation = sqrt((xcom1(1)-xcom2(1))**2+(xcom1(2)-xcom2(2))**2+(xcom1(3)-xcom2(3))**2)


 !
 !-- Loop over all particles to check whether they are at a distance larger 
 !   than L1 distance from the primary and keep count. This way, by comparing
 !   the number of particles inside that region  and the original number of 
 !   particles of the primary we can determine if Roche Lobe overflow has 
 !   happened for the secondary.
 !
 ncount = 0
 do i=Nstar1+1,npart
    dr(:)=xyzh(1:3,i)-xcom1(:)
    if (sqrt(dr(1)**2+dr(2)**2+dr(3)**2) < L1_point(mstar1/mstar2)*separation) then
       ncount=ncount+1
    endif
 enddo
 
 !print *, ncount !ADD A STOPING CONDITION
 
 if (ncount > 0) call fatal('corot_binary_relaxation','Roche Lobe overflow achieved')
 !
 !-- Compute new separation
 !
 !-- dynfac is the factor relating the shrinking timescale with the dynamical timescale dyn/shrink
 dynfac = 20.
 !-- damp = 1/tff  and tff~=tdyn2
 separationdiff = separation * damp * dt / dynfac 
 newseparation  = separation - separationdiff
 
 !-- Adjust positions of both stars, setting the center of mass in the x axis again
 call reset_centreofmass(Nstar1,xyzh(:,1:Nstar1),vxyzu(:,1:Nstar1))
 call reset_centreofmass(Nstar2,xyzh(:,Nstar1+1:npart),vxyzu(:,Nstar1+1:npart))

 xyzh(1,1:Nstar1)       = xyzh(1,1:Nstar1)       - newseparation*mstar2/(mstar1+mstar2)
 xyzh(1,Nstar1+1:npart) = xyzh(1,Nstar1+1:npart) + newseparation*mstar1/(mstar1+mstar2)

end subroutine reduce_separation



subroutine compute_omega(xyzh,vxyzu,fxyzu,npart)
 use io,              only: fatal
 use centreofmass,    only: reset_centreofmass,get_centreofmass
 use extern_corotate, only: omega_corotate 
 use part,            only: igas,massoftype
 
 integer,          intent(in) :: npart
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:)

 integer                      :: i,ierr,Nstar1,Nstar2
 real                         :: mstar1,mstar2,omega1,omega2

 real                         :: xcom1(3),xcom2(3),vcom1(3),vcom2(3),fcm1(3),fcm2(3)
 logical                      :: iexist
 character(len=120)           :: setupfile

 !
 !-- Initialization
 !
 !--------------

! -- Read setup file for Nstar1 and Nstar2
 setupfile = 'wd.setup'
 inquire(file=trim(setupfile),exist=iexist)
 if (iexist) then
    call read_Nstars_from_setup(setupfile,ierr,Nstar1,Nstar2)
    if (ierr==1) call fatal('relax_binary_corotate','problem obtaining data from setup file')
 else
     call fatal('relax_binary_corotate','setup file does not exist')
 endif

 mstar1 = Nstar1 * massoftype(igas)
 mstar2 = Nstar2 * massoftype(igas)

 !-- Get positions of the center of mass of both stars 
 call get_centreofmass(xcom1,vcom1,Nstar1,xyzh(:,1:Nstar1),vxyzu(:,1:Nstar1))

 call get_centreofmass(xcom2,vcom2,Nstar2,xyzh(:,Nstar1+1:npart),vxyzu(:,Nstar1+1:npart))

 !
 !-- Loop over each star to obtain fcm1 and fcm2
 !
 fcm1(:) = 0.
 fcm2(:) = 0.
 do i=1,Nstar1
    fcm1(:) = fcm1(:) + massoftype(igas)*fxyzu(1:3,i)
 enddo

 do i=Nstar1+1,npart
    fcm2(:) = fcm2(:) + massoftype(igas)*fxyzu(1:3,i)
 enddo
 
 fcm1(:)=fcm1(:)/mstar1
 fcm2(:)=fcm2(:)/mstar2

 omega1=sqrt(sqrt(fcm1(1)**2+fcm1(2)**2+fcm1(3)**2)/(mstar1*sqrt(xcom1(1)**2)))
 omega2=sqrt(sqrt(fcm2(1)**2+fcm2(2)**2+fcm2(3)**2)/(mstar2*sqrt(xcom2(1)**2)))

 
 print *, '-----------------------------------------'
 print *, omega1
 print *, omega2
 print *, '-----------------------------------------'

 omega_corotate = (omega1+omega2)/2  !ATTENTION WILL THIS BE REFLECTED IN THE .in file?

end subroutine compute_omega




subroutine read_Nstars_from_setup(filename,ierr,Nstar1,Nstar2)
 use dim,               only: maxvxyzu
 use infile_utils,      only: open_db_from_file,inopts,close_db,read_inopt
 use io,                only: error

 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 integer,          intent(out) :: Nstar1,Nstar2
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)
 real                          :: dummyr
 integer                       :: dummyint
 logical                       :: binary
 !
 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 !write(*, '(1x,2a)') 'Setup_wd: Reading setup mstar option from ',trim(filename)
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
 if (.not.binary) then
    write(*, '(a)') 'Error, single star instead of a binary system'
    ierr=1
    return
 endif
 nerr= nerr + ierr
 call read_inopt(dummyint,'ntotal',db,ierr)
 nerr= nerr + ierr
 call read_inopt(dummyr,'mstar1',db,ierr)
 nerr= nerr + ierr
 call read_inopt(dummyr,'mstar2',db,ierr)
 nerr= nerr + ierr
 call read_inopt(Nstar1,'Nstar1',db,ierr)
 nerr= nerr + ierr
 call read_inopt(Nstar2,'Nstar2',db,ierr)
 nerr= nerr + ierr

 if (nerr > 0) then
    call error ('setup_wd','there were some errors in reading the setup file')
    print "(1x,a,i2,a)",'Setup_wd: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_Nstars_from_setup


!-----------------------------------------------------------------------
!
end module
