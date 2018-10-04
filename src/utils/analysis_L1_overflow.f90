!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine check whether ddensity profile corresponds to that of 
!  a polytrope
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel Blanco
!
!  $Id: $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: 
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'Roche Lobe Overflow'
 logical :: firstcall = .true.

 public  :: do_analysis

 private :: read_mstar_from_setup,lanembden,derivs,rk4

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use io,           only: fatal
 use dim,          only: maxp,maxvxyzu
 use centreofmass, only: reset_centreofmass,get_centreofmass
 use part,         only: igas,iamtype,iphase,maxphase,rhoh
 use prompting,    only: prompt
 use setbinary,    only: L1_point
 use units,        only: umass,udist,utime,unit_density,unit_pressure
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer                      :: i,j,ierr,ncount,Nstar1,Nstar2,iloc
 real                         :: mstar1,mstar2,separation
 real                         :: xcom1(3),xcom2(3),vcom1(3),vcom2(3),dr(3)
 logical                      :: iexist
 character(len=120)           :: setupfile,filename,fileprefix
 character(len=200)           :: fileout

 !
 !-- Initialization
 !
 !--------------
 !iloc = index(dumpfile,'_0')
 !if (iloc > 1) then
 !   fileprefix = trim(dumpfile(1:iloc-1))
 !else
 !   fileprefix = trim(dumpfile)
 !endif
 !setupfile = trim(fileprefix)//'.setup'
 !print *, fileprefix
 !print *, setupfile
 !call prompt('Enter file name for the setup: file', setupfile)
 setupfile = 'wd.setup'
 inquire(file=trim(setupfile),exist=iexist)
 if (iexist) then
    call read_Nstars_from_setup(setupfile,ierr,Nstar1,Nstar2)
    if (ierr==1) call fatal('analysis','problem obtaining data from setup file')
 else
     call fatal('analysis','setup file does not exist')
 endif


 mstar1 = Nstar1 * particlemass
 mstar2 = Nstar2 * particlemass

 call reset_centreofmass(npart,xyzh(:,:),vxyzu(:,:))

 call get_centreofmass(xcom1,vcom1,Nstar1,xyzh(:,1:Nstar1),vxyzu(:,1:Nstar1))

 call get_centreofmass(xcom2,vcom2,Nstar2,xyzh(:,Nstar1+1:npart),vxyzu(:,Nstar1+1:npart))

 separation = sqrt((xcom1(1)-xcom2(1))**2+(xcom1(2)-xcom2(2))**2+(xcom1(3)-xcom2(3))**2)
 
 print *, L1_point(mstar1/mstar2)*separation
 !
 !-- Loop over all particles to check whether they are at a distance larger 
 !   than L1 distance from the primary and keep count. This way, by comparing
 !   the number of particles inside that region  and the original number of 
 !   particles of the primary we can determine if Roche Lobe overflow has 
 !   happened for the secondary.
 !
 ncount = 0
 do i=1,npart
    dr(:)=xyzh(1:3,i)-xcom1(:)
    if (sqrt(dr(1)**2+dr(2)**2+dr(3)**2) < L1_point(mstar1/mstar2)*separation) then
       ncount=ncount+1
    endif
 enddo
 
 print *, Nstar1
 print *, ncount


end subroutine do_analysis

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
 write(*, '(1x,2a)') 'Setup_wd: Reading setup mstar option from ',trim(filename)
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

subroutine read_mstar_from_setup(filename,ierr,mstar)
 use dim,               only: maxvxyzu
 use infile_utils,      only: open_db_from_file,inopts,close_db,read_inopt
 use io,                only: error

 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 real,             intent(out) :: mstar
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)
 real                          :: dummyr
 integer                       :: dummyint
 logical                       :: binary
 !
 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(*, '(1x,2a)') 'Setup_wd: Reading setup mstar option from ',trim(filename)
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
 if (binary) then
    write(*, '(a)') 'Error, binary system instead of a single star'
    ierr=1
    return
 endif

 nerr= nerr + ierr
 call read_inopt(dummyint,'ntotal',db,ierr)
 nerr= nerr + ierr
 call read_inopt(mstar,'mstar',db,ierr)
 nerr= nerr + ierr

if (nerr > 0) then
    call error ('setup_wd','there were some errors in reading the setup file')
    print "(1x,a,i2,a)",'Setup_wd: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_mstar_from_setup


      SUBROUTINE lanembden(maxit,rhoc,endit,masa,r,rho)
      
      use units,       only: udist, umass, utime, unit_density
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      DOUBLE PRECISION, DIMENSION(maxit), INTENT(OUT) :: rho, r
      DOUBLE PRECISION, INTENT(IN) :: rhoc
      INTEGER, INTENT(IN)  :: maxit
      INTEGER, INTENT(OUT) :: endit
!
!--Local variables
!
      DOUBLE PRECISION, DIMENSION(2) :: y, yout
      DOUBLE PRECISION :: x, h, t, alpha, gama, K, masa, pi,      &
                          alpha0    !, pres                                    !REVISE If it works without it just remove   
      INTEGER :: oup, flag, i, j
!
!--Physical parameters in cgs
!
      REAL(8), PARAMETER :: plank=6.6261d-27, me=9.1094d-28,            &
                            mp=1.6726d-24, G=6.6726d-8, msol=1.9891d33, &
                            rsol=6.955d10, mue=2.0d0, npol=1.5d0
!
!--External subroutines
!
!      EXTERNAL derivs               !REVISE NOT NEEDED?
!
!--Initial variables
!
      gama   = (npol+1.)/npol
      pi     = 4.0d0*datan(1.0d0)
      K      = ((3.)**(2./3.)*(plank**2))/(20.*(pi**(2./3.))*             &
               me*((mp*mue)**(gama)))
      alpha0 = (((npol+1.)*K)/(4.*pi*G))**(1./2.)
      alpha  = alpha0*rhoc**((1.-npol)/(2.*npol))
      y(1)   = 1.0d0
      y(2)   = 0.0d0
!
!--Start up the intergration
!
      t      = 0.0d0
      h      = 1.0d-6 
      r(1)   = alpha*t/udist
      rho(1) = rhoc*(y(1)**npol)/unit_density
!      pres   = K*(rho(1)**gama)                                               !REVISE If it works without it just remove 
      masa   = -4*pi*(alpha0**3.)*(rhoc**((3.-npol)/(2*npol)))*         &
               (t**2.)*(y(2))
      masa   = masa/umass
!
!--Integrate lane-embden equation
!
      j = 0
      i = 0
      DO WHILE (y(1) > 0)
         j = j + 1
         CALL rk4(npol,y,2,t,h,yout)
         t = t + h
         y = yout
         IF (MOD(j,10000) == 0) THEN
           i      = i + 1
           r(i)   = alpha*t/udist
           rho(i) = rhoc*(y(1)**npol)/unit_density
!           pres   = K*(rho(i)**gama)                                          !REVISE If it works without it just remove 
           masa   = -4*pi*(alpha0**3.)*(rhoc**((3.-npol)/(2*npol)))*    &
                    (t**2.)*(y(2))
           masa   = masa/umass
         ENDIF
      ENDDO
      endit = i
!
      END SUBROUTINE lanembden
!-----------------------------------------------------------------------
!+
!  Derivs
!+
!-----------------------------------------------------------------------
      SUBROUTINE derivs(npol,x,y,dydx,n)
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      INTEGER, INTENT(IN)                         :: n
      DOUBLE PRECISION, DIMENSION(n), INTENT(IN)  :: y
      DOUBLE PRECISION, INTENT(IN) :: npol, x
      DOUBLE PRECISION, DIMENSION(n), INTENT(OUT)  :: dydx
!
      dydx(1) = y(2)

      IF (DABS(x) < 1.0d-9) THEN
         dydx(2) = -(y(1)**npol)/3.0d0
      ELSE
         dydx(2) = -y(1)**npol-2.0d0*y(2)/x
      ENDIF
!
      END SUBROUTINE derivs
!-----------------------------------------------------------------------
!+
!  rk4
!+
!-----------------------------------------------------------------------
      SUBROUTINE rk4(npol,y,n,x,h,yout)
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      REAL(8), DIMENSION(n), INTENT(IN) :: y
      REAL(8), INTENT(IN) :: x, h, npol
      REAL(8), DIMENSION(n), INTENT(OUT) :: yout
!
!--Local variables
!
      REAL(8), DIMENSION(SIZE(y)) :: k1, k2, k3, k4, yt
      REAL(8) :: h6, hh, xh
      INTEGER :: n
!
      hh = h*0.5d0
      h6 = h/6.d0
!
      CALL derivs(npol,x,y,k1,n)
      xh = x + hh
      yt = y + hh*k1
!
      CALL derivs(npol,xh,yt,k2,n)
      yt = y + hh*k2
!
      CALL derivs(npol,xh,yt,k3,n)
      yt = y + h*k3
!
      if (yt(1) < 0.) then
         yout = yt
      else
         CALL derivs(npol,x+h,yt,k4,n)
         yout = y + h6*(k1 + k4 + 2.d0*(k2 + k3))
      endif
!
      END SUBROUTINE rk4
!-----------------------------------------------------------------------
!
end module
