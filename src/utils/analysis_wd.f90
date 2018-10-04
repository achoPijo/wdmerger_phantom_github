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
 character(len=20), parameter, public :: analysistype = 'wdrelaxation'
 logical :: firstcall = .true.

 public  :: do_analysis

 private :: read_mstar_from_setup,lanembden,derivs,rk4

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use io,           only: fatal
 use dim,          only: maxp,maxvxyzu
 use centreofmass, only: reset_centreofmass
 use physcon,      only: pi,planckh,mass_electron_cgs,mass_proton_cgs
 use part,         only: igas,iamtype,iphase,maxphase,rhoh
 use prompting,    only: prompt
 use units,        only: umass,udist,utime,unit_density,unit_pressure
 use rho_profile,  only: rho_polytrope
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer                      :: i,j,nj,npts,ierr
 integer,           parameter :: nmaxpoints=1000000
 real, dimension(nmaxpoints)  :: rtab, rhotabpoly, rhotabdump,radle,densle
 real,              parameter :: mue=2.0d0
 real                         :: dr,r, dens, deni, denj
 real                         :: gamma, polyk, Mstar, rhocentre
 real                         :: K_cgs, K
 real                         :: rhoc_min,rhoc_max,rhoc,mass,mass_min,mass_max,rel,masstot
 real,              parameter :: tolerance=1.e-3
 logical                      :: done
 real                         :: npartInH,xInH,yInH,zInH,dummyr,dummyr1,dummyr2,dummyr3,dummyr4
 real                         :: dummyr5,dummyr6,dummyr7,dummyr8,dummyr9,dummyr10,dummyr11,dummyr12,dummyr13,dummyr14
 integer                      :: dummyint,dummyint1
 real, dimension(nmaxpoints)  :: rInH,rhoInH,rhotabInH
 logical                      :: iexist
 character(len=120)           :: setupfile,filename
 character(len=200)           :: fileout

 !
 !-- Initialization
 !
 gamma = 5./ 3.
 K_cgs = ((3.)**(2./3.)*(planckh**2))/(20.*(pi**(2./3.))*             &
               mass_electron_cgs*((mass_proton_cgs*mue)**(gamma)))
 K     = K_cgs*(unit_density**gamma)/(unit_pressure)       
 polyk = K

 call prompt('Enter file name for the setup: ', setupfile)
 call read_mstar_from_setup(setupfile,ierr,Mstar)
 if (ierr==1) call fatal('analysis','problem obtaining mstar from setup file')
 print *, Mstar
 !
 !-- Get density profile for polytrope
 !

 call rho_polytrope(gamma,polyk,Mstar,rtab,rhotabpoly,npts,rhocentre)

 !
 !-- Get density profile from dump file
 !
 dr = rtab(2)-rtab(1)
 do j=1,npts
    denj = 0.
    nj = 0
    do i=1,npart
       r = sqrt(xyzh(1,i)**2+xyzh(2,i)**2+xyzh(3,i)**2)
       deni = rhoh(xyzh(4,i),particlemass)
       if (r < rtab(j)+dr/2 .and. r > rtab(j)-dr/2) then
          denj = denj + deni
          nj = nj + 1
       endif
    enddo
    if (nj /= 0) denj = denj/nj
    rhotabdump(j) = denj
 enddo
 !
 !
 ! Write results to file
 !
 fileout = trim(dumpfile)//'_densityprofile.dat'
 open(iunit,file=fileout,status='replace')
 write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'r',&
        2,'rho_le',   &
        3,'rho_sph'

 do j=1,npts
      write(iunit,'(3(1pe18.10,1x))') rtab(j),rhotabpoly(j),rhotabdump(j)
 enddo
 close(iunit)

 !
 !-- Compute lane emden solution numerically with INHouse code 
 !
 
 !--Use a Bisection strategy in order to solve Lane-Embden 
!  equation for the given mass
!
      WRITE(*,*) ' '
      WRITE(*,*) '************************************************** '
      WRITE(*,*) 'Using a bisection scheme to solve the lane-embden'
      WRITE(*,*) 'equation and generate a first density vs radius'
      WRITE(*,*) 'profile of the WD...'
      WRITE(*,*) '************************************************** '
      WRITE(*,*) ' '
      done     = .FALSE.
      rhoc_min = 1.0d4        !cgs
      rhoc_max = 1.0d8        !cgs
      whileloop1 : DO WHILE (done.EQV..FALSE.)
!
         CALL lanembden(nmaxpoints,rhoc_min,npts,mass,radle,densle)
         mass_min = mass
!
         CALL lanembden(nmaxpoints,rhoc_max,npts,mass,radle,densle)
         mass_max = mass
!
         rhoc = 0.5*(rhoc_min + rhoc_max)
         CALL lanembden(nmaxpoints,rhoc,npts,mass,radle,densle)
!
         IF (Mstar > mass) THEN
            rhoc_min = rhoc 
         ELSEIF (Mstar < mass) THEN
            rhoc_max = rhoc
         ENDIF 
!               
         REL = ABS(Mstar-mass)/Mstar
         IF (REL <= tolerance) THEN
            done = .TRUE.
         ENDIF
         WRITE(*,'(A5,1(1F7.4),A8,1F7.4)') 'M_WD=',Mstar,' M_Bisec=',mass
      ENDDO whileloop1
!
!--Call lanembden subroutine a last time, saving this time the
!  values for density, radius, pressure, etc
!
      CALL lanembden(nmaxpoints,rhoc,npts,mass,radle,densle)
      masstot = mass




 !
 !-- READ Result from Inhouse code
 !
 filename= 'star08relaxedInHouse.dat'
 OPEN (UNIT=1, FILE=filename)
!
!--Read time header
!
 READ(1,*)
 READ(1,*) dummyr
!
!--Read Gamma 
!
 READ(1,*)
 READ(1,*) dummyr
!
!--Read number of dust and gas particles
!
 READ(1,*)
 READ(1,*) npartInH
!
!--do loop
!
 partloop : DO j = 1,npartInH
!
!--SPH particle data
!
    READ(1,*) xInH,yInH,zInH,dummyr3,dummyr4,rhoInH(j),dummyr5,dummyr6,dummyr7,dummyr8,dummyr9, &
    dummyr10,dummyr11,dummyr12,dummyr13,dummyr14,dummyint,dummyint1
    rInH(j)=sqrt(xInH**2+yInH**2+zInH**2)*0.1
    rhoInH(j)=rhoInH(j)*6000/unit_density
 ENDDO partloop
 close(UNIT=1)

 do j=1,npts
    dr = radle(j+1)-radle(j-1)
    denj = 0.
    nj = 0
    do i=1,npartInH
       if (rInH(i) < radle(j)+dr/2 .and. rInH(i) > radle(j)-dr/2) then
          denj = denj + rhoInH(i)
          nj = nj + 1
       endif
    enddo
    if (nj /= 0) denj = denj/nj
    rhotabInH(j) = denj
 enddo

!
!-- Write results to file
!

 fileout = trim(dumpfile)//'_densityprofileINHOUSE.dat'
 open(iunit,file=fileout,status='replace')
 write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'r',&
        2,'rho_le_INHOUSE', &
        3,'rho_sph_INHOUSE'

 do j=1,npts
      write(iunit,'(3(1pe18.10,1x))') radle(j),densle(j),rhotabInH(j)
 enddo
 close(iunit)


end subroutine do_analysis


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
