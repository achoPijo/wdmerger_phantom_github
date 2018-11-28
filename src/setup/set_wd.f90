!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: spherical
!
!  DESCRIPTION:
!   This module sets up white dwarf particle distributions
!   
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel Blanco
!
!  $Id: 431c9f0e08d875218af7283673e42735cc61f474 $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: io, physcon, prompting, stretchmap, unifdis
!+
!--------------------------------------------------------------------------
module white_dwarf
 use physcon,    only:pi,twopi,fourpi,piontwo
! use unifdis,    only:set_unifdis
 use io,         only:fatal,warning
 !
 implicit none
 !
 public  :: set_wd
 !
 private :: lanembden, derivs, rk4
 !
contains
!
!-----------------------------------------------------------------------
!+
!  This subroutine positions particles on a sphere with the required 
!  density profile of a white dwarf (complying with lane embden equation)
!  for a given mass
!+
!-----------------------------------------------------------------------
subroutine set_wd(nbody, hfact, masstot, xyzh)
! use io,         only:iprint
 use eos_helmholtz, only: cgsrhomaxhelmeos, cgsrhominhelmeos
 use units,   only:udist, umass, utime, unit_density
 use IFPORT  ! Necessary for the random number generator

 integer,          intent(inout) :: nbody
 real,             intent(in)    :: hfact
 real,             intent(inout) :: masstot                                    
 real,             intent(out)   :: xyzh(:,:)

 integer,          parameter     :: maxit = 1000000
 real,             parameter     :: tolerance = 1.e-3   ! mod_parameters.f90
 integer                         :: p, i, j, l, endit, iseed !k,
 logical                         :: done
! integer,          dimension(nbody) :: nn

 double precision, dimension(maxit) :: radle, densle, rho 
 double precision                   :: xc, yc, zc, dh, dh2, rx, ry, rz, rnk, radmax, &
                                    dr2, rad2, radmin, temp, massp, sum1, sum2,   &
                                    rhoc_min, rhoc_max, rhoc, mass_min, mass_max, &
                                    mass, REL
!-----------------------------------Needed for T and compostion init
! integer                         :: NNNZ
! character(len=13)               :: ONC
!----------------------------------------------------------------------------- 
!
!--Load modules
!
!      USE mod_essentials
!      USE mod_parameters, ONLY : nbody, maxit, unm, uen, uden, hfact,   &
!                                 ka_min, ka_max
!      USE mod_commons
!      USE IFPORT  ! Necessary for the random number generator
!
!     --Include EOS definitions
!
!      INCLUDE 'vector_eos.dek'   REVISE
!

!=====================================================================
!  SUBROUTINE TO GENERATE THE INITIAL STATE OF THE PARTICLES THAT WILL
!  BE USED IN THE SPH PROGRAM
!
!  Last revision: 21/April/2016
!=====================================================================

!
!--Read initial model parameters
!
!      OPEN (UNIT=1,FILE='genparams.txt',STATUS='old')
!      READ(1,*) nbod
!      READ(1,*) masstot
!      READ(1,*) rnk          REVISE (no se si es necesario en phantom)
!      READ(1,*) temp         REVISE donde inicializar la temperatura?
!      READ(1,*) tnow         
!      CLOSE(1)
!
!--Read and setup WD composition                                               ! VAMOS A TRATAR DE IMPLEMENTAR SOLO LA POSICION POR AHORA
!
!      OPEN(UNIT=2,FILE='cdata1.txt',STATUS='old')
!      DO k=1,nel
!         READ(2,*) xss(k,1)
!         xss(k,2:nbody) = xss(k,1)
!      ENDDO
!      CLOSE(2)
!
!--Use a Bisection strategy in order to solve Lane-Embden 
!  equation for the given mass
!
      print *, 11
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
         CALL lanembden(maxit,rhoc_min,endit,mass,radle,densle)
         mass_min = mass
!
         CALL lanembden(maxit,rhoc_max,endit,mass,radle,densle)
         mass_max = mass
!
         rhoc = 0.5*(rhoc_min + rhoc_max)
         CALL lanembden(maxit,rhoc,endit,mass,radle,densle)
!
         IF (masstot > mass) THEN
            rhoc_min = rhoc 
         ELSEIF (masstot < mass) THEN
            rhoc_max = rhoc
         ENDIF 
!               
         REL = ABS(masstot-mass)/masstot
         IF (REL <= tolerance) THEN
            done = .TRUE.
         ENDIF
         WRITE(*,'(A5,1(1F7.4),A8,1F7.4)') 'M_WD=',masstot,' M_Bisec=',mass
      ENDDO whileloop1
      print *, 12
!
!--Call lanembden subroutine a last time, saving this time the
!  values for density, radius, pressure, etc
!
      CALL lanembden(maxit,rhoc,endit,mass,radle,densle)
      masstot = mass
!
!--Generate random particle distribution
!
      radmax = radle(endit)
      rad2   = radmax*radmax
      massp  = mass/FLOAT(nbody)
!!$OMP PARALLEL DEFAULT(none) shared(xyzhm,vxyzut,radle,densle,rho)          &
!!$OMP shared(massp,rad2,radmax,temp,endit,ka1,fh) private(p,dr2,xc,yc,zc,i) &
!!$OMP private(done,iseed)
!!$OMP DO SCHEDULE(runtime)
      print *, 13
      partloop1 :DO p=1,nbody
         print *, 131
         dr2 = huge(xc)
         CALL SEED(p*time())         
         DO WHILE (dr2 > rad2)
            xc  = 2.*radmax*RAND()-radmax
            yc  = 2.*radmax*RAND()-radmax
            zc  = 2.*radmax*RAND()-radmax
            dr2 = xc*xc+yc*yc+zc*zc
         ENDDO
         xyzh(1,p)  = xc
         xyzh(2,p)  = yc
         xyzh(3,p)  = zc
!         vxyzu(1,p) = 0.
!         vxyzu(2,p) = 0. REVISE dont know if needed to initialize here
!         vxyzu(3,p) = 0.
!         ka1(p)      = ka_min
!         fh(p)       = 1.        !       fh is not used in PHANTOM as this
!
!--Compute density using Lane-Embden equation solution
!
         print *, 132
         i = 0
         done = .FALSE.
         whileloop2 : DO WHILE (done.EQV..FALSE.)
            i = i + 1
!
            IF ((radle(i) <= DSQRT(dr2)).AND.                       &
               (DSQRT(dr2) <= radle(i+1))) then
               rho(p)      = densle(i)                                        !REVISE CHECK UNITS OF DENSLE AND RHO
               if (rho(p) < cgsrhominhelmeos/unit_density) then
                  rho(p) = cgsrhominhelmeos/unit_density
                  !print *, "low limit"
               endif
               if (rho(p) > cgsrhomaxhelmeos/unit_density) then
                  rho(p) = cgsrhomaxhelmeos/unit_density
                  !print *, "high limit"
               endif
               xyzh(4,p)   = hfact*((massp/rho(p))**(1./3.))  !hfact*((xyzhm(5,p)/rho(p))**(1./3.))
               done        = .TRUE.
            ENDIF
!
            IF (i == endit) THEN
               PRINT*, 'Particle',p,' not found'
               STOP
            ENDIF
         ENDDO whileloop2
         print *, 133
      ENDDO partloop1
      print *, 14
!!$OMP END DO
!!$OMP END PARALLEL
!
!--Sort particles. Since now factor is a dynamical quantity we need to
!  call sort in order to define it !!
!
!      CALL sort_part()            REVISE
!
!--Now, use the tree to self-consistently calculate the rho-h relation 
!-----This should be computed on first run in PHANTOM
!
!      PRINT*,'Trying to compute h. May take a long time'
!      nstep = 1
!      nout  = 1
!      CALL iter_rhoh
!      PRINT*, 'h calculation finished'
!
!--Read EOS tables
!
!      CALL read_helm_table
!
!      OPEN(UNIT=4,FILE='ZHELI.10b',status='old')
!      aion(1) = 1.0d0
!      zion(1) = 1.0d0
!      DO k=1,nel-2
!         READ(4,'(I4,1X,A5,1X,F4.0,1X,F4.0)') NNNZ,ONC,aion(k+1),zion(k+1)
!      ENDDO
!      CLOSE(UNIT=4)
!
!--Calculate thermal energy Probably add heere an if depending on ISOTHERMAL or not initialize internal energy or not
!      left for later to add T and thermal energy
!
!      DO p = 1, nbody
!         CALL degenerate(p)
!      ENDDO
!      PRINT*,'EOS finished'
!
end subroutine set_wd
!-----------------------------------------------------------------------
!+
!  LANE EMBDEN routine 
!  Please note:
!  rhoc [cgs]
!  masa [code units]
!  r    [code units]
!  rho  [code units]
!+
!-----------------------------------------------------------------------
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


!--------------------------------------------------------------------------
end module white_dwarf
