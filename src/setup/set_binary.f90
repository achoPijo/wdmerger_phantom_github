!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setbinary
!
!  DESCRIPTION:
!   This module is contains utilities for setting up binaries
!
!  REFERENCES:
!   Eggleton (1983) ApJ 268, 368-369 (ref:eggleton83)
!   Lucy (2014), A&A 563, A126
!   https://en.wikipedia.org/wiki/Orbital_elements
!
!  OWNER: Daniel Price
!
!  $Id: 1663dc134ad33587b337a4c2c3e17d5f9220e8eb $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part, physcon
!+
!--------------------------------------------------------------------------
module setbinary
 use physcon, only:pi
 implicit none
 public :: set_binary,Rochelobe_estimate,L1_point,get_a_from_period

 private

contains

!----------------------------------------------------------------
!+
!  setup for a binary
!+
!----------------------------------------------------------------
subroutine set_binary(mprimary,massratio,semimajoraxis,eccentricity, &
                      accretion_radius1,accretion_radius2, &
                      xyzmh_ptmass,vxyz_ptmass,nptmass,omega_corotate,&
                      posang_ascnode,arg_peri,incl,verbose)
 use part,    only:ihacc,ihsoft
 real,    intent(in)    :: mprimary,massratio
 real,    intent(in)    :: semimajoraxis,eccentricity
 real,    intent(in)    :: accretion_radius1,accretion_radius2
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass
 real,    intent(in),  optional :: posang_ascnode,arg_peri,incl
 real,    intent(out), optional :: omega_corotate
 logical, intent(in),  optional :: verbose
 integer :: i1,i2,i
 real    :: m1,m2,mtot,dx(3),dv(3),Rochelobe,Rochelobe2,period
 real    :: x1(3),x2(3),v1(3),v2(3),omega0,cosi,sini,xangle,reducedmass,angmbin
 real    :: a,E,E_dot,P(3),Q(3),omega,big_omega,inc,ecc,tperi,term1,term2
 logical :: do_verbose,flip_x

 do_verbose = .true.
 if (present(verbose)) do_verbose = verbose

 i1 = nptmass + 1
 i2 = nptmass + 2
 nptmass = nptmass + 2

 ! masses
 m1 = mprimary
 m2 = mprimary*massratio
 mtot = m1 + m2

 Rochelobe = Rochelobe_estimate(m1,m2,semimajoraxis)
 Rochelobe2 = Rochelobe_estimate(m2,m1,semimajoraxis)
 period = sqrt(4.*pi**2*semimajoraxis**3/mtot)
 reducedmass = m1*m2/mtot
 angmbin = reducedmass*sqrt(mtot*semimajoraxis*(1. - eccentricity**2))

 if (do_verbose) then
    print "(/,2x,a)",'---------- binary parameters ----------- '
    print "(8(2x,a,g12.3,/),2x,a,g12.3)", &
        'primary mass     :',mprimary, &
        'secondary mass   :',massratio*mprimary, &
        'mass ratio       :',massratio, &
        'reduced mass     :',reducedmass, &
        'semi-major axis  :',semimajoraxis, &
        'period           :',period, &
        'eccentricity     :',eccentricity, &
        'pericentre       :',semimajoraxis*(1. - eccentricity), &
        'apocentre        :',semimajoraxis*(1. + eccentricity)
 endif

 if (accretion_radius1 >  Rochelobe2) then
    print "(/,a,/)",'*** WARNING: accretion radius of primary > Roche lobe'
 endif
 if (accretion_radius2 >  Rochelobe) then
    print "(/,a,/)",'*** WARNING: accretion radius of primary > Roche lobe'
 endif
!
!--check for stupid parameter choices
!
 if (mprimary <= 0.)      stop 'ERROR: primary mass <= 0'
 if (massratio < 0.)      stop 'ERROR: binary mass ratio < 0'
 if (semimajoraxis <= 0.) stop 'ERROR: semi-major axis <= 0'
 if (eccentricity > 1. .or. eccentricity < 0.) &
    stop 'ERROR: eccentricity must be between 0 and 1'

 dx = 0.
 dv = 0.
 if (present(posang_ascnode) .and. present(arg_peri) .and. present(incl)) then
    ! Campbell elements
    a = semimajoraxis
    ecc = eccentricity
    omega     = arg_peri*pi/180. !(arg_peri + 180.)*pi/180.
    flip_x    = .true.
    if (flip_x) omega = omega + pi/2.
    big_omega = posang_ascnode*pi/180.
    inc       = incl*pi/180.
    tperi     = 0.5*period ! time since periastron: use half period to set binary initially at apastron

    ! Solve Kepler equation for eccentric anomaly
    E = get_E(period,eccentricity,tperi)

    ! Positions in plane (Thiele-Innes elements)
    P(1) = cos(omega)*cos(big_omega) - sin(omega)*cos(inc)*sin(big_omega)
    P(2) = cos(omega)*sin(big_omega) + sin(omega)*cos(inc)*cos(big_omega)
    P(3) = sin(omega)*sin(inc)
    Q(1) = -sin(omega)*cos(big_omega) - cos(omega)*cos(inc)*sin(big_omega)
    Q(2) = -sin(omega)*sin(big_omega) + cos(omega)*cos(inc)*cos(big_omega)
    Q(3) = sin(inc)*cos(omega)

    term1 = cos(E)-eccentricity
    term2 = sqrt(1.-(eccentricity*eccentricity))*sin(E)
    E_dot = sqrt((m1 + m2)/(a**3))/(1.-eccentricity*cos(E))

    if (do_verbose) then
       print "(4(2x,a,g12.4,/),2x,a,g12.4)", &
             'Eccentric anomaly:',E, &
             'E_dot            :',E_dot, &
             'inclination (deg):',incl, &
             'arg. pericentre  :',arg_peri, &
             'angle asc. node  :',posang_ascnode
    endif

    ! Rotating everything
    ! Set the positions for the primary and the central secondary
    dx(:) = a*(term1*P(:) + term2*Q(:)) ! + xyzmh_ptmass(1,1)

    ! Set the velocities
    dv(:) = -a*sin(E)*E_dot*P(:) + a*sqrt(1.-(ecc*ecc))*cos(E)*E_dot*Q(:)

    if (flip_x) then
       ! flip x axis (because observers convention is for x axis increasing to the left)
       dx(1) = -dx(1)
       dv(1) = -dv(1)
       ! orbit in the other direction
       dv = -dv
    endif
 else
    ! set binary at apastron
    dx = (/semimajoraxis*(1. + eccentricity),0.,0./)
    dv = (/0.,sqrt(semimajoraxis*(1.-eccentricity**2)*mtot)/dx(1),0./)
 endif

 ! positions of each star so centre of mass is at zero
 x1 = -dx*m2/mtot
 x2 =  dx*m1/mtot

 ! velocities
 v1 = -dv*m2/mtot !(/0.,-m2/mtot*vmag,0./)
 v2 =  dv*m1/mtot !(/0.,m1/mtot*vmag,0./)

 omega0 = dv(2)/semimajoraxis

 ! print info about positions and velocities
 if (do_verbose) then
    print "(7(2x,a,g12.4,/),2x,a,g12.4)", &
        'angular momentum :',angmbin, &
        'mean ang. speed  :',omega0, &
        'Omega_0 (prim)   :',v1(2)/x1(1), &
        'Omega_0 (second) :',v1(2)/x1(1), &
        'R_accretion (1)  :',accretion_radius1, &
        'R_accretion (2)  :',accretion_radius2, &
        'Roche lobe  (1)  :',Rochelobe2, &
        'Roche lobe  (2)  :',Rochelobe
 endif

 if (present(omega_corotate)) then
    if (do_verbose) print "(a)",' SETTING VELOCITIES FOR COROTATING FRAME: '
    omega_corotate = omega0
    v1(2) = v1(2) - omega0*x1(1)
    v2(2) = v2(2) - omega0*x2(1)
    if (do_verbose) print "(2(2x,a,g12.4,/))", &
     'Omega_0 (primary)     :',v1(2)/x1(1), &
     'Omega_0 (secondary)   :',v2(2)/x2(1)
 endif

 ! conclude printout
 if (do_verbose) print "(2x,40('-'),/)"

!
!--positions and accretion radii
!
 xyzmh_ptmass(:,i1:i2) = 0.
 xyzmh_ptmass(1:3,i1) = x1
 xyzmh_ptmass(1:3,i2) = x2
 xyzmh_ptmass(4,i1) = m1
 xyzmh_ptmass(4,i2) = m2
 xyzmh_ptmass(ihacc,i1) = accretion_radius1
 xyzmh_ptmass(ihacc,i2) = accretion_radius2
 xyzmh_ptmass(ihsoft,i1) = 0.0
 xyzmh_ptmass(ihsoft,i2) = 0.0
!
!--velocities
!
 vxyz_ptmass(:,i1) = v1
 vxyz_ptmass(:,i2) = v2
!
! rotate if inclination is non-zero
!
 if (present(incl) .and. .not.(present(arg_peri) .and. present(posang_ascnode))) then
    xangle = incl*pi/180.
    cosi = cos(xangle)
    sini = sin(xangle)
    do i=i1,i2
       call rotate(xyzmh_ptmass(1:3,i),cosi,sini)
       call rotate(vxyz_ptmass(1:3,i),cosi,sini)
    enddo
 endif

end subroutine set_binary

pure subroutine rotate(xyz,cosi,sini)
 real, intent(inout) :: xyz(3)
 real, intent(in)    :: cosi,sini
 real :: xi,yi,zi

 xi = xyz(1)
 yi = xyz(2)
 zi = xyz(3)
 xyz(1) =  xi*cosi + zi*sini
 xyz(2) =  yi
 xyz(3) = -xi*sini + zi*cosi

end subroutine rotate

!--------------------------------
! Solve Kepler's equation for the
! Eccentric anomaly by iteration
!--------------------------------
real function get_E(period,ecc,deltat)
 real, intent(in) :: period,ecc,deltat
 real :: mu,M,E0,M0,E
 real, parameter :: tol = 1.e-10

 mu = 2.*pi/period
 M = mu*deltat ! mean anomaly
 ! first guess
 if (M > tiny(M)) then
    E = M + ecc*sin(M) + ecc**2/M*sin(2.*M)
 else
    E = M
 endif
 E0 = E + 2.*tol

 do while (abs(E - E0) > tol)
    E0 = E
    M0 = M
    M = E - ecc*sin(E)
    E = E + (M - M0)/(1. - ecc*cos(E))
 enddo

 get_E = E

end function get_E

!------------------------------------
! Compute estimate of the Roche Lobe
! Eggleton (1983) ApJ 268, 368-369
!------------------------------------
real function Rochelobe_estimate(m1,m2,sep)
 real, intent(in) :: m1,m2,sep
 real :: q,q13,q23

 q = m2/m1
 q13 = q**(1./3.)
 q23 = q13*q13
 Rochelobe_estimate = sep * 0.49*q23/(0.6*q23 + log(1. + q13))

end function Rochelobe_estimate

!---------------------------------------------
! Find first Lagrange point (L1)
! via Newton-Raphson solution of quintic
!
! INPUT: mass ratio of binary
! OUTPUT: L1 point, as distance from primary
!---------------------------------------------
real function L1_point(qinv)
 real, intent(in) :: qinv
 real :: fL, dfL, dL, L, q11

 q11 = 1./(1.+qinv)
 L = 0.5 + 0.2222222*log10(qinv)

 dL = 1.e7
 do while (abs(dL)>1.e-6)
   fL = qinv/L**2- 1./(1.-L)**2 - (1.+qinv)*L + 1.
   dfL=-2*qinv/L**3 - 2./(1.-L)**3 - (1.+qinv)
   dL = -fL/(dfL*L)
   L = L*(1.+dL)
 enddo

 L1_point = L

end function L1_point

!-------------------------------------------------------------
! Function to determine the semi-major axis given the period
!-------------------------------------------------------------
function get_a_from_period(m1,m2,period) result(a)
 real, intent(in) :: m1,m2,period
 real :: a

 a = ((m1 + m2)*(period/(2.*pi))**2)**(1./3.)

end function get_a_from_period

end module setbinary
