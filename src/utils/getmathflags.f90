!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: getmathflags
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id: 403415f93a767b0cb7d4214dd524a0b4b9d20b24 $
!
!  USAGE: getmathflags [no arguments]
!
!  DEPENDENCIES: testmath
!+
!--------------------------------------------------------------------------
program getmathflags
 use testmath, only:test_math
 implicit none
 integer :: ntests,npass
 logical :: usefsqrt,usefinvsqrt

 ntests = 0
 npass  = 0
 call test_math(ntests,npass,usefsqrt,usefinvsqrt)
 if (usefinvsqrt) then
    print "(a)",'yes'
 endif

end program getmathflags
