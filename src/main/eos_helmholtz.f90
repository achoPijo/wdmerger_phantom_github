!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: eos_mesa
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: 
!
!  $Id: c2ae599291c57d71e81c2c48a0f1bef5956c6b4b $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: helmholtz_microphysiss
!+
!-----------------------------------------------------------------
module eos_helmholtz

 
 use helmholtz_microphysics

 implicit none

 real(8), parameter :: tmaxhelmeos=5.0d9, tminhelmeos=1.0d3
 real(8), parameter :: cgsrhominhelmeos=1.0d0*5899.39, cgsrhomaxhelmeos=1.0d8*5899.39 !in cgs units

 public  :: init_eos_helmholtz, get_eos_press_sound_cv_dPdT_helmholtz, helmholtz_energytemperature_switch
 public  :: tmaxhelmeos, tminhelmeos, cgsrhomaxhelmeos, cgsrhominhelmeos
 private

contains

!----------------------------------------------------------------
!+
!  subroutine initialises the helmholtz eos tables
!+
!----------------------------------------------------------------
subroutine init_eos_helmholtz(ierr)
 
 integer, intent(out) :: ierr
 ierr = 0

 call read_eos_helmholtz(ierr)

end subroutine init_eos_helmholtz

!----------------------------------------------------------------
!+
!  subroutine returns pressure and sound speed as a function of
!  temperature and density
!+
!----------------------------------------------------------------
subroutine get_eos_press_sound_cv_dPdT_helmholtz(temp,den,pres,sound,cv,dPdT)!,abar,zbar) !REVISE ADD other variables obtained by helmeos needed by hydro_rs like dPdt internal energy etc as optional inputs
 use io,         only:fatal
 real, intent(in)    :: temp, den !, abar, zbar
 real, intent(out)   :: pres, sound, cv, dPdT
 integer             :: ierr
 real                :: abar, zbar
 character(len=40), parameter  :: label = 'get_eos_press_sound_cv_dPdT_helmholtz'

 ierr = 0
 abar = 1./(17./240.) !VALUES FOR 40% C 60% O !REVISE when xss is included
 zbar = 0.5/(17./240.)

 call helmeos(temp,den,abar,zbar,ierr,pres,sound,cv_opt=cv,dpresdt_opt=dPdT)

 !print *, sound
 !print *, abar


 if (ierr /= 0 ) call fatal(label,'Temperature or density values outside helmholtz free energy tables range')

end subroutine get_eos_press_sound_cv_dPdT_helmholtz



!----------------------------------------------------------------
!+
!  THIS SUBROUTINE SELECTS BETWEEN INTERNAL ENERGY AND TEMPERATURE
!  EVOLUTION DEPENDING ON HOW BIG ARE TEMPERATURE JUMPS DUE TO 
!  INTERNAL ENERGY CHANGES
!+
!----------------------------------------------------------------
 subroutine helmholtz_energytemperature_switch(temp,ener,den,relflag)!abar,zbar) !REVISE 
!========================================================================

!=========================================================================
!
!--Load modules
!
 use units,          only: umass, unit_density,unit_ergg
! use part,           only: xss, aion, zion, nel
! use burn only:xss, aion, zion !REVISE for nuclear burning
                                 
!
!--Force to declare EVERYTHING
!
 implicit none

!--I/O variables
!
 real,             intent(inout):: temp, ener
 real,             intent(in)   :: den!, abar, zbar !REVISE abar,zbar set manually for now until burn is implemented
 logical,          intent(in)   :: relflag
!
!--Local variables
!
 integer, parameter :: max_newton=50
 real(8), parameter :: eos_tol=1.0d-8
 real(8)  :: abar, zbar, ewant, temp_iter, ener_iter, tnew, errorp, rel, denerdt!, asum,zsum
 integer  :: k, newton, eosflag, ierr
 real     :: rho
!
!--Set initial data for the EOS
!
 if (temp > tmaxhelmeos) temp = tmaxhelmeos
 if (temp < tminhelmeos) temp = tminhelmeos
!
 ewant = ener*unit_ergg ! EOS works in cgs units!!!
 temp_iter  = temp
 rho   = den*unit_density    ! EOS works in cgs units!!!
!
!------------------------------ Timmes & Swesty apj 1999 page 1 (501)
! asum = 0.0
! zsum = 0.0
! do k=1,nel-1
!    asum = asum + xss(k,p)/aion(k)
!    zsum = zsum + xss(k,p)*zion(k)/aion(k)
! enddo
! abar = 1.0/asum
! zbar = zsum/asum
! 
 abar = 1/0.25 !VALUES FOR 100% He 
 zbar = 0.5/0.25
!-------------------------------
!--Loop over particles doing the Newton-Raphson iteration
!
   newton = 0
   errorp = 2.*eos_tol
   do while ((errorp > eos_tol).AND.(newton < max_newton))
      newton = newton + 1
!
!--Call EOS
!
      call helmeos(temp_iter, rho, abar, zbar, ierr, ener_opt=ener_iter, denerdt_opt=denerdt) !REVISE implement error warning mechanism
      
!
!--New temperature
!
      tnew = temp_iter - (ener_iter-ewant)/denerdt
!
!--Do no allow temperature to change more than an order of magnitude per
!  iteration
!
      if (tnew > 10.d0*temp_iter) tnew = 10.0d0*temp_iter
      if (tnew < 0.1d0*temp_iter) tnew = 0.1d0*temp_iter
!
!--Calculate the error
!
      errorp = ABS(tnew-temp_iter)/temp_iter
! 
!--Store new temperature
!
      temp_iter = tnew
!
!--If freezing, keep temperature inside EOS tables
!
      if (temp_iter < tminhelmeos) then
          temp_iter = tminhelmeos
          errorp    = 0.1d0*eos_tol
      endif
!
!--If too hot, keep temperature inside EOS tables
!
      if (temp_iter > tmaxhelmeos) then
          temp_iter = tmaxhelmeos
          errorp    = 0.1d0*eos_tol
      endif
   enddo    ! End Newton-Raphson loop
!
!--If the Newton-Rapshon fails to find a valid temperature, keep it 
!  constant. Check also if temperature and the temperature predicted
!  from the internal energy differe more than a 5%
 rel = dabs(temp_iter-temp)/temp
 if (relflag) then
    eosflag = 2
 else
    if ((newton >= max_newton).or.(rel > 0.05)) then
       eosflag = 2
    else
       eosflag = 1
    end if
 endif
!
!--If eosflag=1 we store temp_iter, whereas for eosflag=2 we store 
!  etot_row
!

 if (eosflag == 2) then  ! Temperature as input
     call helmeos(temp, rho, abar, zbar, ierr, ener_opt=ener_iter)

     ener = ener_iter/unit_ergg
 else                    ! Internal energy as input
     temp = temp_iter
 end if

 end subroutine helmholtz_energytemperature_switch


!----------------------------------------------------------------
!+
!  ADD ANY ADDITIONAL SUBROUTINE
!+
!----------------------------------------------------------------


end module eos_helmholtz
