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
!  Analysis routine to obtain temperature density angular velocity and keplerian velocity 
!  from a given dumpfile
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
 character(len=20), parameter, public :: analysistype = 'Density,Temperature and other radial profiles and relevant data'
 logical :: firstcall = .true.

 public  :: do_analysis

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use io,           only: fatal
 use dim,          only: maxp,maxvxyzu
 use centreofmass, only: reset_centreofmass,get_centreofmass
 use part,         only: igas,iamtype,iphase,maxphase,rhoh,massoftype
 use prompting,    only: prompt
! use setbinary,    only: L1_point
 use units,        only: umass,udist,utime,unit_density,unit_pressure
 use vectorutils,  only: cross_product3D
! use eos_helmholtz,only: xmass,speciesmax
 use eos,           only: ieos, equationofstate, init_eos,finish_eos,relflag
 use eos_helmholtz, only: tmaxhelmeos,tminhelmeos,xmass,speciesmax

 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer, parameter  :: nrpoints = 2000
 integer             :: i,j,ieos,ierr
 real                :: rmax,rTmax,Tmax,dr,mtot,r,ponrhoi,densi,dummyspsoundi
 real                :: press(npart),rtab(npart)
 real                :: vec(3),xcom(3),vcom(3),xdens(3)
 character(len=200)  :: fileout

 !
 !-- Initialization
 !
 press(:)          = 0.
 rtab(:)           = 0.

 ieos = 15
 ierr = 0
 call init_eos(ieos,ierr)

 !
 !-- 
 !

 do i=1,npart
    rtab(i) = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2)
    densi = rhoh(xyzh(4,i),massoftype(iphase(i)))
    call equationofstate(ieos,ponrhoi,dummyspsoundi,densi,xyzh(1,i),xyzh(2,i),xyzh(3,i), &
                            tempi=vxyzu(5,i),xmassi=xmass(:,i))
    press(i)=ponrhoi*densi
 enddo

 !
 !-- Write informtaion to file
 !

 fileout = trim(dumpfile)//'_press.dat'
 open(iunit,file=fileout,status='replace')
 write(iunit,"('#',2(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'r',  &
        2,'Pressure'

 do j=1,nrpoints
      write(iunit,'(2(1pe18.10,1x))') rtab(j),press(j)
 enddo
 close(iunit)


 !
 !-- INCLUDE A FILE OUTPUT WITH TMAX omegamax locations
 ! 


end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
