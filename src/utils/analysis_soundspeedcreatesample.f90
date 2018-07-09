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
!  Analysis routine to obtain soundspeed of  a particular particle
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
 character(len=20), parameter, public :: analysistype = 'SoundspeedcreateSample'
 logical :: firstcall = .true.

 public  :: do_analysis

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use io,           only: fatal
 use dim,          only: maxp,maxvxyzu
 use centreofmass, only: reset_centreofmass,get_centreofmass
 use part,         only: igas,iamtype,iphase,maxphase,rhoh
 use prompting,    only: prompt
 use units,        only: umass,udist,utime,unit_density,unit_velocity,unit_pressure
 use eos,          only: equationofstate,ieos
 use eos_helmholtz,only: xmass
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer                      :: i,iloc,j,nej,i1,i2,nsteps
 real                         :: rTmax,Tmax,Tmin,r, ri,rj,mej,rhoejm,Tejm,vejinfm,phi,bi,deltaT
 real                         :: rhoi1,rhoi2,ponrhoi1,spsoundi1,ponrhoi2,spsoundi2,rhhomax,rhomin,Tin,rhoin,deltarho
 logical                      :: iexist
 character(len=120)           :: fileprefix,fileout
 real                         :: dummyr1,dummyr2,dummyr3



 nsteps = 1000
 rhomax = 110000/unit_density
 rhomin = 1000/unit_density
 deltarho = (rhomax-rhomin)/nsteps
 Tmax   = 10000000000
 Tmin   = 100000
 deltaT = (Tmax-Tmin)/nsteps
 !
 !-- Initialization

 !--------------

 !INITIAlIZE EOS
 call init_eos(ieos,ierr)



 Tin=10000000

 fileout = 'soundspeedconstantT1e7K.dat'
 open(iunit,file=fileout,status='new')
 write(iunit,"('#',2(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'rho[g/cm3]',    &
        2,'soundspeed[cm/s]'


 do i=1,nsteps
    rhoin=rhomin+deltarho*(i-1)

    call equationofstate(ieos,ponrhoi1,spsoundi1,rhoin,0,0,0, &
                           0,Tin,xmass(:,1))

    write(iunit,'(2(1pe18.10,1x))') rhoin*unit_density,spsoundi1*unit_velocity

 enddo

 close(iunit)



 rhoin = 90000/unit_density

 fileout = 'soundspeedconstantrho9e4gcm3.dat'
 open(iunit,file=fileout,status='new')
 write(iunit,"('#',2(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'T[K]',    &
        2,'soundspeed[cm/s]'


 do i=1,nsteps
    Tin=Tmin+deltaT*(i-1)

    call equationofstate(ieos,ponrhoi1,spsoundi1,rhoin,0,0,0, &
                           0,Tin,xmass(:,1))

    write(iunit,'(2(1pe18.10,1x))') Tin,spsoundi1*unit_velocity

 enddo

 close(iunit)
 

end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
