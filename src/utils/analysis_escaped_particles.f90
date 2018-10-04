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
 character(len=20), parameter, public :: analysistype = 'Maximum Temperature'
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
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer                      :: i,iloc,j,nej
 real                         :: rTmax,Tmax,r, ri,rj,mej,rhoejm,Tejm,vejinfm,phi,bi
 logical                      :: iexist
 character(len=120)           :: fileprefix,fileout




 
 !
 !-- Mass Ejecta Analysis

 !--------------

 nej = 0
 mej = 0.
 rhoejm  = 0.
 Tejm    = 0.
 vejinfm = 0.
 do i=1,npart
    ri = sqrt(xyzh(1,i)**2+xyzh(2,i)**2+xyzh(3,i)**2)
    phi= 0.
    do j=1,npart
       if (j /= i) then 
          rj = sqrt(xyzh(1,j)**2+xyzh(2,j)**2+xyzh(3,j)**2)
          phi= phi - particlemass/(abs(ri-rj)) !remember units chosen so that G=1
       endif
    enddo
    bi = phi + 0.5*(sqrt(vxyzu(1,i)**2 + vxyzu(1,i)**2 + vxyzu(1,i)**2)**2) !specific orbital energy
    if (bi > 0.) then
       nej = nej + 1
       mej = mej + particlemass
       Tejm      = Tejm + vxyzu(5,i)
       rhoejm    = rhoejm + rhoh(xyzh(4,i),particlemass)
       vejinfm   = vejinfm + sqrt(2*bi)
    endif
 enddo

 if (nej/=0) then
     mej     = mej
     Tejm    = Tejm/nej
     rhoejm  = rhoejm/nej  * unit_density
     vejinfm = vejinfm/nej * unit_velocity
 endif

 !
 !-- Dump information to file
 !
 

 iloc = index(dumpfile,'_')
 if (iloc > 1) then
    fileprefix = trim(dumpfile(1:iloc-1))
 else
    fileprefix = trim(dumpfile)
 endif

 fileout = trim(fileprefix)//'ejecta.dat'
 inquire(file=trim(fileout),exist=iexist)
 if (iexist) then
    open(iunit,file=fileout,status='old',position='append')
    
    write(iunit,'(5(1pe18.10,1x))') time,mej,Tejm,rhoejm,vejinfm

    close(iunit)
 else
    open(iunit,file=fileout,status='new')
    write(iunit,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'time[code units]',    &
        2,'mej[Msun]', &
        3,'Tejm[K]', &
        4,'rhoejm[g/cm3]', &
        5,'vejinfm[cm/s]'
    
    write(iunit,'(5(1pe18.10,1x))') time,mej,Tejm,rhoejm,vejinfm


    close(iunit)
 endif
end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
