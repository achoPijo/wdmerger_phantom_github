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
 use setbinary,    only: L1_point
 use units,        only: umass,udist,utime,unit_density,unit_pressure
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer                      :: i,iloc
 real                         :: rTmax,Tmax,r
 logical                      :: iexist
 character(len=120)           :: fileprefix,fileout





 !
 !-- Initialization

 !--------------

 rTmax = 0.
 Tmax  = 0.
 do i=1,npart
    r = sqrt(xyzh(1,i)**2+xyzh(2,i)**2+xyzh(3,i)**2)
    Tmax = max(Tmax,vxyzu(5,i))
    if (Tmax == vxyzu(5,i)) then 
       rTmax = r
    endif
 enddo
 

 !
 !-- Dump information to file
 !
 iloc = index(dumpfile,'_')
 if (iloc > 1) then
    fileprefix = trim(dumpfile(1:iloc-1))
 else
    fileprefix = trim(dumpfile)
 endif
 
 fileout = trim(fileprefix)//'maxTemperature.dat'
 inquire(file=trim(fileout),exist=iexist)
 if (iexist) then
    open(iunit,file=fileout,status='old',position='append')
    
    write(iunit,'(3(1pe18.10,1x))') time,Tmax,rTmax

    close(iunit)
 else
    open(iunit,file=fileout,status='new')
    write(iunit,"('#',3(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'r',    &
        2,'Tmax', &
        3,'rTmax'
    
    write(iunit,'(3(1pe18.10,1x))') time,Tmax,rTmax

     close(iunit)
 endif

end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
