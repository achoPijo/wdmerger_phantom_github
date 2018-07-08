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
 character(len=20), parameter, public :: analysistype = 'Soundspeed'
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

 integer                      :: i,iloc,j,nej
 real                         :: rTmax,Tmax,r, ri,rj,mej,rhoejm,Tejm,vejinfm,phi,bi
 real                         :: rhoi1,rhoi2,ponrhoi1,spsoundi1,ponrhoi2,spsoundi2
 logical                      :: iexist
 character(len=120)           :: fileprefix,fileout
 real                         :: dummyr1,dummyr2,dummyr3





 !
 !-- Initialization

 !--------------

! INITIAlIZE EOS
 rTmax = 0.
 Tmax  = 0.
 
 i1 = 7695
 rhoi1 = rhoh(xyzh(4,i1),particlemass)
 i2 = 38828
 rhoi2 = rhoh(xyzh(4,i2),particlemass)

 call equationofstate(ieos,ponrhoi1,spsoundi1,rhoi1,xyzh(1,i1),xyzh(2,i1),xyzh(3,i1), &
                           vxyzu(4,i1),vxyzu(5,i1),xmass(:,i1))
    
 call equationofstate(ieos,ponrhoi2,spsoundi2,rhoi2,xyzh(1,i2),xyzh(2,i2),xyzh(3,i2), &
                           vxyzu(4,i2),vxyzu(5,i2),xmass(:,i2))


 !
 !-- Dump information to file
 !
 iloc = index(dumpfile,'_')
 if (iloc > 1) then
    fileprefix = trim(dumpfile(1:iloc-1))
 else
    fileprefix = trim(dumpfile)
 endif
 
 fileout = trim(fileprefix)//'soundspeed.dat'
 inquire(file=trim(fileout),exist=iexist)
 if (iexist) then
    open(iunit,file=fileout,status='old',position='append')
    
    write(iunit,'(7(1pe18.10,1x))') time,rhoi1*unit_density,vxyzu(5,i1),spsoundi1*unit_velocity,rhoi2*unit_density,vxyzu(5,i2),spsoundi2*unit_velocity

    close(iunit)
 else
    open(iunit,file=fileout,status='new')
    write(iunit,"('#',7(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'time',    &
        2,'rhoi1',    &
        3,'Ti1',    &
        4,'spsoundi1[cm/s]',    &
        5,'rhoi2',    &
        6,'Ti2', &
        7,'spsoundi2[cm/s]'
    
    write(iunit,'(7(1pe18.10,1x))') time,rhoi1*unit_density,vxyzu(5,i1),spsoundi1*unit_velocity,rhoi2*unit_density,vxyzu(5,i2),spsoundi2*unit_velocity

     close(iunit)
 endif


end subroutine do_analysis

!-----------------------------------------------------------------------
!
end module
