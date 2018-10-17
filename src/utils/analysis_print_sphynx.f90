!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis print
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
 character(len=20), parameter, public :: analysistype = 'print sphynx format'
 logical :: firstcall = .true.

 public  :: do_analysis

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use io,           only: fatal
 use dim,          only: maxp,maxvxyzu
 use centreofmass, only: reset_centreofmass
 use physcon,      only: pi,planckh,mass_electron_cgs,mass_proton_cgs
 use part,         only: igas,iamtype,iphase,maxphase,rhoh
 use prompting,    only: prompt
 use units,        only: umass,udist,utime,unit_velocity,unit_density,unit_pressure,unit_energ
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

 !We need to write the values of: x,y,z,vx,vy,vz,h,u,temperature,rho? to a dat file

 !
 !
 ! Write results to file
 !
 fileout = trim(dumpfile)//'_sphynx.dat'
 open(iunit,file=fileout,status='replace')
 write(iunit,"('#',8(1x,'[',i2.2,1x,a11,']',2x))") &
        1,'mass[cgs]',&
        2,'x[cgs]',&
        3,'y[cgs]',   &
        4,'z[cgs]', &
        5,'h[cgs]',&
        6,'rho[cgs]',&
        7,'vx[cgs]',&
        8,'vy[cgs]',   &
        9,'vz[cgs]', &
        10,'T[K]'

!!!!! UNITS UNITS UNIIIIIIIIIIIIIIIIIIIIIIIIIIIITS !!!!!!!!!!!!!!!!!
 do j=1,npart
      write(iunit,'(10(1pe18.10,1x))') particlemass*umass,xyzh(1,j)*udist,xyzh(2,j)*udist,xyzh(3,j)*udist,xyzh(4,j)*udist,rhoh(xyzh(4,j),particlemass)*unit_density,vxyzu(1,j)*unit_velocity,vxyzu(2,j)*unit_velocity,vxyzu(3,j)*unit_velocity,vxyzu(5,j)
 enddo
 close(iunit)

end subroutine do_analysis

end module
