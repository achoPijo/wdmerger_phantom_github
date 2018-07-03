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
!  Analysis routine that computes the fallback luminosities based on
!  Rosswog 2007
!
!  REFERENCES: Rosswog 2007 "Fallback accretion in the aftermath of a
!  compact binary merger"
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
 use units,        only: umass,udist,utime,unit_density,unit_velocity,unit_pressure,unit_energ
 use vectorutils,  only: cross_product3D
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:) !due to reset center of mass
 real,             intent(in)    :: particlemass,time

 integer                      :: lengthtimearray = 3600000 !Miliseconds resolution
 integer                      :: i,iloc,j,nej,
 real                         :: rTmax,Tmax,r, ri,rj,mej,rhoejm,Tejm,vejinfm,phi,bi,Ji(3),ji,ei,Ei,mtot,ai,ti
 real                         :: timearray(lengthtimearray),deltaenergyarray(lengthtimearray),deltat
 logical                      :: iexist
 character(len=120)           :: fileprefix,fileout





 !
 !-- Initialization
 !

 deltat =0.001 ! in seconds
 do i=1,lengthtimearray
    timearray(i)        = deltat/2+(i-1)*deltat
    deltaenergyarray(i) = 0.0
 enddo 


 call reset_centreofmass(npart,xyzh(:,:),vxyzu(:,:))
 
 !   Rosswog 2007
 !
 !-- Fallback luminosities analysis
 !
 Rdisc = 0.1
 nej = 0
 mej = 0.
 rhoejm  = 0.
 Tejm    = 0.
 vejinfm = 0.
 mtot = npart*particlemass
 do i=1,npart
    ri = sqrt(xyzh(1,i)**2+xyzh(2,i)**2+xyzh(3,i)**2)
    bi = -mtot/ri + 0.5*(sqrt(vxyzu(1,i)**2 + vxyzu(1,i)**2 + vxyzu(1,i)**2)**2) !specific orbital energy
    Ei = bi * particlemass
    call cross_product3D(xyzh(1:3,i),vxyzu(1:3,i),Ji)
    ji = particlemass*sqrt(Ji(1)**2 + Ji(2)**2 + Ji(3)**2)
    ei = sqrt(1+2*Ei*(Ji**2)/((particlemass**3)*(mtot**2)))
    ai = -mtot*particlemass/(2*Ei)
    rmaxi = ai*(1+ei)
    rmini = ai*(1-ei)
    if (bi < 0. .and. ri > Rdisc .and. rmini < Rdisc) then
       A = 2*Ei/particlemass
       B = 2*mtot
       C = (ji**2)/(particlemass**2)
       D = 4*A*C-B**2
       if (vxyzu(1,i)*xyzh(1,i) + vxyzu(2,i)*xyzh(2,i) +vxyzu(3,i)*xyzh(3,i) > 0.) then
          call integralRosswog(Iri,A,B,C,D,ri)
          call integralRosswog(Irmax,A,B,C,D,rmaxi)
          call integralRosswog(IRdisc,A,B,C,D,Rdisc)
          ti = (Irmaxi - Iri) + (IRdisc - Irmaxi) 
       else
          call integralRosswog(Iri,A,B,C,D,ri)
          call integralRosswog(IRdisc,A,B,C,D,Rdisc)
          ti = IRdisc - Iri
       endif
       deltae = (Ei - mtot*particlemass/Rdisc)*unit_energ
    do j=1,lengthtimearray
       if (ti * utime < j*deltat) then
          deltaenergyarray(j) = deltaenergyarray(j) + deltae/deltat
          exit
       endif
    enddo
       
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

 fileout = trim(fileprefix)//'fallback.dat'

 open(iunit,file=fileout,status='new')
 write(iunit,"('#',5(1x,'[',i2.2,1x,a11,']',2x))") &
     1,'time[s]', &
     2,'fallback[erg/s]'
 
 do i=1,lengthtimearray

    write(iunit,'(2(1pe18.10,1x))') timearray(i),deltaenergyarray(i)

 enddo

close(iunit)






end subroutine do_analysis

subroutine integralRosswog(I,a,b,c,d,r)
   real, intent(out) :: I
   real, intent(in)  :: a,b,c,d,r

   I = sqrt(a*(r**2)+b*r+c)/a + B*asin((2*a*r+b)/sqrt(-d))/(2*a*sqrt(-a))
end subroutine
!-----------------------------------------------------------------------
!
end module
