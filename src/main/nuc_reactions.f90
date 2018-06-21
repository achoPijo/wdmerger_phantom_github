module nuc_reactions
   implicit none

   public :: nuclear_burning, nuc_burn
   public :: write_options_nuc_burning, read_options_nuc_burning
   public :: init_nuc_burning

   private

   integer :: nuc_burn = 0
   integer :: burn_opt = 2
   real    :: enuctot

   contains

!
! TO DO
! 
! - Create aion,zion, xss, cvs, enuc, luminuc, dtnuc in module part
! - nel = 16 need to define somewhere part? 
! - integrate ZHELI.10b into a subroutine?
! - need cv of each particle, is it worth to store it? maybe if burn is called from force.F90 that value is already computed
! - need to call RPARAM or define the value of those variables anyhow: create them in part e inicializar en setup? buena opción
! - dtmp_min equivalence in PHANTOM?? - dt DONE


      SUBROUTINE nuclear_burning(xyzh,vxyzu,fxyzu,npart,dt)
!=========================================================================
! This subroutine calculates nuclear burning
!
! Last revision: 15/March/2015
!=========================================================================
!
!--Load modules
!
      !USE mod_essentials
      USE units,         ONLY : umass, unit_density, utime, unit_energ, unit_ergg
      USE part,          ONLY : rhoh, massoftype, igas
      use eos_helmholtz, only : xmass,speciesmax
      use eos,           only : ieos,equationofstate
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Helmholtz EOS definitions
!
      !INCLUDE 'vector_eos.dek'

      real, intent(in)    :: xyzh(:,:)
      real, intent(inout) :: vxyzu(:,:)
      real, intent(inout) :: fxyzu(:,:)
      integer, intent(in) :: npart
      real, intent(in)    :: dt
!
!--Local variables
!     
      REAL(8), DIMENSION(speciesmax-1) :: xss2                 
      REAL(8) :: rhopcgs, rhop, temp, cvp, enucp, luminucp, sumAE2, sumdt
      REAL(8), PARAMETER :: rhotiny=5.0d0
      INTEGER :: i, p, JK, k, m, iread, iteration, itermax
      INTEGER, PARAMETER :: iter_max=1000, NIS=speciesmax-2, NSP=NIS+1, NRE=29,&
                            ngrid=60
!
!--Nuclear network variables
!
      REAL(8), DIMENSION(ngrid,NRE,11) :: vgrid 
      REAL(8), DIMENSION(ngrid,11) :: aNegrid
      REAL(8), DIMENSION(ngrid) :: tgrid
      REAL(8), DIMENSION(NRE)   :: V2, flag, ndata, n1, n2, n3, n4, a1, &
                                   a2, a3, a4, Qrad, Qnu, QVAL, AE2
      REAL(8), DIMENSION(NSP,1) :: XXTOT
      REAL(8), DIMENSION(NSP)   :: AN, ZN, BE, YYB, XXXT
      REAL(8) :: DTMOLDP, SUMYY, DTMNEWP, BE0
      INTEGER, DIMENSION(NRE)   :: K1, K2, K3, K4, K5, K6, K7, K8
      INTEGER :: ireac, kgrid, NSNUC, ICO
      CHARACTER(5),  DIMENSION(NSP) :: ONC
      CHARACTER(6),  DIMENSION(NRE) :: z1, z2, z3, z4
      CHARACTER(37), DIMENSION(NRE) :: reaction
!      integer                       :: burn_opt
      real                          :: dummyt,dummyr1,dummyr2,enuctot
!
      COMMON /cvit/V2,tgrid,aNegrid
      COMMON /cvgrid/vgrid,flag,ndata
      COMMON /network/ireac,kgrid,n1,n2,n3,n4,z1,z2,z3,z4,a1,a2,a3,a4,  &
                      reaction
      COMMON /Q/Qrad,Qnu
      COMMON /CM/K1,K2,K3,K4,K5,K6,K7,K8
      COMMON /CNAME/ONC
      COMMON /CNETW/AN,ZN,QVAL,BE,BE0
      COMMON /ABUND/XXTOT
      COMMON /CSTEP/NSNUC,ICO
      COMMON /CSNC2/YYB
      COMMON /NUCLEOS/XXXT
      COMMON /CODER/AE2
!
!--Initializations
!
      iread   = 0
      itermax = 0 
      ICO     = 0 !This is an integer number that should be related to the timestep number or dumpfile number

      print *, "----------------------------------"
      print *, "Cumulative nuclear energy released"
      print *, enuctot*unit_ergg

!
!!$OMP PARALLEL DEFAULT(none) shared(vxyzut,cvs,rho,xss,iread)         &       
!!$OMP shared(dt,enuc,luminuc,aion,zion,tnow) private(m,rhop,temp,cvp) &
!!$OMP private(enucp,luminucp,iteration,sumdt,DTMNEWP,DTMOLDP,sumAE2)        &
!!$OMP private(abar,zbar,temp_row,den_row,cv_row,abar_row,zbar_row,jlo_eos)  &
!!$OMP private(jhi_eos,xss2,XXTOT,XXXT,SUMYY,i,p) reduction(MAX:itermax)
!!$OMP DO SCHEDULE(runtime)
      DO 40 p=1,npart
!
!--Density, temperature and composition initialisation (to cgs units)
!  Skip hydrogen: xss(1) == He ... xss(speciesmax-1) == photons
!
         rhop    = rhoh(xyzh(4,p),massoftype(igas))
         rhopcgs = rhop*unit_density
         temp = vxyzu(5,p)
         call equationofstate(ieos,dummyr1,dummyr2,rhop,xyzh(1,p),xyzh(2,p),xyzh(3,p), &
                           vxyzu(4,p),vxyzu(5,p),xmass(:,p),cvp)
         cvp = cvp * unit_ergg
         DO i=1,speciesmax-1
            xss2(i) = xmass(i+1,p)
         ENDDO 
!
!--Iteration over the nuclear time-steps
!
         iteration = 0
         sumdt     = 0.0
         enucp     = 0.0
         luminucp  = 0.0
         dummyt    = 0.0
         DO WHILE(sumdt.LT.dt*utime)
!
!--Activate nuclear burning only for T>1.0d7
!
            IF (temp.LT.1.0d7) GOTO 541
!
!--Low density particle burning
!
            IF ((iteration.GT.iter_max).AND.(rhopcgs).LT.rhotiny) THEN 
               PRINT*,'Careful, incinerated particle !!!',p,temp
               GOTO 541
            END IF
!
!--Store chemical abundances composition
!
            DO 50 k=1,NSP
               XXTOT(k,1) = xss2(k)
50          ENDDO

!
!--Avoid excessively small time steps
!
            IF (sumdt.EQ.0.0d0) THEN
               DTMNEWP = dt*utime
            ELSE
               DTMNEWP = MIN(DTMNEWP,dt*utime-sumdt)
            ENDIF
!
!--Abundances calculation
!
            DTMOLDP = 0.0d0
            SUMYY   = 0.0d0
            CALL SNUC(dummyt,DTMNEWP,rhopcgs,temp,SUMYY,DTMOLDP,1,iread)
            IF (iread.EQ.0) iread=1
!
!--If nuclear is bigger than SPH time-step, adopt SPH time-step 
!  (shouldnt happen since dtmold <= dtmnew)
!
            IF (DTMOLDP.GT.dt*utime-sumdt) THEN   
               DTMOLDP = dt*utime-sumdt
            ENDIF
            sumdt = sumdt + DTMOLDP
!
!--Abundances change
!
            DO 60 k=1,NSP
               xss2(k) = XXXT(k)
60          ENDDO
!
!--Total nuclear luminosity released. AE2 is in erg g-1 s-1, so 
!  change to code units
!
            sumAE2 = sum(AE2)*DTMOLDP*(umass/unit_energ)  
            enucp  = enucp + sumAE2

            ! Three different options
            ! 1) T is not evolved during the nuclear burning
            ! 2) T is only updated due to photodesintegration,
            ! and only positiove energy contributions are updated to 
            ! luminucp
            ! 3) T is updated every nuclear timestep
            select case (burn_opt)
               case (1)
                  !
                  !--Temperature does not change
                  !
                  luminucp = luminucp + sumAE2/cvp

               case (2)   
                  !
                  !--Change temperature only if photodesintegration is present
                  !  Recalculate specific heats if temperature has changed
                  !
                  if (sumAE2.LT.0.0d0) then
                     temp = temp + sumAE2/cvp
                     call equationofstate(ieos,dummyr1,dummyr2,rhop,xyzh(1,p),xyzh(2,p),xyzh(3,p), &
                                 vxyzu(4,p),temp,xmass(:,p),cvp)
                     cvp = cvp * unit_ergg
                  else
                     luminucp = luminucp + sumAE2/cvp
                  endif

               case (3)
                  !
                  !--Always update temperature
                  ! 
                  temp = temp + sumAE2/cvp
                  call equationofstate(ieos,dummyr1,dummyr2,rhop,xyzh(1,p),xyzh(2,p),xyzh(3,p), &
                              vxyzu(4,p),temp,xmass(:,p),cvp)
                  cvp = cvp * unit_ergg

               case default
                  !
                  !--Temperature does not change
                  !
                  luminucp = luminucp + sumAE2/cvp

            end select
!
            iteration = iteration + 1 
         ENDDO ! finish burn iteration
541      CONTINUE
!
!--Statistics
!
         itermax=max(iteration,itermax)
!         
!--Store changes
!
         vxyzu(5,p) = temp
         fxyzu(4,p) = fxyzu(4,p) + enucp
         fxyzu(5,p) = fxyzu(5,p) + luminucp
         DO k=1,speciesmax-1
            xmass(k+1,p) = xss2(k)
         END DO
!
!--Accumulate total energy released
!
         enuctot = enuctot + enucp

45       CONTINUE
40    ENDDO ! finish particle iteration
!!$OMP END DO 
!!$OMP END PARALLEL
!
      print *, "----------------------------------"
      print *, "Cumulative nuclear energy released"
      print *, enuctot*unit_ergg
      print *, "----------------------------------"      
!
      !IF (rank.EQ.MASTER) PRINT*,'burn: max. num. iteraciones:',itermax
!
      END SUBROUTINE nuclear_burning

!-----------------------------------------------------------------------
!+
!  initialise equation of state (read tables etc.)
!+
!-----------------------------------------------------------------------
subroutine init_nuc_burning(ieos,ierr)
 use io,            only:error
 integer, intent(in)  :: ieos
 integer, intent(out) :: ierr
 !
 ierr = 0
 !
 if (ieos/=15) then
    ierr = 1
    return
 endif

 call RPARAM
 call RNETWORK

 enuctot = 0.

 
 !TODO ATTENTION, NEED TO FIND A WAY TO DETECT PROBLEM IN INITIALIZING 

end subroutine init_nuc_burning

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_nuc_burning(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling nuclear burning'
 call write_inopt(nuc_burn,'nuc_burn','0=nuclear burning off, 1=nuclear burning on',iunit)
 call write_inopt(burn_opt,'burn_opt','1=T constant through nuc_timestep, 2=T only updated due to photodisintegration, 3=Temperature updated in nuc_timestep',iunit)

end subroutine write_options_nuc_burning

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_nuc_burning(name,valstring,imatch,igotall,ierr)
 use io,         only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot  = 0
 character(len=30), parameter  :: label = 'read_options_nuc_burning'

 imatch  = .true.
 select case(trim(name))
 case('nuc_burn')
    read(valstring,*,iostat=ierr) nuc_burn
    ngot = ngot + 1
    if (nuc_burn < 0 .or. nuc_burn > 1) call fatal(label,'value of nuc_burn must be 0 or 1')
 case('burn_opt')
    read(valstring,*,iostat=ierr) burn_opt
    ngot = ngot + 1
    if (burn_opt < 1 .or. burn_opt > 3) call fatal(label,'value of burn_opt must be between 1 and 3')
  
 case default
    imatch = .false.
 end select

 !--make sure we have got all compulsory options (otherwise, rewrite input file)
    igotall = (ngot == 2)

end subroutine read_options_nuc_burning


end module nuc_reactions
