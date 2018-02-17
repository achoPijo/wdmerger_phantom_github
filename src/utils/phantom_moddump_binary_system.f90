!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantommoddump_binary_system
!
!  DESCRIPTION: This program is a simple utility for modifying a dump file
!
!  REFERENCES: None
!
!  OWNER: Jose Miguel Blanco
!
!  $Id: -
!
!  USAGE: moddump dumpfilein1 dumpfilein2 dumpfileout [time] [outformat]
!
!  DEPENDENCIES: checksetup, dim, eos, io, moddump, options, part,
!    prompting, readwrite_dumps, readwrite_infile, setBfield, setup_params
!+
!--------------------------------------------------------------------------
program phantommoddump_binary_system !CHNGCODE
 use dim,             only:maxp,tagline
 use eos,             only:polyk
 use part,            only:xyzh,hfact,massoftype,vxyzu,npart,npartoftype, &
                           Bevol,Bextx,Bexty,Bextz,mhd,maxvecp
 use io,              only:set_io_unit_numbers,iprint,idisk1,warning,fatal,iwritein,id,master
 use readwrite_dumps, only:read_dump,write_fulldump,is_not_mhd
 use setBfield,       only:set_Bfield
 use moddump,         only:modify_dump
 use readwrite_infile,only:write_infile,read_infile
 use options,         only:set_default_options
 use setup_params,    only:ihavesetupB
 use prompting,       only:prompt
 use checksetup,      only:check_setup
 implicit none
 integer :: nargs
 character(len=120) :: dumpfilein1,dumpfilein2,dumpfileout                     !CHNGCODE
 character(len=10) :: string
 real    :: time,timeout,xyzh1(4,npart),xyzh2(4,npart)                         !CHNGCODE
 integer :: i,ierr,nerr,nwarn,iloc,nstar1,nstar2                               !CHNGCODE
 logical :: idumpsphNG,iexist,ians
 integer, parameter          :: lenprefix = 120
 character(len=lenprefix)    :: fileprefix
 character(len=lenprefix+10) :: dumpfile,infile,evfile,logfile

 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 3 .or. nargs > 5) then  !CHNGCODE
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: moddump dumpfilein1 dumpfilein2 dumpfileout [time] [outformat]'
    stop
 endif
 call get_command_argument(1,dumpfilein1)
 call get_command_argument(2,dumpfilein2)
 call get_command_argument(3,dumpfileout)

 print "(/,a,/)",' Phantom moddump: pimp my dumpfiles'
!
!--look for an existing input file with name corresponding to the INPUT dump file
!  read this if it exists
!
 call set_default_options
 iloc = index(dumpfilein1,'_0')
 if (iloc > 1) then
    fileprefix = trim(dumpfilein1(1:iloc-1))
 else
    fileprefix = trim(dumpfilein1)
 endif
 infile = trim(fileprefix)//'.in'
 inquire(file=trim(infile),exist=iexist)
 if (iexist) then
    print "(/,2a,/)",' Reading default values from ', trim(infile)
    call read_infile(infile,logfile,evfile,dumpfile)
 endif
!
!--look for an existing input file with name corresponding to the OUTPUT dump file
!  read this if it exists; this will overwrite the values from the input .in file
!
 iloc = index(dumpfileout,'_0')
 if (iloc > 1) then
    fileprefix = trim(dumpfileout(1:iloc-1))
 else
    fileprefix = trim(dumpfileout)
    dumpfileout = trim(dumpfileout)//'_00000'
 endif
 infile = trim(fileprefix)//'.in'
 inquire(file=trim(infile),exist=iexist)
 if (iexist) then
    print "(2a,/)",' Reading revised default values from ', trim(infile)
    call read_infile(infile,logfile,evfile,dumpfile)
 endif
!
!--reset logfile name
!
 logfile = trim(fileprefix)//'01.log'
 evfile  = trim(fileprefix)//'01.ev'
 if (mhd) then
    ihavesetupB = .true.
 else
    ihavesetupB = .false.
 endif
!
!--read particle setup from dumpfilein1
!
 call read_dump(trim(dumpfilein1),time,hfact,idisk1,iprint,0,1,ierr)
 if (mhd .and. ierr==is_not_mhd) then
    ihavesetupB = .false.
 elseif (ierr /= 0) then
    stop 'error reading dumpfile'
 endif
                                                                               !CHNGCODE
 nstar1 = npart                                                                !CHNGCODE
 xyzh1  = xyzh                                                                 !CHNGCODE
!
!--read particle setup from dumpfilein2                                        !CHNGCODE
!
 call read_dump(trim(dumpfilein2),time,hfact,idisk1,iprint,0,1,ierr)           !CHNGCODE
 if (mhd .and. ierr==is_not_mhd) then                                          !CHNGCODE
    ihavesetupB = .false.                                                      !CHNGCODE
 elseif (ierr /= 0) then                                                       !CHNGCODE
    stop 'error reading dumpfile'                                              !CHNGCODE
 endif                                                                         !CHNGCODE
 
 nstar2 = npart                                                                 !CHNGCODE
 xyzh2  = xyzh                                                                  !CHNGCODE
!
!-- Merge both dumpfiles                                                       !CHNGCODE
!
 npart = nstar1+nstar2                                                          !CHNGCODE
 npartoftype(igas) = npart
 do i = 1,nstar1                                                                !CHNGCODE
    xyzh(:,i) = xyzh1(:,i)                                                     !CHNGCODE
 end do                                                                         !CHNGCODE

 do i = nstar1+1, npart                                                         !CHNGCODE
    xyzh(:,i) = xyzh2(:,i-nstar1)                                              !CHNGCODE
 end do                                                                         !CHNGCODE

 do i = 1,npart                                                                 !CHNGCODE
    vxyzu(1:3,i) = 0                                                           !CHNGCODE
 end do                                                                         !CHNGCODE


!
 call modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu,nstar1,nstar2)       !CHNGCODE BAD METHODOLOGY!!! WE SHOULD BE PASSING nstar1 and nstar2 as a use clause 
!
!--perform sanity checks on the output of modify_dump routine
!
 call check_setup(nerr,nwarn,restart=.true.)
 if (nwarn > 0) call warning('moddump','warnings from modified setup',var='warnings',ival=nwarn)
 if (nerr > 0)  call fatal('moddump','errors in modified setup',var='errors',ival=nerr)

 idumpsphNG = .false.
 if (nargs >= 4) then                                      !CHNGCODE
    call get_command_argument(4,string)                    !CHNGCODE  
    read(string,*,iostat=ierr) timeout
    if (ierr /= 0) then
       if (nargs==4) then
          timeout = time
          print*,' using time = ',timeout,' from dump file'
          if (index(string,'sphNG') /= 0) then
             idumpsphNG = .true.
          endif
       else
          stop 'error reading output time from command line'
       endif
    else
       print*,' setting time = ',timeout
    endif
 else
    timeout = time
 endif
 if (nargs==5) then                                         !CHNGCODE
    call get_command_argument(5,string)                     !CHNGCODE
    if (index(string,'sphNG') /= 0) then
       idumpsphNG = .true.
    endif
 endif

 if (mhd) then
    if (ihavesetupB) then
       ians = .false.
       call prompt(' add/reset magnetic fields?',ians)
    else
       ians = .true.
    endif
    if (ians) then
       call set_Bfield(npart,npartoftype(:),xyzh,massoftype(:),vxyzu,polyk, &
                       Bevol,maxvecp,Bextx,Bexty,Bextz)

    endif
 endif

 call write_fulldump(timeout,dumpfileout,sphNG=idumpsphNG)
!
!--write a fresh input file, whether it exists or not
!
 if (id==master .and. .not.idumpsphNG) then
    call write_infile(infile,logfile,evfile,dumpfileout,iwritein,6)
    print "(a,/,/,a)",' To start the calculation, use: ',' ./phantom '//trim(infile)
 endif

 print "(/,a,/)",' Phantom moddump: another happy customer'

end program phantommoddump_binary_system

