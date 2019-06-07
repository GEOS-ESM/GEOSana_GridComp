#include <preproc.h>
      subroutine restwrt (nstep, beglatspmd, endlatspmd)

#if (defined SPMD )
      use mod_comm, only : mp_add1d, mp_gath_3d
#endif 
#include <vartyp.h> 
#include <lsmpar.h>
#include <lsmspmd.h>
#include <basinpar.h>

! ------------------------ input variables ------------------------
      integer nstep          !current time step
      integer beglatspmd     !beg lat on proc, degenerates to 1 for nonspmd
      integer endlatspmd     !end lat on proc, degenerates to jm for nonspmd
! -----------------------------------------------------------------
#include <dynmem.h>
#include <lsmtcdyn.h>
#include <lsmhis.h>
#include <lsmhisdyn.h>
#include <lsmctl.h>
#include <lsmtvdyn.h>
#include <lsmtc.h>
#include <cplavg.h>
#include <basin.h>

! ------------------------ local variables ------------------------
      integer i              !loop index
      character*80 icom      !dummy mass store comment 
      integer nio            !fortran unit number
      integer getavu         !integer function: next available unit number
      integer lencs          !integer function: length of character string
      integer kloc           !length of character string
      character*80 hname     !temporary file name
      character*80 rname     !temporary file name
      character*80 tmpfn     !temporary file name
      character*80 mssfnr    !mass store restart file name
      character*80 mssfnh    !mass store history file name
#if (defined SPMD )
      integer gid_io, k, n
      integer k1, k2
      gid_io = 0
#endif
* -----------------------------------------------------------------

      icom = ' '

#if (defined SPMD )

      k1 = begkptspmd(beglatspmd)
      k2 = endkptspmd(endlatspmd)

! type-1
      call mp_gath_3d(h2osno, kpt, 1, 1, k1, k2, 
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(h2ocan, kpt, 1, 1, k1, k2, 
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(tv,     kpt, 1, 1, k1, k2, 
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(tg,     kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(moz,    kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(eah,    kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(soot,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(hsno,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(fsno,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(fwet,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(htop,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(tlai,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(tsai,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(elai,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(esai,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(foln,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(stemb,  kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(rootb,  kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(soilc,  kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(igs,    kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)
      call mp_gath_3d(fsun,   kpt, 1, 1, k1, k2,
     &                1, 1, 1, 1, gid_io)

! type-2
      call mp_gath_3d(h2osoi, msl,   1, kpt, 1, msl,   1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(tsoi,   msl, 1,   kpt, 1, msl,   1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(albd,   mband, 1, kpt, 1, mband, 1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(albi,   mband, 1, kpt, 1, mband, 1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(albgrd, mband, 1, kpt, 1, mband, 1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(albgri, mband, 1, kpt, 1, mband, 1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(fabd,   mband, 1, kpt, 1, mband, 1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(fabi,   mband, 1, kpt, 1, mband, 1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(ftdd,   mband, 1, kpt, 1, mband, 1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(ftid,   mband, 1, kpt, 1, mband, 1, 1,
     &                k1, k2, gid_io)
      call mp_gath_3d(ftii,   mband, 1, kpt, 1, mband, 1, 1,
     &                k1, k2, gid_io)

! Type-3
      call mp_gath_3d(slfval, kpt, 1, mslflds, k1, k2,
     &                1, 1, 1, mslflds, gid_io)
      call mp_gath_3d(mlfval, kpt, 1, msl*mmlflds, k1, k2,
     &                1, 1, 1, msl*mmlflds, gid_io)

      call mp_add1d(ndrn, bashist)
#endif

      if (masterproc) then

! open main restart file. write data. close. dispose to mass store

         nio = getavu()
         call opnfil (locfnr, nio, 'u')

! write land surface variables needed for restart and that can be used 
! as initial data for new run

         write (nio) numlon           ! [lsmlat]
         write (nio) latixy           ! [numlon,lsmlat]
         write (nio) longxy
         write (nio) surf2d
         write (nio) soic2d
         write (nio) sand2d
         write (nio) silt2d
         write (nio) clay2d
         write (nio) pctlak
         write (nio) pctwet
         write (nio) drngarea
         write (nio) drngfrac
         write (nio) drngbasn
         write (nio) h2osno
         write (nio) h2ocan
         write (nio) h2osoi
         write (nio) tv    
         write (nio) tg    
         write (nio) tsoi  

! write other main land surface variables needed for restart

         write (nio) moz, eah, soot, hsno, fsno, fwet

! write vegetation for next time step

         write (nio) htop, tlai  ,tsai,  elai , esai,  
     $               foln, stemb ,rootb, soilc, igs

! write albedo for next time step

         write (nio ) albd ,albi ,albgrd ,albgri ,fabd ,
     $                fabi ,ftdd ,ftid   ,ftii   ,fsun    

! write history file related info

         write (nio) nstep,nfil,ntim
         write (nio) nmon,nyr
         write (nio) mcdate_f,mcsec_f,mdcur_f,mscur_f
         write (nio) mcdate_t,mcsec_t,mdcur_t,mscur_t  
         write (nio) slfval,mlfval
         write (nio) slfcnt,mlfcnt
         write (nio) bashist
         write (nio) bascnt
         write (nio) flondat,flnddat, finidat

#ifndef COUP_CCM
! write flux average info
         write (nio) flxave
!        if (flxave) then
!           write (nio) dosend
!           write (nio) icnt
!           write (nio) tauxa, tauya, lhflxa, shflxa, lwupa, qflxa,
!    $                  drnveca
!        endif
#endif

         call relavu (nio)

         write (6,*) 
         write (6,*) '(RESTWRT): successfully wrote local restart ',
     $        'file ',locfnr(1:lencs(locfnr))
         call putfil (locfnr, msspnr, irt, nswrps, icom)

! determine restart file name. these are the mass store names 
! except if irt=0. (no mass store files written) 

         if (irt .eq. 0) then
            rname = locfnr
         else
            kloc = lencs(locfnr)
            do i = kloc, 1, -1
               if (locfnr(i:i).eq.'/') go to 100
            end do
            i = 0
  100       rname = msspnr(1:lencs(msspnr))//'/'// locfnr(i+1:kloc)
         endif

! determine history file name. these are the mass store names 
! except if irt=0. If current history file was closed then set
! ncopnfil to false. (no mass store files written)

         if (ncopnfil) then     
            if (irt .eq. 0)  then
               hname = locfnh(1:lencs(locfnh))
            else
               kloc = lencs(locfnh)
               do i = kloc, 1, -1
                  if (locfnh(i:i).eq.'/') go to 200
               end do
               i = 0
  200          hname = msspnh(1:lencs(msspnh))//'/'// locfnh(i+1:kloc)
            endif
            ncopnfil = .false.
         else
            hname = 'history_file_not_open'
         endif

! open local restart pointer file on local disk. write names of 
! current restart and history files.

         nio = getavu()
         tmpfn= rpnt_dir(1:lencs(rpnt_dir)) //'/'
     &          //rpnt_fil(1:lencs(rpnt_fil))
         call opnfil (tmpfn, nio, 'f')
         write (nio,'(a80)') rname
         write (nio,'(a80)') hname
         call relavu (nio)
         write (6,*) 
         write (6,*) '(RESTWRT): successfully wrote local restart ',
     $        'pointer file ',tmpfn(1:lencs(tmpfn))

      end if
      return
      end
