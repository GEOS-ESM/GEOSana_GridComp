program pcscaling_main

  use satinfo_util, only: DP
  use satinfo_util, only: STDIN
  use satinfo_util, only: die,perr,tell

  use m_pcscaling, only: pcscaling_nmlread
  implicit none

  character(len=*),parameter:: myname='pcscaling_main'
  integer,parameter:: NPRED_=12

  integer:: luinp,luout,ios
  logical:: changed

  character(len=512):: pcinp          ! satbias_pc file for input
  character(len=512):: pcout          ! satbias_pc file for output
  character(len=512):: dbname         ! where pcscaling_nml is stored
  integer:: nymd        ! YYYYMMDD
  integer:: nhms        ! HHMMSS
  integer:: npred       ! number of predictor coefficients

  namelist/setup/ &
    nymd,nhms,npred, &
    pcinp,pcout, &
    dbname

  npred=NPRED_
  nymd=-1
  nhms=-1
  pcinp="satbias_pc"
  pcout="satbias_pc.scaled"
  dbname="pcscaling.db"

  read(STDIN,setup)          ! read a configuration
  !write(*,setup)

        if(nymd<=0.or.nhms<0) then
          call perr(myname,'invalid/missing setup parameter, nymd =',nymd)
          call perr(myname,'                                 nhms =',nhms)
          call  die(myname)
        endif

  ! Open a satbias_pc file for input
  open(newunit=luinp,file=pcinp,form='formatted',status='old',action='read' ,iostat=ios)
        if(ios/=0) call die(myname,'open("'//trim(pcinp)//'",status="old","read") error, iostat =',ios)

  ! Open a satbias_pc file for output
  open(newunit=luout,file=pcout,form='formatted',status='new',action='write',iostat=ios)
        if(ios/=0) call die(myname,'open("'//trim(pcout)//'",status="new","write") error, iostat =',ios)

  ! Load pcscaling factors from file pcscale
  call pcscaling_nmlread(trim(dbname)//"/pcscaling.nml")

  ! Loop over pcinp entries for scaling, for any nymdh matching record
  call procloop_(luinp,luout,(nymd*100+nhms/10000),npred,changed)

  close(luinp)
  if(changed) then
    call tell(myname,'some entry changed and output kept, file =',trim(pcout))
    close(luout,status='keep')
  else
    call tell(myname,'no entry changed and output removed, file =',trim(pcout))
    close(luout,status='delete')
  endif

contains
subroutine procloop_(lunin,luout,nymdh,npred,changed)
  use m_pcscaling, only: pcscaling_yes, pcscaling_apply
  implicit none
  integer,intent(in):: lunin       ! unit for input satbias_pc
  integer,intent(in):: luout       ! unit for output satbias_pc
  integer,intent(in):: nymdh       ! yyyymmddhh of this run
  integer,intent(in):: npred       ! size of varx(:) array
  logical,intent(out):: changed                 ! .true., if any record been scaled.

  character(len=*),parameter:: myname_=myname//"::procloop_"

  character(len=20):: isis
  integer:: iread,ich,ichan,ip,istat
  real(kind=DP):: ostatsx
  real(kind=DP),dimension(npred):: varx

  call tell(myname_,'entered with nymdh =',nymdh)
  call tell(myname_,'             npred =',npred)

!!! See Loop read3 in GSI/radinfo.f90 for a comparison
          istat=0
          iread=0
          changed=.false.
          read3: do
             iread=iread+1
             read(lunin,'(I5,1x,A20,1x,I5,e15.7/2(4x,10e15.7/))',iostat=istat) &
                  ich,isis,ichan,ostatsx,varx(1:npred)
             if (istat/=0) exit

             if(pcscaling_yes) call pcscaling_apply(nymdh,isis,ichan,varx(1:npred),changed=changed)

             write(luout,'(I5,1x,A20,1x,I5,e15.7/2(4x,10e15.7/))') &
                  ich,isis,ichan,ostatsx,varx(1:npred)
          end do read3

          if (istat>0.and..not.is_iostat_end(istat)) then
            call perr(myname_,'loop read3 error, iostat =',istat)
            call perr(myname_,'                   iread =',iread)
            call perr(myname_,'                     ich =',ich)
            call perr(myname_,'                    isis =',isis)
            call perr(myname_,'                   ichan =',ichan)
            call  die(myname_)
          endif
end subroutine procloop_

end program pcscaling_main
