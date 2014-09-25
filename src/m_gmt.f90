!C***
!C*** MODULE gmt
!C***
      MODULE gmt
      USE precisions
      USE constant
      IMPLICIT NONE

      REAL   (DP)            :: VovrVp(0:NMAX-1)
      REAL   (DP)            ::    lat(0:NMAX-1)
      REAL   (DP)            ::    lon(0:NMAX-1)
      INTEGER(KP)            :: nsav_gmt_write    =0
      INTEGER(KP)            :: nsav_gmt_write_pre=0 ! 20120824 for inc plot

      REAL   (DP)            :: ned_x(0:NMAX-1)
      REAL   (DP)            :: ned_y(0:NMAX-1)
      REAL   (DP)            :: ned_z(0:NMAX-1)

      logical                :: PLOT_SCRIPT1=.false.
      logical                :: PLOT_SCRIPT2=.false.
      logical                :: PLOT_SCRIPT3=.false.

      INTEGER(KP), PARAMETER :: PLOT_BOUND=0!==> V=Vpl

      INTEGER(KP), PARAMETER :: code01    =21 ! post01.sh
      INTEGER(KP), PARAMETER :: code02    =22 ! post02.sh
      INTEGER(KP), PARAMETER :: code03    =23 ! post03.sh
      INTEGER(KP), PARAMETER :: ifile_code=24 ! input_gmt.nml&geomery_ll,ned
      INTEGER(KP), PARAMETER :: vfile_code=25 ! velo_#####

      CONTAINS

!C***
!C*** subroutine gmt__init : read files for GMT plot
!C***
      SUBROUTINE gmt__init (ierr)
!C--  for only rank 0
      INTEGER(KP) :: ierr
      INTEGER(KP) :: i, OpenStatus
      namelist /gmt_control/ PLOT_SCRIPT1, PLOT_SCRIPT2,PLOT_SCRIPT3

      ierr=0
      open (ifile_code, file='input_gmt.nml',status='old',              &
     &                                    iostat=OpenStatus)
      if (OpenStatus/=0) then
        print *,"input_gmt.nml open error"
      endif
      read (ifile_code, nml=gmt_control)
      close(ifile_code)

      open (ifile_code, file='geometryn_ll', status='old',              &
     &                                    iostat=OpenStatus)
      if (OpenStatus/=0) then
        print *,"geometryn_ll open error"
      endif
      do i=0, NMAX - 1
        read (ifile_code,  *) lon(i), lat(i)
      enddo
      close(ifile_code)

      open (ifile_code, file='geometryn_ned', status='old',             &
     &                                    iostat=OpenStatus)!20120726
      if (OpenStatus/=0) then
        print *,"geometryn_ned open error"
      endif
      do i=0, NMAX - 1
        read (ifile_code,  *) ned_x(i), ned_y(i), ned_z(i)
      enddo
      close(ifile_code)

      if(PLOT_SCRIPT1) open(code01, file='./post01.sh',status='unknown')
      if(PLOT_SCRIPT2) open(code02, file='./post02.sh',status='unknown')
      if(PLOT_SCRIPT3) open(code03, file='./post03.sh',status='unknown')
      RETURN
      END SUBROUTINE gmt__init
!C***
!C*** subroutine gmt__finalize : just a file closer
!C***
      SUBROUTINE gmt__finalize (ierr)
      integer (KP), intent(INOUT) :: ierr
      ierr=0
      if (PLOT_SCRIPT1) close (code01)
      if (PLOT_SCRIPT2) close (code02)
      if (PLOT_SCRIPT3) close (code03)
      RETURN
      END SUBROUTINE gmt__finalize

!C***
!C*** subroutine gmt__main : call GMT script if necessaly
!C***
      SUBROUTINE gmt__main (next_call_delay,time_sec,Vp,MASK,Vref,myid)
!C--  for all rank
      integer                  :: next_call_delay
      integer (KP), intent(IN) :: myid
      real    (DP), intent(IN) :: time_sec, Vref
      real    (DP), intent(IN) :: Vp(0:NMAX-1), MASK(0:NMAX-1)
      integer(KP)              :: iyr, iday, ihour, imin, isec
      integer(KP)              :: icou=0, i, imaxv
      character(len=20)        :: fname
      character(len=100)       :: gmt_script
!CO   character(len=100)       :: command_clean
      real   (DP)              :: maxv

      if ( myid /=0 ) RETURN

      icou=icou+1
      call S2YDHMS_dp(time_sec,iyr,iday,ihour,imin,isec)

      VovrVp(:)=VovrVp(:)+Vp(:)*MASK(:) !SLIP velocity of ALL CELLS
      maxv=maxval(VovrVp(:))
      imaxv=int(log10(maxv))
      if (imaxv<=0) then
        next_call_delay=1
      else
        next_call_delay=imaxv*nskip_factor_for_eqk+1 !20120929
      endif

      if ( (.not.PLOT_SCRIPT1) .and. (.not.PLOT_SCRIPT2)                &
                               .and. (.not.PLOT_SCRIPT3) ) then
        next_call_delay=0
        return
      endif

!CO   print *, "(icou, next_call_delay)=", icou, next_call_delay

      write(fname,'(a,i5.5)') 'velo_', icou
      open (vfile_code,file=trim(fname), status='unknown' )
      write(vfile_code,'(a1,2x,i5.5,a2,i3.3,a2,2(i2.2,a2),i2.2,a1)' )   &
     &     '#',iyr,'y_',iday,'d_',iHour,'h_',imin,'m_',isec,'s'
      do i=0, NMAX-1
        if (MASK(i)/=0.0) then
          write(vfile_code,'(8(1pe19.9))')                              &
     &                        lon(i),lat(i),log10(VovrVp(i)),           &
     &              time_sec/60.0d0/60.0d0/24.0d0/365.2425d0,           &
     &    ned_x(i), ned_y(i), ned_z(i), (VovrVp(i)*Vp(i)-Vp(i))*Vref
        endif
      enddo
      close(vfile_code)

      if (PLOT_SCRIPT1) then
        gmt_script='sh gmt_script1 '//trim(fname)
        write(code01,'(a)') trim(gmt_script)
!CO     print *, trim(gmt_script)
!CO     call system(trim(gmt_script))
      endif

      if (PLOT_SCRIPT2) then
        gmt_script='sh gmt_script2 '//trim(fname)
        write(code02,'(a)') trim(gmt_script)
!CO     print *, trim(gmt_script)
!CO     call system(trim(gmt_script))
      endif

      if (PLOT_SCRIPT3) then
        if (imaxv<=PLOT_BOUND) then
!CO       gmt_script='sh gmt_script3 '//trim(fname)//' cd.out'
!COCO     gmt_script='mv '//' cd.out'//' cd_'//trim(fname)
!CO       print *, trim(gmt_script)
!COCO     call system(trim(gmt_script))
          gmt_script='sh gmt_script3 '//trim(fname)//' cd_'//trim(fname)
          write(code03,'(a)') trim(gmt_script)
!CO       command_clean='rm -rf cd.out'
!CO       call system(trim(command_clean))
        endif
      endif

!CO   command_clean='rm -rf '//trim(fname)
!CO   print *, trim(command_clean)
!CO   call system(trim(command_clean))
      RETURN
      END SUBROUTINE gmt__main

!C***
!C*** subroutine S2YDHMS_dp
!C***
      SUBROUTINE S2YDHMS_dp                                             &
     &           (tot_seconds,iyears,idays,ihours,iminutes,iseconds)

      real   (DP) , intent(IN ):: tot_seconds
      INTEGER(KP ), intent(OUT):: iyears,idays,ihours,iminutes,iseconds

      real   (DP)              :: remained_seconds
      REAL    (DP)             :: ryears,rdays,rhours,rminutes,rseconds
      REAL(DP), PARAMETER      :: Uminute =    60.00e0_dp
      REAL(DP), PARAMETER      :: Uhour   =    60.00e0_dp*Uminute
      REAL(DP), PARAMETER      :: Uday    =    24.00e0_dp*Uhour
      REAL(DP), PARAMETER      :: Uyr     = 365.2425e0_dp*Uday

      remained_seconds = tot_seconds

      iyears = int (remained_seconds / Uyr  )
      ryears = real(iyears, DP)
      remained_seconds = remained_seconds - ryears*Uyr

      idays  = int (remained_seconds / Uday )
      rdays  = real(idays , DP)
      remained_seconds = remained_seconds - rdays *Uday

      ihours = int (remained_seconds / Uhour)
      rhours = real(ihours, DP)
      remained_seconds = remained_seconds - rhours*Uhour

      iminutes = int (remained_seconds / Uminute)
      rminutes = real(iminutes, DP)
      remained_seconds = remained_seconds - rminutes*Uminute

      rseconds = remained_seconds
      iseconds = int ( rseconds )
      RETURN
      END SUBROUTINE S2YDHMS_dp
      end module gmt
