!C-----------------------------------------------------
!C ＊入力されたコミュニケータ内での同期時間測定版
!C
!C                     LAST CHANGE:2012/05/23 16:56:20.
!C-----------------------------------------------------
      module timer_util

!      use mpi

      IMPLICIT NONE
      include 'mpif.h'

      include 'precision_pri.inc'

      PRIVATE
      type :: timer_type
        real(DP) :: sec_inacycle
        real(DP) ::      begin_time
        real(DP) ::        now_time
        real(DP) ::     elapse_time
        real(DP) :: max_elapse_time
      end type timer_type

      type(timer_type)        :: job_time
      integer(I4B), parameter :: NML_CODE=33
      logical                 :: TIMER_LIMIT
      logical                 :: TIMER_ON
      real(DP)                :: stopwatch=0.0e0_dp

      PUBLIC :: TIMER_LIMIT, TIMER_ON, stopwatch
      PUBLIC :: timer__init, timer__mpi_check

      CONTAINS
!***
!*** public : get reference time
!***
      SUBROUTINE timer__init (icomm)
      integer(I4B), intent(IN) :: icomm  ! ダミー
      real(DP)     :: max_job_time, safety
      real(DP)     :: system_sec
      real(DP)     :: system_cycle
      integer(I4B) :: OpenStatus
      namelist /input_timer/ max_job_time, safety

      open (NML_CODE, file='input_timer.nml', STATUS='old',  &
     &                                   iostat=OpenStatus )

      if ( OpenStatus > 0 ) then
        TIMER_ON=.FALSE.
        print *, 'timer off'
        RETURN
      endif

      read (NML_CODE,nml=input_timer)
!     write(   *    ,nml=input_timer)
      close(NML_CODE)

      job_time%max_elapse_time = max_job_time*safety
      call CPU_TIME (system_sec, 0 )
      job_time%begin_time   = system_sec
      call SYS_CYCLE(system_cycle)
      job_time%sec_inacycle = system_cycle
      TIMER_ON    = .TRUE.
      TIMER_LIMIT = .FALSE.
      RETURN
      end SUBROUTINE timer__init

!***
!*** public : get current time & check timer
!***
      SUBROUTINE timer__mpi_check (icomm)
      integer(I4B), intent(IN) :: icomm
      real(DP)           :: system_sec
      real(DP),     save :: previous_sec = 0.0e0_dp
      real(DP)           :: rtime, rtime_max
      integer(I4B), save :: ncycle = 0
      integer(I4B)       :: ierr

      if ( .not. TIMER_ON  ) RETURN

      call CPU_TIME (system_sec, ncycle)
      job_time%now_time    = system_sec
      if ( system_sec < previous_sec ) then
        ncycle = ncycle + 1
        job_time%now_time  = system_sec         &
     &                      +job_time%sec_inacycle*real(ncycle,DP)
      endif

      job_time%elapse_time = job_time%now_time - job_time%begin_time

      rtime = job_time%elapse_time / job_time%max_elapse_time

      call MPI_allREDUCE(rtime,rtime_max,1,MPI_DOUBLE_PRECISION,   &
     &                   MPI_MAX,          icomm, ierr) ! new
!    &                   MPI_MAX, MPI_COMM_WORLD, ierr) ! old

      if ( rtime_max > 1.0e0_dp ) then
        TIMER_LIMIT = .TRUE.
        stopwatch   = rtime_max*job_time%max_elapse_time
      else
        TIMER_LIMIT = .FALSE.
      endif
!     previous_sec = job_time%elapse_time ! miss 20090323 ?
      previous_sec = job_time%now_time
      RETURN
      end SUBROUTINE timer__mpi_check

!***
!***
!***
      SUBROUTINE SYS_CYCLE(sec)
      real(DP), intent(out) ::sec
      integer(I4B)          :: counter, crate, cmax
      call system_clock(COUNT=counter,COUNT_RATE=crate, COUNT_MAX=cmax)
      sec = real(cmax,DP)/real(crate,DP)
!CO   print *, counter, crate, cmax
      RETURN
      END SUBROUTINE SYS_CYCLE
!***
!***
!***
      SUBROUTINE CPU_TIME(sec, ncyc)
      real(DP)    , intent(out) ::sec
      integer(I4B), intent(in)  :: ncyc
      integer(I4B)              :: counter, crate
      call system_clock(COUNT=counter,COUNT_RATE=crate)
      sec = real(counter,DP)/real(crate,DP)        &
     &     +real(ncyc   ,DP)*job_time%sec_inacycle
      RETURN
      END SUBROUTINE CPU_TIME
      end module timer_util
