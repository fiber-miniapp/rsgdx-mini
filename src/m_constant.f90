!C***
!C*** Module : Constant
!C***
!C***                         LAST CHANGE: 2013/01/17 16:52:06.
!C***
      Module Constant
      USE precisions
      IMPLICIT NONE

      integer(kp), parameter :: NMAX        = 1228800
      integer(kp), parameter :: NPROC       =    32
      integer(kp), parameter :: NPEperGROUP =    32
      integer(kp), parameter :: NSAVMAX     =      10
      integer(kp), parameter :: nskip_factor_for_eqk=5

      integer(kp), parameter :: NBMAX       = NMAX/NPEperGROUP
!C
!C +----------------------------------------------------+
!C | if TRUE : print time information at every timestep |
!C +----------------------------------------------------+
!C===
      logical                :: DEBUG          = .FALSE.
!C===
      logical                :: BLAS           = .FALSE.
      logical                :: OUT_LARGE_FILE = .TRUE.
      logical                :: OUT_TIMESEQ    = .FALSE.

      character(len=20)      :: fn_out, fn_rst
      integer(kp)            :: irest
      real(DP), allocatable  :: cslip(:)
      real(DP)               :: tta    (0:NMAX-1)
      real(DP)               :: ttb    (0:NMAX-1)
      real(DP)               :: ttL    (0:NMAX-1)
      real(DP)               :: ttsigma(0:NMAX-1)
      real(DP)               :: ttmask (0:NMAX-1)
      real(DP)               :: ttVp   (0:NMAX-1)

      End Module Constant