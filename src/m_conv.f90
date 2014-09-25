!C
!C   Created by 兵藤 守 on 11/02/17.
!C   Last Change: 2012/09/13 19:59:34.
!C +------------------------------------------+
!C |                                          |
!C +------------------------------------------+
!C***
!C*** MODULE m_conv : コンボリューション用の変数定義モジュール
!C***
      MODULE m_conv
      IMPLICIT NONE
      include 'precision_pri.inc'
!C
!C +-----------+
!C | file code |
!C +-----------+
!C===
!     integer(I4B), parameter                    ::     gfdcode=99 !ES
      integer(I4B), parameter                    ::     gfdcode=12 !K
      integer(I4B), parameter                    ::     outcode=72
      integer(I4B), parameter                    ::    obsfcode=73
      integer(I4B), parameter                    :: selectfcode=74
      integer(I4B), parameter                    :: inpobsfcode=75
!C===
      integer(I4B)                               ::    NOBS
      integer(I4B)                               :: NSELECT

      integer(I4B),allocatable, dimension(:)     :: INDEX_SELECTOBS

      real   (DP), allocatable, dimension(:,:,:) ::  gfd
!CO   real   (DP), allocatable, dimension(:,:)   ::   us
      real   (DP), allocatable, dimension(:,:)   ::   ud
!CO   real   (DP), allocatable, dimension(:,:)   ::   us_sum
      real   (DP), allocatable, dimension(:,:)   ::   ud_sum
      real   (DP), allocatable, dimension(:)     :: olat,olon
!CO   real   (DP), allocatable, dimension(:)     :: oX  ,oY
!CO   real   (DP), allocatable, dimension(:)     :: slp_array
      real   (DP), allocatable, dimension(:)     :: def_array

      character(len=200)                         ::     gfdname
      character(len=200)                         ::    obsfname
      character(len=200)                         :: selectfname
      character(len= 20)                         ::      cfname

!CO   real   (DP), allocatable, dimension(:,:,:)   ::   tus_sum
      real   (DP), allocatable, dimension(:,:,:)   ::   tud_sum
      real   (DP), allocatable, dimension(:)       ::  timey_ts

      END MODULE m_conv
