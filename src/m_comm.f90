!C***
!C*** module m_comm : make multi-communicators
!C***
!C***           LAST CHANGE: 2012/05/23 16:58:14.
      MODULE m_comm

      IMPLICIT NONE
      include "precision_pri.inc"

      TYPE, public :: comm
        LOGICAL       :: member = .false.
        INTEGER (I4B) :: myrank
        INTEGER (I4B) :: nsize
        INTEGER (I4B) :: comm
        INTEGER (I4B) :: group
      END TYPE comm
      PRIVATE

      TYPE  (comm), POINTER,public :: comms_new(:)
      INTEGER(I4B), ALLOCATABLE    :: ranks    (:)

      INTEGER(I4B)                 :: ngrp
      INTEGER(I4B)                 :: nmember_grp
      INTEGER(I4B), PUBLIC         :: NthG

      public :: comm__init, comm__finalize
      CONTAINS
      SUBROUTINE comm__init(comm_all,group_all,nproc,myrank,NPperGRP)
!CO   USE mpi
      include 'mpif.h'

      INTEGER(I4B), intent(INOUT) ::  comm_all
      INTEGER(I4B), intent(INOUT) :: group_all
      INTEGER(I4B), intent(IN)    ::  nproc, myrank
      INTEGER(I4B), intent(IN)    ::  NPperGRP
      INTEGER(I4B)                :: i, is, ie, ig, irank, ierr
      INTEGER(I4B)                ::  comm_rank
      INTEGER(I4B)                ::  comm_size
      INTEGER(I4B)                ::  comm_div
      INTEGER(I4B)                :: group_div

      nmember_grp=NPperGRP

      call mpi_comm_dup   (MPI_COMM_WORLD, comm_all, ierr)
!C===
      if ( mod (nproc, nmember_grp)/=0 ) STOP 'invalid total processes'
      ngrp = nproc/nmember_grp

      call mpi_comm_group (comm_all, group_all, ierr)

      allocate (ranks(nproc))
      do i=1, nproc
        ranks(i)=i-1
      enddo

      allocate (comms_new (ngrp))
!C
!C +----------------------+
!C | 新しいグループを生成 |
!C +----------------------+
!C===
      do ig = 1, ngrp
        is=(ig-1)*nmember_grp+1-1
        ie=(ig  )*nmember_grp  -1
        call mpi_group_incl                                             &
     &      (group_all, nmember_grp, ranks(is+1:ie+1),group_div,ierr)
        comms_new(ig)%group = group_div
        !C
        !C +----------------------------------+
        !C | 分割グループig番目のメンバー設定 |
        !C +----------------------------------+
        !C===
        do irank= is, ie
          if ( myrank == irank ) then
            comms_new(ig)%member = .TRUE.
          endif
        enddo
        !C===
        call Mpi_barrier    (comm_all,ierr)

        group_div = comms_new(ig)%group
        call mpi_comm_create( comm_all, group_div, comm_div, ierr)
        comms_new(ig)%comm = comm_div

        if (comms_new(ig)%member) then
          comm_div = comms_new(ig)%comm
          call mpi_comm_size ( comm_div, comm_size, ierr )
          call mpi_comm_rank ( comm_div, comm_rank, ierr )
          comms_new(ig)%nsize  = comm_size
          comms_new(ig)%myrank = comm_rank
          call mpi_barrier   ( comm_div,            ierr )
        endif
      enddo
!C===

!C
!C +----------------+
!C | "NthG"をセット |
!C +----------------+
!C===
      do ig=1, ngrp
        if (comms_new(ig)%member) then
          NthG=ig
          exit
        endif
      enddo
!CO   print *,'myrank=',myrank,'comm_rank=',comm_rank,'NthG=',NthG
!C===
      RETURN
      END SUBROUTINE comm__init

      SUBROUTINE comm__finalize(comm_all, group_all)
!CO   USE mpi
      include 'mpif.h'

      INTEGER(I4B), intent(IN) ::  comm_all
      INTEGER(I4B), intent(IN) :: group_all
      INTEGER(I4B)             :: ig, ierr
      INTEGER(I4B)             ::  comm_div
      INTEGER(I4B)             :: group_div

      do ig=1, ngrp
        if ( comms_new(ig)%member ) then
          comm_div  = comms_new(ig)%comm
          call mpi_comm_free  ( comm_div, ierr)
          group_div = comms_new(ig)%group
          call mpi_group_free (group_div, ierr)
        endif
      enddo
      call mpi_comm_free  ( comm_all, ierr)
      call mpi_group_free (group_all, ierr)
      RETURN
      END SUBROUTINE comm__finalize
      END MODULE m_comm
