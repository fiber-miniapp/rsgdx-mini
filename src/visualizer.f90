!C
!C    * add sliponly & deficit etc.  on 20120501
!C                         LAST CHANGE: 2012/09/13 20:05:31.
      MODULE visualizer

      USE   constant , only : NPEperGROUP
      USE m_outputval, only : vsav, vsav_gather,                        &
     &                        ssav, ssav_gather,                        &
     &                       thsav,thsav_gather,                        &
!>mh on 20120501
     &                    sliponly,sliponly_gather,                     &
     &                    deficit , deficit_gather
!>mh on 20120501

      IMPLICIT NONE

      include 'mpif.h'

      private
      public :: merge_data, merge_init, merge_out
      integer, parameter           :: I4B            = 4
      integer, parameter           :: DP             = 8

      INTEGER(I4B), ALLOCATABLE    :: istat (:)
      INTEGER(I4B), POINTER        :: iStack(:), rStack(:)
      INTEGER(I4B), POINTER        :: IBUFs (:), IBUFr (:)
      REAL    (DP), POINTER        :: RBUFs (:), RBUFr (:)

      INTEGER(I4B)            :: ncell
      INTEGER(I4B), pointer :: lindex_item(:)
      INTEGER(I4B), pointer :: stack(:)
      INTEGER(I4B), pointer :: N2O  (:)
      INTEGER(I4B), pointer :: O2N  (:)
      INTEGER(I4B), pointer :: rank_gitem(:)
      INTEGER(I4B), pointer :: rank_litem(:)
      CONTAINS

!C***
!C*** merge_init
!C***
      SUBROUTINE merge_init(myrank)
      integer(I4B), intent(IN) :: myrank
      integer(I4B)             :: NPROCS, i, index_new, index_old
      integer(I4B)             :: icou, icell

      NPROCS=NPEperGROUP

      open (2, file='output_cell_information.dat', status='old' )
      read (2,*) ncell
      allocate(stack(0:NPROCS))
      allocate(N2O(ncell))
      allocate(O2N(ncell))
      allocate(rank_gitem(ncell))
      allocate(rank_litem(ncell))
      read (2, '(10i8)') stack
      read (2, '(10i8)') O2N
      read (2, '(10i8)') N2O
      read (2, '(10i8)') rank_gitem
      read (2, '(10i8)') rank_litem
      close(2)

      icell = stack(myrank+1) - stack(myrank)
      allocate(lindex_item(icell))
      lindex_item=0

      icou=0
      do i=stack(myrank)+1, stack(myrank+1)
        icou=icou+1

        index_new=i
        index_old=N2O(index_new)

!C-- set local id of output cells
         lindex_item(icou)=rank_litem(index_old)
!CO!>check
!CO         print *, *'rank, new, old_index,  gcell, lcell',           &
!CO     & irank, index_new, index_old, rank_gitem (index_old),         &
!CO     &                              rank_litem (index_old)
!CO!<check
      end do
      RETURN
      end subroutine merge_init

!C***
!C*** merge_out
!C***
      SUBROUTINE merge_out (ng, Nt, time)
      integer(I4B), INTENT(IN) :: Nt, ng
      real   (DP )             :: time(1:)
      integer(I4B)             :: it, ip, index_old, index_new, idx
      character(len=100)       :: ofname

      do idx=1, ncell
        write(ofname,'(I3.3,A,I6.6,A)')ng-1,'/slip_history_',idx,'.out'
        open (34, FILE=trim(ofname), status='unknown')

!       do ip=1, ncell
!         index_old=N2O(ip)
!         if (index_old==idx)
!           exit
!         endif
!       enddo

        index_new=O2N(idx)
        ip=index_new
        write(34,*) '#', idx, index_new, N2O(index_new)
        do it=1, Nt
          write(34,'(10(1pe15.7))') time(it),   vsav_gather(it,ip),     &
     &                    ssav_gather(it,ip),  thsav_gather(it,ip),     &
     &                sliponly_gather(it,ip),deficit_gather(it,ip)
        enddo
        close(34)
      enddo
      return
      end subroutine merge_out

!C***
!C*** merge_data : just a test on 2012年  4月 26日 木曜日 11:57:38 JST
!C***
      SUBROUTINE merge_data (myrank, comm, nproc, Ncol, nroot)
!C
!C   +----------------------------------------------+
!C   |     merge nodal group data for TIMESEQUENSE  |
!C   +----------------------------------------------+
!C
      INTEGER(I4B), INTENT(IN)      :: myrank, comm, nproc, Ncol, nroot
      INTEGER(I4B)                  :: Ns, ierr

      Ns=stack(myrank+1)-stack(myrank)
!CO   print *, 'myrank, Ncells', myrank, Ns

!C
!C-- INIT
      allocate ( istat  (MPI_STATUS_SIZE))       !..for MPI (dummy)
      allocate ( iStack ( 0 : Nproc)) ; iStack=0 !..STACK table for INT. ITEM
      allocate ( rStack ( 0 : Nproc)) ; rStack=0 !..STACK table for REAL ITEM

!C
!C   +------------+
!C   | merge vsav |
!C   +------------+
!C===
!     Ncol = (時系列数)                             !vsav
      call set_r2dSend(lindex_item, Ncol, Ns, RBUFs, vsav)
      call set_rStackG(     Ns,   RBUFr,    Ncol,    Nproc,           &
     &                  rStack,   nroot,  myrank,    comm     )
      call mpiGVreal  (     Ns,   RBUFs,  RBUFr,     Ncol,            &
     &                   Nproc,  rStack,  nroot,   myrank,  comm,     &
     &                                             vsav_gather  )
!C===

!C
!C   +------------+
!C   | merge ssav |
!C   +------------+
!C===
!     Ncol = (時系列数)                             !ssav
      call set_r2dSend(lindex_item, Ncol, Ns, RBUFs, ssav)
      call set_rStackG(     Ns,   RBUFr,    Ncol,    Nproc,           &
     &                  rStack,   nroot,  myrank,    comm     )
      call mpiGVreal  (     Ns,   RBUFs,  RBUFr,     Ncol,            &
     &                   Nproc,  rStack,  nroot,   myrank,  comm,     &
     &                                             ssav_gather  )
!C===

!C
!C   +-------------+
!C   | merge thsav |
!C   +-------------+
!C===
!     Ncol = (時系列数)                            !thsav
      call set_r2dSend(lindex_item, Ncol, Ns, RBUFs,thsav)
      call set_rStackG(     Ns,   RBUFr,    Ncol,    Nproc,           &
     &                  rStack,   nroot,  myrank,    comm     )
      call mpiGVreal  (     Ns,   RBUFs,  RBUFr,     Ncol,            &
     &                   Nproc,  rStack,  nroot,   myrank,  comm,     &
     &                                            thsav_gather  )
!C===

!C
!C   +----------------+
!C   | merge sliponly |
!C   +----------------+
!C===
!     Ncol = (時系列数)
      call set_r2dSend(lindex_item, Ncol, Ns, RBUFs,sliponly)
      call set_rStackG(     Ns,   RBUFr,    Ncol,    Nproc,           &
     &                  rStack,   nroot,  myrank,    comm     )
      call mpiGVreal  (     Ns,   RBUFs,  RBUFr,     Ncol,            &
     &                   Nproc,  rStack,  nroot,   myrank,  comm,     &
     &                                         sliponly_gather  )
!C===

!C
!C   +---------------+
!C   | merge deficit |
!C   +---------------+
!C===
!     Ncol = (時系列数)
      call set_r2dSend(lindex_item, Ncol, Ns, RBUFs, deficit)
      call set_rStackG(     Ns,   RBUFr,    Ncol,    Nproc,           &
     &                  rStack,   nroot,  myrank,    comm     )
      call mpiGVreal  (     Ns,   RBUFs,  RBUFr,     Ncol,            &
     &                   Nproc,  rStack,  nroot,   myrank,  comm,     &
     &                                          deficit_gather  )
!C===
      deallocate ( istat  ) !... comment out by MH on 2012/04/25
      deallocate ( iStack ) !... CLEAR stack table
      deallocate ( rStack ) !... CLEAR stack table
      RETURN
      END SUBROUTINE merge_data

!C***
!C***                NEW version 2007/04/12 (MPI_GATHER)
!C***
      SUBROUTINE    set_iStackG  (    Ns,   IRECV,     Nc,   Npe,     &
     &                            ISTACK,    root, myrank,  comm   )
!C
!C   +--------------------------------------------------+
!C   | make Send-Receive TABLE of INTEGER item for root |
!C   +--------------------------------------------------+
!C
      INTEGER(I4B), INTENT(IN)   ::  Ns, Nc, Npe, root, myrank, comm
      INTEGER(I4B), POINTER      :: IRECV(:), ISTACK(:)
      INTEGER(I4B)               :: i,  ierr
      INTEGER(I4B), DIMENSION(1) :: Ns_dum

      Ns_dum(1) = Ns
      call MPI_GATHER ( Ns       , 1, MPI_INTEGER,                   &
                        ISTACK(1), 1, MPI_INTEGER, root, comm, ierr )
      if (myrank == root ) then
        do i = 0, Npe - 1
          ISTACK (i+1) = ISTACK (i+1) + ISTACK (i)
        enddo
      endif
      if ( associated ( IRECV ) ) deallocate(  IRECV  )
      allocate ( IRECV ( ISTACK(Npe)*Nc ) )  !... POINTS * VecLEN
      IRECV = 0
      RETURN
      END SUBROUTINE set_iStackG

!C***
!C***  set_rStackG
!C***
      SUBROUTINE     set_rStackG(     Ns,  RRECV,     Nc,  Npe,     &
     &                            ISTACK,   root, myrank, comm   )
!C
!C   +-------------------------------------------------+
!C   | make Send-Receive TABLE of REAL item for rank 0 |
!C   +-------------------------------------------------+
!C
      INTEGER(I4B), INTENT(IN)   ::  Ns, Nc, Npe, root, myrank, comm
      INTEGER(I4B), POINTER      :: ISTACK(:)
      REAL    (DP), POINTER      :: RRECV(:)
      INTEGER(I4B)               :: i, ierr
      INTEGER(I4B), DIMENSION(1) :: Ns_dum

      Ns_dum(1) = Ns
      call MPI_GATHER ( Ns       , 1, MPI_INTEGER,                    &
                        ISTACK(1), 1, MPI_INTEGER, root, comm,  ierr )
      if (myrank == root ) then
        do i = 0, Npe - 1 !... rank
          istack (i+1) = istack (i+1) + istack (i)
        enddo
      endif
      if ( associated ( RRECV ) ) deallocate(  RRECV  )
      allocate ( RRECV ( ISTACK (Npe)*Nc ) )
      RRECV = 0.0_dp
      RETURN
      END SUBROUTINE set_rStackG

!C***
!C***                NEW version 2007/04/12 (MPI_GATHER)
!C***
      SUBROUTINE  mpiGVinteger (    Ns,  ISEND,  IRECV,   Nc,   Npe,  &
     &                          ISTACK,   root, myrank, comm, iAVS )
!C
!C    +----------------------------+
!C    | Send-Receive INTEGER items |
!C    +----------------------------+
!C
      INTEGER(I4B), INTENT(IN) :: Ns, Nc, Npe, root, myrank, comm
      INTEGER(I4B), POINTER    :: ISEND(:), IRECV(:)
      INTEGER(I4B), POINTER    :: iAVS(:,:)
      INTEGER(I4B), POINTER    :: ISTACK(:)
      INTEGER(I4B)             :: ierr, irk, ic, i
      INTEGER(I4B), DIMENSION(Npe) :: IRCNT, IDISP

      if (myrank== root) then
        do irk = 1, Npe
          IRCNT(irk) = ( ISTACK(irk)-ISTACK(irk-1) )*Nc
          IDISP(irk) =   ISTACK(irk-1)*Nc
        enddo
      endif
      call MPI_GATHERV (  ISEND   , Ns*Nc,        MPI_INTEGER,    &
     &                    IRECV   , IRCNT, IDISP, MPI_INTEGER,    &
     &                               root,  comm, ierr         )
      if (myrank== root) then
        if ( Nc == 1 ) go to 22
!C
!C     +-----------------------+
!C     | for connectivity only |
!C     +-----------------------+
!C===
       do i  = 1, Npe - 1
         do ic = ISTACK(i)*Nc + 1, ISTACK(i+1)*Nc
           IRECV ( ic ) = IRECV( ic ) + IRECV (ISTACK(i)*Nc)
         enddo
       enddo
22     continue
       if ( associated(iAVS) ) deallocate (iAVS)
       iAVS => pallocate_iv (IRECV, Nc, ISTACK(Npe))
!C===
      end if
      if ( associated(ISEND) ) deallocate ( ISEND )
      RETURN
      END SUBROUTINE mpiGVinteger

!C***
!C*** mpiGVreal
!C***
      SUBROUTINE    mpiGVreal (    Ns, RSEND,  RRECV,   Nc,  Npe,   &
     &                         ISTACK,  root, myrank, comm, rAVS   )
!C
!C   +---------------------------+
!C   | Send-Receive REAL   items |
!C   +---------------------------+
!C
      INTEGER(I4B), INTENT(IN) :: Ns, Nc, Npe, root, myrank, comm
      REAL    (DP), POINTER    :: RSEND(:), RRECV(:)
      REAL    (DP), POINTER    :: rAVS(:,:)
      INTEGER(I4B), POINTER    :: ISTACK(:)
      INTEGER(I4B)             :: ierr, irk
      INTEGER(I4B), DIMENSION(Npe) :: IRCNT, IDISP

      if (myrank== root) then
        do irk = 1, Npe
          IRCNT(irk) = ( ISTACK(irk)-ISTACK(irk-1) )*Nc
          IDISP(irk) =   ISTACK(irk-1)*Nc
        enddo
      endif
      call MPI_GATHERV (  RSEND   , Ns*Nc,        MPI_DOUBLE_PRECISION, &
                          RRECV   , IRCNT, IDISP, MPI_DOUBLE_PRECISION, &
     &                                      root, comm, ierr   )
      if (myrank== root) then
        if ( associated(rAVS) ) deallocate (rAVS)
        rAVS => pallocate_rv (RRECV, Nc, ISTACK(Npe))
      end if
      if ( associated(RSEND) ) deallocate ( RSEND )
      RETURN
      END SUBROUTINE mpiGVreal


!C***
!C*** set_i1dsend
!C***
      SUBROUTINE set_i1dsend (Gitem,Ncol,Ns,ISEND,Array1D)

      INTEGER(I4B), INTENT(IN) :: Ncol, Ns
      INTEGER(I4B), POINTER    :: Gitem (:)
      INTEGER(I4B), POINTER    :: ISEND (:), ISEND_DUM(:)
      INTEGER(I4B), INTENT(IN) :: Array1D(:)
      INTEGER(I4B)             :: i, j
      INTEGER(I4B)             :: icou, nnum

      if ( associated(ISEND) ) deallocate ( ISEND )
      allocate (ISEND (Ns*Ncol))

      icou = 0
      do i = 1, Ns !..... グリッド
        icou = icou + 1
        nnum = Gitem(i)
        do j = 1, Ncol  !..... 時系列
          ISEND( (icou-1)*Ncol + j ) = Array1d ((nnum-1)*Ncol + j )
        end do
      end do
      RETURN
      END SUBROUTINE set_i1dsend

!C***
!C*** set_r2dsend
!C***
      SUBROUTINE set_r2dsend (Gitem,Ncol,Ns,RSEND,Array2D)

      INTEGER(I4B), INTENT(IN) :: Ncol, Ns
      REAL    (DP), POINTER    :: RSEND (:)
      INTEGER(I4B), POINTER    :: Gitem (:)
      REAL    (DP), intent(IN) :: Array2D(1:,1:)
      INTEGER(I4B)             :: i, icol
      INTEGER(I4B)             :: icou, nnum
      integer(I4B)             :: na, nb, ndum


      na=Ncol;nb=size(Array2d,2) !時系列サイズ
      ndum = assert_eq2(na, nb,'set_r2dsend') !.. check

      if (associated(RSEND)) deallocate ( RSEND )
      allocate (RSEND (Ns*Ncol))

      icou = 0
      do i = 1, Ns
        icou = icou + 1
        nnum = Gitem(i)       !--- グリッド番号
        do icol = 1, Ncol      !--- 時系列
          RSEND((icou-1)*Ncol+icol)=Array2d(nnum,icol)
        end do
      end do
      RETURN
      END SUBROUTINE set_r2dsend

!C***
!C*** set_r1dsend
!C***
      SUBROUTINE set_r1dsend (Gitem,Ncol,Ns,RSEND,Array1D)

      INTEGER(I4B), INTENT(IN) :: Ncol, Ns
      REAL    (DP), POINTER    :: RSEND(:)
      INTEGER(I4B), pointer    :: Gitem (:)
      REAL    (DP), intent(IN) :: Array1D(:)
      INTEGER(I4B)             :: i, icol
      INTEGER(I4B)             :: icou, nnum

      if (associated(RSEND)) deallocate ( RSEND )
      allocate (RSEND (Ns*Ncol))

      icou = 0
      do i = 1, Ns
        icou = icou + 1
        nnum = Gitem(i)
        do icol = 1, Ncol      !---- 時系列
          RSEND((icou-1)*Ncol+icol) = Array1d((nnum-1)*Ncol+icol)
        end do
      end do
      RETURN
      END SUBROUTINE set_r1dsend











!C***
!C*** SelectionSort old ( miss : Item is broken )
!C***
      SUBROUTINE SelectionSort  (Item, O2N, N2O)

      INTEGER,dimension(:),INTENT(INOUT)::Item,O2N,N2O
      INTEGER                           ::NumItems,SmallestItem,I
      INTEGER                           ::LocationSmallest, idummy
      INTEGER,DIMENSION(1)              ::MINLOC_ARRAY

      NumItems = size ( Item)

      DO I = 1, NumItems
        N2O(I) = I
      ENDDO


      DO I = 1, NumItems -1

        SmallestItem = MINVAL(Item(I:NumItems))
        MINLOC_array = MINLOC(Item(I:NumItems))
        LocationSmallest = (I-1) + MINLOC_array (1)

        idummy = N2O(LocationSmallest)

        Item(LocationSmallest) = Item(I)
        N2O (LocationSmallest) = N2O (I)
        Item(I)                = SmallestItem
        N2O (I)                = idummy
      ENDDO

      do i= 1, NumItems
        O2N(N2O(I)) = I
      enddo
      RETURN
      END SUBROUTINE SelectionSort

!C***
!C*** SelectionSorti new (only make O2N, N2O)
!C***
      SUBROUTINE SelectionSorti  (Item, O2N, N2O)

      INTEGER,dimension(:),INTENT(IN)   ::Item
      INTEGER,dimension(:),INTENT(INOUT)::O2N,N2O
      INTEGER                           ::NumItems,SmallestItem,I
      INTEGER                           ::LocationSmallest, idummy
      INTEGER,DIMENSION(1)              ::MINLOC_ARRAY
      INTEGER,dimension(size(Item))     :: Item_cp

      Item_cp = Item

      NumItems = size ( Item_cp)

      DO I = 1, NumItems
        N2O(I) = I
      ENDDO


      DO I = 1, NumItems -1

        SmallestItem = MINVAL(Item_cp(I:NumItems))
        MINLOC_array = MINLOC(Item_cp(I:NumItems))
        LocationSmallest = (I-1) + MINLOC_array (1)

        idummy = N2O(LocationSmallest)

        Item_cp(LocationSmallest) = Item_cp(I)
        N2O (LocationSmallest)    = N2O (I)
        Item_cp(I)                = SmallestItem
        N2O (I)                   = idummy
      ENDDO

      do i= 1, NumItems
        O2N(N2O(I)) = I
      enddo
      RETURN
      END SUBROUTINE SelectionSorti

!C***
!C*** SelectionSortr (only make O2N, N2O)
!C***
      SUBROUTINE SelectionSortr  (rItem, O2N, N2O)

      REAL(DP),dimension(:),INTENT(IN)   ::rItem
      INTEGER ,dimension(:),INTENT(INOUT)::O2N,N2O
      REAL(DP)                           ::SmallestrItem
      INTEGER                            ::NumItems,I
      INTEGER                            ::LocationSmallest, idummy
      INTEGER,DIMENSION(1)               ::MINLOC_ARRAY
      REAL(DP), dimension(size(rItem))   :: rItem_cp

      rItem_cp = rItem

      NumItems = size ( rItem_cp )

      DO I = 1, NumItems
        N2O(I) = I
      ENDDO


      DO I = 1, NumItems -1

        SmallestrItem = MINVAL(rItem_cp(I:NumItems))
         MINLOC_array = MINLOC(rItem_cp(I:NumItems))
        LocationSmallest = (I-1) + MINLOC_array (1)

        idummy = N2O(LocationSmallest)

        rItem_cp(LocationSmallest) = rItem_cp(I)
        N2O     (LocationSmallest) = N2O (I)
        rItem_cp(I)                = SmallestrItem
        N2O (I)                    = idummy
      ENDDO

      do i= 1, NumItems
        O2N(N2O(I)) = I
      enddo
      RETURN
      END SUBROUTINE SelectionSortr
!C***
!C*** assert_eq2 : PUBLIC
!C***
      FUNCTION assert_eq2(n1,n2,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2
      INTEGER :: assert_eq2
      if (n1==n2) then
        assert_eq2=n1
      else
        write(*,*) 'nrerror: assert_eq failed with this tag:',&
                  string
       STOP 'program terminated by assert_eq2'
      end if
      END FUNCTION assert_eq2

!C***
!C*** pallocate_iv : PUBLIC
!C***
      FUNCTION pallocate_iv(p,n1,n2)
      INTEGER(I4B), DIMENSION(:),POINTER :: p
      INTEGER(I4B), DIMENSION(:,:),POINTER :: pallocate_iv
      INTEGER(I4B), INTENT(IN) :: n1, n2
      INTEGer(I4B) :: ierr, i1, i2, n, na, nb
      allocate(pallocate_iv(n1,n2), stat=ierr)
      if( ierr/=0) then
      write(0,*) n1, n2
      call &
      nrerror_abort('pallocate_iv:problem in attempt to allocate memory')
      endif
      if (.not. associated (p) ) RETURN
      na=size(p)
      nb=n1*n2
!     n=assert_eq2(size(p),n1*n2,'pallocate_iv')
      n=assert_eq2(na,nb,'pallocate_iv')
      do i2 = 1, n2
        do i1 = 1, n1
          pallocate_iv(i1,i2)=p( i1+(i2-1)*n1 )
        enddo
      enddo
      deallocate(p)
      END FUNCTION pallocate_iv

!C***
!C*** pallocate_rv : PUBLIC
!C***
      FUNCTION pallocate_rv(p,n1,n2)
      REAL(DP), DIMENSION(:),POINTER :: p
      REAL(DP), DIMENSION(:,:),POINTER :: pallocate_rv
      INTEGER(I4B), INTENT(IN) :: n1, n2
      INTEGer(I4B) :: ierr, i1, i2, n, na, nb
      allocate(pallocate_rv(n1,n2), stat=ierr)
      if( ierr/=0) then
      write(0,*) n1, n2
      call nrerror_abort                                           &
     &       ('pallocate_rv:problem in attempt to allocate memory')
      endif
      if (.not. associated (p) ) RETURN
      na=size(p)
      nb=n1*n2
!     n=assert_eq2(size(p),n1*n2,'pallocate_rv')
      n=assert_eq2(na,nb,'pallocate_rv')
      do i2 = 1, n2
        do i1 = 1, n1
          pallocate_rv(i1,i2)=p( i1+(i2-1)*n1 )
        enddo
      enddo
      deallocate(p)
      END FUNCTION pallocate_rv
!C***
!C*** nrerror_abort : PUBLIC
!C***
      SUBROUTINE nrerror_abort(string)
      !USE mpi, only : mpi_abort    ! for mpif90
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER(I4B)                 :: ierr
      write(0,*) 'nrerror: ',string
      write(0,*) 'program will terminate by nrerror_abort'
      call mpi_abort (MPI_comm_world, 100, ierr)
      END SUBROUTINE nrerror_abort
!C***
!C*** nrerror : PUBLIC
!C***
      SUBROUTINE nrerror(string)
      !USE mpi, only : mpi_finalize ! for mpif90
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER(I4B)                 :: ierr
      write(0,*) 'nrerror: ',string
      call mpi_finalize (ierr)
      STOP 'program terminated by nrerror'
      END SUBROUTINE nrerror
!C***
!C*** reallocate_rv : PUBLIC
!C***
      FUNCTION reallocate_rv(p,n)
      REAL(DP), DIMENSION(:),POINTER :: p, reallocate_rv
      INTEGER(I4B), INTENT(IN) :: n
      INTEGEr(I4B) :: nold, ierr
      INTEGER(I4B) :: num1
      allocate(reallocate_rv(n), stat=ierr)
      if( ierr/=0) call &
      nrerror('reallocate_rv: problem in attempt to allocate memory')
      if (.not. associated (p) ) RETURN
      nold=size(p)
      num1=min(nold,n)
      reallocate_rv(1:num1)=p(1:num1)
      deallocate(p)
      END FUNCTION reallocate_rv

!C***
!C*** reallocate_iv : PUBLIC
!C***
      FUNCTION reallocate_iv(p,n)
      INTEGER(I4B), DIMENSION(:),POINTER :: p, reallocate_iv
      INTEGER(I4B), INTENT(IN) :: n
      INTEGer(I4B) :: nold, ierr
      INTEGER(I4B) :: num1
      allocate(reallocate_iv(n), stat=ierr)
      if( ierr/=0) call &
      nrerror('reallocate_iv: problem in attempt to allocate memory')
      if (.not. associated (p) ) RETURN
      nold=size(p)
      num1=min(nold,n)
      reallocate_iv(1:num1)=p(1:num1)
      deallocate(p)
      END FUNCTION reallocate_iv
      END MODULE visualizer
