!C +-------------------------------------------------------------------+
!C | PROGRAM:partition_supermatrix(simple 1D division in row direction)|
!C |                                                                   |
!C |  divide "SUPERMATRIX.DAT" to "dfullmat.####" & "rkmat.###"        |
!C |                                                                   |
!C |  write all FULLmatrices & RKmatrices                              |
!C |  for large version                   on 2012/03/07 on assim       |
!C |                                                                   |
!C |                            LAST CHANGE: 2013/06/20 13:22:20.      |
!C +-------------------------------------------------------------------+
!C****
!C**** module : precisions
!C****
      MODULE     precisions
      implicit none
      integer, parameter :: KP=4
      integer, parameter :: DP=8
      END MODULE precisions

!C***
!C*** Module : Constant
!C***
      Module Constant
      USE precisions
      implicit none
!*****************************************************************
      integer(KP), parameter :: NMAX    = 1228800
      integer(KP), parameter :: NPROC   =      64 
!*****************************************************************
      integer(KP), parameter :: NBMAX   = NMAX/NPROC

      integer(KP), parameter :: nrk_inc_buf = 20
      integer(KP), parameter ::  nf_inc_buf = 20

      logical                :: DEBUG   = .FALSE.
      End Module Constant

!C***
!C*** Module : flops
!C***
      Module flops
      USE precisions
      USE constant
      IMPLICIT NONE
      integer(KP) :: nfflop = 0
      integer(KP) :: nrflop = 0
      integer(KP) :: nfflop_pe (NPROC) = 0
      integer(KP) :: nrflop_pe (NPROC) = 0
      End Module flops

!C***
!C*** Module : Hmatrices
!C***
      Module Hmatrices
      USE precisions
      implicit none
      Type fullmatrix
        Integer(KP)          :: row_offset_f
        Integer(KP)          :: col_offset_f
        Integer(KP)          :: rows_f
        Integer(KP)          :: cols_f
        Real   (DP), Pointer :: el_f(:)
      End Type fullmatrix

      Type rkmatrix
        Integer(KP)          :: row_offset_r
        Integer(KP)          :: col_offset_r
        Integer(KP)          :: rows_r
        Integer(KP)          :: cols_r
        Integer(KP)          :: k_r
        Real   (DP), Pointer :: el_a(:)
        Real   (DP), Pointer :: el_b(:)
      End Type rkmatrix

      Type(fullmatrix), pointer :: fulm    (:) => null()
      Type(rkmatrix)  , pointer ::  rkm    (:) => null()
      Type(fullmatrix), pointer :: fulm_dum(:) => null()
      Type(rkmatrix)  , pointer ::  rkm_dum(:) => null()

      Integer(KP)               :: krmax, kfmax

      Type div_data
        integer :: krmax=0
        integer :: kfmax=0
        integer :: N2O_full(100000) !TMP
        integer :: N2O_rk  (204083) !TMP
        type (fullmatrix), pointer :: fulm (:) => null()
        type (  rkmatrix), pointer ::  rkm (:) => null()
      end type div_data

      TYPE(div_data), pointer :: rank_data(:) => null()

      End Module Hmatrices

!C****
!C**** module : util
!C****
      MODULE util
      USE precisions
      USE Hmatrices
      implicit none
      INTERFACE reallocate
        MODULE PROCEDURE reallocate_fulm, reallocate_rkm
      END INTERFACE
      CONTAINS
!C****
        FUNCTION reallocate_fulm(fulm, nsize)
        type(fullmatrix), pointer ::            fulm(:)
        type(fullmatrix), pointer :: reallocate_fulm(:)
        integer(kp)               :: nsize, nsize_old, nsize_new, i

        if ( .not.associated(fulm) ) then
          allocate (reallocate_fulm(nsize))
          RETURN
        endif
        nsize_old=size(fulm, dim=1)
        if ( nsize_old < nsize ) then
          nsize_new=nsize
          nsize_old=nsize_old
        else
          nsize_new=nsize
          nsize_old=nsize
        endif
        allocate(reallocate_fulm(nsize))
        do i=1, nsize_old
          reallocate_fulm(i)%row_offset_f =  fulm(i)%row_offset_f
          reallocate_fulm(i)%col_offset_f =  fulm(i)%col_offset_f
          reallocate_fulm(i)%rows_f       =  fulm(i)%rows_f
          reallocate_fulm(i)%cols_f       =  fulm(i)%cols_f
          reallocate_fulm(i)%el_f         => fulm(i)%el_f
        enddo
        deallocate(fulm)
        END FUNCTION reallocate_fulm
!C****
        FUNCTION reallocate_rkm(rkm, nsize)
        type(rkmatrix), pointer ::            rkm(:)
        type(rkmatrix), pointer :: reallocate_rkm(:)
        integer(kp)             :: nsize, nsize_old, nsize_new, i

        if ( .not.associated(rkm) ) then
          allocate (reallocate_rkm(nsize))
          return
        endif
        nsize_old=size( rkm, dim=1)
        if ( nsize_old<nsize ) then
          nsize_new=nsize
          nsize_old=nsize_old
        else
          nsize_new=nsize
          nsize_old=nsize
        endif
        allocate(reallocate_rkm(nsize))
        do i=1, nsize_old
          reallocate_rkm(i)%row_offset_r =  rkm(i)%row_offset_r
          reallocate_rkm(i)%col_offset_r =  rkm(i)%col_offset_r
          reallocate_rkm(i)%rows_r       =  rkm(i)%rows_r
          reallocate_rkm(i)%cols_r       =  rkm(i)%cols_r
          reallocate_rkm(i)%k_r          =  rkm(i)%k_r
          reallocate_rkm(i)%el_a         => rkm(i)%el_a
          reallocate_rkm(i)%el_b         => rkm(i)%el_b
        enddo
        deallocate(rkm)
        END FUNCTION reallocate_rkm
      END MODULE util

!C****
!C**** program main
!C****
      program partition_supermatrix
      USE precisions
      USE Constant
      USE Hmatrices
      USE flops
      implicit none
      character(LEN=30) :: fname
      integer(KP)       :: kf, kr
      integer(KP)       :: istart,iend

      integer(KP)       :: i, j, icou, index
      integer(KP), pointer :: rows_f(:)=>null()
      integer(KP), pointer :: cols_f(:)=>null()
      integer(KP), pointer :: row_offset_f(:)=>null()
      integer(KP), pointer :: col_offset_f(:)=>null()
      integer(KP), pointer :: rows_r(:)=>null()
      integer(KP), pointer :: cols_r(:)=>null()
      integer(KP), pointer :: row_offset_r(:)=>null()
      integer(KP), pointer :: col_offset_r(:)=>null()
      integer(KP), pointer :: k_r(:)=>null()

      istart = 0
      iend   = NMAX-1

      kf=0
      kr=0

      allocate(rank_data(NPROC))

      Write(fname,'(A)') "supermatrix.dat"
      open(12,file=fname, status='old')
      Read(12,*)
      Call Read_supermatrix ( istart, iend, kf, kr )
      Close(12)

      kfmax=kf
      krmax=kr

      print *,"end read supermatrix !"
      print *,"independent full-matrices =",kfmax
      print *,"independent   Rk-matrices =",krmax
!DBG  print *, 'total FLOPS', nfflop, nrflop, nfflop+nrflop

!C
!C +------------------------------------------+
!C | write operation numbers for each process |
!C +------------------------------------------+
!C====
      j=0
      open(3, file='operation_conunt_etc.out', status='unknown')
      do i = 1, NPROC
          write (3,'(6i10)')                                            &
     &     i, nfflop_pe(i), nrflop_pe(i), nfflop_pe(i)+nrflop_pe(i),    &
     &     rank_data(i)%kfmax, rank_data(i)%krmax

          j=j+nfflop_pe(i)+nrflop_pe(i)
      enddo
      close(3)
!C====

!C
!C +-----------------------------------------+
!C | set FULL & RK matrices for each process |
!C +-----------------------------------------+
!C====
      do i = 1, NPROC
        !C-- set full matrix for i-th process
        icou=rank_data(i)%kfmax
        allocate (rank_data(i)%fulm(icou))
        do j = 1, icou
          index=rank_data(i)%N2O_full(j)
          rank_data(i)%fulm(j)=fulm(index)
        enddo

        !C-- set rk matrix for i-th process
        icou=rank_data(i)%krmax
        allocate (rank_data(i)%rkm(icou))
        do j=1, icou
          index=rank_data(i)%N2O_rk(j)
          rank_data(i)%rkm(j)=rkm(index)
        enddo
      enddo
!DBG  print *, rank_data(1)%fulm(1)%el_f(1)
!DBG  print *, rank_data(1)%rkm (1)%el_a(1)
!DBG  print *, rank_data(1)%rkm (1)%el_b(1)
!C====

      icou=0
      do i=1,NPROC
        do j=1, rank_data(i)%krmax
          icou=icou+rank_data(i)%rkm(j)%k_r
        enddo
        print *, 'NPROC=', i, 'iloop=',icou
      enddo
     !STOP

!C
!C +-----------------------------------------+
!C | set FULL & RK matrices for each process |
!C +-----------------------------------------+
!C====
      do i=1, NPROC
        write(*,'(A)') '    '
        write(*,'(A)') '===='
        !C-- FULL_MATRICES
        write(*,'(A,i4,A)') 'DATA for',i-1,'-th MPI process'
        icou=rank_data(i)%kfmax
        allocate(rows_f      (icou))
        allocate(cols_f      (icou))
        allocate(row_offset_f(icou))
        allocate(col_offset_f(icou))

        !C-- set integer data associated with "rank_data(i)%fulmat"
        !C-- to 1D arrays("rows_f","cols_f","row_offset_f","col_offset_f")
        do j = 1, icou
          rows_f      (j)=rank_data(i)%fulm(j)%rows_f
          cols_f      (j)=rank_data(i)%fulm(j)%cols_f
          row_offset_f(j)=rank_data(i)%fulm(j)%row_offset_f
          col_offset_f(j)=rank_data(i)%fulm(j)%col_offset_f
        enddo

        !C-- write binary file for full matrices
        write(fname, '(a9,i4.4)') 'dfullmat.', i-1
        open (88, file=fname, form='unformatted',                       &
     &            convert="little_endian", status='unknown')
        write(88) icou
        write(88) rows_f
        write(88) cols_f
        write(88) row_offset_f
        write(88) col_offset_f
        do j=1, icou
          write(88) rank_data(i)%fulm(j)%el_f
        enddo
        close(88)
        deallocate(rows_f)
        deallocate(cols_f)
        deallocate(row_offset_f)
        deallocate(col_offset_f)
        write(*,'(A,i6)') "            Full-matrices =",icou
!       write(*,'(A)') ' End Full matrices data'

        !C-- RK_MATRICES
        icou=rank_data(i)%krmax
        allocate(rows_r(icou))
        allocate(cols_r(icou))
        allocate(row_offset_r(icou))
        allocate(col_offset_r(icou))
        allocate(k_r(icou))

        !C-- set integer data associated with "rank_data(i)%rkmat"
        !C-- to 1D arrays("rows_r","cols_r","row_offset_r","col_offset_r","k_r")
        do j = 1, icou
          rows_r      (j)=rank_data(i)%rkm(j)%rows_r
          cols_r      (j)=rank_data(i)%rkm(j)%cols_r
          row_offset_r(j)=rank_data(i)%rkm(j)%row_offset_r
          col_offset_r(j)=rank_data(i)%rkm(j)%col_offset_r
          k_r         (j)=rank_data(i)%rkm(j)%k_r
        enddo

        !C-- write binary file for rk matrices
        write(fname, '(a7,i4.4)') 'drkmat.', i-1
        open (88, file=fname, form='unformatted',                       &
     &            convert="little_endian", status='unknown')
        write(88) icou
        write(88) rows_r
        write(88) cols_r
        write(88) row_offset_r
        write(88) col_offset_r
        write(88) k_r
        do j=1, icou
          write(88) rank_data(i)%rkm(j)%el_a
          write(88) rank_data(i)%rkm(j)%el_b
        enddo
        close(88)
        deallocate(rows_r)
        deallocate(cols_r)
        deallocate(row_offset_r)
        deallocate(col_offset_r)
        deallocate(k_r)
        write(*,'(A,i6)') "              Rk-matrices =",icou
!       write(*,'(A)') ' End RK   matrices data'
      enddo
      STOP
      END PROGRAM partition_supermatrix

!C***
!C*** subroutine : read_supermatrix
!C***
      Recursive Subroutine Read_supermatrix                             &
     & ( istart, iend, kf, kr )
      USE precisions
      Use Hmatrices
      Use util
      USE Constant, only: nrk_inc_buf, nf_inc_buf, NBMAX
      USE flops

      implicit none

      character (len= 1)    :: c1dam
      character (len= 2)    :: c2dam
      character (len= 4)    :: c4dam
      character (len=10)    :: c10dam
      Integer   (KP)        :: row_offset,col_offset
      Integer   (KP)        :: rows,cols
      Integer   (KP)        :: block_rows,block_cols
      Integer   (KP)        :: istart,iend,kf,kr,rs,re,ityp,num
      Integer   (KP)        :: i,j,k,kt
      Real(DP)              :: dam
      Real(DP), allocatable :: eltmp(:,:)
      Integer               :: ninc_rk, ninc_f
      character(LEN=15)     :: fn, rkn
      integer(KP)           :: npe, icou
      integer(KP)           :: ipe,istart_pe,iend_pe,nrs,nre,ncs,nce

      npe=size(nfflop_pe)
!C
!C +--------------------------------------------------+
!C | dummy size for reducing call-count of reallocate |
!C +--------------------------------------------------+
!C===
      ninc_rk=nrk_inc_buf
      ninc_f = nf_inc_buf
!C===
      Read (12,*) c4dam,ityp

      If (ityp  == 3) then
        Read(12,*) c10dam,c1dam,row_offset
        Read(12,*) c10dam,c1dam,col_offset
        Read(12,*) c4dam,c1dam,rows
        Read(12,*) c4dam,c1dam,cols
        Read(12,*) c10dam,c1dam,block_rows
        Read(12,*) c10dam,c1dam,block_cols
        Do i = 1,block_rows
          Do j = 1,block_cols
            Call Read_supermatrix (istart,iend,kf,kr)
          EndDo
        EndDo
      ElseIf (ityp == 2)then
        Read(12,*) c10dam,c1dam,row_offset
        Read(12,*) c10dam,c1dam,col_offset
        Read(12,*) c4dam,c1dam,rows
        Read(12,*) c4dam,c1dam,cols

        rs=row_offset
        re=row_offset+rows-1
        If ((re < istart) .or. (rs > iend)) then
          Do i=1,rows*cols
            Read(12,*) dam
          EndDo
        Else
          !C--- remember the full-matrix
          kf=kf+1
          if (.not.associated (fulm).or.size(fulm, dim=1) < kf) then
            fulm_dum => reallocate(fulm    , kf+ninc_f) !MH
            fulm     => reallocate(fulm_dum, kf+ninc_f) !MH
          endif
          fulm(kf)%row_offset_f=row_offset
          fulm(kf)%col_offset_f=col_offset
          fulm(kf)%rows_f=rows
          fulm(kf)%cols_f=cols

          !C--- get the size of matrix
          num=rows*cols
          Allocate(eltmp(0:rows-1,0:cols-1))
          Allocate(fulm(kf)%el_f(0:num-1))
          Do i=0,cols-1
            Do j=0,rows-1
              Read(12,*) eltmp(j,i)
            EndDo
          EndDo

          num=0
          Do i=0,rows-1
            Do j=0,cols-1
              fulm(kf)%el_f(num)=eltmp(i,j)
              num=num+1
            EndDO
          EndDo
          Deallocate(eltmp)

          nfflop=nfflop+rows*cols

          !C--- set index "rank_data(ipe)%N2O_full"
          do ipe = 1, npe
            istart_pe= (ipe-1)*NBMAX
            iend_pe  = (ipe  )*NBMAX - 1
            if (iend_pe   < row_offset       ) cycle
            if (istart_pe > row_offset+rows-1) cycle

            rank_data(ipe)%kfmax=rank_data(ipe)%kfmax+1
            icou=rank_data(ipe)%kfmax
            rank_data(ipe)%N2O_full(icou)=kf

            nrs=max(istart_pe, row_offset)
            nre=min(  iend_pe, row_offset+rows-1)
            ncs=col_offset
            nce=col_offset+cols-1
            nfflop_pe(ipe)=nfflop_pe(ipe)+(nre-nrs+1)*(nce-ncs+1)
!DBG        print *, 'ful',nrs, nre, ncs, nce
          enddo

!DBG      write(22,*) kf, row_offset, col_offset, rows, cols
        EndIf
      ElseIf(ityp == 1)then
        Read(12,*) c10dam,c1dam,row_offset
        Read(12,*) c10dam,c1dam,col_offset
        Read(12,*) c4dam,c1dam,rows
        Read(12,*) c4dam,c1dam,cols
        Read(12,*) c1dam,c1dam,k
        Read(12,*) c2dam,c1dam,kt
        rs=row_offset
        re=row_offset+rows-1
        If ((re < istart) .or. (rs > iend)) then
          Do i=0,k*rows-1
            Read(12,*) dam
          EndDo
          Do i=0,k*cols-1
            Read(12,*) dam
          EndDo
        Else
          !C--- remember the Rk-matrix
          kr=kr+1
          if (.not.associated (rkm).or.size(rkm, dim=1) < kr) then
            rkm_dum => reallocate(rkm    , kr+ninc_rk) !MH
            rkm     => reallocate(rkm_dum, kr+ninc_rk) !MH
          endif
          rkm(kr)%row_offset_r=row_offset
          rkm(kr)%col_offset_r=col_offset
          rkm(kr)%rows_r=rows
          rkm(kr)%cols_r=cols
          rkm(kr)%k_r=k
!DBG      write(23,*) kr, row_offset, col_offset, rows, cols, k

          !C--- get the size of vector Ai(i=1,--,k)
          num=k*rows
          Allocate(rkm(kr)%el_a(0:num-1))
          Do i=0,num-1
            Read(12,*) rkm(kr)%el_a(i)
          EndDo

          !C--- get the size of vector Bi(i=1,--,k)
          num=k*cols
          Allocate(rkm(kr)%el_b(0:num-1))
          Do i=0,num-1
            Read(12,*) rkm(kr)%el_b(i)
          EndDo

          nrflop=nrflop+(rows+cols)*k

          !C--- set index "rank_data(ipe)%N2O_rk"
          do ipe = 1, npe
            istart_pe= (ipe-1)*NBMAX
            iend_pe  = (ipe  )*NBMAX - 1
            if (iend_pe   < row_offset       ) cycle
            if (istart_pe > row_offset+rows-1) cycle

            rank_data(ipe)%krmax=rank_data(ipe)%krmax+1
            icou=rank_data(ipe)%krmax
            rank_data(ipe)%N2O_rk(icou) =kr

            nrs=max(istart_pe, row_offset)
            nre=min(  iend_pe, row_offset+rows-1)
            ncs=col_offset
            nce=col_offset+cols-1
            nrflop_pe(ipe)=nrflop_pe(ipe)+((nre-nrs+1)+(nce-ncs+1))*k
!DBG        print *, 'rkm',nrs, nre, ncs, nce
          enddo
        EndIf
      EndIf
      End Subroutine Read_supermatrix
