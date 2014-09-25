!C
!C   Created by 兵藤 守 on 09/04/23.
!C   Copyright 2009 IFREE. All rights reserved.
!C
!C +---------------------------------------------------------+
!C | 地震性イベント時の累積辷り量                            |
!C |   ■　V=1cm/yr, V=Vplより速い滑りをともに出力 on 07/02   |
!C |   ■　対応するモーメントマグニチュードも出力  on 07/02   |
!C |   ■　イベントの前後との間隔、継続時間を出力  on 07/03   |
!C |   ■　X-深さグリッドへの対応(微小断層面積)    on 10/02   |
!C |                                                         |
!C | cell面積をgeometryn_areaより読み込む         on 12/1001 |
!C +---------------------------------------------------------+
!C

!C***
!C*** MODULE nrutil
!C***
      MODULE nrutil_cos
      IMPLICIT NONE
!     private
!     public :: reallocate
      include 'precision_pri.inc'
      INTERFACE reallocate
        MODULE PROCEDURE reallocate_iv, reallocate_sv, reallocate_dv
      END INTERFACE
      CONTAINS
      SUBROUTINE nrerror(string)
      CHARACTER(LEN=*), INTENT(IN) :: string
!CO   INTEGER(I4B)                 :: ierr
      write(0,*) 'nrerror: ',string
      write(6,*) 'nrerror: ',string
      STOP 'program terminated by nrerror'
      END SUBROUTINE nrerror

      FUNCTION reallocate_dv(p,n)
      REAL   (DP) , DIMENSION(:),POINTER :: p, reallocate_dv
      INTEGER(I4B), INTENT(IN)           :: n
      INTEGER(I4B)                       :: nold, ierr
      allocate(reallocate_dv(n), stat=ierr)
      if( ierr/=0) call                                               &
     & nrerror('reallocate_dv: problem in attempt to allocate memory')
      if (.not. associated (p) ) RETURN
      nold=size(p)
      reallocate_dv(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
      END FUNCTION reallocate_dv

      FUNCTION reallocate_sv(p,n)
      REAL   (SP) , DIMENSION(:),POINTER :: p, reallocate_sv
      INTEGER(I4B), INTENT(IN)           :: n
      INTEGER(I4B)                       :: nold, ierr
      allocate(reallocate_sv(n), stat=ierr)
      if( ierr/=0) call                                               &
     & nrerror('reallocate_sv: problem in attempt to allocate memory')
      if (.not. associated (p) ) RETURN
      nold=size(p)
      reallocate_sv(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
      END FUNCTION reallocate_sv

      FUNCTION reallocate_iv(p,n)
      INTEGER(I4B), DIMENSION(:),POINTER :: p, reallocate_iv
      INTEGER(I4B), INTENT(IN)           :: n
      INTEGER(I4B)                       :: nold, ierr
      allocate(reallocate_iv(n), stat=ierr)
      if( ierr/=0) call                                               &
     & nrerror('reallocate_iv: problem in attempt to allocate memory')
      if (.not. associated (p) ) RETURN
      nold=size(p)
      reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
      END FUNCTION reallocate_iv
      END MODULE nrutil_cos

!C***
!C*** MODULE myutil
!C***
      MODULE myutil_cos
      IMPLICIT NONE
!     private
!     public :: S2YDHMS
      include 'precision_pri.inc'
      INTERFACE S2YDHMS
        MODULE PROCEDURE S2YDHMS_sp, S2YDHMS_dp
      END INTERFACE

      CONTAINS

      SUBROUTINE S2YDHMS_sp                                             &
     &           (tot_seconds,iyears,idays,ihours,iminutes,iseconds)

      real   (SP) , intent(IN ):: tot_seconds
      INTEGER(I4B), intent(OUT):: iyears,idays,ihours,iminutes,iseconds

      real   (SP)              :: remained_seconds
      REAL   (SP)              :: ryears,rdays,rhours,rminutes,rseconds
      REAL(SP), PARAMETER      :: Uminute =    60.00e0_sp
      REAL(SP), PARAMETER      :: Uhour   =    60.00e0_sp*Uminute
      REAL(SP), PARAMETER      :: Uday    =    24.00e0_sp*Uhour
      REAL(SP), PARAMETER      :: Uyr     = 365.2425e0_sp*Uday

      remained_seconds = tot_seconds

      iyears = int (remained_seconds / Uyr  )
      ryears = real(iyears, SP)
      remained_seconds = remained_seconds - ryears*Uyr

      idays  = int (remained_seconds / Uday )
      rdays  = real(idays , SP)
      remained_seconds = remained_seconds - rdays *Uday

      ihours = int (remained_seconds / Uhour)
      rhours = real(ihours, SP)
      remained_seconds = remained_seconds - rhours*Uhour

      iminutes = int (remained_seconds / Uminute)
      rminutes = real(iminutes, SP)
      remained_seconds = remained_seconds - rminutes*Uminute

      rseconds = remained_seconds
      iseconds = int ( rseconds )
      RETURN
      END SUBROUTINE S2YDHMS_sp

      SUBROUTINE S2YDHMS_dp                                             &
     &           (tot_seconds,iyears,idays,ihours,iminutes,iseconds)

      real   (DP) , intent(IN ):: tot_seconds
      INTEGER(I4B), intent(OUT):: iyears,idays,ihours,iminutes,iseconds

      real   (DP)              :: remained_seconds
      REAL   (DP)              :: ryears,rdays,rhours,rminutes,rseconds
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
      END MODULE myutil_cos


!C***
!C*** MODULE coseismic_slip
!C***
      MODULE coseismic_slip

      USE constant
      USE m_comm

      USE nrutil_cos
      USE myutil_cos

      IMPLICIT NONE
      include 'mpif.h'
      private
      public :: eqk_slip_calc

      integer     , parameter ::  I4B  = 4
      integer(I4B), parameter ::   SP  = 4

      integer(I4B), parameter ::   STDIN   = 05
      integer(I4B), parameter ::   STDOUT  = 06
      integer(I4B), parameter :: NROOT_PROC= 0

      real   (SP) , parameter :: year  =60.e0*60.e0*24.e0*365.2425e0
      real   (SP) , parameter :: Vp0   = 6.5e0_sp/year  ![cm/s]
      real   (SP) , parameter :: rigid = 3.0e10_sp

      real(SP), allocatable, dimension(:,:) :: dtallslip
      real(SP), allocatable, dimension(:,:) ::  seisslip
      real(SP), allocatable, dimension(:,:) :: qseisslip

      real(SP), allocatable, dimension(:)   :: dtallslip_all
      real(SP), allocatable, dimension(:)   ::  seisslip_all
      real(SP), allocatable, dimension(:)   :: qseisslip_all

      real(SP), allocatable, dimension(:)   ::  slip_all
      real(SP), allocatable, dimension(:)   ::  seis_all
      real(SP), allocatable, dimension(:)   :: qseis_all

      real(SP), allocatable, dimension(:)   :: event_interval

      real(SP)    , pointer :: event_duration(:) => null()
      integer(I4B), pointer :: nstepB        (:) => null()
      integer(I4B), pointer :: nstepE        (:) => null()

!------------------mod on 2012Oct30---------------!
      real(SP)    , pointer :: event_duration_dum(:) => null()
      integer(I4B), pointer :: nstepB_dum        (:) => null()
      integer(I4B), pointer :: nstepE_dum        (:) => null()
!------------------mod on 2012Oct30---------------!

      real(SP)              :: dt_eps
      real(SP)              :: patch_len

      integer (I4B)         :: nstep_same
      integer (I4B)         :: NSTEP, NEVENT


!C=== Moment
      real(SP)::  seismic_moment, m_magnitude
      real(SP):: qseismic_moment,qm_magnitude
      real(SP),allocatable, dimension(:)  :: patch_square
!C===

      CONTAINS
!C***
!C*** input_eqk.nml     : 装置番号11 --> オープンエラー114
!C*** geometrynファイル : 装置番号13 --> オープンエラー134
!C*** nvpcos1000#       : 装置番号1
!C*** nvpcos1000#.txt   : 装置番号2
!C***
!!    call eqk_slip_calc(nsav,tsav(0:nsav-1),dtsav(0:nsav-1),       &
!!   &     vsav(0:NBMAX-1,0:nsav-1),ssav(0:NBMAX-1,0:nsav-1),ierror)
      SUBROUTINE eqk_slip_calc(nsav,tsav,dtsav,vsav,ssav,ierror)
      integer(I4B), INTENT(IN) :: nsav
      integer(I4B)             :: ierror, ierr
!C for SP version
!     real    (SP), INTENT(IN) :: tsav(:)  ,dtsav(:)
!     real    (SP), INTENT(IN) :: vsav(1:,:), ssav(1:,:)
!C for DP version
      real    (DP), INTENT(IN) :: tsav(1:nsav)  ,dtsav(1:nsav)
      real    (DP), INTENT(IN) :: vsav(1:NBMAX,1:nsav), ssav(1:NBMAX,1:nsav)

      real    (SP)             :: time_val(1:nsav)
      real    (SP)             :: time_inc(1:nsav)

      integer(I4B)             :: icount, icou
      integer(I4B)             :: nbegin, nend
      integer(I4B)             :: ibegin, i, j, k, kk
      integer(I4B)             :: ievn, it
      integer(I4B)             :: idum, openstatus
      integer (I4B)            :: i1,i2,i3,i4 !CO,i5,i6
      real    (SP)             :: r1,r2 !CO,r3,r4,r5,r6
      real    (SP)             :: dt_inc
      real    (SP)             ::  dum
      real    (SP)             :: Vseis
      real    (SP)             :: Vconv
      real    (SP)             :: normalized_Vseis
      real    (SP)             :: normalized_Vconv

      character(len=11)        :: ofname
      character(len=99)        :: geofile

      LOGICAL                  :: no_more_event_cand

      namelist /eqk/ dt_eps, nstep_same,          &
     &                   geofile         !="geometryn_area"

      Vseis            = 1.0e0_sp  ![cm/s]
      Vconv            = Vp0
      normalized_Vseis = Vseis/Vp0 ![    ]
      normalized_Vconv = Vconv/Vp0 ![    ] == 1

      do i=1, nsav
        time_val(i)=real( tsav(i)/year,SP)
        time_inc(i)=real(dtsav(i)/year,SP)
      enddo
      NSTEP = nsav
      open (11, file='input_eqk.nml', status='old', iostat=OpenStatus)
      if ( OpenStatus /= 0 ) then
        ierror=114
        RETURN
      endif
      read (11, nml=eqk )
      if ( comms_new(NthG)%myrank == NROOT_PROC ) then
        write(6,nml=eqk)
      endif
      close(11)
!C
!C +---------------+
!C | search events |
!C +---------------+
!C====
      icount = 0
      nend   = 1
      do while ( nend /= NSTEP )
        ibegin=nend
        no_more_event_cand = .true.
!C-- search begin step of 'next_event' from the end of 'current_event'
        do k=ibegin, NSTEP - 1
          if ( time_inc(k+1) < dt_eps ) then
            no_more_event_cand = .false.
            nbegin=k
            dt_inc=0.0e0_sp
!C-- search the end step of 'current_event'
            do kk=nbegin, NSTEP-1
              dt_inc = dt_inc + time_inc(kk+1)
              if ( dt_inc > dt_eps ) then
                dt_inc = dt_inc - time_inc(kk+1)
                nend    =  kk
                goto 999
              endif
            enddo
            nend=NSTEP
            go to 999
          endif
        enddo
        if ( no_more_event_cand ) go to 333
 999    continue
!C-- if a candidate consists of many steps(>nstep_same),
!C                          it is added to a member of events
        if ( nend - nbegin + 1 > nstep_same ) then
          icount = icount + 1
          if (.not.associated(nstepB).and..not.associated(nstepE)) then
            allocate (nstepB(icount),nstepE(icount))
            allocate (event_duration(icount))
          else
            nstepB_dum=>reallocate(nstepB,icount)
            nstepE_dum=>reallocate(nstepE,icount)
            event_duration_dum=>reallocate(event_duration,icount)

!------------------mod on 2012Oct30---------------!
            nstepB=>nstepB_dum
            nstepE=>nstepE_dum
            event_duration=>event_duration_dum
            nullify(nstepB_dum)
            nullify(nstepE_dum)
            nullify(event_duration_dum)
!------------------mod on 2012Oct30---------------!
          endif
          nstepB        (icount) = nbegin
          nstepE        (icount) = nend
          event_duration(icount) = dt_inc
        endif
      enddo
 333  continue
      NEVENT = icount

      if ( NEVENT == 0) then
        ierror=49
        RETURN
      endif

      allocate ( event_interval(NEVENT-1) )
      event_interval = 0.0e0

      do ievn= 1, NEVENT-1
        do k = nstepE(ievn)+1, nstepB(ievn+1)
          event_interval(ievn)=event_interval(ievn)+time_inc(k)
        enddo
      enddo

!C>-----------block only
      allocate( dtallslip (NBMAX,NEVENT),                               &
     &           seisslip (NBMAX,NEVENT),                               &
     &          qseisslip (NBMAX,NEVENT) )
                dtallslip=0.0e0_sp
                 seisslip=0.0e0_sp
                qseisslip=0.0e0_sp

      do ievn = 1, NEVENT
        do it = nstepB(ievn)+1, nstepE(ievn)
          do j = 1, NBMAX
            i=j
            dtallslip(i,ievn) = dtallslip(i,ievn) + real(ssav(j,it),SP)
            if ( real(vsav(j,it),SP) >= normalized_Vseis ) then
              seisslip(i,ievn) = seisslip(i,ievn) + real(ssav(j,it),SP)
            endif
            if ( real(vsav(j,it),SP) >= normalized_Vconv ) then
             qseisslip(i,ievn) =qseisslip(i,ievn) + real(ssav(j,it),SP)
            endif
          enddo
        enddo
      enddo
!C<-----------block only

      if ( comms_new(NthG)%myrank == NROOT_PROC ) then
        allocate(dtallslip_all(NMAX*NEVENT))
        allocate( seisslip_all(NMAX*NEVENT))
        allocate(qseisslip_all(NMAX*NEVENT))
      endif

      call MPI_GATHER(dtallslip    ,NBMAX*NEVENT, MPI_REAL,           &
     &                dtallslip_all,NBMAX*NEVENT, MPI_REAL,           &
     &                NROOT_PROC, comms_new(NthG)%comm, ierr )
      call MPI_GATHER( seisslip    ,NBMAX*NEVENT, MPI_REAL,           &
     &                 seisslip_all,NBMAX*NEVENT, MPI_REAL,           &
     &                NROOT_PROC, comms_new(NthG)%comm, ierr )
      call MPI_GATHER(qseisslip    ,NBMAX*NEVENT, MPI_REAL,           &
     &                qseisslip_all,NBMAX*NEVENT, MPI_REAL,           &
     &                NROOT_PROC, comms_new(NthG)%comm, ierr )

      deallocate(dtallslip,seisslip,qseisslip)

      if ( comms_new(NthG)%myrank == NROOT_PROC ) then
        open(13, file=trim(geofile), status='old', iostat=openstatus)
        if ( openstatus/= 0 ) then
          ierror=134
          RETURN
        endif

        allocate (patch_square(NMAX))
        do j=1, NMAX
          read(13, *) idum, dum,dum,dum, patch_square(j), patch_len
          patch_square(j)=patch_square(j)*1000.0e0_sp*1000.0e0_sp !**in METER
        enddo
        close(13)

        allocate( slip_all(NMAX))
        allocate( seis_all(NMAX))
        allocate(qseis_all(NMAX))

        do j=1,NEVENT
          icou=0
          slip_all=0.0e0_sp
          seis_all=0.0e0_sp
         qseis_all=0.0e0_sp
          do i=1, NPEperGROUP
            do k=1, NBMAX
              kk=(i-1)*NBMAX*NEVENT+(j-1)*NBMAX+k
              icou=icou+1
              slip_all(icou)=dtallslip_all(kk)
              seis_all(icou)= seisslip_all(kk)
             qseis_all(icou)=qseisslip_all(kk)
            enddo
          enddo
          if (icou/=NMAX) STOP

          write(ofname,'(a6,i5.5)')'nvpcos',10000+j
          open(1, file=trim(ofname), status='unknown')
          do i=1,NMAX
            write(1,1002) i, slip_all(i), seis_all(i), qseis_all(i)
          end do
          close(1)

          seismic_moment=0.0e0_sp
         qseismic_moment=0.0e0_sp
          do i=1,NMAX
            seismic_moment= seismic_moment+patch_square(i)* seis_all(i)*rigid
           qseismic_moment=qseismic_moment+patch_square(i)*qseis_all(i)*rigid
          enddo
          m_magnitude=(log10( seismic_moment)-9.1e0_sp)/1.5e0_sp
         qm_magnitude=(log10(qseismic_moment)-9.1e0_sp)/1.5e0_sp

!C-- for GMT pstext
          open(2, file=trim(ofname)//'.txt', status='unknown')
!C--    moment
          i1=int( m_magnitude)
          i3=int(qm_magnitude)
          r1=( m_magnitude-i1)*100.0e0
          r2=(qm_magnitude-i3)*100.0e0
          i2=int(r1)
          i4=int(r2)
          write(2,'(a,i1.1,a,i2.2,a,i1.1,a,i2.2,a)')                    &
       &        'M@-w_seis@-=',i1,'.',i2,'(M@-w_aseis@-=',i3,'.',i4,')'
!C-- interval from previous event
          if ( j==1 ) then
            r1=0.0e0_sp
          else
            r1=event_interval(j-1)
          endif
          write(2,'(a,1pe10.4)' ) '@~D@~T@-pre@-=',r1
!C-- event_begin_time, event_duration
          r1=time_val(nstepB(j))
          r2=event_duration(j)
          write(2,'(a,i3.3,a,1pe10.4,a,i3.3,a,1pe10.4,a)' )           &
       &        'T@-',j,'@-=',r1,'(@~D@~t@-',j,'@-=', r2 , ')'
!C-- inverval to next event
          if ( j==NEVENT ) then
            r1=0.0e0_sp
          else
             r1=event_interval(j)
          endif
          write(2,'(a,1pe10.4)' ) '@~D@~T@-post@-=',r1

          close(2)
!C===
        end do
      else
        RETURN
      endif
 1002 format(I7,3(1x,E13.6))
      RETURN
      END SUBROUTINE eqk_slip_calc

      end MODULE coseismic_slip
