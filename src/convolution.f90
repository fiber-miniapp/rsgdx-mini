!C
!C   Created by ʼƣ �� on 11/02/17.
!C   Last Change: 2012/09/13 20:25:21.
!C
!C +-----------------------------------------+
!C | ��deallocation��˺��Ƥ������������    |
!C |   deallocate ( gfd)                     |
!C |   deallocate (olat, olon)               |
!C |   deallocate (INDEX_SELECTOBS)          |
!C |   deallocate (slp_array, def_array )    |
!C |   deallocate (us, ud)                   |
!C |                   on 2011ǯ  2�� 23��   |
!C | ��������˥Ρ����ֹ�(NthG)���դ���      |
!C |                   on 2011ǯ  2�� 25��   |
!C | ��ñ�η׻��Τ���iii,NthG�����          |
!C |                   on 2012ǯ  6�� 11��   |
!C | ��������˥Ρ����ֹ�(NthG)���դ���(��)  |
!C |   ������     "deformation_###.out"      |
!C |                   on 2012ǯ  6�� 12��   |
!C | �������ʬ���Ф����ѷ���׻�����褦��  |
!C |                   on 2012ǯ  8�� 24��   |
!C | ��NaN���򤱤뤿��convolution__conv��    |
!C |   tinc=0�ΰ�������on 2012ǯ  9�� 12��   |
!C | ����¬���Ѱ̤�ñ���󥯤˽����Ф����Τ�|
!C |  MPI_ALLREDUCE�ǤϤʤ�MPI_REDUCE       |
!C |                   on 2012ǯ  9�� 12��   |
!C |                                         |
!C | STAGE-IN�Υ��顼�ǥե����뤬���Ĥ���ʤ�|
!C | �Ȥ�����륨�顼��ȯ��������������    |
!C | �� input_deform.nml=>all                |
!C | �� GFD_SRF         =>all                |
!C | �� ��¬����ɸ�ꥹ��=>nroot_rank         |
!C | �� ��¬������ꥹ��=>all                |
!C |                   on 2012ǯ  9�� 13��   |
!C +-----------------------------------------+
!C

!C***
!C*** MODULE convolution
!C***
      MODULE convolution

      USE m_conv
      USE m_comm !... comms_new, NthG
      IMPLICIT NONE
      include 'precision_pri.inc'
      private
      public :: convolution__init
      public :: convolution__conv
      public :: convolution__writeTS

      !C-- called count of "convolution__conv"
      integer(I4B), private   :: icou       = 0
      integer(I4B), parameter :: nroot_rank = 0 !20120912

      CONTAINS
!C***
!C*** initialization
!C***
      SUBROUTINE convolution__init(NBMAX, NTIME)!NTIME is added on 20120824
      integer(I4B), intent(IN) :: NBMAX, NTIME
      integer(I4B)             :: nobs2, myid
      integer(I4B)             :: ksrc,kobs
      integer(I4B)             :: islct
      character(len=100)       :: gfdnameall
      integer(I4B)             :: OpenStatus ! 12Sep12

      namelist /obsparams/ NOBS
      namelist /fnames   / gfdname, obsfname, selectfname

      myid=comms_new(NthG)%myrank
!C
!C +------------------------------------------+
!C | "input_deform.nml"�򤹤٤ƤΥ�󥯤��ɹ� |
!C +------------------------------------------+
!C===
      open ( inpobsfcode, file='input_deform.nml', status='old',        &
     &                                       iostat=OpenStatus  )
      if (OpenStatus/=0) then
        print *,"input_deform.nml open error :",NthG, myid
      endif
      read ( inpobsfcode, nml = obsparams )
      read ( inpobsfcode, nml = fnames    )
      close( inpobsfcode )
!C===

!C
!C +-------------+
!C | ALLOCATAION |
!C +-------------+
!C===
        allocate ( gfd(3,NOBS, NBMAX)) ; gfd =0.0e0_dp

      if ( myid == nroot_rank ) then
        allocate (olat(  NOBS       )) ; olat=0.0e0_dp
        allocate (olon(  NOBS       )) ; olon=0.0e0_dp
!CO     allocate (oX  (  NOBS       )) ; oX  =0.0e0_dp
!CO     allocate (oY  (  NOBS       )) ; oY  =0.0e0_dp
      endif
!C===

!C
!C +----------------------------+
!C | �Ѱ�GF���ɤ߹���(�����) |
!C +----------------------------+
!C===
      write (gfdnameall,'(i3.3,a,i4.4)') NthG-1, trim(gfdname), myid
      open  ( gfdcode, FILE=trim(gfdnameall ),                          &
     &           status='old', form='unformatted', iostat=Openstatus )
      if (OpenStatus/=0) then
        print *,"gfd_file open error :", NthG, myid
      endif
      do ksrc=1, NBMAX
         read(gfdcode) gfd(:,:,ksrc)
      end do
      close ( gfdcode )
!C===

!C
!C +----------------------------------------------+
!C | ��¬����ɸ���ɤ߹���(GFD�б�:nroot_rank�Τ�) |
!C +----------------------------------------------+
!C===
      if ( myid == nroot_rank ) then
        open  ( obsfcode, FILE=trim(obsfname), status='old',            &
     &                                   iostat=OpenStatus )
        if (OpenStatus /= 0) then
          print *,"obsfile open error :nroot_rank", NthG, myid
        endif
        read  ( obsfcode, *) nobs2
        if (nobs2/=NOBS) STOP "nobs2 must be equal to NOBS"
        do kobs = 1, NOBS
          read( obsfcode, *) olat(kobs), olon(kobs)
        end do
        close ( obsfcode )
      endif
!C===

!C
!C +----------------------------------------------+
!C | ���Ϥ����¬����ɸ����ǥå���(�����)�ɹ� |
!C +----------------------------------------------+
!C===
      open  ( selectfcode,FILE=trim(selectfname), status='old',         &
     &                                      iostat=OpenStatus )
      if (OpenStatus/=0) then
        print *,"select_file open error :", NthG, myid
      endif
      read  ( selectfcode, * ) NSELECT
      allocate (INDEX_SELECTOBS(NSELECT))
      do islct = 1, NSELECT
        read( selectfcode, *) INDEX_SELECTOBS(islct)
      end do
      close(selectfcode)
!C===

!C
!C +-------------+
!C | ALLOCATAION |
!C +-------------+
!C===
!CO   allocate (slp_array(NBMAX) )
      allocate (def_array(NBMAX) )

!CO   allocate (us    (3, NSELECT)) ; us     = 0.0e0_dp
      allocate (ud    (3, NSELECT)) ; ud     = 0.0e0_dp

      if ( myid == nroot_rank ) then
        !C----for MPI_REDUCE
!CO     allocate (us_sum(3, NSELECT)) ; us_sum = 0.0e0_dp
        allocate (ud_sum(3, NSELECT)) ; ud_sum = 0.0e0_dp

        !C----for ��������ϥǡ��������� on 20120824
!CO     allocate (tus_sum(3,NSELECT,NTIME)) ; tus_sum = 0.0e0_dp
        allocate (tud_sum(3,NSELECT,NTIME)) ; tud_sum = 0.0e0_dp
        allocate (timey_ts         (NTIME)) ; timey_ts= 0.0e0_dp
      endif
!C===
      RETURN
      END SUBROUTINE convolution__init
!C***
!C***
!C***
!CO   SUBROUTINE convolution__conv(NBMAX,V, V_Vpl,     timey) !miss
      SUBROUTINE convolution__conv(NBMAX,   V_Vpl,tinc,timey)
      include 'mpif.h'
      integer(I4B), intent(IN) :: NBMAX
      real   (8)  , intent(IN) :: V_Vpl(1:NBMAX) !m/yr
!CO   real   (8)  , intent(IN) :: V    (1:NBMAX) !m/yr
      real   (8)  , intent(IN) :: tinc,timey     !year
      integer(I4B)             :: ksrc,kobs, ierr
      integer(I4B)             :: ith_select
      icou = icou+1
!CO   slp_array(1:NBMAX)=V    (1:NBMAX) ! [m/yr]
      def_array(1:NBMAX)=V_Vpl(1:NBMAX) ! [m/yr]

!C
!C +----------------------+
!C | calc (slip deficits) |
!C +----------------------+
!C===
!CO   us=0.0e0_dp
      ud=0.0e0_dp
      do ith_select = 1, NSELECT
        kobs=INDEX_SELECTOBS(ith_select)
        do ksrc=1, NBMAX
!CO       us(:,kobs)=us(:,kobs)+gfd(:,kobs,ksrc)*slp_array(ksrc)
          ud(:,kobs)=ud(:,kobs)+gfd(:,kobs,ksrc)*def_array(ksrc)
        enddo
      enddo

!CO   call MPI_REDUCE   (us , us_sum, 3*NSELECT, MPI_REAL8,             &
!CO  &                   MPI_SUM, nroot_rank,comms_new(NthG)%comm, ierr)
      call MPI_REDUCE   (ud , ud_sum, 3*NSELECT, MPI_REAL8,             &
     &                   MPI_SUM, nroot_rank,comms_new(NthG)%comm, ierr)

      if (comms_new(NthG)%myrank==nroot_rank) then
        write(cfname  ,'(A,I5.5)') 'cd_velo_',icou !--mod cfname on 120825
        open (outcode,FILE=trim(cfname), status='unknown' )
        do kobs = 1, NSELECT
          if (tinc /= 0.0) then !avoid NaN  20120912
            write(outcode,'(8(1pe20.12))')                              &
     &    olat(kobs),olon(kobs),ud_sum(:,kobs)/tinc!,us_sum(:,kobs)/tinc
          else
            write(outcode,'(8(1pe20.12))')                              &
     &    olat(kobs),olon(kobs),ud_sum(:,kobs)     !,us_sum(:,kobs)
          endif
        enddo
        close(outcode)

        if (icou==1) then
!CO       tus_sum(:,:,icou)=us_sum(:,:) !save as t.s.
          tud_sum(:,:,icou)=ud_sum(:,:) !save as t.s.
        else
!CO       tus_sum(:,:,icou)=tus_sum(:,:,icou-1)+us_sum(:,:) !save as t.s.
          tud_sum(:,:,icou)=tud_sum(:,:,icou-1)+ud_sum(:,:) !save as t.s.
        endif
        timey_ts(icou)=timey            !save as t.s.
      endif
!C===
      RETURN
      END SUBROUTINE convolution__conv

!C***
!C*** GPS��¬��������ǡ����񤭽Ф��롼����
!C***
      SUBROUTINE convolution__writeTS
      integer(I4B) :: kobs
      integer(I4B) :: isnp
      integer(I4B) :: ith_select

      if ( comms_new(NthG)%myrank == nroot_rank ) then
        write(cfname,'(A,I3.3,A)') 'deformation_',NthG-1,'.out'
        open (outcode,FILE=trim(cfname), status='unknown' )
        do ith_select=1, NSELECT
          kobs=INDEX_SELECTOBS(ith_select)
          write(outcode,'(a1,2(1x,E13.6))') "#", olat(kobs), olon(kobs)
          do isnp=1, icou
            write(outcode ,'(7(1pe20.12))')                             &
     &        timey_ts(isnp),tud_sum(:,kobs,isnp)!,tus_sum(:,kobs,isnp)
          enddo
        enddo
        close(outcode)
      endif
!C===
      RETURN
      END SUBROUTINE convolution__writeTS

      END module convolution
