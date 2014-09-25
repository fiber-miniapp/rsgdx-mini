!C
!C    for calculation with various frictional parameters on 20120602(miss)
!C    < slip response function is compressed by h-matrices >
!C    < save step is controled by dt >
!C      save step data if dslip >= dslipsav
!C
!C    * add sliponly & deficit as save variables on 20120501
!C    * add "fn_rst" to namelist "input9"        on 20120531
!C
!C    * output only recent data(irest=1 for restart) on 20120808
!C
!C    * コンボリューションを速度でなく変位増分に対し実施 on 20120824
!C
!C                         LAST CHANGE: 2013/01/17 16:37:23.
!C
      program RSGDX

      USE precisions
      USE Constant
      USE timer_util
      USE m_comm
      USE m_outputval, only : vsav, ssav, thsav, sliponly, deficit
      USE gmt,         only : VovrVp

!K    USE gmt
!K    USE convolution
!K    USE coseismic_slip
!K    USE visualizer
!K    USE mpi_comm_timer

      IMPLICIT NONE
      INCLUDE 'mpif.h'
!C
!C +-----------+
!C |  for MPI  |
!C +-----------+
!C===
      integer(kp) :: np, my_rank, ierr
!C===

!C
!C +-----------------+
!C |  for time step  |
!C +-----------------+
!C===
!C--  variables
      integer(kp) :: nstep
      real   (dp) :: t,dt
      real   (dp) :: h,hdone
      real   (dp) :: V     (0:NBMAX-1)
      real   (dp) :: Theta (0:NBMAX-1)
      real   (dp) :: slip  (0:NBMAX-1)

!C--  for save data
      integer(kp) :: nsav
      integer(kp) :: nfmax,nfile
      real   (dp) :: dslpmax
      real   (dp) ::   tsav(0:NSAVMAX-1)
      real   (dp) ::  dtsav(0:NSAVMAX-1)
      real   (dp) :: hprint(0:NSAVMAX-1)
      character(len=30) :: fname
!C--  for control
      integer(kp) :: kmax
      real   (dp) :: eps,dslipsav,t2
!C--  else
      integer(kp) :: nb
      real   (dp) :: year

      integer(kp) :: nmvc
      integer(kp) :: comm_all, group_all
      integer(kp) :: ninc_gmt_write
!C===
      real(DP)    :: t_inc ! for convolution of disp. increment
      real(DP)    ::del_deficit(0:NBMAX-1)

!C
!C +--------------+
!C |  for common  |
!C +--------------+
!C===
!c--  common
      common /forMPI/ my_rank
      common /control/ eps,dslipsav,t2,kmax
      common /variables/ t,V,Theta,slip,dt
      common /consts/ year
      common /saves/ tsav
      common /mvcount/nmvc

      real(DP) :: mask(0:NBMAX-1)
      real(DP) :: Vp  (0:NBMAX-1)
      real(DP) :: Vp0,eta
      common /for_mask/ mask
      common /Vpara/ Vp0,Vp,eta
!C===
      nmvc=0
!C
!C +------------------------+
!C | initialize MPI setting |
!C +------------------------+
!C===
!C--  set MPI_COMM_WORLD
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, np, ierr )
      if ( np /= NPROC ) then
        print *, "wrong np", np, NPROC
        stop
      endif
      call MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr )

!C--  set communicators "NthG"
      call comm__init(comm_all, group_all, np, my_rank, NPEperGROUP)
!C===

!C--  set timer routine
      call timer__init ( comms_new(NthG)%comm )

!K    if (DEBUG) call mpi_timer_start( 1 )
      call init_set(hdone, h)
!K    if (DEBUG) call mpi_timer_start( 1 )

      nfmax=11

!c--- for save data
!c    hdone=0.0_dp
      dt   =0.0_dp
      nfile=10
      nsav = 0
      tsav (nsav)= t
      dtsav(nsav)=dt
      do nb=0,NBMAX-1
        vsav (nb,nsav)=V    (nb)
        thsav(nb,nsav)=Theta(nb)
        ssav (nb,nsav)=slip (nb)
        if ( irest /= 0 ) then !for RESTART CALC
          sliponly(nb,nsav)= cslip (nb)
          deficit (nb,nsav)=                                            &
     &            (Vp0*Vp(nb)*tsav(nsav)-sliponly(nb,nsav))*mask(nb)
        endif
      end do

!C--Initialize Crustal Deformation Routine
!K    call convolution__init (NBMAX, NSAVMAX)
!K    if ( comms_new(NthG)%myrank == 0 ) then
!K      call gmt__init (ierr)
!K    endif

      nstep=1
10001 continue
      call deriv_etc(h, hdone, dslpmax, VovrVp)
!K    if ( nsav == nsav_gmt_write ) then
!C----- << disp. increment version on 20120824
!K      if (nsav==0) then
!K        del_deficit(:)=deficit(:,nsav)-deficit(:,nsav)
!K        t_inc=1.0d0
!K        call convolution__conv                                        &
!K   &                 (NBMAX, del_deficit, t_inc, tsav(nsav)/year)
!K        nsav_gmt_write_pre=nsav
!K      else
!K        del_deficit(:)=-deficit(:,nsav)+deficit(:,nsav_gmt_write_pre)
!K        t_inc=(tsav(nsav)-tsav(nsav_gmt_write_pre))/year
!K        call convolution__conv                                        &
!K   &                 (NBMAX, del_deficit, t_inc, tsav(nsav)/year)
!K        nsav_gmt_write_pre=nsav
!K      endif
!C----- >> disp. increment version on 20120824
!K    endif

!K    if ( nsav==nsav_gmt_write ) then
!K      call gmt__main(ninc_gmt_write, t, ttVp, ttmask,                 &
!K   &                          Vp0*year,comms_new(NthG)%myrank )
!K      nsav_gmt_write=nsav+ninc_gmt_write
!K      call MPI_Bcast( nsav_gmt_write, 1, MPI_INTEGER, 0,              &
!K   &                               comms_new(NthG)%comm, ierr )
!K    endif

      if ( DEBUG .and. comms_new(NthG)%myrank == 0 ) then
        print *, nstep, nsav, t/year, hdone, dslpmax, h
      end if

!c--  save & write data
      if (dslpmax >= dslipsav) then
        nsav=nsav+1
        if (nsav < NSAVMAX) then
          call isave_vals(nsav)
          if ( comms_new(NthG)%myrank ==0 .and. mod(nsav,100) == 0 ) then
            print *, nsav, tsav(nsav)/year, dtsav(nsav)
          end if
        else
          call iwrite_monitor(NSAVMAX, 10)
          nfile=nfile+1
          if (OUT_LARGE_FILE) then
            call iwrite_binary_outfile(        1,NSAVMAX, 12  )
          endif
!K        if (DEBUG) call mpi_timer_end( 1 )
          if (nfile == nfmax) goto 10000
          call iprep_next
        end if
      end if
!C
!C +-------------+
!C | timer check |
!C +-------------+
!C===
      call timer__mpi_check ( comms_new(NthG)%comm )
      if ( TIMER_ON .and. TIMER_LIMIT ) then
        if (comms_new(NthG)%myrank == 0) then
        print *,nsav,stopwatch,'timer limit! stopping your calculation!'
        endif
        call iwrite_monitor(nsav, 10)
        nfile=nfile+1
        if (OUT_LARGE_FILE) then
          call iwrite_binary_outfile(     1,nsav, 12 )
        endif
!K      if (DEBUG) call mpi_timer_end( 1 )
        call timer__mpi_check ( comms_new(NthG)%comm )
        if (comms_new(NthG)%myrank == 0) then
          print *, nsav, 'end write files', stopwatch
        endif
        goto 10000
      end if
!C===

!c--   time integration
      if ((t+h-t2) > 0) h=t2-t
      call time_integ(h,hdone,eps)
      if (t > t2) then
        call iwrite_monitor(NSAVMAX, 10)
        if (OUT_LARGE_FILE) then
          call iwrite_binary_outfile(        1,NSAVMAX, 12   )
        endif
        if (comms_new(NthG)%myrank == 0) then
          print *, 'reach tend'
        endif
        goto 10000
      end if

      nstep=nstep+1
      goto 10001

 1001 format(a2,i5,2(1x,e13.6))
 1011 format(a2,i5,3(1x,e13.6))
 1002 format(3(e13.6,1x))

10000 continue

!K    call convolution__writeTS ! write Time Sequence for OBS
!K    call eqk_slip_calc(nsav,tsav(0:nsav-1),dtsav(0:nsav-1),           &
!K   &     vsav(0:NBMAX-1,0:nsav-1),ssav(0:NBMAX-1,0:nsav-1),ierr)
!K    if (ierr/=0) then
!K      print *, 'eqk_slip ierr=', ierr, comms_new(NthG)%myrank
!K    endif
!K    call MPI_BARRIER(comms_new(NthG)%comm,ierr)
!K
!K    if ( comms_new(NthG)%myrank == 0 ) then
!K      call gmt__finalize (ierr)
!K    endif
!K
!K    if (DEBUG) call mpi_timer_end  ( 1 )
      call MPI_BARRIER(comms_new(NthG)%comm,ierr)
!C
!C +---------------------------------------+
!C | output cell slip-history if necessary |
!C +---------------------------------------+
!C===
!K      if (DEBUG) call mpi_timer_start( 2 )
!K      if (OUT_TIMESEQ) then
!K        call merge_init(comms_new(NthG)%myrank)
!K        call merge_data(comms_new(NthG)%myrank, comms_new(NthG)%comm ,  &
!K     &                                       NPEperGROUP, NSAVMAX, 0)
!K      endif
!K      if (DEBUG) call mpi_timer_end  ( 2 )
!K      if (DEBUG) call mpi_timer_start( 3 )
!K      if (OUT_TIMESEQ) then
!K      if (comms_new(NthG)%myrank==0) then
!K        call merge_out(NthG, NSAV, tsav(:)/year)
!K      endif
!K      endif
!K      if (DEBUG) call mpi_timer_end  ( 3 )
!C===

!K    if (DEBUG) call mpi_timer_result                                  &
!K   &          (NthG,comms_new(NthG)%myrank,comms_new(NthG)%comm )
      if (comms_new(NthG)%myrank==0) write(6,*) '***** nmvc=', nmvc
      call comm__finalize(comm_all, group_all)
      call MPI_FINALIZE(ierr)

      STOP

      CONTAINS
!C
!C +---------------------------+
!C | save variables for output |
!C +---------------------------+
!C===
      SUBROUTINE isave_vals(nsav)
      integer (kp), intent(IN) :: nsav
      integer (kp)             :: nb
      tsav  (nsav) =  t
      hprint(nsav) =  h
       dtsav(nsav) = dt
       dt          = 0.0_dp
      do nb=0,NBMAX-1
        vsav (nb,nsav)=V(nb)
        thsav(nb,nsav)=Theta(nb)
        ssav (nb,nsav)=slip(nb)
!       stsav(nb,nsav)=stress(nb)  ! for STRESS
        slip (nb)     =0.0_dp
      end do
      sliponly(0:,nsav)=ssav(0:,nsav)+sliponly(0:,nsav-1)

      deficit (0:,nsav)=(Vp0*Vp(0:)*tsav(nsav)-sliponly(0:,nsav))*mask(0:)
      RETURN
      END SUBROUTINE isave_vals
!C===

!C
!C +---------------+
!C | write monitor |
!C +---------------+
!C===
      SUBROUTINE iwrite_monitor (NSAVS, fcode)
      integer (kp), intent(IN) :: NSAVS, fcode
      integer (kp)             :: isav
      if (comms_new(NthG)%myrank==0) then
        write(fname,'(A,I3.3)') 'monitor_', NthG-1
        open (fcode,FILE=fname, status='unknown')
        do isav=0,NSAVS-1
          write(fcode,1001)                                             &
     &     'k=',isav, tsav(isav)/year, dtsav(isav), hprint(isav)
        end do
        close(fcode)
      end if
 1001 format(a2,i5,4(1x,e13.6))
      RETURN
      END SUBROUTINE iwrite_monitor
!C===

!C
!C +-------------------------+
!C | write binary outputfile |
!C +-------------------------+
!C===
      SUBROUTINE iwrite_binary_outfile(NSAVS,NSAVE, fcode)
      integer(kp), intent(IN) :: NSAVS, NSAVE, fcode
      integer(kp)             :: isav, myid
      myid=comms_new(NthG)%myrank
      write(fname,'(I3.3,A,A,I4.4)') NthG-1,trim(fn_out),'.p_',myid
      open (fcode,FILE=trim(fname),form='unformatted',status='unknown')
      write(fcode) NSAVE-NSAVS+1
      write(fcode)                                                      &
     &   tsav(NSAVS-1:NSAVE-1)/year, dtsav(NSAVS-1:NSAVE-1),            &
     & hprint(NSAVS-1:NSAVE-1)
      do isav=NSAVS-1, NSAVE-1
        write(fcode) vsav(0:NBMAX-1,isav),thsav(0:NBMAX-1,isav),        &
     &               ssav(0:NBMAX-1,isav)!,stsav(0:NBMAX-1,isav)
      end do
      close(fcode)
      END SUBROUTINE iwrite_binary_outfile
!C===

!C
!C +---------------+
!C | dummy routine |
!C +---------------+
!C===
      SUBROUTINE iprep_next
      integer (kp) :: nb
      nsav=0
       tsav(nsav)=  t
      dtsav(nsav)= dt
      dt         = 0.0_dp
      do nb=0,NBMAX-1
        vsav (nb,nsav)=    V (nb)
        thsav(nb,nsav)=Theta (nb)
        ssav (nb,nsav)=slip  (nb)
!       stsav(nb,nsav)=stress(nb)
        slip (nb)     =0.0_dp
      end do
      END SUBROUTINE iprep_next
!C===
      END PROGRAM RSGDX

!C**
!C**   time integration : 5th-order Runge-Kutta
!C**
      subroutine time_integ(h,hdone,eps)
      Use precisions
      Use Constant
      USE m_comm
      IMPLICIT NONE
      INCLUDE 'mpif.h'

      real(dp)::A2,A3,A4,A5,A6
      real(dp)::B21,B31,B32,B41,B42,B43,B51,B52,B53,B54,B61,B62,B63,B64,B65
      real(dp)::C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,             &
     &B21=.2d0,B31=3.d0/40.d0,B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,        &
     &B43=1.2d0,B51=-11.d0/54.d0,B52=2.5d0,B53=-70.d0/27.d0,            &
     &B54=35.d0/27.d0,B61=1631.d0/55296.d0,B62=175.d0/512.d0,           &
     &B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,    &
     &C1=37.d0/378.d0,C3=250.d0/621.d0,C4=125.d0/594.d0,                &
     &C6=512.d0/1771.d0,DC1=C1-2825.d0/27648.d0,                        &
     &DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,                &
     &DC5=-277.d0/14336.d0,DC6=C6-.25d0)
      real(dp)::SAFETY,PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)
!c    constant
      real(dp)::mask(0:NBMAX-1)
!c    variable
      real(dp)::t,dt
      real(dp)::V    (0:NBMAX-1)
      real(dp)::Theta(0:NBMAX-1)
      real(dp)::slip (0:NBMAX-1)
!c
      real(dp)::Vtmp  (0:NBMAX-1)
      real(dp)::Thtmp (0:NBMAX-1)
      real(dp)::dVdt  (0:NBMAX-1)
      real(dp)::dThdt (0:NBMAX-1)
      real(dp)::dVdt2 (0:NBMAX-1)
      real(dp)::dThdt2(0:NBMAX-1)
      real(dp)::dVdt3 (0:NBMAX-1)
      real(dp)::dThdt3(0:NBMAX-1)
      real(dp)::dVdt4 (0:NBMAX-1)
      real(dp)::dThdt4(0:NBMAX-1)
      real(dp)::dVdt5 (0:NBMAX-1)
      real(dp)::dThdt5(0:NBMAX-1)
      real(dp)::dVdt6 (0:NBMAX-1)
      real(dp)::dThdt6(0:NBMAX-1)
      real(dp)::Verr  (0:NBMAX-1)
      real(dp)::Therr (0:NBMAX-1)
!c    for MPI
      integer(kp)::ierr
      integer(kp)::my_rank
!c    define new h
      real(dp)::h,hdone
      real(dp)::eps
      real(dp)::err
      real(dp)::errmax
      real(dp)::errmax_gb
      real(dp)::Vscal (0:NBMAX-1)
      real(dp)::Thscal(0:NBMAX-1)
      real(dp)::htemp
!c    else
      integer(kp)::nbobs
      real(dp)   :: year

!c--- common
      common /forMPI/ my_rank
      common /variables/ t,V,Theta,slip,dt
      common /for_derivs/ dVdt,dThdt
      common /scale/ Vscal,Thscal
      common /for_mask/ mask
      common /consts/ year

 1    do nbobs=0,NBMAX-1
          Vtmp(nbobs)=    V(nbobs)+B21*h* dVdt(nbobs)
         Thtmp(nbobs)=Theta(nbobs)+B21*h*dThdt(nbobs)
      end do

      call derivs(Vtmp,Thtmp,dVdt2,dThdt2)

      do nbobs=0,NBMAX-1
         Vtmp (nbobs)=V    (nbobs)                                      &
     &        +h*(B31*dVdt (nbobs)+B32*dVdt2 (nbobs))
         Thtmp(nbobs)=Theta(nbobs)                                      &
     &        +h*(B31*dThdt(nbobs)+B32*dThdt2(nbobs))
      end do

      call derivs(Vtmp,Thtmp,dVdt3,dThdt3)

      do nbobs=0,NBMAX-1
         Vtmp(nbobs) =V    (nbobs)                                      &
     &        +h*(B41*dVdt (nbobs)+B42*dVdt2 (nbobs)+B43*dVdt3 (nbobs))
         Thtmp(nbobs)=Theta(nbobs)                                      &
     &        +h*(B41*dThdt(nbobs)+B42*dThdt2(nbobs)+B43*dThdt3(nbobs))
      end do

      call derivs(Vtmp,Thtmp,dVdt4,dThdt4)

      do nbobs=0,NBMAX-1
         Vtmp(nbobs) =V    (nbobs)                                      &
     &        +h*(B51*dVdt (nbobs)+B52*dVdt2 (nbobs)+B53*dVdt3 (nbobs)  &
     &        +B54*dVdt4 (nbobs))
         Thtmp(nbobs)=Theta(nbobs)                                      &
     &        +h*(B51*dThdt(nbobs)+B52*dThdt2(nbobs)+B53*dThdt3(nbobs)  &
     &        +B54*dThdt4(nbobs))
      end do

      call derivs(Vtmp,Thtmp,dVdt5,dThdt5)

      do nbobs=0,NBMAX-1
         Vtmp(nbobs) =V    (nbobs)                                      &
     &        +h*(B61*dVdt (nbobs)+B62*dVdt2 (nbobs)+B63*dVdt3 (nbobs)  &
     &        +B64*dVdt4 (nbobs)+B65*dVdt5(nbobs))
         Thtmp(nbobs)=Theta(nbobs)                                      &
     &        +h*(B61*dThdt(nbobs)+B62*dThdt2(nbobs)+B63*dThdt3(nbobs)  &
     &        +B64*dThdt4(nbobs)+B65*dThdt5(nbobs))
      end do

      call derivs(Vtmp,Thtmp,dVdt6,dThdt6)

      do nbobs=0,NBMAX-1
         Vtmp(nbobs) =V    (nbobs)                                      &
     &        +h*(C1*dVdt (nbobs)+C3*dVdt3 (nbobs)                      &
     &        +C4*dVdt4 (nbobs)+C6*dVdt6 (nbobs))
         Thtmp(nbobs)=Theta(nbobs)                                      &
     &        +h*(C1*dThdt(nbobs)+C3*dThdt3(nbobs)                      &
     &        +C4*dThdt4(nbobs)+C6*dThdt6(nbobs))
         Verr(nbobs)=                                                   &
     &        h*(DC1*dVdt(nbobs)+DC3*dVdt3(nbobs)                       &
     &        +DC4*dVdt4(nbobs)+DC5*dVdt5(nbobs)+DC6*dVdt6(nbobs))
         Therr(nbobs)=                                                  &
     &        h*(DC1*dThdt(nbobs)+DC3*dThdt3(nbobs)                     &
     &        +DC4*dThdt4(nbobs)+DC5*dThdt5(nbobs)+DC6*dThdt6(nbobs))
      end do

!c---  define new h
      errmax=0.0_dp
      do  nbobs=0,NBMAX-1
         err=max(                                                       &
     &        abs(Verr(nbobs)/Vscal(nbobs)*mask(nbobs)),                &
     &        abs(Therr(nbobs)/Thscal(nbobs)*mask(nbobs))               &
     &        )
         errmax=max(err,errmax)
      end do

      call MPI_BARRIER  (comms_new(NthG)%comm,ierr)
      call MPI_ALLREDUCE(errmax , errmax_gb, 1, MPI_REAL8,              &
     &                   MPI_MAX, comms_new(NthG)%comm, ierr)

      errmax=errmax_gb/eps

      if (errmax > 1.) then
!c--    retry
        htemp=SAFETY*h*(errmax**PSHRNK)
        h = sign(max(abs(htemp), 0.1*abs(h)), h)

        if ( t+h == t .or. h < 1d-10 ) then
          write(6,*) 'stepsize underflow'
          write(6,*) 'tnew,t,h',t+h,t,h
        end if
        goto 1
      else
!c--    done
        hdone=h
        t=t+hdone
        do nbobs=0,NBMAX-1
          V    (nbobs)= Vtmp(nbobs)
          Theta(nbobs)=Thtmp(nbobs)
        end do
        if (errmax > ERRCON)then
          h=SAFETY*h*(errmax**PGROW)
        else
          h=5.*h
        endif
      end if
      RETURN
      end SUBROUTINE time_integ

!C***
!C*** subroutine : derivs
!C***
      subroutine derivs(Vtmp,Thtmp,dVdt,dThdt)
      USE precisions
      Use Constant

      USE m_comm

      implicit none
      INCLUDE 'mpif.h'

!c---  constants
      real(dp) :: Vp0
      real(dp) :: eta
      real(dp) :: mask (0:NBMAX-1)
      real(dp) :: Vp   (0:NBMAX-1)
!c--   for fric
      real(dp) :: mu0  (0:NBMAX-1)
      real(dp) :: a    (0:NBMAX-1)
      real(dp) :: b    (0:NBMAX-1)
      real(dp) :: L    (0:NBMAX-1)
      real(dp) :: sigma(0:NBMAX-1)
!c--- variable
      real(dp) :: Vtmp (0:NBMAX-1)
      real(dp) :: Thtmp(0:NBMAX-1)
      real(dp) :: dVdt (0:NBMAX-1)
      real(dp) :: dThdt(0:NBMAX-1)
!c
      real(dp) :: sum  (0:NBMAX-1)
      real(dp) :: cfac
      integer(kp) :: nftime, ierr
      real(dp) :: dVall(0:NMAX-1)
      real(dp) :: dV   (0:NBMAX-1)
!c--- for MPI
      integer(kp) :: my_rank

      integer(kp) :: nmvc
!c    else
      integer(kp) ::  nbobs,nbsrc,istart,iend

!c--- common
      common /forMPI/ my_rank
      common /for_mask/ mask
      common /Vpara/ Vp0,Vp,eta
      common /fricpara/ mu0,a,b,L,sigma

      common /mvcount/nmvc

      istart=NBMAX* comms_new(NthG)%myrank
      iend  =NBMAX*(comms_new(NthG)%myrank+1)-1

      nftime=0
      do nbsrc=0,NBMAX-1
        dV(nbsrc)=(Vtmp(nbsrc)-Vp(nbsrc))*mask(nbsrc)
      end do
      !CDIR NODEP
      call MPI_BARRIER  ( comms_new(NthG)%comm, ierr )
      call MPI_ALLGATHER( dV, NBMAX, MPI_REAL8, dVall, NBMAX,           &
     &                    MPI_REAL8, comms_new(NthG)%comm, ierr )
      call MPI_BARRIER  ( comms_new(NthG)%comm, ierr )

      nbsrc=NMAX

      if (comms_new(NthG)%myrank==0) nmvc=nmvc+1

      call multipl ( istart, iend, dVall, sum )

!CDIR NODEP
      do nbobs=0,NBMAX-1

        cfac=1.0_dp
!c      cfac=1.d0+Vtmp(nbobs)*Vp0/(0.1d-3)

!c--   dTheta/dt for composite law (Vc=1.0e-8_dp)
        if (b(nbobs) == 0 .or. mask(nbobs) == 0) then
          dThdt(nbobs)=0.0_dp
        else
          dThdt(nbobs)=Vp0*(                                            &
     &        b(nbobs)/L(nbobs)                                         &
     &         *exp(-Thtmp(nbobs)/b(nbobs)-Vtmp(nbobs)/(1.0e-8_dp/Vp0)) &
     &        -Vtmp(nbobs)/L(nbobs)                                     &
     &         *(Thtmp(nbobs)-b(nbobs)*log(cfac*Vp(nbobs)/Vtmp(nbobs))) &
     &           )

!c slip law
!c           dThdt(nbobs)=Vp0*(
!c    +           -Vtmp(nbobs)/L(nbobs)*(
!c    +           Thtmp(nbobs)-b(nbobs)*
!c    +            log(cfac*Vp(nbobs)/Vtmp(nbobs))
!c    +           )
!c    +           )
        end if

!c--   dV/dt
        dVdt(nbobs)=-1.0_dp/(                                           &
     &        a(nbobs)*sigma(nbobs)/Vtmp(nbobs)/cfac+eta*Vp0            &
     &        )*(                                                       &
     &        sigma(nbobs)*dThdt(nbobs)+Vp0*sum(nbobs)                  &
     &        )*mask(nbobs)
      end do

      RETURN
      END SUBROUTINE derivs

!C***
!C*** subroutine : init_set
!C***
      subroutine init_set(hdone,h)
      USE precisions
      Use Constant
      Use Hmatrices

      Use Hmatrices_matrixxvector_index
      USE m_comm
      IMPLICIT NONE
      INCLUDE 'mpif.h'

      integer(kp) :: iroot, ierr
!c--- for MPI
      integer(kp) :: my_rank
      real   (dp) :: lmd,mu,dip,LL,WW,dL,depth0,rake
      integer(kp) :: icalc

!c--  for fric
!c    input
      real(dp) :: fmu,aa,bb,bb2,L0,rho,rhow,g
      real(dp) :: d10,d20,d30,b_e1,b_e2,d_mask1,d_mask2
!c    calc
      real(dp) :: mu0  (0:NBMAX-1)
      real(dp) :: a    (0:NBMAX-1)
      real(dp) :: b    (0:NBMAX-1)
      real(dp) :: L    (0:NBMAX-1)
      real(dp) :: sigma(0:NBMAX-1)
      real(dp) :: mask (0:NBMAX-1)
!c
!c--  for plate velocity
!c    input
      real(dp) :: Vp0
      real(dp) :: Vp_c1,Vp_c2,Vp_x
!c    calc
      real(dp) :: Vp   (0:NBMAX-1)

!c--  for control
      real   (dp) :: eps,dslipsav,t2,etafac
      real   (dp) :: h0,V0
      integer(kp) :: kmax!,irest
!     character(len=9) :: fn
!     character(len=9) :: fn_rst
!c    calc
      real(dp) :: Vs
      real(dp) :: eta
!c--- variables
      real(dp) :: t,dt
      real(dp) :: h,hdone
      real(dp) :: V    (0:NBMAX-1)
      real(dp) :: Vpast(0:NBMAX-1)
      real(dp) :: Theta(0:NBMAX-1)
      real(dp) :: slip (0:NBMAX-1)
!c--- else
!     integer nb,nall
      integer(kp) :: nb
      character(len=30) :: fname
      real(dp) :: year
!c
!      real(dp) :: tta    (0:NMAX-1)
!      real(dp) :: ttb    (0:NMAX-1)
!      real(dp) :: ttL    (0:NMAX-1)
!      real(dp) :: ttsigma(0:NMAX-1)
!      real(dp) :: ttmask (0:NMAX-1)
!      real(dp) :: ttVp   (0:NMAX-1)
!c
      integer(kp) :: im

      integer(kp) :: k, kk, kkk, i, num_b
      integer(kp) :: icols, icol_offset, ncs, nce
      character(len=20) :: inp_namelist

!c--- common
      common /forMPI/ my_rank
      common /for_mask/ mask
      common /Vpara/ Vp0,Vp,eta
      common /fricpara/ mu0,a,b,L,sigma
      common /variables/ t,V,Theta,slip,dt
      common /control/ eps,dslipsav,t2,kmax
      common /consts/ year
      common /past/ Vpast
      integer(KP), pointer :: rows      (:)=>null()
      integer(KP), pointer :: cols      (:)=>null()
      integer(KP), pointer :: row_offset(:)=>null()
      integer(KP), pointer :: col_offset(:)=>null()
      integer(KP), pointer :: k_r       (:)=>null()

      namelist /input1/ lmd,mu,dip,LL,WW,dL,depth0,rake
      namelist /input2/ icalc
      namelist /input3/ fmu,aa,bb,bb2,L0,rho,rhow,g
      namelist /input4/ d10,d20,d30,b_e1,b_e2,d_mask1,d_mask2
      namelist /input5/ Vp0,Vp_c1,Vp_c2,Vp_x
      namelist /input6/ eps,dslipsav,t2,etafac
      namelist /input7/ h0,V0
      namelist /input8/ kmax,irest
      namelist /input9/ fn_rst, fn_out

      write(fname,'(i3.3,a,i4.4)')                                      &
     &             NthG-1, '/dfullmat.', comms_new(NthG)%myrank
      Open (12,file=trim(fname),status='unknown',form='unformatted')
      read (12) kfmax
      allocate(rows      (kfmax))
      allocate(cols      (kfmax))
      allocate(row_offset(kfmax))
      allocate(col_offset(kfmax))
      read (12) rows
      read (12) cols
      read (12) row_offset
      read (12) col_offset

      allocate(fulm(kfmax))
      do im=1, kfmax
        fulm(im)%rows_f=rows(im)
        fulm(im)%cols_f=cols(im)
        fulm(im)%row_offset_f=row_offset(im)
        fulm(im)%col_offset_f=col_offset(im)
        allocate(fulm(im)%el_f(0:rows(im)*cols(im)-1))
        read(12) fulm(im)%el_f
      enddo
      deallocate(rows      )
      deallocate(cols      )
      deallocate(row_offset)
      deallocate(col_offset)
      close(12)
!CO   Print *,"my_rank=",my_rank,"end read d_fullmatrix !"
!CO   Print *,"my_rank=",my_rank,"full-matrices =",kfmax

      write(fname,'(i3.3,a,i4.4)')                                      &
     &             NthG-1,'/drkmat.',comms_new(NthG)%myrank
      Open (12,file=trim(fname),status='unknown',form='unformatted')
      read (12) krmax
      allocate(rows      (krmax))
      allocate(cols      (krmax))
      allocate(row_offset(krmax))
      allocate(col_offset(krmax))
      allocate(k_r       (krmax))
      read (12) rows
      read (12) cols
      read (12) row_offset
      read (12) col_offset
      read (12) k_r

      allocate(rkm(krmax))
      do im=1, krmax
        rkm(im)%rows_r=rows(im)
        rkm(im)%cols_r=cols(im)
        rkm(im)%row_offset_r=row_offset(im)
        rkm(im)%col_offset_r=col_offset(im)
        rkm(im)%k_r=k_r(im)
        allocate(rkm(im)%el_a(0:rows(im)*k_r(im)-1))
        allocate(rkm(im)%el_b(0:cols(im)*k_r(im)-1))
        read(12) rkm(im)%el_a
        read(12) rkm(im)%el_b
      enddo
      close(12)
!CO   Print *,"my_rank=",comms_new(NthG)%myrank,"end read d_rkmatrix !"
!CO   Print *,"my_rank=",comms_new(NthG)%myrank,"Rk-matrices =",krmax

!fj>MH
      !C-- create index for rkmatrixxvector()
      kmax_tot_rkm = 0
      do k=1, krmax
        kmax_tot_rkm = kmax_tot_rkm + rkm(k)%k_r
      enddo

      allocate(  k_rkm     (kmax_tot_rkm) )
      allocate( kk_rkm     (kmax_tot_rkm) )
      allocate( num_bst_rkm(kmax_tot_rkm) )

!C-- k_rkm:通しをつけたrankが何番の圧縮行列に属するか
!C--kk_rkm:通しをつけたrankがその圧縮行列の何番のrankに対応するか
!C--kkkに対応したnum_bst_rkm el_bのスタートアドレス
      kkk=0
      do k=1, krmax
        icols=rkm(k)%cols_r
        icol_offset=rkm(k)%col_offset_r
        ncs=icol_offset
        nce=icol_offset+icols-1
        kmax=rkm(k)%k_r
        num_b=0
        do kk=1, kmax
          kkk=kkk+1
           k_rkm     (kkk)=k
          kk_rkm     (kkk)=kk
          num_bst_rkm(kkk)=num_b
          do i=ncs, nce
            num_b=num_b+1
          enddo
        enddo
      enddo
!fj<MH

!c---  open files
      write(inp_namelist,'(A,i3.3,A)') 'inputpara.g_',NthG-1,'.nml'
      open(10,FILE=trim(inp_namelist),status='old')
!c---  input parameters
      read (10,nml=input1)
      read (10,nml=input2)
      read (10,nml=input3)
      read (10,nml=input4)
      read (10,nml=input5)
      read (10,nml=input6)
      read (10,nml=input7)
      read (10,nml=input8)
      read (10,nml=input9)

      close(10)

!c---  calc parameters

      mu0( 1)=fmu
      Vs=sqrt(mu/rho)
      eta=mu*.5/Vs*etafac
!*[cm/yr]->[m/sec]
      year=60.d0*60.d0*24.d0*365.2425d0
      Vp0=Vp0*1.d-2/year
      if ( comms_new(NthG)%myrank == 0 ) then
        write(fname,'(i3.3,a)') NthG-1,'/hlib_a_b.p_all'
        open(3,FILE=fname, status='old')

!C     own block only
        do nb=0,NMAX-1
           read(3,*) tta(nb),ttb(nb),ttL(nb),ttsigma(nb),ttVp(nb),      &
     &               ttmask(nb)
        end do
        close(3)
      end if

      call MPI_BARRIER(comms_new(NthG)%comm,ierr)
      iroot=0
      call MPI_SCATTER(tta,NBMAX,MPI_REAL8,a,NBMAX,MPI_REAL8,           &
     &                 iroot,comms_new(NthG)%comm,ierr)
      call MPI_SCATTER(ttb,NBMAX,MPI_REAL8,b,NBMAX,MPI_REAL8,           &
     &                 iroot,comms_new(NthG)%comm,ierr)
      call MPI_SCATTER(ttL,NBMAX,MPI_REAL8,L,NBMAX,MPI_REAL8,           &
     &                 iroot,comms_new(NthG)%comm,ierr)
      call MPI_SCATTER(ttsigma,NBMAX,MPI_REAL8,sigma,NBMAX,MPI_REAL8,   &
     &                 iroot,comms_new(NthG)%comm,ierr)
      call MPI_SCATTER(ttVp,NBMAX,MPI_REAL8,Vp,NBMAX,MPI_REAL8,         &
     &                 iroot,comms_new(NthG)%comm,ierr)
      call MPI_SCATTER(ttmask,NBMAX,MPI_REAL8,mask,NBMAX,MPI_REAL8,     &
     &                 iroot,comms_new(NthG)%comm,ierr)

!c--  irest=1:restart, =0:normal start
      if (irest == 0) then
        t=0.
        h=h0
        hdone=0.
        do nb=0,NBMAX-1
          V(nb)=V0
          Vpast(nb)=V0
          Theta(nb)=-b(nb)*log(V(nb)/Vp(nb))
          slip(nb)=0.d0
        end do
      else
        write(fname,'(I3.3,A,A,I4.4)')                                  &
     &           NthG-1,trim(fn_rst),'.p_',comms_new(NthG)%myrank
        call iread_binary_output(12, fname)
      endif
 1001 format(a2,i5,2(1x,e13.6))
 1002 format(3(e13.6,1x))
      RETURN
      CONTAINS
!C
!C  +----------------------------+
!C  | read binary 11_output file |
!C  +----------------------------+
!C====
      SUBROUTINE iread_binary_output (fcode, fname)
      integer(kp )    , intent(IN) :: fcode
      character(LEN=*), intent(IN) :: fname
      integer(kp )                 :: nsav, nsav_rst, nb
      integer(kp )                 :: openstatus
      real(DP),  allocatable       ::   t_tmp(:)
      real(DP),  allocatable       ::  dt_tmp(:)
      real(DP),  allocatable       ::   h_tmp(:)

      real(DP),  allocatable       ::  vsav_tmp(:)
      real(DP),  allocatable       :: thsav_tmp(:)
      real(DP),  allocatable       ::  ssav_tmp(:)
!CO   real(DP),  allocatable       :: stsav_tmp(:)

!CO   print *, trim(fname)
      open(fcode, FILE=trim(fname), form='unformatted', status='old',   &
     &                     iostat=openstatus )
      if (openstatus/=0) then
        print *,"restart file open error   :", trim(fname), NthG
      else
        print *,"restart file open success :", trim(fname), NthG
      endif

      read(fcode) nsav_rst
      if (nsav_rst-1<irest) STOP'error : nsav_rst is smaller than irest'
      allocate( t_tmp(0:nsav_rst-1));  t_tmp = 0.0_dp
      allocate(dt_tmp(0:nsav_rst-1)); dt_tmp = 0.0_dp
      allocate( h_tmp(0:nsav_rst-1));  h_tmp = 0.0_dp
      read(fcode) t_tmp, dt_tmp, h_tmp
      t=t_tmp(irest)*year
      h=h_tmp(irest)
      deallocate(t_tmp,dt_tmp,h_tmp)
      allocate( vsav_tmp(0:NBMAX-1));  vsav_tmp = 0.0_dp
      allocate(thsav_tmp(0:NBMAX-1)); thsav_tmp = 0.0_dp
      allocate( ssav_tmp(0:NBMAX-1));  ssav_tmp = 0.0_dp
!     allocate(stsav_tmp(0:NBMAX-1)); stsav_tmp = 0.0_dp
      allocate( cslip   (0:NBMAX-1));     cslip = 0.0_dp
      do nsav=0, irest
        read(fcode) vsav_tmp,thsav_tmp,ssav_tmp!,stsav_tmp
        cslip=cslip+ssav_tmp
      enddo
      close(fcode)

      do nb=0, NBMAX-1
        Theta(nb)= thsav_tmp(nb)
        V    (nb)=  vsav_tmp(nb)
        Vpast(nb)=         V(nb)
        slip(nb) =  ssav_tmp(nb)
      enddo
      deallocate(vsav_tmp ,thsav_tmp , ssav_tmp)
!     deallocate(vsav_tmp ,thsav_tmp , ssav_tmp ,stsav_tmp)
      RETURN
      END SUBROUTINE iread_binary_output
!C====
      END SUBROUTINE init_set

!C***
!C*** SUBROUTINE : deriv_etc
!C*** calc derivative & scale for V & Theta also calc slip, stress & dsmax
      subroutine deriv_etc(h, hdone, dslpmax, VoV)
      USE precisions
      USE m_comm
      Use Constant

      IMPLICIT NONE
      INCLUDE 'mpif.h'
      real(dp) :: tiny
      parameter(TINY=1.0d-30)

      real(dp) :: VoV   (0:NMAX-1)

!c--- for MPI
      integer(kp) :: my_rank,ierr
      integer(kp) :: istart, iend
!c--- constants
      real(dp) :: Vp0
      real(dp) :: eta
      real(dp) :: mask  (0:NBMAX-1)
      real(dp) :: Vp    (0:NBMAX-1)
!c--  for fric
      real(dp) :: mu0   (0:NBMAX-1)
      real(dp) :: a     (0:NBMAX-1)
      real(dp) :: b     (0:NBMAX-1)
      real(dp) :: L     (0:NBMAX-1)
      real(dp) :: sigma (0:NBMAX-1)
!c--- variables
      real(dp) :: t,dt
      real(dp) :: h,hdone
      real(dp) :: V     (0:NBMAX-1)
      real(dp) :: Vpast (0:NBMAX-1)
      real(dp) :: Theta (0:NBMAX-1)
      real(dp) :: slip  (0:NBMAX-1)

      real(dp) :: sum_V (0:NBMAX-1)
      real(dp) :: dVdt  (0:NBMAX-1)
      real(dp) :: dThdt (0:NBMAX-1)
      real(dp) :: Vscal (0:NBMAX-1)
      real(dp) :: Thscal(0:NBMAX-1)
      real(dp) :: cfac

      real(dp) :: dslpmax_gb,dslpmax
      real(dp) :: dVall (0:NMAX -1)
      real(dp) :: dV    (0:NBMAX-1)
      real(dp) :: dslip
!c--- else
      integer(kp) :: nbobs,nbsrc

      integer(kp) :: nmvc
!c--- common
      common /forMPI/ my_rank
      common /for_mask/ mask
      common /Vpara/ Vp0,Vp,eta
      common /fricpara/ mu0,a,b,L,sigma
      common /variables/ t,V,Theta,slip,dt
      common /for_derivs/ dVdt,dThdt
      common /scale/ Vscal,Thscal
      common /past/ Vpast

!c--- add common
      common /mvcount/nmvc
!c
!c--- gather V & calc slip
      dt=dt+hdone
      do nbsrc=0,NBMAX-1
        dslip=(Vpast(nbsrc)+V(nbsrc))*Vp0*0.5*hdone*mask(nbsrc)
        Vpast(nbsrc)=V(nbsrc)
        slip(nbsrc)=slip(nbsrc)+dslip
      end do

      do nbsrc=0,NBMAX-1
        dV(nbsrc)=(V(nbsrc)-Vp(nbsrc))*mask(nbsrc)
      end do

      call MPI_ALLGATHER(dV,NBMAX,MPI_REAL8,dVall,NBMAX,                &
     &     MPI_REAL8,comms_new(NthG)%comm,ierr)

      VoV=dVall

      nbsrc=NMAX
      if (my_rank == 0) nmvc=nmvc+1

      istart= comms_new(NthG)%myrank   *NBMAX
      iend  =(comms_new(NthG)%myrank+1)*NBMAX-1
      Call multipl(istart,iend,dVall,sum_V)

!c--- calc derivatives etc.
      do nbobs=0,NBMAX-1
        cfac=1.0_dp
!c      cfac=1.+V(nbobs)*Vp0/(0.1d-3)

!c--  dTheta/dt for composite law (Vc=1.0e-8)
        if (b(nbobs) == 0 .or. mask(nbobs) == 0) then
          dThdt(nbobs)=0.0_dp
        else
          dThdt(nbobs)=Vp0*(                                            &
     &        b(nbobs)/L(nbobs)                                         &
     &         *exp(-Theta(nbobs)/b(nbobs)-V(nbobs)/(1.0e-8_dp/Vp0))    &
     &        -V(nbobs)/L(nbobs)                                        &
     &         *(Theta(nbobs)-b(nbobs)*log(cfac*Vp(nbobs)/V(nbobs)))    &
     &           )
!c slip law
!c           dThdt(nbobs)=Vp0*(
!c    +           -V(nbobs)/L(nbobs)*(
!c    +           Theta(nbobs)-b(nbobs)*
!c    +           log(cfac*Vp(nbobs)/V(nbobs))
!c    +           )
!c    +           )
        end if
!c--    dV/dt
        dVdt(nbobs)=-1.e0_dp/(                                          &
     &        a(nbobs)*sigma(nbobs)/V(nbobs)/cfac+eta*Vp0               &
     &        )*(                                                       &
     &        sigma(nbobs)*dThdt(nbobs)+Vp0*sum_V(nbobs)                &
     &        )*mask(nbobs)

!c--   scale for error estimation
        Vscal(nbobs)=( abs(V(nbobs))+abs(h*dVdt(nbobs)) )               &
     &        *mask(nbobs)+TINY

        Thscal(nbobs)=( abs(Theta(nbobs))+abs(h*dThdt(nbobs)) )         &
     &        *mask(nbobs)+TINY
      end do
!c--- mean dslip

      dslpmax=0._dp
      do nbobs=0,NBMAX-1
        dslpmax=max(slip(nbobs),dslpmax)
      end do

      call MPI_ALLREDUCE(dslpmax,dslpmax_gb,1,MPI_REAL8,                &
     &     MPI_MAX,comms_new(NthG)%comm,ierr)
      dslpmax=dslpmax_gb

      RETURN
      END SUBROUTINE deriv_etc

!C***
!C*** subroutine : multipl
!C***
      Subroutine multipl(istart,iend,vec,res_vec)
      USE precisions
      Use Constant
      Use Hmatrices
!fj>MH
      use Hmatrices_matrixxvector_index
!fj<MH
!     use mpi_timer
      IMPLICIT NONE

      Integer(kp)                   :: k,istart,iend,i
      Real(dp),dimension(0:NMAX -1) :: vec
      Real(dp),Dimension(0:NBMAX-1) :: res_vec
!fj>MH
      Real(dp), allocatable, save :: res_vec_wk(:,:)
      Integer(kp)          , save :: iflg_allocated = 0
      Integer(kp)          , save :: num_threads
      Integer(kp)                 :: ithread_no
      Integer(kp)                 :: kk, kkk, num_b_st
      Integer                     :: omp_get_max_threads
      integer                     :: omp_get_thread_num
!fj<MH

      Do i=0,NBMAX-1
        res_vec(i)=0._dp
      EndDo

!fj>MH
      if ( iflg_allocated == 0 ) then
        iflg_allocated = 1
        num_threads = omp_get_max_threads()
!MISS   allocate ( res_vec_wk(0:NMAX -1,0:num_threads-1) ) !by KURATA
        allocate ( res_vec_wk(0:NBMAX-1,0:num_threads-1) )
      endif
!fj<MH

!fj>MH
!$OMP parallel default(none)   &
!$OMP shared(istart, iend,  krmax, kfmax) &
!$OMP shared(vec, res_vec, res_vec_wk) &
!$OMP shared(num_threads) &
!$OMP shared(kmax_tot_rkm, k_rkm,kk_rkm, num_bst_rkm) &
!$OMP private(i,k,kk,kkk,num_b_st,ithread_no)

      ithread_no = omp_get_thread_num()

!MISS do k=0, NMAX-1 !by KURATA
      do k=0, NBMAX-1
        res_vec_wk(k,ithread_no) = 0.0
      enddo

!$OMP do SCHEDULE(static,1)
      Do kkk=1,kmax_tot_rkm
        k = k_rkm(kkk)
        kk=kk_rkm(kkk)
        num_b_st=num_bst_rkm(kkk)
        Call rkmatrixxvector(istart,iend,k,kk,num_b_st,vec,             &
     &                       res_vec_wk(0,ithread_no))
      EndDo
!$OMP end do

!C    type2:fullmatrix x vector
!$OMP do SCHEDULE(static,1)
      Do k=1,kfmax
        Call fullmatrixxvector(istart,iend,k,vec,res_vec_wk(0,ithread_no))
      EndDo
!$OMP end do

!$OMP do
      do k=0, NBMAX-1
        do i=0, num_threads-1
          res_vec(k)=res_vec(k)+res_vec_wk(k,i)
        enddo
      enddo
!$OMP end do

!$OMP end parallel
      RETURN
      End Subroutine multipl

!C***
!C*** subroutine : rkmatrixxvector
!C***
      Subroutine rkmatrixxvector(istart,iend,k,kk,num_b_st,vec,res_vec)
      Use precisions
      Use Constant
      Use Hmatrices
      IMPLICIT NONE    !MH
!fj>MH
      Integer(kp), intent(IN) :: istart, iend, k, kk, num_b_st
      Real (dp), dimension(0:NMAX -1), intent(IN)    :: vec
      Real (dp), dimension(0:NBMAX-1), intent(INOUT) :: res_vec

      Integer(kp) :: row_offset,col_offset,rows,cols,kmax
      Integer(kp) :: nrs,nre,ncs,nce,num_a1,num_b,num_a,i
      Real(dp)    :: coef
!fj<MH
      real(dp)    :: ddot

      rows=rkm(k)%rows_r
      cols=rkm(k)%cols_r
      row_offset=rkm(k)%row_offset_r
      col_offset=rkm(k)%col_offset_r
      kmax=rkm(k)%k_r

      nrs=max(istart,row_offset)
      nre=min(iend,row_offset+rows-1)
      ncs=col_offset
      nce=col_offset+cols-1
      num_a1=nrs-row_offset

!fj>MH

!C
!C +---------------+
!C | inner-product |
!C +---------------+
!C====
      num_b=num_b_st
      coef=0.0_dp
      if ( BLAS ) then
!C--     BLAS
        coef=ddot(cols, rkm(k)%el_b(num_b:),1,vec(ncs:),1)
      else
!C--  no BLAS
        Do i=ncs,nce
          coef=coef+rkm(k)%el_b(num_b)*vec(i)
          num_b=num_b+1
        EndDo
      endif
!C====

!C
!C +-------+
!C | DAXPY |
!C +-------+
!C====
      num_a=num_a1+(kk-1)*rows
      if ( BLAS ) then
!C--     BLAS
        call daxpy ( nre-nrs+1, coef, rkm(k)%el_a(num_a:),              &
     &                       1,  res_vec(nrs-istart), 1   )
      else
!C--  no BLAS
        Do i=nrs,nre
          res_vec(i-istart)=res_vec(i-istart)+coef*rkm(k)%el_a(num_a)
          num_a=num_a+1
        EndDo
!C====
      endif

!fj<MH
      RETURN
      End Subroutine rkmatrixxvector

!C***
!C*** subroutine : fullmatrixxvector
!C***
      Subroutine fullmatrixxvector(istart,iend,k,vec,res_vec)
      USE precisions
      Use Constant
      Use Hmatrices
      IMPLICIT NONE

      Integer(kp) :: row_offset,col_offset,rows,cols
      Integer(kp) :: nrs,nre,ncs,nce,num1,num,i,j
      Integer(kp) ::istart,iend,k
      integer(kp) :: ncol,nrow
      Real(dp),dimension(0:NMAX -1) ::vec
      Real(dp),dimension(0:NBMAX-1) ::res_vec

      rows=fulm(k)%rows_f
      cols=fulm(k)%cols_f
      row_offset=fulm(k)%row_offset_f
      col_offset=fulm(k)%col_offset_f

      nrs=max(istart,row_offset)
      nre=min(iend,row_offset+rows-1)
      ncs=col_offset
      nce=col_offset+cols-1
      num1=(nrs-row_offset)*cols
!C
!C +-------+
!C | DGEMV |
!C +-------+
!C====
      if (BLAS) then
!C--     BLAS
        ncol=cols
        nrow=nre-nrs+1
        call dgemv("T",ncol,nrow,1.0d0,fulm(k)%el_f(num1),cols,         &
     &                       vec(ncs),1,1.0d0,res_vec(nrs-istart),1)
      else
!C--  no BLAS
        Do i=nrs,nre
          num=num1+cols*(i-nrs)
          Do j=ncs,nce
            res_vec(i-istart)=res_vec(i-istart)+fulm(k)%el_f(num)*vec(j)
            num=num+1
          EndDo
        EndDo
!C====
      endif
      RETURN
      End Subroutine fullmatrixxvector
