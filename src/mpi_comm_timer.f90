!C +------------------------------------------------+
!C |  add character "comment"         on 2012/04/13 |
!C |                                                |
!C | * for multi calc                               |
!C |                                                |
!C |            LAST CHANGE: 2012/05/23 16:56:56.   |
!C +------------------------------------------------+
module mpi_comm_timer
use mpi
  implicit none
! include 'mpif.h'

  integer, private :: i, myrank, ierr
  integer, private, parameter :: size = 100
  real(8), private, save :: start(1:size) = 0.0_8, end(1:size) = 0.0_8, &
 &     elapse(1:size) = 0.0_8, max_elapse(1:size)
  character(len=100), private, save :: comment(1:size)

       private
       public :: mpi_timer_start, mpi_timer_end, mpi_timer_result
contains

! subroutine mpi_timer_start(num,parag)
! integer, intent(in) :: num
! character(len=*) :: parag
! start(num) = mpi_wtime()
! comment(num)= parag
! end subroutine mpi_timer_start

  subroutine mpi_timer_start(num)
  integer, intent(in) :: num
  start(num) = mpi_wtime()
  end subroutine mpi_timer_start


  subroutine mpi_timer_end(num)
  integer, intent(in) :: num
  end(num) = mpi_wtime()
  elapse(num) = elapse(num) + (end(num) - start(num))
! write(*, *) 'mpi_timer number = ', num, 'elapse = ', elapse(num)
  end subroutine mpi_timer_end

  subroutine mpi_timer_result(ith,myrank, comm)
  integer(4), intent(IN) :: ith, myrank, comm

!  comment( 1)='derivs      : main      '
!  comment( 2)='derivs      : visual merge'
!  comment( 3)='derivs      : visual write'
  write(comment(1),'(A,i3.3)') 'derivs      :main except IO group=',ith
  write(comment(2),'(A,i3.3)') 'derivs      : visual merge  group=',ith
  write(comment(3),'(A,i3.3)') 'derivs      : visual write  group=',ith

!CO  comment( 1)='derivs      : total     '
!CO  comment( 2)='derivs      : set dV    '
!CO  comment( 3)='derivs      : MPI_BARRIER before MPI_ALLGATHER '
!CO  comment( 4)='derivs      : MPI_ALLGATHER, set dVall         '
!CO  comment( 5)='derivs      : MPI_BARRIER after  MPI_ALLGATHER '
!CO  comment( 6)='derivs      : count nmvc                       '
!CO  comment( 7)='derivs      : subroutine multipl               '
!CO  comment( 8)='derivs      : set dThdt, dVdt                  '
!CO
!CO  comment( 9)='multipl     : except MPI_REDUCE_SCATTER        '
!CO  comment(10)='multipl     : only   MPI_REDUCE_SCATTER        '
!CO  comment(11)='multipl     : total                            '
!CO
!CO  comment(12)='derivs_etc  : set slip   '
!CO  comment(13)='derivs_etc  : MPI_BARRIER before MPI_ALLGATHER '
!CO  comment(14)='derivs_etc  : MPI_ALLGATHER, set dVall         '
!CO  comment(15)='derivs_etc  : count nmvc                       '
!CO  comment(16)='derivs_etc  : subroutine multipl               '
!CO  comment(17)='derivs_etc  : set dThdt, dVdt, Vscal, Thscal   '
!CO  comment(18)='derivs_etc  : calc dslpmax                     '
!CO  comment(19)='derivs_etc  : MPI_BARRIER before MPI_ALLREDUCE '
!CO  comment(20)='derivs_etc  : MPI_ALLREDUCE for max slip       '
!CO  comment(21)='derivs_etc  : total                            '
!CO
!CO  comment(31)='RSGDX:write binary output for ALL data'
!CO  comment(32)='RSGDX:write ascii  output for ALL data '
!CO  comment(33)='RSGDX:write binary output for timer stoped case'
!CO  comment(34)='RSGDX:write ascii  output for timer stoped case'
! call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)

  do i = 1, size
!    call mpi_reduce(elapse(i), max_elapse(i), 1, MPI_DOUBLE_PRECISION, &
!         MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call mpi_reduce(elapse(i), max_elapse(i), 1, MPI_DOUBLE_PRECISION, &
          MPI_MAX, 0,     comm      , ierr)

     if(myrank == 0) then
        if(max_elapse(i) /= 0.0_8) then
           write(*, *) '##>', trim(comment(i))
           write(*, *) 'mpi_timer number = ', i, 'max elapse = ', max_elapse(i)
        end if
     end if
!    call mpi_reduce(elapse(i), max_elapse(i), 1, MPI_DOUBLE_PRECISION, &
!         MPI_MIN, 0, MPI_COMM_WORLD, ierr)
     call mpi_reduce(elapse(i), max_elapse(i), 1, MPI_DOUBLE_PRECISION, &
          MPI_MIN, 0,     comm      , ierr)

     if(myrank == 0) then
        if(max_elapse(i) /= 0.0_8) then
           write(*, *) 'mpi_timer number = ', i, 'min elapse = ', max_elapse(i)
           write(*, *) '<##'
        end if
     end if
  end do
  end subroutine mpi_timer_result

end module mpi_comm_timer
