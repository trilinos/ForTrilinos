program main
  use ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING
  use simpleinterface
  use forteuchos
  use mpi
  implicit none

  integer(c_int) :: i
  integer(c_int) :: n, nnz;
  integer(c_int) :: my_rank, num_procs

  integer(c_int), dimension(:), allocatable :: row_inds, col_inds, row_ptrs
  real(c_double), dimension(:), allocatable :: lhs, rhs, values

  integer(c_int) :: cur_pos, offset

  type(ParameterList) :: plist

  integer :: ierr

  n = 50
  nnz = 3*n

  my_rank = 0
  num_procs = 1

  ! Initialize MPI subsystem
  call MPI_INIT(ierr)

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank,   ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

  ! Read in the parameterList
  call plist%create("Stratimikos")
  call load_from_xml(plist, "stratimikos.xml")

  ! Step 0: Construct tri-diagonal matrix, and rhs
  allocate(row_inds(n))
  allocate(row_ptrs(n+1))
  allocate(col_inds(nnz))
  allocate(values(nnz))
  row_ptrs(1) = 0
  cur_pos = 1
  offset  = n * my_rank
  do i = 1, n
    ! Is logic evaluated the same way in fortran?
    if (i .ne. 1 .or. my_rank > 0) then
      col_inds(cur_pos) = offset + i-1
      values  (cur_pos) = -1.0
      cur_pos = cur_pos + 1
    end if
    col_inds(cur_pos) = offset + i
    values  (cur_pos) = 2.0
    cur_pos = cur_pos + 1
    if (i .ne. n .or. my_rank .ne. num_procs-1) then
      col_inds(cur_pos) = offset + i+1
      values  (cur_pos) = -1.0
      cur_pos = cur_pos + 1
    end if
    row_ptrs(i+1) = cur_pos-1;

    row_inds(i) = offset + i
  end do
  nnz = cur_pos-1

  allocate(lhs(n))
  allocate(rhs(n))
  do i = 1, n
    rhs(i) = 1.0
  end do

  ! Step 1: initialize a handle
  call init()
  ! call init(mpi_comm)

  ! Step 2: setup the problem
  call setup_matrix(n, row_inds, row_ptrs, nnz, col_inds, values)

  ! // Step 3: setup the solver
  call setup_solver(plist)

  ! // Step 4: solve the system
  call solve(n, rhs, lhs)

  ! Step 5: clean up
  call finalize()

  call MPI_FINALIZE(ierr)

end program
