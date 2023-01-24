! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program main
! --------------------------------------------------------------------------- !
! Power method for estimating the eigenvalue of maximum magnitude of
! a matrix.
!
! We don't intend for you to write your own eigensolvers; the Anasazi
! package provides them.  You should instead see this class as a surrogate
! for a ForTrilinos interface to the Tpetra package.
! --------------------------------------------------------------------------- !
#include "ForTrilinos_config.h"

use iso_fortran_env
use, intrinsic :: iso_c_binding

#include "ForTrilinos.h"
use forteuchos
use fortpetra

#if FORTRILINOS_USE_MPI
use mpi
#endif

implicit none

! -- ForTrilinos objects
type(TeuchosComm) :: comm
type(TpetraMap) :: map
type(TpetraCrsMatrix) :: A

! -- Scalars
integer :: ierr
real(scalar_type) :: lambda
integer(global_size_type) :: num_gbl_indices
integer :: my_rank
integer(size_type) :: num_entries_in_row, max_entries_per_row, i
integer :: lcl_row, row_nnz, n
integer :: num_my_elements, iconv
integer(global_ordinal_type) gbl_row, id_of_first_row

! -- Arrays
integer(global_ordinal_type), allocatable :: cols(:)
real(scalar_type), allocatable :: vals(:)
! --------------------------------------------------------------------------- !


! Initialize MPI subsystem, if applicable
#if FORTRILINOS_USE_MPI
call MPI_INIT(ierr)
if (ierr /= 0) then
  stop "MPI failed to init"
endif
comm = TeuchosComm(MPI_COMM_WORLD)
#else
comm = TeuchosComm()
#endif

my_rank = comm%getRank()

! The number of rows and columns in the matrix.
num_gbl_indices = 50

map = TpetraMap(num_gbl_indices, comm)
FORTRILINOS_CHECK_IERR()

! Check that the map was created with the appropriate number of elements
if (map%getGlobalNumElements() /= num_gbl_indices) then
  write(error_unit, '(A,I3,A)') 'Expected ', num_gbl_indices, ' global indices'
  stop 1
end if


if (my_rank == 0) &
  write(*, *) "Creating the sparse matrix"

! Create a Tpetra sparse matrix whose rows have distribution given by the Map.

! TpetraOperator implements a function from one Tpetra(Multi)Vector to another
! Tpetra(Multi)Vector.  TpetraCrsMatrix implements TpetraOperator; its apply()
! method computes a sparse matrix-(multi)vector multiply.  It's typical for
! numerical algorithms that use Tpetra objects to be templated on the type of
! the TpetraOperator specialization.
max_entries_per_row = 3
A = TpetraCrsMatrix(map, max_entries_per_row)

! Fill the sparse matrix, one row at a time.
allocate(vals(3))
allocate(cols(3))
num_my_elements = int(map%getLocalNumElements(), kind=kind(num_my_elements))
fill: do lcl_row = 1, num_my_elements
  gbl_row = map%getGlobalElement(lcl_row)
  if (gbl_row == 1) then
    ! A(1, 1:2) = [2, -1]
    row_nnz = 2
    cols(1:2) = [gbl_row, gbl_row+1]
    vals(1:2) = [2d0, -1d0]
  else if (gbl_row == num_gbl_indices) then
    ! A(N, N-1:N) = [-1, 2]
    row_nnz = 2
    cols(1:2) = [gbl_row-1, gbl_row]
    vals(1:2) = [-1d0, 2d0]
  else
    ! A(i, i-1:i+1) = [-1, 2, -1]
    row_nnz = 3
    cols(1:3) = [gbl_row-1, gbl_row, gbl_row+1]
    vals(1:3) = [-1d0, 2d0, -1d0]
  end if
  call A%insertGlobalValues(gbl_row, cols(1:row_nnz), vals(1:row_nnz))
end do fill
deallocate(vals)
deallocate(cols)

! Tell the sparse matrix that we are done adding entries to it.
call A%fillComplete()

! Run the power method and report the result.
call PowerMethod(A, lambda, iconv, comm)
if (iconv > 0 .and. my_rank == 0) then
  write(*, *) "Estimated max eigenvalue: ", lambda
end if
if (abs(lambda-3.99)/3.99 >= .005) then
  write(error_unit, '(A)') 'Largest estimated eigenvalue is not correct!'
  stop 1
end if

! Now we're going to change values in the sparse matrix and run the power method
! again.  We'll increase the value of the (1,1) component of the matrix to make
! the matrix more diagonally dominant.  It should decrease the number of
! iterations required for the power method to converge. In Tpetra, if
! fillComplete() has been called, you have to call resumeFill() before you may
! change the matrix(either its values or its structure).

! Increase diagonal dominance
if (my_rank == 0) &
  write(*, *) "Increasing magnitude of A(1,1), solving again"

! Must call resumeFill() before changing the matrix, even its values.
call A%resumeFill()

id_of_first_row = 1
if (map%isNodeGlobalElement(id_of_first_row)) then
  ! Get a copy of the row with with global index 1.  Modify the diagonal entry
  ! of that row.  Submit the modified values to the matrix.
  num_entries_in_row = A%getNumEntriesInGlobalRow(id_of_first_row);
  n = int(num_entries_in_row, kind=kind(n))
  allocate(vals(n))
  allocate(cols(n))

  ! Fill vals and cols with the values resp.(global) column indices of the
  ! sparse matrix entries owned by the calling process.
  !
  ! Note that it's legal(though we don't exercise it in this example) for the
  ! row Map of the sparse matrix not to be one to one.  This means that more
  ! than one process might own entries in the first row.  In general, multiple
  ! processes might own the (1,1) entry, so that the global A(1,1) value is
  ! really the sum of all processes' values for that entry.  However, scaling
  ! the entry by a constant factor distributes across that sum, so it's OK to do
  ! so.
  call A%getGlobalRowCopy(id_of_first_row, cols, vals, num_entries_in_row)
  do i = 1, n
    if (cols(i) == id_of_first_row) then
      vals(i) = vals(i) * 10.
    end if
  end do

  ! "Replace global values" means modify the values, but not the structure of
  ! the sparse matrix.  If the specified columns aren't already populated in
  ! this row on this process, then this method throws an exception.  If you want
  ! to modify the structure(by adding new entries), you'll need to call
  ! insertGlobalValues().
  i = A%replaceGlobalValues(id_of_first_row, cols, vals)

  deallocate(vals)
  deallocate(cols)

end if

! Call fillComplete() again to signal that we are done changing the
! matrix.
call A%fillComplete()

! Run the power method again.
call PowerMethod(A, lambda, iconv, comm)
if (iconv > 0 .and. my_rank == 0) then
  write(*, *) "Estimated max eigenvalue: ", lambda
end if
if (abs(lambda-20.05555)/20.05555 >= .0001) then
  write(error_unit, '(A)') 'Largest estimated eigenvalue is not correct!'
  stop 1
end if

call A%release()
call map%release()
call comm%release()

#if FORTRILINOS_USE_MPI
! Finalize MPI must be called after releasing all handles
call MPI_FINALIZE(ierr)
#endif


contains

  subroutine PowerMethod(A, lambda, iconv, comm, verbose_in)
    ! ----------------------------------------------------------------------- !
    ! Power method for estimating the eigenvalue of maximum magnitude of
    ! a matrix.  This function returns the eigenvalue estimate.
    ! ----------------------------------------------------------------------- !
    use fortpetra
    implicit None
    integer(int_type), intent(out) :: iconv
    real(scalar_type), intent(out) :: lambda
    logical, intent(in), optional :: verbose_in
    type(TpetraCrsMatrix), intent(in) :: A
    type(TeuchosComm), intent(in) :: comm
    ! ----------------------------------------------------------------------- !
    ! Arguments
    ! ---------
    ! A          (input) The sparse matrix
    !
    ! lambda     (output) The estimated eigenvalue of maximum magnitude
    !
    ! iconv      (output) Flag indicated whether or not the procedure
    !    converged:
    !    iconv = 0  did not converge
    !          = 1  converged based on tol1 (more stringent)
    !          = 2  converged based on tol2 (less stringent)
    !
    ! verbose_in (input, optional) Print information diagnositics
    ! ----------------------------------------------------------------------- !
    ! -- Local Parameters
    real(scalar_type), parameter :: one=1., zero=0.
    integer(size_type), parameter :: num_vecs=1
    integer(int_type), parameter :: maxit1=500, maxit2=750
    real(scalar_type), parameter :: tol1=1.e-5, tol2=1.e-2

    ! -- Local Scalars
    real(scalar_type) :: normz, residual
    integer(int_type) :: report_frequency, it
    integer(size_type) :: my_rank
    logical :: verbose

    ! -- Local Arrays
    real(scalar_type) :: norms(num_vecs), dots(num_vecs)

    ! -- Local ForTrilinos Objects
    type(TpetraMultiVector) :: q, z, resid
    type(TpetraMap) :: domain_map, range_map
    ! ----------------------------------------------------------------------- !

    verbose = .false.
    if (present(verbose_in)) verbose = verbose_in
    my_rank = comm%getRank()

    ! Set output only arguments
    lambda = 0.
    iconv = 0

    ! Create three vectors for iterating the power method.  Since the power
    ! method computes z = A*q, q should be in the domain of A and z should be in
    ! the range. (Obviously the power method requires that the domain and the
    ! range are equal, but it's a good idea to get into the habit of thinking
    ! whether a particular vector "belongs" in the domain or range of the
    ! matrix.)  The residual vector "resid" is of course in the range of A.
    domain_map = A%getDomainMap()
    range_map = A%getRangeMap()
    q = TpetraMultiVector(domain_map, num_vecs)
    z = TpetraMultiVector(range_map, num_vecs)
    resid = TpetraMultiVector(range_map, num_vecs)

    ! Fill the iteration vector z with random numbers to start.  Don't have
    ! grand expectations about the quality of our pseudorandom number generator,
    ! but it is usually good enough for eigensolvers.
    call z%randomize()

    ! lambda: Current approximation of the eigenvalue of maximum magnitude.
    ! normz: 2-norm of the current iteration vector z.
    ! residual: 2-norm of the current residual vector 'resid'.
    lambda = zero
    normz = zero
    residual = zero

    ! How often to report progress in the power method.  Reporting progress
    ! requires computing a residual, which can be expensive.  However, if you
    ! don't compute the residual often enough, you might keep iterating even
    ! after you've converged.
    report_frequency = 10

    ! Do the power method, until the method has converged or the
    ! maximum iteration count has been reached.
    iconv = 0
    iters: do it = 1, maxit2
      call z%norm2(norms)
      normz = norms(1)
      call q%scale(one / normz, z)  !  q := z / normz
      call A%apply(q, z)  !  z := A * q
      call q%dot(z, dots)
      lambda = dots(1)  !  Approx. max eigenvalue

      ! Compute and report the residual norm every report_frequency
      ! iterations, or if we've reached the maximum iteration count.
      if (mod(it-1, report_frequency) == 0 .or. it == maxit2) then
        call resid%update(one, z, -lambda, q, zero)  ! z := A*q - lambda*q
        call resid%norm2(norms)
        residual = norms(1)  ! 2-norm of the residual vector
        if (verbose .and. my_rank == 0) then
          write(*, '(A,I3,A)') "Iteration ", it, ":"
          write(*, '(A,F8.4)') "- lambda = ", lambda
          write(*, '(A,E8.2)') "- ||A*q - lambda*q||_2 = ", residual
        end if
      end if

      if (it <= maxit1) then
        if(residual < tol1) then
          iconv = 1
          exit iters
        end if
      else
        if (residual < tol2) then
          iconv = 2
          exit iters
        end if
      end if

    end do iters

    if (iconv == 0 .and. my_rank == 0) then
      write(*, *) "PowerMethod failed to converge after ", maxit2, " iterations"
    end if

    call q%release()
    call z%release()
    call resid%release()

    return

  end subroutine PowerMethod

end program main
