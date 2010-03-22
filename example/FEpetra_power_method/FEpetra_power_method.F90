program main
  ! This file is the object-oriented equivalent of verySimple.F90.  In Trilinos 10.0,
  ! this is a snapshot of an unstable (evolving) file expected to become stable in a
  ! subsequent release.  This file exercises the derived types defined in 
  ! ForTrilinos/src/epetra/Epetra*.F90, which wrap the interface bodies in 
  ! ForTrilinos/src/epetra/forepetra.F90.   
    
  ! This file represents the preferred style for using ForTrilinos and is recommended for 
  ! Fortran users whose compilers support the object-oriented features of Fortran 2003.
  ! As of the Trilinos 10.0 release date, the latest versions of the IBM and Cray compilers 
  ! nominally support the required features.  The Numerical Algorithms Group (NAG) and Intel 
  ! compilers support all features but one: final subroutines.  (In each case, the support is
  ! somewhat immature and buggy.)  ForTrilinos/src/ForTrilinos_hermetic.F90 contains utilities 
  ! that help users work around the lack of final subroutines.

#include "ForTrilinos_config.h"
#ifdef HAVE_MPI
  use FEpetra_MpiComm      ,only : epetra_mpicomm
#else
  use FEpetra_SerialComm   ,only : epetra_serialcomm
#endif
  use FEpetra_Map          ,only : epetra_map
  use FEpetra_Vector       ,only : epetra_vector
  use ForTrilinos_utils    ,only : valid_kind_parameters
  use iso_c_binding        ,only : c_int,c_double
  implicit none
  interface 
   subroutine power_method(A,lambda,niters,tolerance,verbose)
    use iso_c_binding        ,only : c_int,c_double
    use FEpetra_Vector, only : epetra_Vector
    !use FEpetra_CrsMatrix, only : epetra_CrsMatrix
    !type(epetra_CrsMatrix), intent(inout) :: A
    type(epetra_Vector) :: q,z,resid,A
    real(c_double) :: lambda
    integer(c_int) :: niters
    real(c_double) :: tolerance
    logical        :: verbose
   end subroutine
  end interface
  ! Data declarations 
#ifdef HAVE_MPI  
  type(epetra_mpicomm) :: communicator
  integer(c_int)       :: ierror
#else
  type(epetra_serialcomm) :: communicator
#endif
  type(epetra_map)    :: map
  type(epetra_vector) :: x, b
  integer(c_int) :: NumGlobalElements, NumGlobalElements_return
  integer(c_int),dimension(:),allocatable :: MyGlobalElements
  integer(c_int) :: NumMyElements
  integer(c_int) :: Index_Base=1
  real(c_double) :: bnorm(1), xnorm(1)
  real(c_double) :: err_tol,expected_bnorm,expected_xnorm,bnorm_err,xnorm_err 
  real(c_double) :: two = 2.0, zero = 0.0
  logical        :: success = .true.,zero_initial=.true.
  integer,parameter:: MPI_COMM_WORLD=91
  integer        :: MyPID, NumProc
  logical        :: verbose

  if (.not. valid_kind_parameters()) stop 'C interoperability not supported on this platform.'
  
  ! Executable code
  
! Create a comm
#ifdef HAVE_MPI
  call MPI_INIT(ierror)
  communicator= epetra_mpicomm(MPI_COMM_WORLD) 
# else
  communicator= epetra_serialcomm() 
#endif

  MyPID   = communicator%MyPID()
  NumProc = communicator%NumProc()
  verbose = MyPID==0
  print *,MyPID,NumProc
! Create a map 
  NumGlobalElements = 100 
  if (NumGlobalElements < NumProc) stop 'Number of global elements cannot be less that number of processors'
  map = epetra_map(NumGlobalElements,Index_Base,communicator)
  NumGlobalElements_return = map%NumGlobalElements()
  NumMyElements = map%NumMyElements()
  print *,'NumGlobalElements = ', numGlobalElements_return
  print *,'NumMyElements=', map%NumMyElements()
  if ( NumGlobalElements /= NumGlobalElements_return ) &
    stop 'In ForTrilinos (verySimpleObjectOriented.F90: return mismatch'
  allocate(MyGlobalElements(NumMyElements))
  MyGlobalElements = map%MyGlobalElements()
   
  ! Create vectors
  x = epetra_vector(map,zero_initial)
  b = epetra_vector(map,zero_initial)
 
  ! Do some vector operations
  call b%PutScalar(two)
  call x%Update(two, b, zero) ! /* x = 2*b */
 
  bnorm = b%Norm2()
  xnorm = x%Norm2()
 
  print *, "2 norm of x = ", xnorm(1) 
  print *, "2 norm of b = ", bnorm(1) 

! Test the expected value 
  err_tol = 1.0e-14
  expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements_return )
  expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements_return )
  bnorm_err = abs( expected_bnorm - bnorm(1) ) / expected_bnorm
  xnorm_err = abs( expected_xnorm - xnorm(1) ) / expected_xnorm
  print *, "error in 2 norm of x = ",bnorm_err
  print *, "error in 2 norm of b = ",xnorm_err
  if (bnorm_err > err_tol) success = .false.
  if (xnorm_err > err_tol) success = .false.
 
  ! Clean up memory (in reverse order).  This step is not required
  ! with compilers that support Fortran 2003 type finalization:
  call b%force_finalization()
  call x%force_finalization()
  call map%force_finalization()
  call communicator%force_finalization()
 
  if (success) then
    print *  
    print *, "End Result: TEST PASSED" 
  else
    print *  
    print *, "End Result: TEST FAILED"
  end if
#ifdef HAVE_MPI
  call MPI_Finalize(ierror)
#endif
end program main

subroutine power_method(A,lambda,niters,tolerance,verbose)
 use iso_c_binding        ,only : c_int,c_double
 use FEpetra_Vector, only : epetra_Vector
 !use FEpetra_CrsMatrix, only : epetra_CrsMatrix
 !type(epetra_CrsMatrix), intent(inout) :: A
 type(epetra_Vector) :: q,z,resid,A
 real(c_double) :: lambda
 integer(c_int) :: niters
 real(c_double) :: tolerance 
 logical        :: verbose
 

 !q = epetra_vector(A%RowMap()) 
 !z = epetra_vector(A%RowMap()) 
 !resid = epetra_vector(A%RowMap()) 
end subroutine
