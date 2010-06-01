program main
  ! This file is the object-oriented fortran equivalent of Epetra_power_method.cpp. In Trilinos 10.0,
  ! this is a snapshot of an unstable (evolving) file expected to become stable in a
  ! subsequent release.  This file exercises the derived types defined in 
  ! ForTrilinos/src/epetra/FEpetra*.F90, which wrap the interface bodies in 
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
  use mpi
  use FEpetra_MpiComm      ,only : Epetra_MpiComm
#else
  use FEpetra_SerialComm   ,only : Epetra_SerialComm
#endif
  use FEpetra_Map          ,only : Epetra_Map
  use FEpetra_Vector       ,only : Epetra_Vector
  use FEpetra_CrsMatrix    ,only : Epetra_CrsMatrix
  use ForTrilinos_utils    ,only : valid_kind_parameters
  use ForTrilinos_enum_wrappers
  use ForTrilinos_error
  use iso_c_binding        ,only : c_int,c_double
  implicit none
! Data declarations 
#ifdef HAVE_MPI  
  type(Epetra_MpiComm) :: communicator
#else
  type(Epetra_SerialComm) :: communicator
#endif
  type(Epetra_Map)    :: map
  type(Epetra_CrsMatrix) :: A
  type(error)         :: err
  integer(c_int)      :: NumGlobalElements, NumGlobalElements_return
  integer(c_int),dimension(:),allocatable :: MyGlobalElements
  integer(c_int),dimension(:),allocatable :: NumNz
  integer(c_int)      :: NumMyElements,i
  integer(c_int)      :: Index_Base=1
  integer(c_int)      :: MyPID, NumProc
  logical             :: verbose
  integer(c_int)      :: indices(2), NumEntries
  real(c_double)      :: two = 2.0,values(2)
  integer             :: rc, ierr,ierr_pm=0 
  real(c_double)      :: lambda=0.0,tolerance=1.0E-2
  integer(c_int)      :: niters, numvals
  integer(c_int),dimension(:),allocatable::Rowinds
  real(c_double),dimension(:),allocatable::Rowvals

  if (.not. valid_kind_parameters()) stop 'C interoperability not supported on this platform.'
  
  ! Executable code
  
! Create a comm
#ifdef HAVE_MPI
  call MPI_INIT(ierr)
  communicator= Epetra_MpiComm(MPI_COMM_WORLD) 
# else
  communicator= Epetra_SerialComm() 
#endif

  MyPID   = communicator%MyPID()
  NumProc = communicator%NumProc()
  verbose = MyPID==0
  print *,verbose
  print *,MyPID,NumProc

! Create a map 
  NumGlobalElements = 100 
  if (NumGlobalElements < NumProc) stop 'Number of global elements cannot be less that number of processors'
  map = Epetra_Map(NumGlobalElements,Index_Base,communicator)
  NumGlobalElements_return = map%NumGlobalElements()   ! test line

! Get update list and number of local equations from newly created Map
  NumMyElements = map%NumMyElements()
  print *,'NumGlobalElements = ', numGlobalElements_return  ! test line
  print *,'NumMyElements=', map%NumMyElements()             ! test line
  if ( NumGlobalElements /= NumGlobalElements_return ) &
    stop 'In ForTrilinos (verySimpleObjectOriented.F90: return mismatch'
  allocate(MyGlobalElements(NumMyElements))
  MyGlobalElements = map%MyGlobalElements()

! Create an integer vector NumNz tat is used to build the Epetra Matrix
! NumNz(i) is the number of OFF-DIAGONAL term for the ith global equation
! on this processor
  allocate(NumNz(NumMyElements))

! We are building a tridiagonal matrix where each row has (-1 2 -1)
! So we need 2 off-diagonal terms (except for the first and last equation)
  do i=1,NumMyElements
   if(MyGlobalElements(i)==1.or.MyGlobalElements(i)==NumGlobalElements) then
     NumNz(i) = 2_c_int
   else
     NumNz(i) = 3_c_int
   end if
  end do
  
! Create a Epetra_Matrix
  A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

! Add rows one at a time
! Need some vectors to help
! off diagonal values will always be -1
  values(1) = -1.0
  values(2) = -1.0
  do i=1,NumMyElements
    if (MyGlobalElements(i)==1) then
      indices(1) = 1
      NumEntries = 1
    else if(MyGlobalElements(i)==NumGlobalElements) then
      indices(1) = NumGlobalElements-1
      NumEntries = 1
    else
      indices(1) = MyGlobalElements(i)-1
      indices(2) = MyGlobalElements(i)+1
      NumEntries = 2
    end if
     call A%InsertGlobalValues(MyGlobalElements(i),NumEntries,values,indices,err)
     if (err%error_code()/=0) stop 'A%InsertGlobalValues: failed'
  !Put in the diaogonal entry
     call A%InsertGlobalValues(MyGlobalElements(i),1,[two],MyGlobalElements,err)
     if (err%error_code()/=0) stop 'A%InsertGlobalValues: failed'
  end do
 
  !Finish up
  call A%FillComplete()
 
  !Create vectors for power methods
  !variable needed for interation
  niters=NumGlobalElements*10 
  call power_method(A,lambda,niters,tolerance,verbose,ierr_pm)
  ierr_pm=ierr_pm+1

  !Iterate
  if (A%MyGlobalRow(0_c_int)) then
    numvals=A%NumGlobalEntries(1_c_int)
    allocate(Rowvals(numvals),Rowinds(numvals))
  end if
  call A%ExtractGlobalRowCopy(0_c_int,numvals,numvals,Rowvals,Rowinds) ! get A(0,0)
  do i=1,numvals
   call A%ReplaceGlobalValues(0_c_int,numvals,Rowvals,Rowinds)
  enddo  

  !Iterate again
  call power_method(A,lambda,niters,tolerance,verbose,ierr_pm)
  ierr_pm=ierr_pm+1

 
  ! Clean up memory (in reverse order).  This step is not required
  ! with compilers that support Fortran 2003 type finalization:
  call A%force_finalize()
  call map%force_finalize()
  call communicator%force_finalize()
 
#ifdef HAVE_MPI
  call MPI_FINALIZE(rc)
#endif

contains

subroutine power_method(A,lambda,niters,tolerance,verbose,ierr_pm)
 use iso_c_binding        ,only : c_int,c_double
 use FEpetra_Vector, only : Epetra_Vector
 use FEpetra_CrsMatrix, only : Epetra_CrsMatrix
 implicit none
 type(Epetra_CrsMatrix), intent(inout) :: A
 integer(c_int),    intent(inout):: ierr_pm
 real(c_double) :: lambda
 integer(c_int), intent(in) :: niters
 real(c_double), intent(in) :: tolerance 
 logical, intent(in)        :: verbose
 real(c_double), allocatable,dimension(:) :: normz,residual
 type(Epetra_Vector) :: q,z,resid
 integer(c_int) :: iter
 print *,'Inside power method'
 ierr_pm=1
 q = Epetra_Vector(A%RowMap()) 
 z = Epetra_Vector(A%RowMap()) 
 resid = Epetra_Vector(A%RowMap()) 

 !Fill z with random numbers
 call z%Random()

 do iter=1,niters
  normz=z%Norm2()              !Compute 2-norm of z
  call q%Scale(1.0/normz(1),z)  
  call A%Multiply(.false.,q,z) ! Compute z=A*q
  call q%Dot(z,[lambda])       ! Approximate maximum eignvalue
  if (mod(iter,100)==0.or.iter==niters) then
   call resid%Update(1.0_c_double,z,-lambda,q,0.0_c_double) ! Compute A*q-lambda*q
   residual=resid%Norm2()
   if (verbose) then
    print *,'Iter=',iter,'lambda=',lambda,'Resisual of A*q-lambda*q=',residual(1)
   endif
  endif
  if (residual(1)<tolerance) then
   ierr_pm=0
   exit
  endif
 enddo
end subroutine

end 
