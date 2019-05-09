#include "ForTrilinosTpetra_config.hpp"
module test_Tpetra_multivector_helper
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  implicit none

contains
  ! ----------------------- Create A Vector ----------------------------- !
  !   For uniform contiguous distribution
  subroutine Tpetra_MV_Create(comm, num_local, num_vecs, vec, init_val)
    integer, parameter :: dp = kind(0.d0)
    type(TeuchosComm), intent(in) :: comm
    integer, intent(in) :: num_local
    integer(size_type) :: num_vecs
    type(TpetraMultiVector), intent(out) :: vec
    type(TpetraMap) :: src
    real(dp), intent(in), optional :: init_val

    src = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm)
    vec = TpetraMultiVector(src, num_vecs)
    if (present(init_val)) then
       call vec%putScalar(init_val)

    end if

    call src%release()
  end subroutine Tpetra_MV_Create

end module test_Tpetra_multivector_helper
