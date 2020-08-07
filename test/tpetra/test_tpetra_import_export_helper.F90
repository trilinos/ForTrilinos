#include "ForTrilinosTpetra_config.hpp"
module test_Tpetra_import_export_helper
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  implicit none

contains
  ! ----------------------- Create An Import ----------------------------- !
  !   For uniform contiguous distribution
  subroutine Tpetra_IE_MakeImport(comm, srcDist, tgtDist, imp)
    type(TeuchosComm), intent(in) :: comm
    integer, intent(in) :: srcDist, tgtDist
    type(TpetraImport), intent(out) :: imp
    type(TpetraMap) :: src, tgt

    src = TpetraMap(TPETRA_GLOBAL_INVALID, srcDist, comm)
    tgt = TpetraMap(TPETRA_GLOBAL_INVALID, tgtDist, comm)

    imp = TpetraImport(src, tgt)

    call src%release()
    call tgt%release()
    return
  end subroutine Tpetra_IE_MakeImport

  !  For uniform distribution to nonuniform
  subroutine Tpetra_IE_MakeImport_NU(comm, srcDist, tgtDist, imp)
    type(TeuchosComm), intent(in) :: comm
    integer(size_type), intent(in) :: tgtDist(:)
    integer :: srcDist
    type(TpetraImport), intent(out) :: imp
    type(TpetraMap) :: src, tgt

    src = TpetraMap(TPETRA_GLOBAL_INVALID, srcDist, comm)
    tgt = TpetraMap(TPETRA_GLOBAL_INVALID, tgtDist, comm)

    imp = TpetraImport(src, tgt)

    call src%release()
    call tgt%release()
    return
  end subroutine Tpetra_IE_MakeImport_NU

  ! ----------------------- Create An Export ----------------------------- !
  !   For uniform contiguous distribution
  subroutine Tpetra_IE_MakeExport(comm, srcDist, tgtDist, exp)
    type(TeuchosComm), intent(in) :: comm
    integer, intent(in) :: srcDist, tgtDist
    type(TpetraExport), intent(out) :: exp
    type(TpetraMap) :: src, tgt

    src = TpetraMap(TPETRA_GLOBAL_INVALID, srcDist, comm)
    tgt = TpetraMap(TPETRA_GLOBAL_INVALID, tgtDist, comm)

    exp = TpetraExport(src, tgt)

    call src%release()
    call tgt%release()
    return
  end subroutine Tpetra_IE_MakeExport

  !  For nonuniform distribution to uniform
  subroutine Tpetra_IE_MakeExport_NU(comm, srcDist, tgtDist, exp)
    type(TeuchosComm), intent(in) :: comm
    integer(size_type), intent(in) :: srcDist(:)
    integer :: tgtDist
    type(TpetraExport), intent(out) :: exp
    type(TpetraMap) :: src, tgt

    src = TpetraMap(TPETRA_GLOBAL_INVALID, srcDist, comm)
    tgt = TpetraMap(TPETRA_GLOBAL_INVALID, tgtDist, comm)

    exp = TpetraExport(src, tgt)

    call src%release()
    call tgt%release()
    return
  end subroutine Tpetra_IE_MakeExport_NU
end module test_Tpetra_import_export_helper
