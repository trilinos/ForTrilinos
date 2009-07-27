module ForTrilinos_enums
  use iso_c_binding ,only : c_int ! Kind parameter (precision specifier)
  implicit none                   ! Prevent implicit typing

  ! This file is the Fortran companion to CTrilinos/src/CTrilinos_enums.h.
  ! The definitions here are the nearest Fortran equivalents to the contents of that file.
  ! The Fortran 2003 standard guarantees that the types (e.g., 'integer(c_int)') and 
  ! enumeration values are interoperable with those employed by a companion C compiler.

  enum ,bind(C)
    enumerator ::                     &
      FT_Invalid_ID,                  &
      FT_Epetra_Distributor_ID,       &
      FT_Epetra_SerialComm_ID,        &
      FT_Epetra_BLAS_ID,              &
      FT_Epetra_Comm_ID,              &
      FT_Epetra_Operator_ID,          &
      FT_Epetra_MultiVector_ID,       &
      FT_Epetra_OffsetIndex_ID,       &
      FT_Epetra_Object_ID,            &
      FT_Epetra_RowMatrix_ID,         &
      FT_Epetra_CompObject_ID,        &
      FT_Epetra_Directory_ID,         &
      FT_Epetra_Flops_ID,             &
      FT_Epetra_SrcDistObject_ID,     &
      FT_Epetra_MpiComm_ID,           &
      FT_Epetra_CrsMatrix_ID,         &
      FT_Epetra_CrsGraph_ID,          &
      FT_Epetra_DistObject_ID,        &
      FT_Epetra_Vector_ID,            &
      FT_Epetra_Export_ID,            &
      FT_Epetra_Map_ID,               &
      FT_Epetra_BlockMap_ID,          &
      FT_Epetra_Import_ID
  end enum

  ! Since the Fortran 2003 standard guarantees that enum values correspond to C int values, we can create
  ! the alias below for c_int with certainty that it can be used as the Fortran kind parameter that
  ! makes Fortran integer values interoperable with C enumeration values. This alias is the Fortran
  ! counterpart to CTrilinos_Type_ID_t in CTrilinos/src/CTrilinos_enums.h:

  integer(kind(c_int)) ,parameter :: ForTrilinos_Type_ID_t = c_int

  ! The type below is interoperable with CTrilinos_Object_ID_t in CTrilinos/src/CTrilinos_enums.h:

  type ,bind(C) :: ForTrilinos_Object_ID_t
    integer(ForTrilinos_Type_ID_t) :: type      ! Object data type (interoperable with CTrilinos_Type_ID_t)
    integer(c_int)                 :: index     ! Array index of the object
    logical(c_bool)                :: is_const  ! Whether or not object is declared const
  end type

  ! Each type definition below is identical in form to the ForTrilinos_Object_ID_t definition but with a name corresponding
  ! to the CT_Epetra_*_ID_t type (where *=Distributor,SerialComm,...) with which it is designed to be interoperable.

  type ,bind(C) :: FT_Epetra_Distributor_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_SerialComm_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_BLAS_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Comm_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Operator_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_MultiVector_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_OffsetIndex_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Object_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_RowMatrix_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_CompObject_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Directory_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Flops_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_SrcDistObject_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_MpiComm_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_CrsMatrix_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_CrsGraph_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_DistObject_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Vector_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Export_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Map_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_BlockMap_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Import_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; logical(c_bool) :: is_const
  end type

  ! Epetra_DataAcces

  integer(kind(c_int)) ,parameter :: FT_Epetra_DataAccess_E_t = c_int

  enum ,bind(C)
    enumerator ::                  &
      FT_Epetra_DataAccess_E_Copy, &
      FT_Epetra_DataAccess_E_View
  end enum

  ! Epetra_CombineMode

  integer(kind(c_int)) ,parameter :: FT_Epetra_CombineMode_E_t = c_int

  enum ,bind(C)
    enumerator ::                        &
      FT_Epetra_CombineMode_E_Add,       &
      FT_Epetra_CombineMode_E_Zero,      &
      FT_Epetra_CombineMode_E_Insert,    &
      FT_Epetra_CombineMode_E_InsertAdd, &
      FT_Epetra_CombineMode_E_Average,   &
      FT_Epetra_CombineMode_E_AbsMax
  end enum

end module ForTrilinos_enums
