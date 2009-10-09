module ForTrilinos_enums
  use iso_c_binding ,only : c_int        ! Kind parameter (precision specifier)
  implicit none                          ! Prevent implicit typing

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
      FT_Epetra_Import_ID,            &
      FT_Epetra_Time_ID,              &
      FT_Epetra_JadMatrix_ID,         &
      FT_Epetra_LinearProblem_ID,     &
      FT_Epetra_LAPACK_ID,            &
      FT_Teuchos_CommandLineProcessor_ID,& 
      FT_Teuchos_ParameterList_ID,    &
      FT_Teuchos_ParameterEntry_ID,   &
      FT_Teuchos_any_ID,              &
      FT_Amesos_BaseSolver_ID,        &
      FT_Amesos_ID,                   &
      FT_Epetra_FECrsMatrix_ID
  end enum

  ! Since the Fortran 2003 standard guarantees that enum values correspond to C int values, we can create
  ! the alias below for c_int with certainty that it can be used as the Fortran kind parameter that
  ! makes Fortran integer values interoperable with C enumeration values. This alias is the Fortran
  ! counterpart to CTrilinos_Type_ID_t in CTrilinos/src/CTrilinos_enums.h:

  integer(kind(c_int)) ,parameter :: ForTrilinos_Type_ID_t = c_int

  ! CTrilinos uses integers to pass C++ bools back and forth, so we should never use type c_bool when
  ! interfacing with CTrilinos. Use FT_boolean_t instead:

  integer(kind(c_int)) ,parameter :: FT_boolean_t = c_int

  integer(FT_boolean_t) ,parameter :: FT_FALSE = 0
  integer(FT_boolean_t) ,parameter :: FT_TRUE  = 1

  ! The type below is interoperable with CTrilinos_Object_ID_t in CTrilinos/src/CTrilinos_enums.h:

  type ,bind(C) :: ForTrilinos_Object_ID_t
    integer(ForTrilinos_Type_ID_t) :: type      ! Object data type (interoperable with CTrilinos_Type_ID_t)
    integer(c_int)                 :: index     ! Array index of the object
    integer(FT_boolean_t)          :: is_const  ! Whether or not object is declared const
  end type

  ! Each type definition below is identical in form to the ForTrilinos_Object_ID_t definition but with a name corresponding
  ! to the CT_Epetra_*_ID_t type (where *=Distributor,SerialComm,...) with which it is designed to be interoperable.

  type ,bind(C) :: FT_Epetra_Distributor_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_SerialComm_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_BLAS_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Comm_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Operator_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_MultiVector_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_OffsetIndex_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Object_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_RowMatrix_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_CompObject_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Directory_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Flops_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_SrcDistObject_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
#ifdef HAVE_MPI
  type ,bind(C) :: FT_Epetra_MpiComm_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
#endif /* HAVE_MPI */
  type ,bind(C) :: FT_Epetra_CrsMatrix_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_CrsGraph_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_DistObject_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Vector_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Export_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Map_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_BlockMap_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Import_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_Time_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_JadMatrix_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_LinearProblem_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Epetra_LAPACK_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Teuchos_CommandLineProcessor_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Teuchos_ParameterList_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Teuchos_ParameterEntry_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
  type ,bind(C) :: FT_Teuchos_any_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
#ifdef HAVE_FORTRILINOS_AMESOS
  type ,bind(C) :: FT_Amesos_BaseSolver_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
#endif /* HAVE_FORTRILINOS_AMESOS */
#ifdef HAVE_FORTRILINOS_AMESOS
  type ,bind(C) :: FT_Amesos_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
  end type
#endif /* HAVE_FORTRILINOS_AMESOS */
  type ,bind(C) :: FT_Epetra_FECrsMatrix_ID_t
    integer(ForTrilinos_Type_ID_t) :: type; integer(c_int) :: index; integer(FT_boolean_t) :: is_const
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

  ! ProblemDifficultyLevel

  integer(kind(c_int)) ,parameter :: FT_ProblemDifficultyLevel_E_t = c_int

  enum ,bind(C)
    enumerator ::                           &
      FT_ProblemDifficultyLevel_E_easy,     &
      FT_ProblemDifficultyLevel_E_moderate, &
      FT_ProblemDifficultyLevel_E_hard,     &
      FT_ProblemDifficultyLevel_E_unsure
  end enum

  ! EParseCommandLineReturn

  integer(kind(c_int)) ,parameter :: FT_EParseCommandLineReturn_E_t = c_int

  enum ,bind(C)
    enumerator ::             &
      FT_EParseCommandLineReturn_E_PARSE_SUCCESSFUL = 0,   &
      FT_EParseCommandLineReturn_E_PARSE_HELP_PRINTED = 1, &
      FT_EParseCommandLineReturn_E_UNRECOGNIZED_OPTION = 2
  end enum

  ! EOptType

  integer(kind(c_int)) ,parameter :: FT_EOptType_E_t = c_int

  enum ,bind(C)
    enumerator ::     &
      FT_EOptType_E_OPT_NONE,       &
      FT_EOptType_E_OPT_BOOL_TRUE,  &
      FT_EOptType_E_OPT_BOOL_FALSE, &
      FT_EOptType_E_OPT_INT,        &
      FT_EOptType_E_OPT_DOUBLE,     &
      FT_EOptType_E_OPT_STRING,     &
      FT_EOptType_E_OPT_ENUM_INT
  end enum

  ! EValidateUsed

  integer(kind(c_int)) ,parameter :: FT_EValidateUsed_E_t = c_int

  enum ,bind(C)
    enumerator ::            &
      FT_EValidateUsed_E_VALIDATE_USED_ENABLED, &
      FT_EValidateUsed_E_VALIDATE_USED_DISABLED
  end enum

  ! EValidateDefaults

  integer(kind(c_int)) ,parameter :: FT_EValidateDefaults_E_t = c_int

  enum ,bind(C)
    enumerator ::                &
      FT_EValidateDefaults_E_VALIDATE_DEFAULTS_ENABLED, &
      FT_EValidateDefaults_E_VALIDATE_DEFAULTS_DISABLED
  end enum

end module ForTrilinos_enums
