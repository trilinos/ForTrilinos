#ifdef HAVE_FORTRILINOS_AZTECOO

module foraztecoo
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  implicit none   ! Prevent implicit typing

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/aztecoo/CAztecoo*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

  ! _________________ AztecOO interface bodies _________________


  ! Original C++ prototype:
  ! AztecOO(Epetra_Operator * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
  ! CTrilinos prototype:
  ! CT_AztecOO_ID_t AztecOO_Create_FromOperator ( CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID );

  type(FT_AztecOO_ID_t) function AztecOO_Create_FromOperator ( AID, XID, BID ) &
        bind(C,name='AztecOO_Create_FromOperator')
    import :: FT_AztecOO_ID_t ,FT_Epetra_Operator_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
  end function


  ! Original C++ prototype:
  ! AztecOO(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
  ! CTrilinos prototype:
  ! CT_AztecOO_ID_t AztecOO_Create_FromRowMatrix ( CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID );

  type(FT_AztecOO_ID_t) function AztecOO_Create_FromRowMatrix ( AID, XID, BID ) &
        bind(C,name='AztecOO_Create_FromRowMatrix')
    import :: FT_AztecOO_ID_t ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
  end function


  ! Original C++ prototype:
  ! AztecOO(const Epetra_LinearProblem& LinearProblem);
  ! CTrilinos prototype:
  ! CT_AztecOO_ID_t AztecOO_Create_FromLinearProblem ( CT_Epetra_LinearProblem_ID_t LinearProblemID );

  type(FT_AztecOO_ID_t) function AztecOO_Create_FromLinearProblem ( LinearProblemID ) &
        bind(C,name='AztecOO_Create_FromLinearProblem')
    import :: FT_AztecOO_ID_t ,FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: LinearProblemID
  end function


  ! Original C++ prototype:
  ! AztecOO();
  ! CTrilinos prototype:
  ! CT_AztecOO_ID_t AztecOO_Create (  );

  type(FT_AztecOO_ID_t) function AztecOO_Create (  ) bind(C,name='AztecOO_Create')
    import :: FT_AztecOO_ID_t
    
  end function


  ! Original C++ prototype:
  ! AztecOO(const AztecOO& Solver);
  ! CTrilinos prototype:
  ! CT_AztecOO_ID_t AztecOO_Duplicate ( CT_AztecOO_ID_t SolverID );

  type(FT_AztecOO_ID_t) function AztecOO_Duplicate ( SolverID ) &
        bind(C,name='AztecOO_Duplicate')
    import :: FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: SolverID
  end function


  ! Original C++ prototype:
  ! virtual ~AztecOO(void);
  ! CTrilinos prototype:
  ! void AztecOO_Destroy ( CT_AztecOO_ID_t * selfID );

  subroutine AztecOO_Destroy ( selfID ) bind(C,name='AztecOO_Destroy')
    import :: FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)                                         :: selfID
  end subroutine


  ! Original C++ prototype:
  ! int SetProblem(const Epetra_LinearProblem& prob, bool call_SetPrecMatrix=false);
  ! CTrilinos prototype:
  ! int AztecOO_SetProblem ( CT_AztecOO_ID_t selfID, CT_Epetra_LinearProblem_ID_t probID, boolean call_SetPrecMatrix );

  integer(c_int) function AztecOO_SetProblem ( selfID, probID, call_SetPrecMatrix ) &
        bind(C,name='AztecOO_SetProblem')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_LinearProblem_ID_t ,FT_boolean_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: probID
    integer(FT_boolean_t)       ,intent(in)   ,value              :: call_SetPrecMatrix
  end function


  ! Original C++ prototype:
  ! int SetUserOperator(Epetra_Operator * UserOperator);
  ! CTrilinos prototype:
  ! int AztecOO_SetUserOperator ( CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t UserOperatorID );

  integer(c_int) function AztecOO_SetUserOperator ( selfID, UserOperatorID ) &
        bind(C,name='AztecOO_SetUserOperator')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_Operator_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: UserOperatorID
  end function


  ! Original C++ prototype:
  ! int SetUserMatrix(Epetra_RowMatrix * UserMatrix, bool call_SetPrecMatrix=false);
  ! CTrilinos prototype:
  ! int AztecOO_SetUserMatrix ( CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t UserMatrixID, boolean call_SetPrecMatrix );

  integer(c_int) function AztecOO_SetUserMatrix ( selfID, UserMatrixID, call_SetPrecMatrix ) &
        bind(C,name='AztecOO_SetUserMatrix')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_RowMatrix_ID_t ,FT_boolean_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: UserMatrixID
    integer(FT_boolean_t)       ,intent(in)   ,value              :: call_SetPrecMatrix
  end function


  ! Original C++ prototype:
  ! int SetLHS(Epetra_MultiVector * X);
  ! CTrilinos prototype:
  ! int AztecOO_SetLHS ( CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t XID );

  integer(c_int) function AztecOO_SetLHS ( selfID, XID ) bind(C,name='AztecOO_SetLHS')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
  end function


  ! Original C++ prototype:
  ! int SetRHS(Epetra_MultiVector * B);
  ! CTrilinos prototype:
  ! int AztecOO_SetRHS ( CT_AztecOO_ID_t selfID, CT_Epetra_MultiVector_ID_t BID );

  integer(c_int) function AztecOO_SetRHS ( selfID, BID ) bind(C,name='AztecOO_SetRHS')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
  end function


  ! Original C++ prototype:
  ! int SetPrecMatrix(Epetra_RowMatrix * PrecMatrix);
  ! CTrilinos prototype:
  ! int AztecOO_SetPrecMatrix ( CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t PrecMatrixID );

  integer(c_int) function AztecOO_SetPrecMatrix ( selfID, PrecMatrixID ) &
        bind(C,name='AztecOO_SetPrecMatrix')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_RowMatrix_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: PrecMatrixID
  end function


  ! Original C++ prototype:
  ! int SetPrecOperator(Epetra_Operator * PrecOperator);
  ! CTrilinos prototype:
  ! int AztecOO_SetPrecOperator ( CT_AztecOO_ID_t selfID, CT_Epetra_Operator_ID_t PrecOperatorID );

  integer(c_int) function AztecOO_SetPrecOperator ( selfID, PrecOperatorID ) &
        bind(C,name='AztecOO_SetPrecOperator')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_Operator_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: PrecOperatorID
  end function


  ! Original C++ prototype:
  ! int ConstructPreconditioner(double & condest);
  ! CTrilinos prototype:
  ! int AztecOO_ConstructPreconditioner ( CT_AztecOO_ID_t selfID, double * condest );

  integer(c_int) function AztecOO_ConstructPreconditioner ( selfID, condest ) &
        bind(C,name='AztecOO_ConstructPreconditioner')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    real(c_double)              ,intent(inout)                    :: condest
  end function


  ! Original C++ prototype:
  ! int DestroyPreconditioner();
  ! CTrilinos prototype:
  ! int AztecOO_DestroyPreconditioner ( CT_AztecOO_ID_t selfID );

  integer(c_int) function AztecOO_DestroyPreconditioner ( selfID ) &
        bind(C,name='AztecOO_DestroyPreconditioner')
    import :: c_int ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! double Condest() const;
  ! CTrilinos prototype:
  ! double AztecOO_Condest ( CT_AztecOO_ID_t selfID );

  real(c_double) function AztecOO_Condest ( selfID ) bind(C,name='AztecOO_Condest')
    import :: c_double ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! int CheckInput() const;
  ! CTrilinos prototype:
  ! int AztecOO_CheckInput ( CT_AztecOO_ID_t selfID );

  integer(c_int) function AztecOO_CheckInput ( selfID ) bind(C,name='AztecOO_CheckInput')
    import :: c_int ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! Epetra_LinearProblem * GetProblem() const;
  ! CTrilinos prototype:
  ! CT_Epetra_LinearProblem_ID_t AztecOO_GetProblem ( CT_AztecOO_ID_t selfID );

  type(FT_Epetra_LinearProblem_ID_t) function AztecOO_GetProblem ( selfID ) &
        bind(C,name='AztecOO_GetProblem')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! Epetra_Operator * GetUserOperator() const;
  ! CTrilinos prototype:
  ! CT_Epetra_Operator_ID_t AztecOO_GetUserOperator ( CT_AztecOO_ID_t selfID );

  type(FT_Epetra_Operator_ID_t) function AztecOO_GetUserOperator ( selfID ) &
        bind(C,name='AztecOO_GetUserOperator')
    import :: FT_Epetra_Operator_ID_t ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! Epetra_RowMatrix * GetUserMatrix() const;
  ! CTrilinos prototype:
  ! CT_Epetra_RowMatrix_ID_t AztecOO_GetUserMatrix ( CT_AztecOO_ID_t selfID );

  type(FT_Epetra_RowMatrix_ID_t) function AztecOO_GetUserMatrix ( selfID ) &
        bind(C,name='AztecOO_GetUserMatrix')
    import :: FT_Epetra_RowMatrix_ID_t ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! Epetra_Operator * GetPrecOperator() const;
  ! CTrilinos prototype:
  ! CT_Epetra_Operator_ID_t AztecOO_GetPrecOperator ( CT_AztecOO_ID_t selfID );

  type(FT_Epetra_Operator_ID_t) function AztecOO_GetPrecOperator ( selfID ) &
        bind(C,name='AztecOO_GetPrecOperator')
    import :: FT_Epetra_Operator_ID_t ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! Epetra_RowMatrix * GetPrecMatrix() const;
  ! CTrilinos prototype:
  ! CT_Epetra_RowMatrix_ID_t AztecOO_GetPrecMatrix ( CT_AztecOO_ID_t selfID );

  type(FT_Epetra_RowMatrix_ID_t) function AztecOO_GetPrecMatrix ( selfID ) &
        bind(C,name='AztecOO_GetPrecMatrix')
    import :: FT_Epetra_RowMatrix_ID_t ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! Epetra_MultiVector * GetLHS() const;
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t AztecOO_GetLHS ( CT_AztecOO_ID_t selfID );

  type(FT_Epetra_MultiVector_ID_t) function AztecOO_GetLHS ( selfID ) &
        bind(C,name='AztecOO_GetLHS')
    import :: FT_Epetra_MultiVector_ID_t ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! Epetra_MultiVector * GetRHS() const;
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t AztecOO_GetRHS ( CT_AztecOO_ID_t selfID );

  type(FT_Epetra_MultiVector_ID_t) function AztecOO_GetRHS ( selfID ) &
        bind(C,name='AztecOO_GetRHS')
    import :: FT_Epetra_MultiVector_ID_t ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! void PrintLinearSystem(const char* name);
  ! CTrilinos prototype:
  ! void AztecOO_PrintLinearSystem ( CT_AztecOO_ID_t selfID, const char * name );

  subroutine AztecOO_PrintLinearSystem ( selfID, name ) &
        bind(C,name='AztecOO_PrintLinearSystem')
    import :: FT_AztecOO_ID_t ,c_char
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: name
  end subroutine


  ! Original C++ prototype:
  ! int SetParameters(Teuchos::ParameterList& parameterlist, bool cerr_warning_if_unused=false);
  ! CTrilinos prototype:
  ! int AztecOO_SetParameters ( CT_AztecOO_ID_t selfID, CT_Teuchos_ParameterList_ID_t parameterlistID, boolean cerr_warning_if_unused );

  integer(c_int) function AztecOO_SetParameters ( selfID, parameterlistID, &
        cerr_warning_if_unused ) bind(C,name='AztecOO_SetParameters')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Teuchos_ParameterList_ID_t ,FT_boolean_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: parameterlistID
    integer(FT_boolean_t)       ,intent(in)   ,value              :: cerr_warning_if_unused
  end function


  ! Original C++ prototype:
  ! int SetAztecDefaults();
  ! CTrilinos prototype:
  ! int AztecOO_SetAztecDefaults ( CT_AztecOO_ID_t selfID );

  integer(c_int) function AztecOO_SetAztecDefaults ( selfID ) &
        bind(C,name='AztecOO_SetAztecDefaults')
    import :: c_int ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! int SetAztecOption(int option, int value);
  ! CTrilinos prototype:
  ! int AztecOO_SetAztecOption ( CT_AztecOO_ID_t selfID, int option, int value );

  integer(c_int) function AztecOO_SetAztecOption ( selfID, option, value ) &
        bind(C,name='AztecOO_SetAztecOption')
    import :: c_int ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: option
    integer(c_int)              ,intent(in)   ,value              :: value
  end function


  ! Original C++ prototype:
  ! int GetAztecOption(int option);
  ! CTrilinos prototype:
  ! int AztecOO_GetAztecOption ( CT_AztecOO_ID_t selfID, int option );

  integer(c_int) function AztecOO_GetAztecOption ( selfID, option ) &
        bind(C,name='AztecOO_GetAztecOption')
    import :: c_int ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: option
  end function


  ! Original C++ prototype:
  ! int SetAztecParam(int param, double value);
  ! CTrilinos prototype:
  ! int AztecOO_SetAztecParam ( CT_AztecOO_ID_t selfID, int param, double value );

  integer(c_int) function AztecOO_SetAztecParam ( selfID, param, value ) &
        bind(C,name='AztecOO_SetAztecParam')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: param
    real(c_double)              ,intent(in)   ,value              :: value
  end function


  ! Original C++ prototype:
  ! const int* GetAllAztecOptions() const;
  ! CTrilinos prototype:
  ! const int * AztecOO_GetAllAztecOptions ( CT_AztecOO_ID_t selfID );

  type(c_ptr) function AztecOO_GetAllAztecOptions ( selfID ) &
        bind(C,name='AztecOO_GetAllAztecOptions')
    import :: c_ptr ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! const double* GetAllAztecParams() const;
  ! CTrilinos prototype:
  ! const double * AztecOO_GetAllAztecParams ( CT_AztecOO_ID_t selfID );

  type(c_ptr) function AztecOO_GetAllAztecParams ( selfID ) &
        bind(C,name='AztecOO_GetAllAztecParams')
    import :: c_ptr ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! int SetAllAztecOptions(const int * options);
  ! CTrilinos prototype:
  ! int AztecOO_SetAllAztecOptions ( CT_AztecOO_ID_t selfID, const int * options );

  integer(c_int) function AztecOO_SetAllAztecOptions ( selfID, options ) &
        bind(C,name='AztecOO_SetAllAztecOptions')
    import :: c_int ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)         ,dimension(*) :: options
  end function


  ! Original C++ prototype:
  ! int SetAllAztecParams(const double * params);
  ! CTrilinos prototype:
  ! int AztecOO_SetAllAztecParams ( CT_AztecOO_ID_t selfID, const double * params );

  integer(c_int) function AztecOO_SetAllAztecParams ( selfID, params ) &
        bind(C,name='AztecOO_SetAllAztecParams')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    real(c_double)              ,intent(in)         ,dimension(*) :: params
  end function


  ! Original C++ prototype:
  ! int Iterate(int MaxIters, double Tolerance);
  ! CTrilinos prototype:
  ! int AztecOO_Iterate_Current ( CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance );

  integer(c_int) function AztecOO_Iterate_Current ( selfID, MaxIters, Tolerance ) &
        bind(C,name='AztecOO_Iterate_Current')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: MaxIters
    real(c_double)              ,intent(in)   ,value              :: Tolerance
  end function


  ! Original C++ prototype:
  ! int Iterate(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B, int MaxIters, double Tolerance);
  ! CTrilinos prototype:
  ! int AztecOO_Iterate ( CT_AztecOO_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID, int MaxIters, double Tolerance );

  integer(c_int) function AztecOO_Iterate ( selfID, AID, XID, BID, MaxIters, Tolerance ) &
        bind(C,name='AztecOO_Iterate')
    import :: c_int ,FT_AztecOO_ID_t ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_MultiVector_ID_t , &
          c_double
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
    integer(c_int)              ,intent(in)   ,value              :: MaxIters
    real(c_double)              ,intent(in)   ,value              :: Tolerance
  end function


  ! Original C++ prototype:
  ! int recursiveIterate(int MaxIters, double Tolerance);
  ! CTrilinos prototype:
  ! int AztecOO_recursiveIterate ( CT_AztecOO_ID_t selfID, int MaxIters, double Tolerance );

  integer(c_int) function AztecOO_recursiveIterate ( selfID, MaxIters, Tolerance ) &
        bind(C,name='AztecOO_recursiveIterate')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: MaxIters
    real(c_double)              ,intent(in)   ,value              :: Tolerance
  end function


  ! Original C++ prototype:
  ! const double *GetAztecStatus() const;
  ! CTrilinos prototype:
  ! const double * AztecOO_GetAztecStatus ( CT_AztecOO_ID_t selfID );

  type(c_ptr) function AztecOO_GetAztecStatus ( selfID ) &
        bind(C,name='AztecOO_GetAztecStatus')
    import :: c_ptr ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! int SetUseAdaptiveDefaultsTrue();
  ! CTrilinos prototype:
  ! int AztecOO_SetUseAdaptiveDefaultsTrue ( CT_AztecOO_ID_t selfID );

  integer(c_int) function AztecOO_SetUseAdaptiveDefaultsTrue ( selfID ) &
        bind(C,name='AztecOO_SetUseAdaptiveDefaultsTrue')
    import :: c_int ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! int SetAdaptiveParams(int NumTrials, double * athresholds, double * rthresholds, double condestThreshold, double maxFill, int maxKspace);
  ! CTrilinos prototype:
  ! int AztecOO_SetAdaptiveParams ( CT_AztecOO_ID_t selfID, int NumTrials, double * athresholds, double * rthresholds, double condestThreshold, double maxFill, int maxKspace );

  integer(c_int) function AztecOO_SetAdaptiveParams ( selfID, NumTrials, athresholds, &
        rthresholds, condestThreshold, maxFill, maxKspace ) &
        bind(C,name='AztecOO_SetAdaptiveParams')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumTrials
    real(c_double)                                  ,dimension(*) :: athresholds
    real(c_double)                                  ,dimension(*) :: rthresholds
    real(c_double)              ,intent(in)   ,value              :: condestThreshold
    real(c_double)              ,intent(in)   ,value              :: maxFill
    integer(c_int)              ,intent(in)   ,value              :: maxKspace
  end function


  ! Original C++ prototype:
  ! int AdaptiveIterate(int MaxIters, int MaxSolveAttempts, double Tolerance);
  ! CTrilinos prototype:
  ! int AztecOO_AdaptiveIterate ( CT_AztecOO_ID_t selfID, int MaxIters, int MaxSolveAttempts, double Tolerance );

  integer(c_int) function AztecOO_AdaptiveIterate ( selfID, MaxIters, MaxSolveAttempts, &
        Tolerance ) bind(C,name='AztecOO_AdaptiveIterate')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: MaxIters
    integer(c_int)              ,intent(in)   ,value              :: MaxSolveAttempts
    real(c_double)              ,intent(in)   ,value              :: Tolerance
  end function


  ! Original C++ prototype:
  ! int NumIters() const;
  ! CTrilinos prototype:
  ! int AztecOO_NumIters ( CT_AztecOO_ID_t selfID );

  integer(c_int) function AztecOO_NumIters ( selfID ) bind(C,name='AztecOO_NumIters')
    import :: c_int ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! double TrueResidual() const;
  ! CTrilinos prototype:
  ! double AztecOO_TrueResidual ( CT_AztecOO_ID_t selfID );

  real(c_double) function AztecOO_TrueResidual ( selfID ) &
        bind(C,name='AztecOO_TrueResidual')
    import :: c_double ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! double ScaledResidual() const;
  ! CTrilinos prototype:
  ! double AztecOO_ScaledResidual ( CT_AztecOO_ID_t selfID );

  real(c_double) function AztecOO_ScaledResidual ( selfID ) &
        bind(C,name='AztecOO_ScaledResidual')
    import :: c_double ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! double RecursiveResidual() const;
  ! CTrilinos prototype:
  ! double AztecOO_RecursiveResidual ( CT_AztecOO_ID_t selfID );

  real(c_double) function AztecOO_RecursiveResidual ( selfID ) &
        bind(C,name='AztecOO_RecursiveResidual')
    import :: c_double ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! double SolveTime() const;
  ! CTrilinos prototype:
  ! double AztecOO_SolveTime ( CT_AztecOO_ID_t selfID );

  real(c_double) function AztecOO_SolveTime ( selfID ) bind(C,name='AztecOO_SolveTime')
    import :: c_double ,FT_AztecOO_ID_t
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! int GetAllAztecStatus(double * status);
  ! CTrilinos prototype:
  ! int AztecOO_GetAllAztecStatus ( CT_AztecOO_ID_t selfID, double * status );

  integer(c_int) function AztecOO_GetAllAztecStatus ( selfID, status ) &
        bind(C,name='AztecOO_GetAllAztecStatus')
    import :: c_int ,FT_AztecOO_ID_t ,c_double
    
    type(FT_AztecOO_ID_t)       ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: status
  end function


  ! _________________ AztecOO_StatusTest interface bodies _________________


  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest_Cast ( CTrilinos_Object_ID_t id );

  type(FT_AztecOO_StatusTest_ID_t) function AztecOO_StatusTest_Cast ( id ) &
        bind(C,name='AztecOO_StatusTest_Cast')
    import :: FT_AztecOO_StatusTest_ID_t ,ForTrilinos_Object_ID_t
    
    type(ForTrilinos_Object_ID_t)   ,intent(in)   ,value              :: id
  end function


  ! CTrilinos prototype:
  ! CTrilinos_Object_ID_t AztecOO_StatusTest_Abstract ( CT_AztecOO_StatusTest_ID_t id );

  type(ForTrilinos_Object_ID_t) function AztecOO_StatusTest_Abstract ( id ) &
        bind(C,name='AztecOO_StatusTest_Abstract')
    import :: ForTrilinos_Object_ID_t ,FT_AztecOO_StatusTest_ID_t
    
    type(FT_AztecOO_StatusTest_ID_t),intent(in)   ,value              :: id
  end function


  ! Original C++ prototype:
  ! virtual ~AztecOO_StatusTest();
  ! CTrilinos prototype:
  ! void AztecOO_StatusTest_Destroy ( CT_AztecOO_StatusTest_ID_t * selfID );

  subroutine AztecOO_StatusTest_Destroy ( selfID ) &
        bind(C,name='AztecOO_StatusTest_Destroy')
    import :: FT_AztecOO_StatusTest_ID_t
    
    type(FT_AztecOO_StatusTest_ID_t)                                  :: selfID
  end subroutine


  ! Original C++ prototype:
  ! virtual bool ResidualVectorRequired() const = 0;
  ! CTrilinos prototype:
  ! boolean AztecOO_StatusTest_ResidualVectorRequired ( CT_AztecOO_StatusTest_ID_t selfID );

  integer(FT_boolean_t) function AztecOO_StatusTest_ResidualVectorRequired ( selfID ) &
        bind(C,name='AztecOO_StatusTest_ResidualVectorRequired')
    import :: FT_boolean_t ,FT_AztecOO_StatusTest_ID_t
    
    type(FT_AztecOO_StatusTest_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst, bool SolutionUpdated) = 0;
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusType_E_t AztecOO_StatusTest_CheckStatus ( CT_AztecOO_StatusTest_ID_t selfID, int CurrentIter, CT_Epetra_MultiVector_ID_t CurrentResVectorID, double CurrentResNormEst, boolean SolutionUpdated );

  integer(FT_AztecOO_StatusType_E_t) function AztecOO_StatusTest_CheckStatus ( selfID, &
        CurrentIter, CurrentResVectorID, CurrentResNormEst, SolutionUpdated ) &
        bind(C,name='AztecOO_StatusTest_CheckStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTest_ID_t ,c_int , &
          FT_Epetra_MultiVector_ID_t ,c_double ,FT_boolean_t
    
    type(FT_AztecOO_StatusTest_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: CurrentIter
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: CurrentResVectorID
    real(c_double)                  ,intent(in)   ,value              :: CurrentResNormEst
    integer(FT_boolean_t)           ,intent(in)   ,value              :: SolutionUpdated
  end function


  ! Original C++ prototype:
  ! virtual AztecOO_StatusType GetStatus() const = 0;
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusType_E_t AztecOO_StatusTest_GetStatus ( CT_AztecOO_StatusTest_ID_t selfID );

  integer(FT_AztecOO_StatusType_E_t) function AztecOO_StatusTest_GetStatus ( selfID ) &
        bind(C,name='AztecOO_StatusTest_GetStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTest_ID_t
    
    type(FT_AztecOO_StatusTest_ID_t),intent(in)   ,value              :: selfID
  end function


  ! _________________ AztecOO_StatusTestCombo interface bodies _________________


  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Cast ( CTrilinos_Object_ID_t id );

  type(FT_AztecOO_StatusTestCombo_ID_t) function AztecOO_StatusTestCombo_Cast ( id ) &
        bind(C,name='AztecOO_StatusTestCombo_Cast')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,ForTrilinos_Object_ID_t
    
    type(ForTrilinos_Object_ID_t)        ,intent(in)   ,value              :: id
  end function


  ! CTrilinos prototype:
  ! CTrilinos_Object_ID_t AztecOO_StatusTestCombo_Abstract ( CT_AztecOO_StatusTestCombo_ID_t id );

  type(ForTrilinos_Object_ID_t) function AztecOO_StatusTestCombo_Abstract ( id ) &
        bind(C,name='AztecOO_StatusTestCombo_Abstract')
    import :: ForTrilinos_Object_ID_t ,FT_AztecOO_StatusTestCombo_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: id
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusTestCombo(ComboType t);
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create ( CT_ComboType_E_t t );

  type(FT_AztecOO_StatusTestCombo_ID_t) function AztecOO_StatusTestCombo_Create ( t ) &
        bind(C,name='AztecOO_StatusTestCombo_Create')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,FT_ComboType_E_t
    
    integer(FT_ComboType_E_t)            ,intent(in)   ,value              :: t
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a);
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_OneTest ( CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID );

  type(FT_AztecOO_StatusTestCombo_ID_t) function AztecOO_StatusTestCombo_Create_OneTest ( t, &
        aID ) bind(C,name='AztecOO_StatusTestCombo_Create_OneTest')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,FT_ComboType_E_t , &
          FT_AztecOO_StatusTest_ID_t
    
    integer(FT_ComboType_E_t)            ,intent(in)   ,value              :: t
    type(FT_AztecOO_StatusTest_ID_t)     ,intent(in)   ,value              :: aID
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a, AztecOO_StatusTest& b);
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_Create_TwoTests ( CT_ComboType_E_t t, CT_AztecOO_StatusTest_ID_t aID, CT_AztecOO_StatusTest_ID_t bID );

  type(FT_AztecOO_StatusTestCombo_ID_t) function AztecOO_StatusTestCombo_Create_TwoTests ( &
        t, aID, bID ) bind(C,name='AztecOO_StatusTestCombo_Create_TwoTests')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,FT_ComboType_E_t , &
          FT_AztecOO_StatusTest_ID_t
    
    integer(FT_ComboType_E_t)            ,intent(in)   ,value              :: t
    type(FT_AztecOO_StatusTest_ID_t)     ,intent(in)   ,value              :: aID
    type(FT_AztecOO_StatusTest_ID_t)     ,intent(in)   ,value              :: bID
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusTestCombo& AddStatusTest(AztecOO_StatusTest& a);
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo_AddStatusTest ( CT_AztecOO_StatusTestCombo_ID_t selfID, CT_AztecOO_StatusTest_ID_t aID );

  type(FT_AztecOO_StatusTestCombo_ID_t) function AztecOO_StatusTestCombo_AddStatusTest ( &
        selfID, aID ) bind(C,name='AztecOO_StatusTestCombo_AddStatusTest')
    import :: FT_AztecOO_StatusTestCombo_ID_t ,FT_AztecOO_StatusTest_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
    type(FT_AztecOO_StatusTest_ID_t)     ,intent(in)   ,value              :: aID
  end function


  ! Original C++ prototype:
  ! virtual ~AztecOO_StatusTestCombo();
  ! CTrilinos prototype:
  ! void AztecOO_StatusTestCombo_Destroy ( CT_AztecOO_StatusTestCombo_ID_t * selfID );

  subroutine AztecOO_StatusTestCombo_Destroy ( selfID ) &
        bind(C,name='AztecOO_StatusTestCombo_Destroy')
    import :: FT_AztecOO_StatusTestCombo_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t)                                  :: selfID
  end subroutine


  ! Original C++ prototype:
  ! bool ResidualVectorRequired() const;
  ! CTrilinos prototype:
  ! boolean AztecOO_StatusTestCombo_ResidualVectorRequired ( CT_AztecOO_StatusTestCombo_ID_t selfID );

  integer(FT_boolean_t) function AztecOO_StatusTestCombo_ResidualVectorRequired ( selfID ) &
        bind(C,name='AztecOO_StatusTestCombo_ResidualVectorRequired')
    import :: FT_boolean_t ,FT_AztecOO_StatusTestCombo_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst, bool SolutionUpdated);
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_CheckStatus ( CT_AztecOO_StatusTestCombo_ID_t selfID, int CurrentIter, CT_Epetra_MultiVector_ID_t CurrentResVectorID, double CurrentResNormEst, boolean SolutionUpdated );

  integer(FT_AztecOO_StatusType_E_t) function AztecOO_StatusTestCombo_CheckStatus ( selfID, &
        CurrentIter, CurrentResVectorID, CurrentResNormEst, SolutionUpdated ) &
        bind(C,name='AztecOO_StatusTestCombo_CheckStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestCombo_ID_t ,c_int , &
          FT_Epetra_MultiVector_ID_t ,c_double ,FT_boolean_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                       ,intent(in)   ,value              :: CurrentIter
    type(FT_Epetra_MultiVector_ID_t)     ,intent(in)   ,value              :: CurrentResVectorID
    real(c_double)                       ,intent(in)   ,value              :: CurrentResNormEst
    integer(FT_boolean_t)                ,intent(in)   ,value              :: SolutionUpdated
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusType GetStatus() const;
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusType_E_t AztecOO_StatusTestCombo_GetStatus ( CT_AztecOO_StatusTestCombo_ID_t selfID );

  integer(FT_AztecOO_StatusType_E_t) function AztecOO_StatusTestCombo_GetStatus ( selfID ) &
        bind(C,name='AztecOO_StatusTestCombo_GetStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestCombo_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! ComboType GetComboType() const;
  ! CTrilinos prototype:
  ! CT_ComboType_E_t AztecOO_StatusTestCombo_GetComboType ( CT_AztecOO_StatusTestCombo_ID_t selfID );

  integer(FT_ComboType_E_t) function AztecOO_StatusTestCombo_GetComboType ( selfID ) &
        bind(C,name='AztecOO_StatusTestCombo_GetComboType')
    import :: FT_ComboType_E_t ,FT_AztecOO_StatusTestCombo_ID_t
    
    type(FT_AztecOO_StatusTestCombo_ID_t),intent(in)   ,value              :: selfID
  end function


  ! _________________ AztecOO_StatusTestMaxIters interface bodies _________________


  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Cast ( CTrilinos_Object_ID_t id );

  type(FT_AztecOO_StatusTestMaxIters_ID_t) function AztecOO_StatusTestMaxIters_Cast ( id ) &
        bind(C,name='AztecOO_StatusTestMaxIters_Cast')
    import :: FT_AztecOO_StatusTestMaxIters_ID_t ,ForTrilinos_Object_ID_t
    
    type(ForTrilinos_Object_ID_t)           ,intent(in)   ,value              :: id
  end function


  ! CTrilinos prototype:
  ! CTrilinos_Object_ID_t AztecOO_StatusTestMaxIters_Abstract ( CT_AztecOO_StatusTestMaxIters_ID_t id );

  type(ForTrilinos_Object_ID_t) function AztecOO_StatusTestMaxIters_Abstract ( id ) &
        bind(C,name='AztecOO_StatusTestMaxIters_Abstract')
    import :: ForTrilinos_Object_ID_t ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: id
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusTestMaxIters(int MaxIters);
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters_Create ( int MaxIters );

  type(FT_AztecOO_StatusTestMaxIters_ID_t) function AztecOO_StatusTestMaxIters_Create ( &
        MaxIters ) bind(C,name='AztecOO_StatusTestMaxIters_Create')
    import :: FT_AztecOO_StatusTestMaxIters_ID_t ,c_int
    
    integer(c_int)                          ,intent(in)   ,value              :: MaxIters
  end function


  ! Original C++ prototype:
  ! virtual ~AztecOO_StatusTestMaxIters();
  ! CTrilinos prototype:
  ! void AztecOO_StatusTestMaxIters_Destroy ( CT_AztecOO_StatusTestMaxIters_ID_t * selfID );

  subroutine AztecOO_StatusTestMaxIters_Destroy ( selfID ) &
        bind(C,name='AztecOO_StatusTestMaxIters_Destroy')
    import :: FT_AztecOO_StatusTestMaxIters_ID_t
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t)                                  :: selfID
  end subroutine


  ! Original C++ prototype:
  ! bool ResidualVectorRequired() const;
  ! CTrilinos prototype:
  ! boolean AztecOO_StatusTestMaxIters_ResidualVectorRequired ( CT_AztecOO_StatusTestMaxIters_ID_t selfID );

  integer(FT_boolean_t) function AztecOO_StatusTestMaxIters_ResidualVectorRequired ( selfID ) &
        bind(C,name='AztecOO_StatusTestMaxIters_ResidualVectorRequired')
    import :: FT_boolean_t ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst, bool SolutionUpdated);
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusType_E_t AztecOO_StatusTestMaxIters_CheckStatus ( CT_AztecOO_StatusTestMaxIters_ID_t selfID, int CurrentIter, CT_Epetra_MultiVector_ID_t CurrentResVectorID, double CurrentResNormEst, boolean SolutionUpdated );

  integer(FT_AztecOO_StatusType_E_t) function AztecOO_StatusTestMaxIters_CheckStatus ( &
        selfID, CurrentIter, CurrentResVectorID, CurrentResNormEst, SolutionUpdated ) &
        bind(C,name='AztecOO_StatusTestMaxIters_CheckStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestMaxIters_ID_t ,c_int , &
          FT_Epetra_MultiVector_ID_t ,c_double ,FT_boolean_t
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                          ,intent(in)   ,value              :: CurrentIter
    type(FT_Epetra_MultiVector_ID_t)        ,intent(in)   ,value              :: CurrentResVectorID
    real(c_double)                          ,intent(in)   ,value              :: CurrentResNormEst
    integer(FT_boolean_t)                   ,intent(in)   ,value              :: SolutionUpdated
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusType GetStatus() const;
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusType_E_t AztecOO_StatusTestMaxIters_GetStatus ( CT_AztecOO_StatusTestMaxIters_ID_t selfID );

  integer(FT_AztecOO_StatusType_E_t) function AztecOO_StatusTestMaxIters_GetStatus ( selfID ) &
        bind(C,name='AztecOO_StatusTestMaxIters_GetStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! int GetMaxIters() const;
  ! CTrilinos prototype:
  ! int AztecOO_StatusTestMaxIters_GetMaxIters ( CT_AztecOO_StatusTestMaxIters_ID_t selfID );

  integer(c_int) function AztecOO_StatusTestMaxIters_GetMaxIters ( selfID ) &
        bind(C,name='AztecOO_StatusTestMaxIters_GetMaxIters')
    import :: c_int ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! int GetNumIters() const;
  ! CTrilinos prototype:
  ! int AztecOO_StatusTestMaxIters_GetNumIters ( CT_AztecOO_StatusTestMaxIters_ID_t selfID );

  integer(c_int) function AztecOO_StatusTestMaxIters_GetNumIters ( selfID ) &
        bind(C,name='AztecOO_StatusTestMaxIters_GetNumIters')
    import :: c_int ,FT_AztecOO_StatusTestMaxIters_ID_t
    
    type(FT_AztecOO_StatusTestMaxIters_ID_t),intent(in)   ,value              :: selfID
  end function


  ! _________________ AztecOO_StatusTestResNorm interface bodies _________________


  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTestResNorm_ID_t AztecOO_StatusTestResNorm_Cast ( CTrilinos_Object_ID_t id );

  type(FT_AztecOO_StatusTestResNorm_ID_t) function AztecOO_StatusTestResNorm_Cast ( id ) &
        bind(C,name='AztecOO_StatusTestResNorm_Cast')
    import :: FT_AztecOO_StatusTestResNorm_ID_t ,ForTrilinos_Object_ID_t
    
    type(ForTrilinos_Object_ID_t)          ,intent(in)   ,value              :: id
  end function


  ! CTrilinos prototype:
  ! CTrilinos_Object_ID_t AztecOO_StatusTestResNorm_Abstract ( CT_AztecOO_StatusTestResNorm_ID_t id );

  type(ForTrilinos_Object_ID_t) function AztecOO_StatusTestResNorm_Abstract ( id ) &
        bind(C,name='AztecOO_StatusTestResNorm_Abstract')
    import :: ForTrilinos_Object_ID_t ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: id
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusTestResNorm(const Epetra_Operator & Operator, const Epetra_Vector & LHS, const Epetra_Vector & RHS,double Tolerance);
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusTestResNorm_ID_t AztecOO_StatusTestResNorm_Create ( CT_Epetra_Operator_ID_t OperatorID, CT_Epetra_Vector_ID_t LHSID, CT_Epetra_Vector_ID_t RHSID, double Tolerance );

  type(FT_AztecOO_StatusTestResNorm_ID_t) function AztecOO_StatusTestResNorm_Create ( &
        OperatorID, LHSID, RHSID, Tolerance ) &
        bind(C,name='AztecOO_StatusTestResNorm_Create')
    import :: FT_AztecOO_StatusTestResNorm_ID_t ,FT_Epetra_Operator_ID_t , &
          FT_Epetra_Vector_ID_t ,c_double
    
    type(FT_Epetra_Operator_ID_t)          ,intent(in)   ,value              :: OperatorID
    type(FT_Epetra_Vector_ID_t)            ,intent(in)   ,value              :: LHSID
    type(FT_Epetra_Vector_ID_t)            ,intent(in)   ,value              :: RHSID
    real(c_double)                         ,intent(in)   ,value              :: Tolerance
  end function


  ! Original C++ prototype:
  ! virtual ~AztecOO_StatusTestResNorm();
  ! CTrilinos prototype:
  ! void AztecOO_StatusTestResNorm_Destroy ( CT_AztecOO_StatusTestResNorm_ID_t * selfID );

  subroutine AztecOO_StatusTestResNorm_Destroy ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_Destroy')
    import :: FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t)                                  :: selfID
  end subroutine


  ! Original C++ prototype:
  ! int DefineResForm( ResType TypeOfResidual, NormType TypeOfNorm, Epetra_Vector * Weights = 0);
  ! CTrilinos prototype:
  ! int AztecOO_StatusTestResNorm_DefineResForm ( CT_AztecOO_StatusTestResNorm_ID_t selfID, CT_ResType_E_t TypeOfResidual, CT_NormType_E_t TypeOfNorm, CT_Epetra_Vector_ID_t WeightsID );

  integer(c_int) function AztecOO_StatusTestResNorm_DefineResForm ( selfID, TypeOfResidual, &
        TypeOfNorm, WeightsID ) bind(C,name='AztecOO_StatusTestResNorm_DefineResForm')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t ,FT_ResType_E_t ,FT_NormType_E_t , &
          FT_Epetra_Vector_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    integer(FT_ResType_E_t)                ,intent(in)   ,value              :: TypeOfResidual
    integer(FT_NormType_E_t)               ,intent(in)   ,value              :: TypeOfNorm
    type(FT_Epetra_Vector_ID_t)            ,intent(in)   ,value              :: WeightsID
  end function


  ! Original C++ prototype:
  ! int DefineScaleForm( ScaleType TypeOfScaling, NormType TypeOfNorm, Epetra_Vector * Weights = 0, double ScaleValue = 1.0);
  ! CTrilinos prototype:
  ! int AztecOO_StatusTestResNorm_DefineScaleForm ( CT_AztecOO_StatusTestResNorm_ID_t selfID, CT_ScaleType_E_t TypeOfScaling, CT_NormType_E_t TypeOfNorm, CT_Epetra_Vector_ID_t WeightsID, double ScaleValue );

  integer(c_int) function AztecOO_StatusTestResNorm_DefineScaleForm ( selfID, TypeOfScaling, &
        TypeOfNorm, WeightsID, ScaleValue ) &
        bind(C,name='AztecOO_StatusTestResNorm_DefineScaleForm')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t ,FT_ScaleType_E_t ,FT_NormType_E_t , &
          FT_Epetra_Vector_ID_t ,c_double
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    integer(FT_ScaleType_E_t)              ,intent(in)   ,value              :: TypeOfScaling
    integer(FT_NormType_E_t)               ,intent(in)   ,value              :: TypeOfNorm
    type(FT_Epetra_Vector_ID_t)            ,intent(in)   ,value              :: WeightsID
    real(c_double)                         ,intent(in)   ,value              :: ScaleValue
  end function


  ! Original C++ prototype:
  ! int ResetTolerance(double Tolerance);
  ! CTrilinos prototype:
  ! int AztecOO_StatusTestResNorm_ResetTolerance ( CT_AztecOO_StatusTestResNorm_ID_t selfID, double Tolerance );

  integer(c_int) function AztecOO_StatusTestResNorm_ResetTolerance ( selfID, Tolerance ) &
        bind(C,name='AztecOO_StatusTestResNorm_ResetTolerance')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t ,c_double
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                         ,intent(in)   ,value              :: Tolerance
  end function


  ! Original C++ prototype:
  ! int SetMaxNumExtraIterations(int maxNumExtraIterations);
  ! CTrilinos prototype:
  ! int AztecOO_StatusTestResNorm_SetMaxNumExtraIterations ( CT_AztecOO_StatusTestResNorm_ID_t selfID, int maxNumExtraIterations );

  integer(c_int) function AztecOO_StatusTestResNorm_SetMaxNumExtraIterations ( selfID, &
        maxNumExtraIterations ) &
        bind(C,name='AztecOO_StatusTestResNorm_SetMaxNumExtraIterations')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                         ,intent(in)   ,value              :: maxNumExtraIterations
  end function


  ! Original C++ prototype:
  ! int GetMaxNumExtraIterations();
  ! CTrilinos prototype:
  ! int AztecOO_StatusTestResNorm_GetMaxNumExtraIterations ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  integer(c_int) function AztecOO_StatusTestResNorm_GetMaxNumExtraIterations ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_GetMaxNumExtraIterations')
    import :: c_int ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! bool ResidualVectorRequired() const;
  ! CTrilinos prototype:
  ! boolean AztecOO_StatusTestResNorm_ResidualVectorRequired ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  integer(FT_boolean_t) function AztecOO_StatusTestResNorm_ResidualVectorRequired ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_ResidualVectorRequired')
    import :: FT_boolean_t ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst, bool SolutionUpdated);
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusType_E_t AztecOO_StatusTestResNorm_CheckStatus ( CT_AztecOO_StatusTestResNorm_ID_t selfID, int CurrentIter, CT_Epetra_MultiVector_ID_t CurrentResVectorID, double CurrentResNormEst, boolean SolutionUpdated );

  integer(FT_AztecOO_StatusType_E_t) function AztecOO_StatusTestResNorm_CheckStatus ( &
        selfID, CurrentIter, CurrentResVectorID, CurrentResNormEst, SolutionUpdated ) &
        bind(C,name='AztecOO_StatusTestResNorm_CheckStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestResNorm_ID_t ,c_int , &
          FT_Epetra_MultiVector_ID_t ,c_double ,FT_boolean_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                         ,intent(in)   ,value              :: CurrentIter
    type(FT_Epetra_MultiVector_ID_t)       ,intent(in)   ,value              :: CurrentResVectorID
    real(c_double)                         ,intent(in)   ,value              :: CurrentResNormEst
    integer(FT_boolean_t)                  ,intent(in)   ,value              :: SolutionUpdated
  end function


  ! Original C++ prototype:
  ! AztecOO_StatusType GetStatus() const;
  ! CTrilinos prototype:
  ! CT_AztecOO_StatusType_E_t AztecOO_StatusTestResNorm_GetStatus ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  integer(FT_AztecOO_StatusType_E_t) function AztecOO_StatusTestResNorm_GetStatus ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_GetStatus')
    import :: FT_AztecOO_StatusType_E_t ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! void ResetStatus();
  ! CTrilinos prototype:
  ! void AztecOO_StatusTestResNorm_ResetStatus ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  subroutine AztecOO_StatusTestResNorm_ResetStatus ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_ResetStatus')
    import :: FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  ! Original C++ prototype:
  ! double GetTolerance() const;
  ! CTrilinos prototype:
  ! double AztecOO_StatusTestResNorm_GetTolerance ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  real(c_double) function AztecOO_StatusTestResNorm_GetTolerance ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_GetTolerance')
    import :: c_double ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! double GetTestValue() const;
  ! CTrilinos prototype:
  ! double AztecOO_StatusTestResNorm_GetTestValue ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  real(c_double) function AztecOO_StatusTestResNorm_GetTestValue ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_GetTestValue')
    import :: c_double ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! double GetResNormValue() const;
  ! CTrilinos prototype:
  ! double AztecOO_StatusTestResNorm_GetResNormValue ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  real(c_double) function AztecOO_StatusTestResNorm_GetResNormValue ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_GetResNormValue')
    import :: c_double ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! double GetScaledNormValue() const;
  ! CTrilinos prototype:
  ! double AztecOO_StatusTestResNorm_GetScaledNormValue ( CT_AztecOO_StatusTestResNorm_ID_t selfID );

  real(c_double) function AztecOO_StatusTestResNorm_GetScaledNormValue ( selfID ) &
        bind(C,name='AztecOO_StatusTestResNorm_GetScaledNormValue')
    import :: c_double ,FT_AztecOO_StatusTestResNorm_ID_t
    
    type(FT_AztecOO_StatusTestResNorm_ID_t),intent(in)   ,value              :: selfID
  end function


  end interface
end module foraztecoo

#endif /* HAVE_FORTRILINOS_AZTECOO */
