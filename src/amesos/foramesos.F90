#ifdef HAVE_FORTRILINOS_AMESOS

module foramesos
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  implicit none   ! Prevent implicit typing

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/amesos/CAmesos*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

  ! _________________ Amesos_BaseSolver interface bodies _________________


  ! CTrilinos prototype:
  ! CT_Amesos_BaseSolver_ID_t Amesos_BaseSolver_Cast ( CTrilinos_Object_ID_t id );

  type(FT_Amesos_BaseSolver_ID_t) function Amesos_BaseSolver_Cast ( id ) &
        bind(C,name='Amesos_BaseSolver_Cast')
    import :: FT_Amesos_BaseSolver_ID_t ,ForTrilinos_Object_ID_t
    
    type(ForTrilinos_Object_ID_t)  ,intent(in)   ,value              :: id
  end function


  ! CTrilinos prototype:
  ! CTrilinos_Object_ID_t Amesos_BaseSolver_Abstract ( CT_Amesos_BaseSolver_ID_t id );

  type(ForTrilinos_Object_ID_t) function Amesos_BaseSolver_Abstract ( id ) &
        bind(C,name='Amesos_BaseSolver_Abstract')
    import :: ForTrilinos_Object_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: id
  end function


  ! Original C++ prototype:
  ! virtual ~Amesos_BaseSolver();
  ! CTrilinos prototype:
  ! void Amesos_BaseSolver_Destroy ( CT_Amesos_BaseSolver_ID_t * selfID );

  subroutine Amesos_BaseSolver_Destroy ( selfID ) bind(C,name='Amesos_BaseSolver_Destroy')
    import :: FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t)                                  :: selfID
  end subroutine


  ! Original C++ prototype:
  ! virtual int SymbolicFactorization() = 0;
  ! CTrilinos prototype:
  ! int Amesos_BaseSolver_SymbolicFactorization ( CT_Amesos_BaseSolver_ID_t selfID );

  integer(c_int) function Amesos_BaseSolver_SymbolicFactorization ( selfID ) &
        bind(C,name='Amesos_BaseSolver_SymbolicFactorization')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual int NumericFactorization() = 0;
  ! CTrilinos prototype:
  ! int Amesos_BaseSolver_NumericFactorization ( CT_Amesos_BaseSolver_ID_t selfID );

  integer(c_int) function Amesos_BaseSolver_NumericFactorization ( selfID ) &
        bind(C,name='Amesos_BaseSolver_NumericFactorization')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual int Solve() = 0;
  ! CTrilinos prototype:
  ! int Amesos_BaseSolver_Solve ( CT_Amesos_BaseSolver_ID_t selfID );

  integer(c_int) function Amesos_BaseSolver_Solve ( selfID ) &
        bind(C,name='Amesos_BaseSolver_Solve')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual int SetUseTranspose(bool UseTranspose) = 0;
  ! CTrilinos prototype:
  ! int Amesos_BaseSolver_SetUseTranspose ( CT_Amesos_BaseSolver_ID_t selfID, boolean UseTranspose );

  integer(c_int) function Amesos_BaseSolver_SetUseTranspose ( selfID, UseTranspose ) &
        bind(C,name='Amesos_BaseSolver_SetUseTranspose')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t ,FT_boolean_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)          ,intent(in)   ,value              :: UseTranspose
  end function


  ! Original C++ prototype:
  ! virtual bool UseTranspose() const = 0;
  ! CTrilinos prototype:
  ! boolean Amesos_BaseSolver_UseTranspose ( CT_Amesos_BaseSolver_ID_t selfID );

  integer(FT_boolean_t) function Amesos_BaseSolver_UseTranspose ( selfID ) &
        bind(C,name='Amesos_BaseSolver_UseTranspose')
    import :: FT_boolean_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual int SetParameters( Teuchos::ParameterList &ParameterList ) = 0;
  ! CTrilinos prototype:
  ! int Amesos_BaseSolver_SetParameters ( CT_Amesos_BaseSolver_ID_t selfID, CT_Teuchos_ParameterList_ID_t ParameterListID );

  integer(c_int) function Amesos_BaseSolver_SetParameters ( selfID, ParameterListID ) &
        bind(C,name='Amesos_BaseSolver_SetParameters')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ParameterListID
  end function


  ! Original C++ prototype:
  ! virtual const Epetra_LinearProblem* GetProblem() const = 0;
  ! CTrilinos prototype:
  ! CT_Epetra_LinearProblem_ID_t Amesos_BaseSolver_GetProblem ( CT_Amesos_BaseSolver_ID_t selfID );

  type(FT_Epetra_LinearProblem_ID_t) function Amesos_BaseSolver_GetProblem ( selfID ) &
        bind(C,name='Amesos_BaseSolver_GetProblem')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual bool MatrixShapeOK() const = 0;
  ! CTrilinos prototype:
  ! boolean Amesos_BaseSolver_MatrixShapeOK ( CT_Amesos_BaseSolver_ID_t selfID );

  integer(FT_boolean_t) function Amesos_BaseSolver_MatrixShapeOK ( selfID ) &
        bind(C,name='Amesos_BaseSolver_MatrixShapeOK')
    import :: FT_boolean_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual const Epetra_Comm & Comm() const = 0;
  ! CTrilinos prototype:
  ! CT_Epetra_Comm_ID_t Amesos_BaseSolver_Comm ( CT_Amesos_BaseSolver_ID_t selfID );

  type(FT_Epetra_Comm_ID_t) function Amesos_BaseSolver_Comm ( selfID ) &
        bind(C,name='Amesos_BaseSolver_Comm')
    import :: FT_Epetra_Comm_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual int NumSymbolicFact() const = 0;
  ! CTrilinos prototype:
  ! int Amesos_BaseSolver_NumSymbolicFact ( CT_Amesos_BaseSolver_ID_t selfID );

  integer(c_int) function Amesos_BaseSolver_NumSymbolicFact ( selfID ) &
        bind(C,name='Amesos_BaseSolver_NumSymbolicFact')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual int NumNumericFact() const = 0;
  ! CTrilinos prototype:
  ! int Amesos_BaseSolver_NumNumericFact ( CT_Amesos_BaseSolver_ID_t selfID );

  integer(c_int) function Amesos_BaseSolver_NumNumericFact ( selfID ) &
        bind(C,name='Amesos_BaseSolver_NumNumericFact')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual int NumSolve() const = 0;
  ! CTrilinos prototype:
  ! int Amesos_BaseSolver_NumSolve ( CT_Amesos_BaseSolver_ID_t selfID );

  integer(c_int) function Amesos_BaseSolver_NumSolve ( selfID ) &
        bind(C,name='Amesos_BaseSolver_NumSolve')
    import :: c_int ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual void PrintStatus() const = 0;
  ! CTrilinos prototype:
  ! void Amesos_BaseSolver_PrintStatus ( CT_Amesos_BaseSolver_ID_t selfID );

  subroutine Amesos_BaseSolver_PrintStatus ( selfID ) &
        bind(C,name='Amesos_BaseSolver_PrintStatus')
    import :: FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  ! Original C++ prototype:
  ! virtual void PrintTiming() const = 0;
  ! CTrilinos prototype:
  ! void Amesos_BaseSolver_PrintTiming ( CT_Amesos_BaseSolver_ID_t selfID );

  subroutine Amesos_BaseSolver_PrintTiming ( selfID ) &
        bind(C,name='Amesos_BaseSolver_PrintTiming')
    import :: FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  ! Original C++ prototype:
  ! virtual void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  ! CTrilinos prototype:
  ! void Amesos_BaseSolver_setParameterList ( CT_Amesos_BaseSolver_ID_t selfID, CT_Teuchos_ParameterList_ID_t paramListID );

  subroutine Amesos_BaseSolver_setParameterList ( selfID, paramListID ) &
        bind(C,name='Amesos_BaseSolver_setParameterList')
    import :: FT_Amesos_BaseSolver_ID_t ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: paramListID
  end subroutine


  ! Original C++ prototype:
  ! virtual Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_getNonconstParameterList ( CT_Amesos_BaseSolver_ID_t selfID );

  type(FT_Teuchos_ParameterList_ID_t) function Amesos_BaseSolver_getNonconstParameterList ( &
        selfID ) bind(C,name='Amesos_BaseSolver_getNonconstParameterList')
    import :: FT_Teuchos_ParameterList_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_unsetParameterList ( CT_Amesos_BaseSolver_ID_t selfID );

  type(FT_Teuchos_ParameterList_ID_t) function Amesos_BaseSolver_unsetParameterList ( &
        selfID ) bind(C,name='Amesos_BaseSolver_unsetParameterList')
    import :: FT_Teuchos_ParameterList_ID_t ,FT_Amesos_BaseSolver_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! virtual void GetTiming( Teuchos::ParameterList &TimingParameterList ) const;
  ! CTrilinos prototype:
  ! void Amesos_BaseSolver_GetTiming ( CT_Amesos_BaseSolver_ID_t selfID, CT_Teuchos_ParameterList_ID_t TimingParameterListID );

  subroutine Amesos_BaseSolver_GetTiming ( selfID, TimingParameterListID ) &
        bind(C,name='Amesos_BaseSolver_GetTiming')
    import :: FT_Amesos_BaseSolver_ID_t ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Amesos_BaseSolver_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: TimingParameterListID
  end subroutine


  ! _________________ Amesos interface bodies _________________


  ! Original C++ prototype:
  ! Amesos();
  ! CTrilinos prototype:
  ! CT_Amesos_ID_t Amesos_Create (  );

  type(FT_Amesos_ID_t) function Amesos_Create (  ) bind(C,name='Amesos_Create')
    import :: FT_Amesos_ID_t
    
  end function


  ! Original C++ prototype:
  ! ~Amesos();
  ! CTrilinos prototype:
  ! void Amesos_Destroy ( CT_Amesos_ID_t * selfID );

  subroutine Amesos_Destroy ( selfID ) bind(C,name='Amesos_Destroy')
    import :: FT_Amesos_ID_t
    
    type(FT_Amesos_ID_t)                                          :: selfID
  end subroutine


  ! Original C++ prototype:
  ! Amesos_BaseSolver *Create(const char *ClassType, const Epetra_LinearProblem& LinearProblem );
  ! CTrilinos prototype:
  ! CT_Amesos_BaseSolver_ID_t Amesos_CreateSolver ( CT_Amesos_ID_t selfID, const char * ClassType, CT_Epetra_LinearProblem_ID_t LinearProblemID );

  type(FT_Amesos_BaseSolver_ID_t) function Amesos_CreateSolver ( selfID, ClassType, &
        LinearProblemID ) bind(C,name='Amesos_CreateSolver')
    import :: FT_Amesos_BaseSolver_ID_t ,FT_Amesos_ID_t ,c_char , &
          FT_Epetra_LinearProblem_ID_t
    
    type(FT_Amesos_ID_t)        ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: ClassType
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: LinearProblemID
  end function


  ! Original C++ prototype:
  ! bool Query(const char * ClassType);
  ! CTrilinos prototype:
  ! boolean Amesos_Query ( CT_Amesos_ID_t selfID, const char * ClassType );

  integer(FT_boolean_t) function Amesos_Query ( selfID, ClassType ) &
        bind(C,name='Amesos_Query')
    import :: FT_boolean_t ,FT_Amesos_ID_t ,c_char
    
    type(FT_Amesos_ID_t)        ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: ClassType
  end function


  ! Original C++ prototype:
  ! static Teuchos::ParameterList GetValidParameters();
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Amesos_GetValidParameters (  );

  type(FT_Teuchos_ParameterList_ID_t) function Amesos_GetValidParameters (  ) &
        bind(C,name='Amesos_GetValidParameters')
    import :: FT_Teuchos_ParameterList_ID_t
    
  end function


  end interface
end module foramesos

#endif /* HAVE_FORTRILINOS_AMESOS */
