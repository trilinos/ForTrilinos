#ifdef HAVE_FORTRILINOS_GALERI

module forgaleri
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  implicit none   ! Prevent implicit typing

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/galeri/CGaleri*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

  ! _________________ Galeri_Utils interface bodies _________________



  ! Original C++ prototype:
  ! Epetra_MultiVector* CreateCartesianCoordinates(const string CoordType, const Epetra_BlockMap* BlockMap, Teuchos::ParameterList& List);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Galeri_Utils_CreateCartesianCoordinates ( const char CoordType[], CT_Epetra_BlockMap_ID_t BlockMapID, CT_Teuchos_ParameterList_ID_t ListID );

  type(FT_Epetra_MultiVector_ID_t) function Galeri_Utils_CreateCartesianCoordinates ( &
        CoordType, BlockMapID, ListID ) &
        bind(C,name='Galeri_Utils_CreateCartesianCoordinates')
    import :: FT_Epetra_MultiVector_ID_t ,c_char ,FT_Epetra_BlockMap_ID_t , &
          FT_Teuchos_ParameterList_ID_t
    
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: CoordType
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: BlockMapID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ListID
  end function


  ! Original C++ prototype:
  ! void Solve(const Epetra_LinearProblem Problem);
  ! CTrilinos prototype:
  ! void Galeri_Utils_Solve_LinearProblem ( CT_Epetra_LinearProblem_ID_t ProblemID );

  subroutine Galeri_Utils_Solve_LinearProblem ( ProblemID ) &
        bind(C,name='Galeri_Utils_Solve_LinearProblem')
    import :: FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: ProblemID
  end subroutine


  ! Original C++ prototype:
  ! void Solve(const Epetra_RowMatrix* Matrix, const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS);
  ! CTrilinos prototype:
  ! void Galeri_Utils_Solve_Matrix ( CT_Epetra_RowMatrix_ID_t MatrixID, CT_Epetra_MultiVector_ID_t LHSID, CT_Epetra_MultiVector_ID_t RHSID );

  subroutine Galeri_Utils_Solve_Matrix ( MatrixID, LHSID, RHSID ) &
        bind(C,name='Galeri_Utils_Solve_Matrix')
    import :: FT_Epetra_RowMatrix_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: MatrixID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: LHSID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: RHSID
  end subroutine


  ! Original C++ prototype:
  ! double ComputeNorm(const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS);
  ! CTrilinos prototype:
  ! double Galeri_Utils_ComputeNorm ( CT_Epetra_MultiVector_ID_t LHSID, CT_Epetra_MultiVector_ID_t RHSID );

  real(c_double) function Galeri_Utils_ComputeNorm ( LHSID, RHSID ) &
        bind(C,name='Galeri_Utils_ComputeNorm')
    import :: c_double ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: LHSID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: RHSID
  end function


  ! Original C++ prototype:
  ! double ComputeNorm(const Epetra_RowMatrix* A, const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS);
  ! CTrilinos prototype:
  ! double Galeri_Utils_ComputeNorm_Matrix ( CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t LHSID, CT_Epetra_MultiVector_ID_t RHSID );

  real(c_double) function Galeri_Utils_ComputeNorm_Matrix ( AID, LHSID, RHSID ) &
        bind(C,name='Galeri_Utils_ComputeNorm_Matrix')
    import :: c_double ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: LHSID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: RHSID
  end function


  ! Original C++ prototype:
  ! string toString(const int& x);
  ! CTrilinos prototype:
  ! const char * Galeri_Utils_toString_Int ( int x );

  character(kind=c_char) function Galeri_Utils_toString_Int ( x ) &
        bind(C,name='Galeri_Utils_toString_Int')
    import :: c_char ,c_int
    
    integer(c_int)              ,intent(in)   ,value              :: x
  end function


  ! Original C++ prototype:
  ! string toString(const unsigned int& x);
  ! CTrilinos prototype:
  ! const char * Galeri_Utils_toString_UInt ( unsigned int x );

  character(kind=c_char) function Galeri_Utils_toString_UInt ( x ) &
        bind(C,name='Galeri_Utils_toString_UInt')
    import :: c_char ,c_int
    
    integer(c_int)              ,intent(in)   ,value              :: x
  end function


  ! Original C++ prototype:
  ! string toString(const double& x);
  ! CTrilinos prototype:
  ! const char * Galeri_Utils_toString_Double ( double x );

  character(kind=c_char) function Galeri_Utils_toString_Double ( x ) &
        bind(C,name='Galeri_Utils_toString_Double')
    import :: c_char ,c_double
    
    real(c_double)              ,intent(in)   ,value              :: x
  end function


  ! Original C++ prototype:
  ! void GetNeighboursCartesian2d(const int i, const int nx, const int ny, int & left, int & right, int & lower, int & upper);
  ! CTrilinos prototype:
  ! void Galeri_Utils_GetNeighboursCartesian2d ( const int i, const int nx, const int ny, int * left, int * right, int * lower, int * upper );

  subroutine Galeri_Utils_GetNeighboursCartesian2d ( i, nx, ny, left, right, lower, upper ) &
        bind(C,name='Galeri_Utils_GetNeighboursCartesian2d')
    import :: c_int
    
    integer(c_int)              ,intent(in)   ,value              :: i
    integer(c_int)              ,intent(in)   ,value              :: nx
    integer(c_int)              ,intent(in)   ,value              :: ny
    integer(c_int)              ,intent(inout)                    :: left
    integer(c_int)              ,intent(inout)                    :: right
    integer(c_int)              ,intent(inout)                    :: lower
    integer(c_int)              ,intent(inout)                    :: upper
  end subroutine


  ! Original C++ prototype:
  ! void GetNeighboursCartesian2d(const int i, const int nx, const int ny, int& left, int& right, int& lower, int& upper, int& left2, int& right2, int& lower2, int& upper2);
  ! CTrilinos prototype:
  ! void Galeri_Utils_GetNeighboursCartesian2d_Both ( const int i, const int nx, const int ny, int * left, int * right, int * lower, int * upper, int * left2, int * right2, int * lower2, int * upper2 );

  subroutine Galeri_Utils_GetNeighboursCartesian2d_Both ( i, nx, ny, left, right, lower, &
        upper, left2, right2, lower2, upper2 ) &
        bind(C,name='Galeri_Utils_GetNeighboursCartesian2d_Both')
    import :: c_int
    
    integer(c_int)              ,intent(in)   ,value              :: i
    integer(c_int)              ,intent(in)   ,value              :: nx
    integer(c_int)              ,intent(in)   ,value              :: ny
    integer(c_int)              ,intent(inout)                    :: left
    integer(c_int)              ,intent(inout)                    :: right
    integer(c_int)              ,intent(inout)                    :: lower
    integer(c_int)              ,intent(inout)                    :: upper
    integer(c_int)              ,intent(inout)                    :: left2
    integer(c_int)              ,intent(inout)                    :: right2
    integer(c_int)              ,intent(inout)                    :: lower2
    integer(c_int)              ,intent(inout)                    :: upper2
  end subroutine


  ! Original C++ prototype:
  ! void GetNeighboursCartesian3d(const int i, const int nx, const int ny, const int nz, int& left, int& right, int& lower, int& upper, int& below, int& above);
  ! CTrilinos prototype:
  ! void Galeri_Utils_GetNeighboursCartesian3d ( const int i, const int nx, const int ny, const int nz, int * left, int * right, int * lower, int * upper, int * below, int * above );

  subroutine Galeri_Utils_GetNeighboursCartesian3d ( i, nx, ny, nz, left, right, lower, &
        upper, below, above ) bind(C,name='Galeri_Utils_GetNeighboursCartesian3d')
    import :: c_int
    
    integer(c_int)              ,intent(in)   ,value              :: i
    integer(c_int)              ,intent(in)   ,value              :: nx
    integer(c_int)              ,intent(in)   ,value              :: ny
    integer(c_int)              ,intent(in)   ,value              :: nz
    integer(c_int)              ,intent(inout)                    :: left
    integer(c_int)              ,intent(inout)                    :: right
    integer(c_int)              ,intent(inout)                    :: lower
    integer(c_int)              ,intent(inout)                    :: upper
    integer(c_int)              ,intent(inout)                    :: below
    integer(c_int)              ,intent(inout)                    :: above
  end subroutine


  ! Original C++ prototype:
  ! void PrintStencil2D(const Epetra_CrsMatrix* Matrix, const int nx, const int ny, int GID = -1);
  ! CTrilinos prototype:
  ! void Galeri_Utils_PrintStencil2D ( CT_Epetra_CrsMatrix_ID_t MatrixID, const int nx, const int ny, int GID );

  subroutine Galeri_Utils_PrintStencil2D ( MatrixID, nx, ny, GID ) &
        bind(C,name='Galeri_Utils_PrintStencil2D')
    import :: FT_Epetra_CrsMatrix_ID_t ,c_int
    
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: MatrixID
    integer(c_int)              ,intent(in)   ,value              :: nx
    integer(c_int)              ,intent(in)   ,value              :: ny
    integer(c_int)              ,intent(in)   ,value              :: GID
  end subroutine


  ! _________________ Galeri_Maps interface bodies _________________


  ! Original C++ prototype:
  ! Epetra_Map* CreateMap(string MapType, Epetra_Comm& Comm, Teuchos::ParameterList& List);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Galeri_Maps_CreateMap ( char MapType[], CT_Epetra_Comm_ID_t CommID, CT_Teuchos_ParameterList_ID_t ListID );

  type(FT_Epetra_Map_ID_t) function Galeri_Maps_CreateMap ( MapType, CommID, ListID ) &
        bind(C,name='Galeri_Maps_CreateMap')
    import :: FT_Epetra_Map_ID_t ,c_char ,FT_Epetra_Comm_ID_t , &
          FT_Teuchos_ParameterList_ID_t
    
    character(kind=c_char)                          ,dimension(*) :: MapType
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: CommID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ListID
  end function


  ! _________________ Galeri_CrsMatrices interface bodies _________________


  ! Original C++ prototype:
  ! Epetra_CrsMatrix* CreateCrsMatrix(string MatrixType, const Epetra_Map* Map, Teuchos::ParameterList& List);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Galeri_CrsMatrices_CreateCrsMatrix ( char MatrixType[], CT_Epetra_Map_ID_t MapID, CT_Teuchos_ParameterList_ID_t ListID );

  type(FT_Epetra_CrsMatrix_ID_t) function Galeri_CrsMatrices_CreateCrsMatrix ( MatrixType, &
        MapID, ListID ) bind(C,name='Galeri_CrsMatrices_CreateCrsMatrix')
    import :: FT_Epetra_CrsMatrix_ID_t ,c_char ,FT_Epetra_Map_ID_t , &
          FT_Teuchos_ParameterList_ID_t
    
    character(kind=c_char)                              ,dimension(*) :: MatrixType
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: MapID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ListID
  end function


  end interface
end module forgaleri

#endif /* HAVE_FORTRILINOS_GALERI */
