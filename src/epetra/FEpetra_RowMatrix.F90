module FEpetra_RowMatrix
  use ForTrilinos_universal ,only : universal
  use ForTrilinos_enums !,only: FT_Epetra_RowMatrix_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use FEpetra_Map, only: Epetra_Map
  use FEpetra_MultiVector, only: Epetra_MultiVector
  use FEpetra_Vector, only: Epetra_Vector
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: Epetra_RowMatrix ! Expose type/methods

  type ,abstract ,extends(universal) :: Epetra_RowMatrix
   private
   type(FT_Epetra_RowMatrix_ID_t)         :: RowMatrix_id 
  contains
    ! Developers only
    procedure                                     :: get_EpetraRowMatrix_ID
    procedure                                     :: set_EpetraRowMatrix_ID
    procedure                 ,nopass             :: alias_EpetraRowMatrix_ID
    procedure ,non_overridable                    :: generalize_EpetraRowMatrix
    procedure(EpetraCrsMatrix_assign)   ,deferred :: CrsMatrix_assign 
    procedure(EpetraRowMatrix_assign)   ,deferred :: RowMatrix_assign 
    procedure                                     :: RowMatrix_assign_ID
    generic :: assignment(=)=>RowMatrix_assign,CrsMatrix_assign,RowMatrix_assign_ID
    !Matrix data extraction routines
    procedure(NumMyRowEntries_interface),deferred :: NumMyRowEntries
    procedure(MaxNumEntries_interface)  ,deferred :: MaxNumEntries
    !procedure                           ,deferred :: ExtractMyRowCopy
    !procedure                           ,deferred :: ExtractDiagonalCopy    
    ! Computational Methods
    procedure(Multiply_interface) ,deferred :: Multiply
  !  procedure(Multiply_Vector_interface) ,deferred :: Multiply_Vector
  !  procedure(Multiply_MultiVector_interface) ,deferred :: Multiply_MultiVector
  !  generic   :: Multiply=> Multiply_Vector,Multiply_MultiVector
    !Atribute access functions
    procedure(RowMatrixRowMap_interface),deferred :: RowMatrixRowMap
    !I/O methods
    !Memory Management
    procedure,non_overridable :: force_finalization_EpetraRowMatrix
  end type
  
  abstract interface
    subroutine EpetraRowMatrix_assign(lhs,rhs)
      import:: Epetra_RowMatrix
      class(Epetra_RowMatrix) ,intent(in)    :: rhs
      class(Epetra_RowMatrix) ,intent(inout) :: lhs
    end subroutine
    subroutine EpetraCrsMatrix_assign(lhs,rhs)
      use ForTrilinos_enums, only: FT_Epetra_CrsMatrix_ID_t
      import:: Epetra_RowMatrix
      type(FT_Epetra_CrsMatrix_ID_t) ,intent(in)    :: rhs
      class(Epetra_RowMatrix) ,intent(inout) :: lhs
    end subroutine
    integer(c_int) function NumMyRowEntries_interface(this,MyRow)
      use iso_c_binding, only : c_int
      import:: Epetra_RowMatrix
      class(Epetra_RowMatrix), intent(in) :: this
      integer(c_int),          intent(in) :: MyRow
    end function
    integer(c_int) function MaxNumEntries_interface(this)
      use iso_c_binding, only : c_int
      import:: Epetra_RowMatrix
      class(Epetra_RowMatrix), intent(in) :: this
    end function
  !  subroutine ExtractMyRowCopy_interface(this,MyRow,length,NumEntries,values,indices,error)
  !    use iso_c_binding, only : c_int, c_double
  !    import:: Epetra_RowMatrix
  !    class(Epetra_RowMatrix), intent(in) :: this
  !    integer(c_int),          intent(in) :: MyRow
  !    integer(c_int),          intent(in) :: NumEntries 
  !    integer(c_int),          intent(in) :: length 
  !    real(c_double),dimension(:), intent(out):: values,indices
  !    integer(c_int),optional,intent(out) :: error
  !  end subroutine
  !  subroutine ExtractDiagonalCopy_interface(this,vector,error)
  !    use iso_c_binding, only : c_int
  !    use FEpetra_Vector, only:Epetra_Vector
  !    import:: Epetra_RowMatrix
  !    class(Epetra_RowMatrix), intent(in) :: this
  !    class(Epetra_Vector), intent(out) :: vector
  !    integer(c_int),optional,intent(out) :: error
  !  end subroutine
  !   subroutine Multiply_Vector_interface(this,TransA,x,y,error)
  !    use iso_c_binding, only: c_int
  !    use ForTrilinos_enums, only: FT_boolean_t
  !    import :: Epetra_RowMatrix,Epetra_Vector
  !    class(Epetra_RowMatrix), intent(in) :: this
  !    integer(FT_boolean_t), intent(in) :: TransA
  !    class(Epetra_Vector), intent(in) :: x
  !    class(Epetra_Vector), intent(in) :: y
  !    integer(c_int), optional,intent(inout) :: error
  !   end subroutine
    ! subroutine Multiply_MultiVector_interface(this,TransA,x,y,error)
     subroutine Multiply_interface(this,TransA,x,y,error)
      use iso_c_binding, only: c_int
      import :: Epetra_RowMatrix,Epetra_MultiVector
      class(Epetra_RowMatrix), intent(in) :: this
      logical, intent(in) :: TransA
      class(Epetra_MultiVector), intent(in) :: x
      class(Epetra_MultiVector), intent(in) :: y
      integer(c_int), optional,intent(inout) :: error
     end subroutine
    function RowMatrixRowMap_interface(this) 
     import:: Epetra_RowMatrix,Epetra_Map
     class(Epetra_RowMatrix), intent(in) :: this
     !class(Epetra_Map), allocatable :: RowMatrixRowMap_interface
     type(Epetra_Map) :: RowMatrixRowMap_interface
    end function 
  end interface

  contains
  
  type(FT_Epetra_RowMatrix_ID_t) function get_EpetraRowMatrix_ID(this)
    class(Epetra_RowMatrix) ,intent(in) :: this
    get_EpetraRowMatrix_ID = this%RowMatrix_id
  end function
  
  subroutine set_EpetraRowMatrix_ID(this,id)
    class(Epetra_RowMatrix)        ,intent(inout) :: this
    type(FT_Epetra_RowMatrix_ID_t) ,intent(in)    :: id 
    this%RowMatrix_id=id
  end subroutine 
  
  type(FT_Epetra_RowMatrix_ID_t) function alias_EpetraRowMatrix_ID(generic_id)
    use iso_c_binding, only : c_loc
    use ForTrilinos_table_man
    use ForTrilinos_enums
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_RowMatrix_ID))
    alias_EpetraRowMatrix_ID=degeneralize_EpetraRowMatrix(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize_EpetraRowMatrix(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_RowMatrix) ,intent(in) ,target :: this
   generalize_EpetraRowMatrix = generalize_all( c_loc(this%RowMatrix_id) )
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_RowMatrix) ,intent(in) ,target :: this
   ! generalize_EpetraRowMatrix = Epetra_RowMatrix_Generalize ( this%RowMatrix_id )
   ! ____ Use for CTrilinos function implementation ______
  end function
  
  type(FT_Epetra_RowMatrix_ID_t) function degeneralize_EpetraRowMatrix(generic_id) bind(C)
    !use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_RowMatrix_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)              ,value  :: generic_id
    type(FT_Epetra_RowMatrix_ID_t),pointer:: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraRowMatrix = local_ptr
  end function
 
  subroutine RowMatrix_assign_ID(lhs,rhs)
    use ForTrilinos_enums
    type(FT_Epetra_RowMatrix_ID_t) ,intent(in)   :: rhs
    class(Epetra_RowMatrix)        ,intent(inout):: lhs
    print *,'RowMatrix_assign_ID'
    lhs%RowMatrix_id=rhs
  end subroutine
 
  subroutine force_finalization_EpetraRowMatrix(this)
    class(Epetra_RowMatrix) ,intent(inout) :: this
    print *,'Destroy_RowMatrix'
    call Epetra_RowMatrix_Destroy( this%RowMatrix_id )
  end subroutine
end module 
