module FEpetra_Comm
  use ForTrilinos_universal ,only : universal
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: epetra_comm ! Expose type/methods

  type ,abstract ,extends(universal) :: epetra_comm
  contains
    procedure(clone_interface) ,deferred :: clone
  end type
  
  abstract interface

    ! Original C++ prototype:
    ! virtual Epetra_Comm * Clone() const = 0;
    ! CTrilinos prototype:
    ! CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID );
  
    function clone_interface(this) 
      import :: epetra_comm
      class(epetra_comm) ,intent(in) :: this
      class(epetra_comm) ,allocatable :: clone_interface
    end function

  end interface
  
end module 
