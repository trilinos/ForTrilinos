module FEpetra_Comm
  use ForTrilinos_universal ,only : universal
  use ForTrilinos_enums !,only: FT_Epetra_Comm_ID_t,ForTrilinos_Object_ID_t
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: epetra_comm ! Expose type/methods

  !type ,abstract, extends(universal)  :: epetra_comm
  type ,abstract  :: epetra_comm
    type(FT_Epetra_Comm_ID_t) :: comm_id
  contains
!    procedure(generalize_interface) ,deferred :: generalize
!    procedure(clone_interface) ,deferred :: clone
    procedure :: get_ID => EpetraComm_ID
  end type
  
!  abstract interface

!    type(ForTrilinos_Object_ID_t) function generalize_interface(this)
!      import :: epetra_comm,ForTrilinos_Object_ID_t
!      class(epetra_comm) ,intent(in) ,target :: this
!    end function

    ! Original C++ prototype:
    ! virtual Epetra_Comm * Clone() const = 0;
    ! CTrilinos prototype:
    ! CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID );
  
 !   function clone_interface(this) 
 !     import :: epetra_comm
 !     class(epetra_comm) ,intent(in) :: this
 !     class(epetra_comm) ,allocatable :: clone_interface
 !  end function

!  end interface

  contains 
    type(FT_Epetra_Comm_ID_t) function EpetraComm_ID(this) 
      class(epetra_comm) ,intent(in) :: this
      EpetraComm_ID = this%comm_id
    end function

end module 
