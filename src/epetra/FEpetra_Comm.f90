module FEpetra_Comm
  use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t ,ForTrilinos_Object_ID_t
  use ForTrilinos_universal ,only : universal
  use forepetra
  implicit none
  private               ! Hide everything by default
  public :: epetra_comm ! Expose type/constructor/methods

  type ,extends(universal) :: epetra_comm
    private
    type(FT_Epetra_Comm_ID_t) :: selfID
  contains
    procedure :: clone
    procedure :: generalize => epetra_comm_generalize
    procedure :: assign_struct
    generic   :: assignment(=) => assign_struct 
    procedure :: force_finalization => call_finalize
    final :: finalize
  end type
  
  interface epetra_comm
    procedure constructor
  end interface

contains
  subroutine call_finalize(this)
    class(epetra_comm) ,intent(inout) :: this
    call finalize(this)
  end subroutine

  subroutine assign_struct(lhs,rhs)
    class(epetra_comm) ,intent(inout) :: lhs
    type(FT_Epetra_Comm_ID_t) ,intent(in) :: rhs
    lhs%selfID = rhs
  end subroutine

  ! CTrilinos_Object_ID_t Epetra_Comm_Abstract (CT_Epetra_Comm_ID_t id );

  type(ForTrilinos_Object_ID_t) function epetra_comm_generalize(this) 
    use iso_c_binding ,only : c_loc
    class(epetra_comm) ,intent(in) ,target :: this

   !interface
   !  type(ForTrilinos_Object_ID_t) function FEpetra_Comm_Abstract(object_id_ptr) bind(C,name="for_linking_only")
   !    use iso_c_binding ,only : c_ptr 
   !    import :: ForTrilinos_Object_ID_t
   !    type(c_ptr) ,value :: object_id_ptr
   !  end function
   !end interface

    epetra_comm_generalize = Epetra_Comm_Abstract ( this%selfID )
   !epetra_comm_generalize = FEpetra_Comm_Abstract( c_loc(this%selfID) )
  end function
  
  type(epetra_comm) function constructor(comm)
    class(universal) ,intent(in) :: comm
   !constructor%self_ID = Epetra_Comm_Cast(comm%generalize())
    select type(comm)
     !class is (epetra_comm)
      type is (epetra_comm)
        constructor%selfID = Epetra_Comm_Cast( comm%generalize() )
      class default
        stop 'epetra_comm constructor: argument type not supported'
    end select
  end function

  ! Original C++ prototype:
  ! virtual Epetra_Comm * Clone() const = 0;
  ! CTrilinos prototype:
  ! CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID );

  type(epetra_comm) function clone(this)
    class(epetra_comm) ,intent(in) :: this
    clone%selfID = Epetra_Comm_Clone(this%selfID)
  end function

  ! Original C++ prototype:
  ! virtual ~Epetra_Comm();
  ! CTrilinos prototype:
  ! void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID );

  subroutine finalize(this)
    type(epetra_comm) ,intent(inout) :: this
    call Epetra_Comm_Destroy ( this%selfID ) 
  end subroutine
end module 
