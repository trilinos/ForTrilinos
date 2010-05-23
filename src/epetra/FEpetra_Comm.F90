module FEpetra_Comm
  use ForTrilinos_universal_new ,only : universal
  use ForTrilinos_enums !,only: FT_Epetra_Comm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use forepetra
#include "ForTrilinos_config.h"
  implicit none
  private               ! Hide everything by default
  public :: Epetra_Comm ! Expose type/methods

  type ,abstract ,extends(universal) :: Epetra_Comm
   private
   type(FT_Epetra_Comm_ID_t)         :: comm_id 
  contains
    !Constructors
    procedure(clone_interface)          ,deferred :: clone
    ! Developers only
    procedure,private :: remote_dealloc_EpetraComm
    procedure                                     :: get_EpetraComm_ID
    procedure                                     :: set_EpetraComm_ID
    procedure                 ,nopass             :: alias_EpetraComm_ID
    procedure ,non_overridable                    :: generalize_EpetraComm
!    procedure(EpetraComm_assign)        ,deferred :: Comm_assign 
!    procedure                                     :: Comm_assign_ID
!    procedure(EpetraSerialComm_assign)  ,deferred :: SerialComm_assign 
!#ifdef HAVE_MPI
!    procedure(EpetraMpiComm_assign)     ,deferred :: MpiComm_assign 
!!    generic :: assignment(=)=>MpiComm_assign,SerialComm_assign,Comm_assign,Comm_assign_ID
!#else
!    generic :: assignment(=)=>SerialComm_assign,Comm_assign,Comm_assign_ID
!#endif /* HAVE_MPI */
    !Barrier Methods
    procedure(barrier_interface)          ,deferred          ::barrier
    !Broadcast Methods
    procedure(broadcast_double_interface) ,deferred  ::broadcast_double
    procedure(broadcast_int_interface)    ,deferred  ::broadcast_int
    procedure(broadcast_long_interface)   ,deferred  ::broadcast_long
    procedure(broadcast_char_interface)   ,deferred  ::broadcast_char
    !generic :: broadcast=>broadcast_double,broadcast_int,broadcast_long,broadcast_char
    generic :: Broadcast=>broadcast_double,broadcast_int,broadcast_char
    !Gather Methods
    !generic :: GatherAll=>
    procedure(gather_double_interface)   ,deferred  ::gather_double
    !generic :: GatherAll=>gather_double,gather_int,gather_long,gather_char
    generic :: GatherAll=>gather_double!,gather_int,gather_char
    !Sum Methods
    !generic :: SumAll=>
    !Max/Min Methods
    !generic :: MaxAll=>
    !generic :: MinAll=>
    !Parallel Prefix Methods
    !Attribute Accessor Methods
    procedure(MyPID_interface)           ,deferred::MyPID
    procedure(NumProc_interface)         ,deferred::NumProc
    !Gather/catter and Directory Constructors
    !I/O methods
  end type
  
  abstract interface

    ! Original C++ prototype:
    ! virtual Epetra_Comm * Clone() const = 0;
    ! CTrilinos prototype:
    ! CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID );
  
    function clone_interface(this) 
      import:: Epetra_Comm
      class(Epetra_Comm) ,intent(in)  :: this
      class(Epetra_Comm) ,allocatable :: clone_interface
    end function
!    subroutine EpetraSerialComm_assign(lhs,rhs)
!      use ForTrilinos_enums
!      import:: Epetra_Comm
!      type(FT_Epetra_SerialComm_ID_t),intent(in)    :: rhs
!      class(Epetra_Comm)             ,intent(inout) :: lhs
!    end subroutine
!#ifdef HAVE_MPI
!    subroutine EpetraMpiComm_assign(lhs,rhs)
!      use ForTrilinos_enums
!      import:: Epetra_Comm
!      type(FT_Epetra_MpiComm_ID_t)   ,intent(in)    :: rhs
!      class(Epetra_Comm)             ,intent(inout) :: lhs
!    end subroutine
!#endif /* HAVE_MPI */
!    subroutine EpetraComm_assign(lhs,rhs)
!      use ForTrilinos_enums
!      import:: Epetra_Comm
!      class(Epetra_Comm) ,intent(in)    :: rhs
!      class(Epetra_Comm) ,intent(inout) :: lhs
!    end subroutine
    subroutine barrier_interface(this) 
      import:: Epetra_Comm
      class(Epetra_Comm) ,intent(in)  :: this
    end subroutine
    subroutine broadcast_double_interface(this,MyVals,count,root) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm
      class(Epetra_Comm)           ,intent(in)    :: this
      real(c_double) ,dimension(:) ,intent(inout) :: MyVals
      integer(c_int)               ,intent(in)    :: count
      integer(c_int)               ,intent(in)    :: root
    end subroutine
    subroutine broadcast_int_interface(this,MyVals,count,root) 
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_int) ,dimension(:) ,intent(inout) :: MyVals
      integer(c_int)               ,intent(in)    :: count
      integer(c_int)               ,intent(in)    :: root
    end subroutine
    subroutine broadcast_long_interface(this,MyVals,count,root) 
      use iso_c_binding ,only: c_int,c_long
      import:: Epetra_Comm
      class(Epetra_Comm)           ,intent(in)    :: this
      integer(c_long),dimension(:) ,intent(inout) :: MyVals
      integer(c_int)               ,intent(in)    :: count
      integer(c_int)               ,intent(in)    :: root
    end subroutine
    subroutine broadcast_char_interface(this,MyVals,count,root) 
      use iso_c_binding ,only: c_int,c_char
      import:: Epetra_Comm
      class(Epetra_Comm)                 ,intent(in)    :: this
      character(kind=c_char),dimension(:),intent(inout) :: MyVals
      integer(c_int)                     ,intent(in)    :: count
      integer(c_int)                     ,intent(in)    :: root
    end subroutine
    subroutine gather_double_interface(this,MyVals,AllVals,count) 
      use iso_c_binding ,only: c_int,c_double
      import:: Epetra_Comm
      class(Epetra_Comm)                 ,intent(in)    :: this
      real(c_double), dimension(:)        :: MyVals
      real(c_double), dimension(:)       :: AllVals
      integer(c_int)                     ,intent(in)    :: count
    end subroutine
    integer(c_int) function MyPID_interface(this)
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm
      class(Epetra_Comm), intent(in) :: this
    end function
    integer(c_int) function NumProc_interface(this)
      use iso_c_binding ,only: c_int
      import:: Epetra_Comm
      class(Epetra_Comm), intent(in) :: this
    end function
  end interface

  contains
  
  type(FT_Epetra_Comm_ID_t) function get_EpetraComm_ID(this)
    class(Epetra_Comm) ,intent(in) :: this
    get_EpetraComm_ID = this%comm_id
  end function
  
  subroutine set_EpetraComm_ID(this,id)
    class(Epetra_Comm)        ,intent(inout) :: this
    type(FT_Epetra_Comm_ID_t) ,intent(in)    :: id 
    this%comm_id=id
  end subroutine 
  
  type(FT_Epetra_Comm_ID_t) function alias_EpetraComm_ID(generic_id)
    use iso_c_binding, only : c_loc
    use ForTrilinos_table_man
    use ForTrilinos_enums
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Comm_ID))
    alias_EpetraComm_ID=degeneralize_EpetraComm(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize_EpetraComm(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_Comm) ,intent(in) ,target :: this
   generalize_EpetraComm = generalize_all( c_loc(this%comm_id) )
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_Comm) ,intent(in) ,target :: this
   ! generalize_EpetraComm = Epetra_Comm_Generalize ( this%comm_id )
   ! ____ Use for CTrilinos function implementation ______
  end function
  
  type(FT_Epetra_Comm_ID_t) function degeneralize_EpetraComm(generic_id) bind(C)
    !use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_Comm_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)              ,value  :: generic_id
    type(FT_Epetra_Comm_ID_t),pointer:: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraComm = local_ptr
  end function
 
  subroutine remote_dealloc_EpetraComm(this)
    class(Epetra_Comm) ,intent(inout) :: this
    call Epetra_Comm_Destroy( this%comm_id )
  end subroutine
end module 
