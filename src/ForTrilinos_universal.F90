module ForTrilinos_universal

  ! This module implements a base type that all ForTrilinos derived types (except 'hermetic') extend.
  ! It provides a universal dummy argument class to which any actual argument can be passed in an 
  ! Epetra type-bound procedure.  The deferred binding "generalize" ensures that each Epetra derived
  ! type implements a type-bound procedure that can be invoked to create an equivalent general
  ! entity of derived type ForTrilinos_Object_ID_t, which can then be converted to any other Epetra
  ! derived type.

  use ForTrilinos_hermetic ,only : hermetic
  use ForTrilinos_enums ,only : ForTrilinos_Object_ID_t
  implicit none
  type ,abstract ,public ,extends(hermetic) :: universal
  contains
    procedure(generalize_interface) ,deferred :: generalize 
  end type

  ! Implementations of this procedure will take in a specific struct and return
  ! a generic struct that can be used to call a casting function.

  abstract interface
    type(ForTrilinos_Object_ID_t) function generalize_interface(this)
      import :: universal,ForTrilinos_Object_ID_t
      class(universal) ,intent(in) ,target :: this
    end function
  end interface
end module
