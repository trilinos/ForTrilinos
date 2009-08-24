module ForTrilinos_universal
  use ForTrilinos_hermetic ,only : hermetic
  use ForTrilinos_enums ,only : ForTrilinos_Object_ID_t
  implicit none
  type ,abstract ,public ,extends(hermetic) :: universal
  contains
    procedure(generalize_interface) ,deferred :: generalize 
   !procedure(specialize_interface) ,deferred :: specialize 
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
