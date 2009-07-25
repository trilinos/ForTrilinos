module universal_module
  use ForTrilinos_hermetic ,only : hermetic
  implicit none
  type ,abstract ,public ,extend(hermetic) :: universal
  contains
    procedure(generalize_interface) ,deferred :: generalize 
  end type

  ! Implementations of this procedure will take in a specific struct and return
  ! a generic struct that can be used to call a casting function.

  abstract interface
    type(ForTrilinos_Object_ID_t) function generalize_interface(this)
      class(universal) ,intent(in) :: this
    end function
  end interface
end module
