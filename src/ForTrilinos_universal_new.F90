module ForTrilinos_universal_new
!#include "ForTrilinos_config.h"

  ! This module implements a base type that all ForTrilinos derived types (except 'hermetic') extend.
  ! It provides a universal dummy argument class to which any actual argument can be passed in an 

  use ForTrilinos_hermetic_new ,only : hermetic
!  use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t
  implicit none
  type ,abstract ,public ,extends(hermetic) :: universal
  end type

end module
