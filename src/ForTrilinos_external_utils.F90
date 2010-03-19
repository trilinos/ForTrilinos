module ForTrilinos_external_utils
#include "ForTrilinos_config.h"
  use iso_c_binding ,only : c_int        ! Kind parameter (precision specifier)
  use ForTrilinos_enums

  implicit none                          ! Prevent implicit typing

  interface

#ifdef HAVE_MPI

  ! /*! Create an Epetra_MpiComm from Fortran */
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Fortran_Create ( int fcomm );

  type(FT_Epetra_MpiComm_ID_t) function Epetra_MpiComm_Fortran_Create( fcomm ) &
        bind(C,name='Epetra_MpiComm_Fortran_Create')
    import :: FT_Epetra_MpiComm_ID_t, c_int

    integer(c_int),intent(in),value :: fcomm
  end function

#endif

  end interface
end module ForTrilinos_external_utils
