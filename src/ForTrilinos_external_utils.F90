module ForTrilinos_external_utils
  use iso_c_binding ,only : c_int        ! Kind parameter (precision specifier)
  use ForTrilinos_enums
#include "ForTrilinos_config.h"

  implicit none                          ! Prevent implicit typing

  interface

#ifdef HAVE_MPI

  ! /*! Create an Epetra_MpiComm from Fortran */
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Fortran_Create ( MPI_Fint fcomm );

  type(FT_Epetra_MpiComm_ID_t) function Epetra_MpiComm_Fortran_Create( fcomm ) &
        bind(C,name='Epetra_MpiComm_Fortran_Create')
    import :: FT_Epetra_MpiComm_ID_t

    integer,intent(in),value :: fcomm
  end function

#endif

  end interface
end module ForTrilinos_external_utils
