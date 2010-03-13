module ForTrilinos_table_man
  use iso_c_binding ,only : c_int        ! Kind parameter (precision specifier)
  use ForTrilinos_enums

  implicit none                          ! Prevent implicit typing

  interface

  ! /*! Copies the RCP from one table into a second table. The new ID
  !  *  will be returned from the function. Both the old and the new
  !  *  IDs will need to be removed from the tables in order to destroy
  !  *  the object. */
  ! CTrilinos_Universal_ID_t CT_Alias(CTrilinos_Universal_ID_t aid,
  ! CTrilinos_Table_ID_t new_table);

  type(ForTrilinos_Universal_ID_t) function CT_Alias( selfID, new_table ) &
        bind(C,name='CT_Alias')
    import :: ForTrilinos_Universal_ID_t, ForTrilinos_Table_ID_t

    type(ForTrilinos_Universal_ID_t),intent(in),value :: selfID
    integer(ForTrilinos_Table_ID_t) ,intent(in),value :: new_table
  end function

  ! /*! Removes the RCP from one table and puts it in another. *aid will
  !  *  hold the new struct value afterward. Only the new RCP will need
  !  *  to be removed in order to destroy the object. */
  ! void CT_Migrate(CTrilinos_Universal_ID_t *aid, CTrilinos_Table_ID_t new_table);

  subroutine CT_Migrate( selfID, new_table ) &
        bind(C,name='CT_Migrate')
    import :: ForTrilinos_Universal_ID_t, ForTrilinos_Table_ID_t

    type(ForTrilinos_Universal_ID_t)                  :: selfID
    integer(ForTrilinos_Table_ID_t) ,intent(in),value :: new_table
  end subroutine

  ! /*! Checks to see if the underlying object referenced by a table
  !  *  entry is dynamic_cast'able to a given type (can be used to
  !  *  distinguish, e.g., an Epetra_SerialComm from an Epetra_MpiComm). */
  ! boolean CT_TypeCheck(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t type);

  subroutine CT_TypeCheck( selfID, typeid ) &
        bind(C,name='CT_TypeCheck')
    import :: ForTrilinos_Universal_ID_t, ForTrilinos_Table_ID_t

    type(ForTrilinos_Universal_ID_t),intent(in),value :: selfID
    integer(ForTrilinos_Table_ID_t) ,intent(in),value :: typeid
  end subroutine

  end interface
end module ForTrilinos_table_man
