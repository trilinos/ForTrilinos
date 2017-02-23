#define DBC_INSIST_FAIL 993299
#define DBC_ASSERT_FAIL 993399
module DBCF_M
!---------------------------------
! Author: Jordan P. Lefebvre, lefebvrejp@ornl.gov
! Description: DBCF Module for fortran to use.
!---------------------------------
     use, intrinsic :: iso_fortran_env, only: error_unit
     use, intrinsic :: iso_fortran_env, only: output_unit
     implicit none
     public
     private :: error_unit
     private :: output_unit
     contains
     subroutine dbcf_dbc_stop(c, file, line)
        implicit none
       character (len = *), intent(in) :: c; ! logic error
       character (len = *), intent(in) :: file; ! file path
       integer, intent(in) :: line; !line number
       flush(error_unit)
       flush(output_unit)
       write(error_unit, '(A,A,A,A,A,I5,A)') "Assertion: ", c, ", failed in ", file, ", line: ", line, "."
       write(output_unit, '(A,A,A,A,A,I5,A)') "Assertion: ", c, ", failed in ", file, ", line: ", line, "."
       flush(error_unit)
       flush(output_unit)
       stop DBC_ASSERT_FAIL;
     end subroutine

     subroutine dbcf_insist_stop(c, m, file, line)
         implicit none
         character (len = *), intent(in) :: c; ! logic error
         character (len = *), intent(in) :: m; ! message passed in
         character (len = *), intent(in) :: file; ! file path
         integer, intent(in) :: line; !line numberi
         write(error_unit, '(A,A,A)') "The following message was provided: '", m, "' "
         write(error_unit, '(A,A,A,A,A,I5,A)') "Insist: ", c, ", failed in ", file, ", line: ", line, "."
         write(output_unit, '(A,A,A)') "The following message was provided: '", m, "' "
         write(output_unit, '(A,A,A,A,A,I5,A)') "Insist: ", c, ", failed in ", file, ", line: ", line, "."
         flush(error_unit)
         flush(output_unit)
         stop DBC_INSIST_FAIL;
     end subroutine
end module DBCF_M
