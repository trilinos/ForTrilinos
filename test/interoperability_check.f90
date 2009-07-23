module interoperability_check
  use ,intrinsic :: iso_c_binding ,only : &
    c_int,c_char,c_double,c_ptr,c_long,c_bool,c_null_char
  use ,intrinsic :: iso_fortran_env ,only : error_unit ,output_unit
  implicit none
  private
  public :: valid_kind_parameters

  character(len=*) ,parameter :: no_fortran_kind= &
    'The companion C processor deﬁnes the corresponding C type, &
    & but there is no interoperating Fortran processor kind.'
  character(len=*) ,parameter :: no_c_type= &
    'The C processor does not deﬁne the corresponding C type.'
  character(len=*) ,parameter :: interoperable = &
    'An interoperating Fortran kind has a corresponding C type.'
  character(len=*) ,parameter :: imprecise= &
    'The C processor’s type does not have a precision equal to the precision of any of the Fortran processor’s real kinds.'
  character(len=*) ,parameter :: limited= &
    'The C processor’s type does not have a range equal to the range of any of the Fortran processor’s real kinds.'
  character(len=*) ,parameter :: limited_and_imprecise= &
    'The C processor’s type has neither the precision nor range of any of the Fortran processor’s real kinds.'
  character(len=*) ,parameter :: not_interoperable_nonspecific = &
    'There is no interoperating Fortran processor kind for unspecified reasons.'

contains
  
  logical function valid_kind_parameters()
    select case(c_long)
      case(-1)
        write(error_unit ,fmt='(2a)') 'c_long error: ',no_fortran_kind
      case(-2)
        write(error_unit ,fmt='(2a)') 'c_long error: ',no_c_type
      case default
        write(output_unit,fmt='(2a)') 'c_long: ',interoperable
    end select
  
    select case(c_double)
      case(-1)
        write(error_unit ,fmt='(2a)') 'c_double error: ',imprecise
      case(-2)
        write(error_unit ,fmt='(2a)') 'c_double error: ',limited
      case(-3)
        write(error_unit ,fmt='(2a)') 'c_double error: ',limited_and_imprecise
      case(-4)
        write(error_unit ,fmt='(2a)') 'c_double error: ',not_interoperable_nonspecific
      case default
        write(output_unit,fmt='(2a)') 'c_double: ',interoperable
    end select
  
    select case(c_bool)
      case(-1)
        write(error_unit ,fmt='(a)') 'c_bool error: invalid value for a logical kind parameter on the processor.'
      case default
        write(output_unit ,fmt='(a)') 'c_bool:  valid value for a logical kind parameter on the processor.'
    end select
  
    select case(c_char)
      case(-1)
        write(error_unit ,fmt='(a)') 'c_char error: invalid value for a character kind parameter on the processor.'
      case default
        write(output_unit ,fmt='(a)') 'c_char:  valid value for a character kind type parameter on the processor.'
    end select
   
    valid_kind_parameters = .true.
  
  end function

end module interoperability_check
