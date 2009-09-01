program main
  use fortrilinos_utils ,only : valid_kind_parameters
  use iso_fortran_env   ,only : error_unit ,output_unit
  implicit none
  if (valid_kind_parameters(verbose=.true.)) then
    write(output_unit,fmt='(a)') 'interoperability_check.F90: test passed. :)'
  else
    write(error_unit,fmt='(a)') 'interoperability_check.F90: test failed. :('
  end if
end program
