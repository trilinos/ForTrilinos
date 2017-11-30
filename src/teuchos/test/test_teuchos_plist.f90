!Copyright 2017, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program main

#include "FortranTestMacros.h"

    use ISO_FORTRAN_ENV
    implicit none

    call test_plist()
contains

subroutine test_plist()
    use ISO_FORTRAN_ENV
    use, intrinsic :: ISO_C_BINDING
    use forteuchos
    implicit none

    type(ParameterList) :: plist, sublist, sublistref
    integer, dimension(6) :: test_int = (/ -1, 1, 3, 3, 5, 7 /)
    real(C_DOUBLE), dimension(4) :: test_dbl = (/ 0.1d0, 1.9d0, -2.0d0, 4.0d0 /)

    integer :: ival
    real(C_DOUBLE) :: dval
    character(kind=C_CHAR, len=16) :: sval

    write(0, *) "Constructing..."
    call plist%create("myname")
    EXPECT_EQ(0, ierr)
    EXPECT_TRUE(c_associated(plist%swigptr))

    ! Test a function that raises an exception
    call load_from_xml(plist, "nonexistent_path.xml")
    EXPECT_EQ(-3, ierr)
    EXPECT_TRUE(c_associated(plist%swigptr))
    ierr = 0

    ! Get and set a vlaue
    call plist%set("myint", 4)
    call plist%get("myint", ival)
    EXPECT_EQ(4, ival)

    call plist%set("mydbl", 1.25d0)
    call plist%get("mydbl", dval)
    EXPECT_EQ(1.25d0, dval)

    call plist%set("intarr", test_int)
    call plist%set("dblarr", test_dbl)

    EXPECT_EQ(6, plist%get_length('intarr'))
    EXPECT_EQ(4, plist%get_length('dblarr'))

    ! Wrong parameter type
    call plist%get("intarr", test_dbl)
    EXPECT_EQ(-3, ierr)
    ierr = 0
    ! Wrong array size
    call plist%get("intarr", test_int(:4))
    EXPECT_EQ(-4, ierr)
    ierr = 0

    call plist%set("deleteme", 123)
    EXPECT_TRUE(plist%is_parameter('deleteme'))
    call plist%remove("deleteme")
    EXPECT_FALSE(plist%is_parameter('deleteme'))

    EXPECT_FALSE(c_associated(sublist%swigptr))
    sublist = plist%sublist("sublist")
    EXPECT_TRUE(c_associated(sublist%swigptr))
    call sublist%set("anotherval", 4.0d0)
    call sublist%set("stringval", "some string!")
    EXPECT_EQ(12, sublist%get_length('stringval'))
    call sublist%get("stringval", sval)
    EXPECT_EQ('some string!', trim(sval))

    ! Set a string that's too long for sval
    call sublist%set("stringval", "the string is too damn long!")
    EXPECT_EQ(0, ierr)
    call sublist%get("stringval", sval)
    EXPECT_EQ(-4, ierr)
    ierr = 0

    call plist%set("sublist2", sublist)

    ! Add a sub-sublist
    sublist = sublist%sublist("subsublist")

    ! Add another parameter on the original sublist
    call sublist%set("late_arrival", 1.0d0)
    call sublist%release()

    ! XXX: make this interface cleaner. Currently getting a reference to a
    ! sublist requires that another parameterlist is allocated first
    call sublistref%create()
    call plist%get("sublist", sublistref)
    call sublistref%set("added_to_copy", 1.0d0)
    write(0, *) "Printing copy of sublist..."
    call sublistref%print()
    call sublistref%release()

    write(0, *) "Printing..."
    call plist%print()
    write(0, *) "Saving to XML file"
    call save_to_xml(plist, "myparams.xml")

    call plist%release()
end subroutine

end program
