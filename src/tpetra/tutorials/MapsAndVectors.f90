! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program main

#include "FortranTestMacros.h"
#include "ForTrilinosTpetra_config.hpp"

use iso_fortran_env
use, intrinsic :: iso_c_binding
use forteuchos
use fortpetra

#ifdef HAVE_MPI
use mpi
#endif

implicit none

! -- Parameters
real(scalar_type), parameter :: one=1.d0, fortytwo=42.d0
integer(size_type), parameter :: num_vecs=1

! -- ForTrilinos objects
type(TeuchosComm) :: comm
type(TpetraMap) :: contig_map, contig_map2, contig_map3, cyclic_map
type(TpetraMultiVector) :: x, y, z

! -- Scalars
integer(size_type) :: my_rank, num_procs, k
integer(size_type) :: num_local_entries, num_elements_per_proc
integer(global_size_type) :: num_global_entries, invalid
integer(global_ordinal_type) :: index_base
logical(bool_type) :: zero_out, x_test
real(scalar_type) :: alpha, beta, gamma
real(mag_type) :: the_norm
real(mag_type), allocatable :: norms(:)
integer(global_ordinal_type), allocatable :: element_list(:)

! --------------------------------------------------------------------------- !


! Initialize MPI subsystem, if applicable
#ifdef HAVE_MPI
call MPI_INIT(ierr)
if (ierr /= 0) then
  stop "MPI failed to init"
endif
call comm%create(MPI_COMM_WORLD)
#else
call comm%create()
#endif

my_rank = comm%getRank()
num_procs = comm%getSize()

! Tpetra objects are templated on several parameters.  ForTrilinos
! provides concrete specializations of Tpetra objects and exposes them to
! Fortran users.  The module fortpetra provides named constants of type
! default integer for Tpetra data types that can be used as the kind type
! parameters.

! Tpetra has local and global Maps.  Local maps describe objects that are
! replicated over all participating MPI processes.  Global maps describe
! distributed objects.  You can do imports and exports between local and global
! maps; this is how you would turn locally replicated objects into distributed
! objects and vice versa.
!
! num_local_entries: The local (on the calling MPI process) number of
! entries (indices) in the first Map that we create.  ForTpetra
! expects a size_type for this value.
num_local_entries = 5

! num_global_entries: The total (global, i.e., over all MPI processes) number of
! entries (indices) in the Map.  ForTpetra expects global_size_t for this value.
! This type is at least 64 bits long on 64-bit machines.
!
! For this example, we scale the global number of entries in the Map with the
! number of MPI processes.  That way, you can run this example with any number
! of MPI processes and every process will still have at least one entry
num_global_entries = num_procs * num_local_entries

! Tpetra can index the entries of a Map starting with 0 (C style),
! 1 (Fortran style), or any base you want.  1-based indexing is
! handy when interfacing with Fortran.  We choose 1-based indexing
! here.
index_base = 1

! Create some Maps.  All Map constructors must be called as a collective over
! the input communicator.  Not all Map constructors necessarily require
! communication, but some do, so it's best to treat them all as collectives.
!
! Create a Map that puts the same number of equations on each processor.  The
! resulting Map is "contiguous and uniform."
call contig_map%create(num_global_entries, index_base, comm)

! contig_map is contiguous by construction.  Test this at run time.
EXPECT_TRUE(contig_map%isContiguous())

! contig_map2: Create a Map which is the same as contig_map, but uses
! a different Map constructor.  This one asks for the number of entries on each
! MPI process.  The resulting Map is "contiguous" but not necessarily uniform,
! since the numbers of entries on different MPI processes may differ.  In this
! case, the number of entries on each MPI process is the same, but that doesn't
! always have to be the case.
call contig_map2%create(num_global_entries, num_local_entries, index_base, comm);

! Since contig_map and contig_map2 have the same communicators, and the same
! number of entries on all MPI processes in their communicators, they are "the
! same."
EXPECT_TRUE(contig_map%isSameAs(contig_map2))

! contig_map3: Use the same Map constructor as contig_map3, but don't specify
! the global number of entries.  This is helpful if you only know how many
! entries each MPI process has, but don't know the global number.  Instead of
! num_global_entries, we use -1
invalid = -1
call contig_map3%create(invalid, num_local_entries, index_base, comm)

! Even though we made contig_map3 without specifying the global number of
! entries, it should still be the same as contig_map2.
EXPECT_TRUE(contig_map2%isSameAs(contig_map3))

! Create a Map which has the same number of global entries per process as
! contig_map, but distributes them differently, in round-robin (1-D cyclic)
! fashion instead of contiguously.

! We'll use the version of the Map constructor that takes, on each MPI process,
! a list of the global entries in the Map belonging to that process.  You can
! use this constructor to construct an overlapping (also called "not 1-to-1")
! Map, in which one or more entries are owned by multiple processes.  We don't
! do that here; we make a nonoverlapping (also called "1-to-1") Map.
num_elements_per_proc = 5
allocate(element_list(num_elements_per_proc))
do k = 1, num_elements_per_proc
  element_list(k) = int(my_rank + k * num_procs, kind=global_ordinal_type)
end do
call cyclic_map%create(num_global_entries, element_list, index_base, comm);
deallocate(element_list)

! If there's more than one MPI process in the communicator, then cyclic_map is
! definitely NOT contiguous.
if (num_procs > 1) then
  x_test = .not. cyclic_map%isContiguous()
  EXPECT_TRUE(x_test)
end if

! contig_map and cyclic_map should always be compatible.  However, if the
! communicator contains more than 1 process, then contig_map and cyclic_map are
! NOT the same.
EXPECT_TRUE(contig_map%isCompatible(cyclic_map))

if (num_procs > 1) then
  x_test = .not. contig_map%isSameAs(cyclic_map)
  EXPECT_TRUE(x_test)
end if

! We have maps now, so we can create vectors.

! Create a Vector with the contiguous Map.  This version of the constructor will
! fill in the vector with zeros.
call x%create(contig_map, num_vecs)

! The two-argument copy constructor with second argument Teuchos::Copy performs
! a deep copy.  x and y have the same Map.  The one-argument copy constructor
! does a _shallow_ copy.
call y%create(x, Copy)

! Create a Vector with the 1-D cyclic Map.  Calling the constructor
! with false for the second argument leaves the data uninitialized,
! so that you can fill it later without paying the cost of
! initially filling it with zeros.
zero_out = .false.
call z%create(cyclic_map, num_vecs, zero_out)

! Set the entries of z to (pseudo)random numbers.  Please don't consider this
! a good parallel pseudorandom number generator.
call z%randomize()

! Set the entries of x to all ones.
call x%putScalar(one)

alpha = 3.14159;
beta = 2.71828;
gamma = -10.0;

! x = beta*x + alpha*z
!
! This is a legal operation!  Even though the Maps of x and z are not the same,
! their Maps are compatible.  Whether it makes sense or not depends on your
! application.
call x%update(alpha, z, beta)

call y%putScalar(fortytwo)

! y = gamma*y + alpha*x + beta*z
call y%update(alpha, x, beta, z, gamma)
call y%update(alpha, x, beta)

! Compute the 2-norm of y.
allocate(norms(num_vecs))
call y%norm2(norms)

the_norm = norms(1)

! Print the norm of y on Proc 0.
if (my_rank == 0) then
    write(*, '(A,f8.3)') 'Norm of y:', the_norm
end if

! Release objects created and no longer in use
call z%release()
call y%release()
call x%release()
call cyclic_map%release()
call contig_map3%release()
call contig_map2%release()
call contig_map%release()
call comm%release()
deallocate(norms)

#ifdef HAVE_MPI
! Finalize MPI must be called after releasing all handles
call MPI_FINALIZE(ierr)
EXPECT_EQ(0, ierr)
#endif

end program main
