!!>
! @file:   DBCF.h
! @author: Jordan P. Lefebvre, lefebvrejp@ornl.gov
! @author: Tom Evans, 9te@ornl.gov
! @date Created on October 4, 2011, 7:53 AM
! @brief Design by contract for Fortran
!!
#ifndef DBCF_H
#define DBCF_H
!!>---------------------------------------------------------------------------//
!!
! \page Nemesis_DBC Using the Nemesis Design-by-Contract Macros
!
! \section ddbc Using the Nemesis Design-by-Contract Macros
!
! The assertion macros are intended to be used for validating preconditions
! which must be true in order for following code to be correct, etc.  For
! example,
!
! \code
! Assert( x > 0. );
! y = sqrt(x);
! \endcode
!
! If the assertion fails, the code should just bomb.  Philosophically, it
! should be used to ferret out bugs in preceding code, making sure that prior
! results are within reasonable bounds before proceeding to use those
! results in further computation, etc.
!
! These macros are provided to support the Design By Contract formalism.
! The activation of each macro is keyed off a bit in the DBC macro which can
! be specified on the command line:
!\verbatim
!     Bit     DBC macro affected
!     ---     ------------------
!      0      Require
!      1      Check
!      2      Ensure
!\endverbatim
!
! So for instance, \c -DDBC=7 turns them all on, \c -DDBC=0 turns them all
! off, and \c -DDBC=1 turns on \c Require but turns off \c Check and \c
! Ensure.  The default is to have them all enabled.
!
! The \c Insist macro is akin to the other macros, but it provides the
! opportunity to specify an instructive message.  The idea here is that you
! should use Insist for checking things which are more or less under user
! control.  If the user makes a poor choice, we "insist" that it be
! corrected, providing a corrective hint.
!
! \note We provide a way to eliminate assertions, but not insistings.  The
! idea is that \c Check is used to perform sanity checks during program
! development, which you might want to eliminate during production runs for
! performance sake.  Insist is used for things which really really must be
! true, such as "the file must've been opened", etc.  So, use \c Check for
! things which you want taken out of production codes (like, the check might
! inhibit inlining or something like that), but use Insist for those things
! you want checked even in a production code.
!!
!!>
! \def Require(condition)
!
! Pre-condition checking macro.  On when DBC & 1 is true.
!/
!!>
! \def Check(condition)
!
! Intra-scope checking macro.  On when DBC & 2 is true.
!/
!!>
!\def Ensure(condition)
!
! Post-condition checking macro.  On when DBC & 4 is true.
!/
!!>
! \def Insist(condition, message)
!
! Inviolate check macro.  Insist is always on.
!/

! #include "dbcf_config.h"
! XXX
! Replace this with actual configured macros
#define DBC 0

#if !defined(DBC)
#define DBC 0
#endif
#if DBC & 1
#define REQUIRE_ON
#if defined(__INTEL_COMPILER)
#define Require(c) if (.not.(c)) call dbcf_dbc_stop( #c, __FILE__, __LINE__ )
#else
#define Require(c) if (.not.(c)) call dbcf_dbc_stop( "c", __FILE__, __LINE__ )
#endif
#else
#define Require(c)
#endif

#if DBC & 2
#define CHECK_ON
#if defined(__INTEL_COMPILER)
#define Check(c) if (.not.(c)) call dbcf_dbc_stop( #c, __FILE__, __LINE__ )
#define Assert(c) if (.not.(c)) call dbcf_dbc_stop( #c, __FILE__, __LINE__ )
#else
#define Check(c) if (.not.(c)) call dbcf_dbc_stop( "c", __FILE__, __LINE__ )
#define Assert(c) if (.not.(c)) call dbcf_dbc_stop( "c", __FILE__, __LINE__ )
#endif
#else
#define Check(c)
#define Assert(c)
#endif

#if DBC & 4
#define REMEMBER_ON
#if defined(__INTEL_COMPILER)
#define Ensure(c) if (.not.(c)) call dbcf_dbc_stop( #c, __FILE__, __LINE__ )
#else
#define Ensure(c) if (.not.(c)) call dbcf_dbc_stop( "c", __FILE__, __LINE__ )
#endif
#define Remember(c) c
#else
#define Ensure(c)
#define Remember(c)
#endif

#if defined(__INTEL_COMPILER)
#define Insist(c,m) if (.not.(c)) call dbcf_insist_stop( #c, m, __FILE__, __LINE__ )
#else
#define Insist(c,m) if (.not.(c)) call dbcf_insist_stop( "c", m, __FILE__, __LINE__ )
#endif

#endif
!! DBCF_H */
