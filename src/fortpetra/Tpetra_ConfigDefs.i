/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

//---------------------------------------------------------------------------//
// in "kokkos-kernels/src/Kokkos_ArithTraits.hpp"
// in "kokkos-kernels/src/Kokkos_InnerProductSpaceTraits.hpp"
// ... from
namespace Kokkos {
  namespace Details {

    template<class T>
    class ArithTraits {
    public:
      typedef T val_type;
      typedef T mag_type;
    };

    template<class T>
    class InnerProductSpaceTraits {
    public:
      typedef T dot_type;
    };

  }

  // KokkosKernels also promotes it to non-Details namespace
  using Details::ArithTraits;
}
%template() Kokkos::Details::ArithTraits<SC>;
%template() Kokkos::Details::InnerProductSpaceTraits<SC>;

//---------------------------------------------------------------------------//
%{
#include <Tpetra_ConfigDefs.hpp>
%}

// ignore these defines from TpetraCore_config.h
#define TPETRA_DEPRECATED
#define KOKKOS_INLINE_FUNCTION

// ignore Details namespace
%ignore Tpetra::Details;

// Ignore some enums
%ignore EPrivateComputeViewConstructor;
%ignore EPrivateHostViewConstructor;

%include "Tpetra_ConfigDefs.hpp"

//---------------------------------------------------------------------------//
%{
#include <Tpetra_CombineMode.hpp>
%}

%include "Tpetra_CombineMode.hpp"

