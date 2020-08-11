//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cmake-modern/ScaleDependencies/CheckTpetraInstantiations.cc
 * \brief  CheckTpetraInstantiations class definitions.
 * \note   Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <TpetraCore_config.h>

#if defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) \
    && !defined(HAVE_TPETRA_INST_INT_LONG_LONG)
#ifdef HAVE_TPETRA_INST_INT_INT
#warning "Template instantiation <int, int> is defined"
#endif
#ifdef HAVE_TPETRA_INST_INT_LONG
#warning "Template instantiation <int, long> is defined"
#endif
#error "SCALE requires the <LO=int, GO=long long> instantiation of Tpetra"
#endif

int main()
{
    return 0;
}

//---------------------------------------------------------------------------//
// end of cmake-modern/ScaleDependencies/CheckTpetraInstantiations.cc
//---------------------------------------------------------------------------//
