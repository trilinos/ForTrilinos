/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%{
#include <Tpetra_Map.hpp>
%}

%constant GO TPETRA_GLOBAL_INVALID = -1;
%constant LO TPETRA_LOCAL_INVALID = 0;

// =======================================================================
// Ignore permanently
// =======================================================================
%ignore Tpetra::Map::Map(const global_size_t numGlobalElements, const GlobalOrdinal indexList[],
         const LocalOrdinal indexListSize, const GlobalOrdinal indexBase,
         const Teuchos::RCP< const Teuchos::Comm< int > > &comm);       // superseded by Teuchos::ArrayView version
%ignore Tpetra::Map::Map(Map<local_ordinal_type, global_ordinal_type, node_type>&&);        // move constructor
%ignore Tpetra::Map::operator=(Map<local_ordinal_type, global_ordinal_type, node_type>&&);  // move assignment

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::Map::getLocalMap;               // no need to expose this yet
%ignore Tpetra::Map::getMyGlobalIndices;        // return type is not exposed externally, requires using `auto`; for now, use getNodeElementList

// =======================================================================
// Typemaps
// =======================================================================

// Function signatures for local quantities are incorrectly declared as size_t
%apply LO { size_t getNodeNumElements,
            size_t numLocalElements };

// Convert from C to/from Fortran indices
%apply int INDEX { int localIndex, int getMinLocalIndex, int getMaxLocalIndex,
int getLocalElement};

%apply const Teuchos::ArrayView<const int>& INDEX { const Teuchos::ArrayView<const LO>& indices }

%apply const Teuchos::ArrayView<int>& INDEX { const Teuchos::ArrayView<LO>& nodeIDList, const Teuchos::ArrayView<LO>& LIDList  }

// Set indexBase to 1 (Fortran) and don't pass it to the user
%typemap(in,numinputs=0) long long indexBase {
	$1 = 1;
}

%include "Tpetra_Map_decl.hpp"

%teuchos_rcp(Tpetra::Map<LO,GO,NO>)
%template(TpetraMap) Tpetra::Map<LO,GO,NO>;
