/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"
%import <std_string.i>
%import <Teuchos_Comm.i>

%{
#include "Teuchos_RCP.hpp"
#include "Tpetra_Map.hpp"
%}

%constant GO TPETRA_GLOBAL_INVALID = -1;
%constant LO TPETRA_LOCAL_INVALID = 0;

// =======================================================================
// Ignore permanently
// =======================================================================
%ignore Tpetra::Map::Map(const global_size_t numGlobalElements, const GlobalOrdinal indexList[],
         const LocalOrdinal indexListSize, const GlobalOrdinal indexBase,
         const Teuchos::RCP< const Teuchos::Comm< int > > &comm);       // superseded by Teuchos::ArrayView version
%ignore Tpetra::Map::getNode;

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::Map::Map(const global_size_t numGlobalElements,
         const Kokkos::View<const GlobalOrdinal*, device_type>& indexList,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm); // needs Kokkos::View
%ignore Tpetra::Map::describe;                  // needs Teuchos::FancyOStream
%ignore Tpetra::Map::getLocalMap;               // no need to expose this yet
%ignore Tpetra::Map::getMyGlobalIndices;        // return type is not exposed externally, requires using `auto`; for now, use getNodeElementList

// Function signatures for local quantities are incorrectly declared as size_t
%apply LO { size_t getNodeNumElements,
            size_t numLocalElements };

// =======================================================================
// Fix ±1 issues
// =======================================================================
%typemap(in)  int localIndex        %{$1 = *$input - 1;%}
%typemap(out) int getMinLocalIndex  %{$result = $1 + 1;%}
%typemap(out) int getMaxLocalIndex  %{$result = $1 + 1;%}
%typemap(out) int getLocalElement   %{$result = $1 + 1;%}
%typemap(argout) const Teuchos::ArrayView<int>& nodeIDList %{
  for (int i = 0; i < $1->size(); i++)
    (*$1)[i]++;
%}
%typemap(argout) const Teuchos::ArrayView<int>& LIDList %{
  for (int i = 0; i < $1->size(); i++)
    (*$1)[i]++;
%}

%extend Tpetra::Map<LO,GO,NO> {
    Map(global_size_t numGlobalElements, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalGlobal lg=GloballyDistributed) {
      return new Tpetra::Map<LO,GO,NO>(numGlobalElements, 1/*indexBase*/, comm, lg);
    }
    Map(global_size_t numGlobalElements, size_t numLocalElements, const Teuchos::RCP<const Teuchos::Comm<int> > &comm) {
      return new Tpetra::Map<LO,GO,NO>(numGlobalElements, numLocalElements, 1/*indexBase*/, comm);
    }
    Map(const global_size_t numGlobalElements, Teuchos::ArrayView<const GO> indexList, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
      return new Tpetra::Map<LO,GO,NO>(numGlobalElements, indexList, 1/*indexBase*/, comm);
    }
}

%ignore Tpetra::Map::Map(global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Node> &node=defaultArgNode<Node>());
%ignore Tpetra::Map::Map(global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node=defaultArgNode<Node>());
%ignore Tpetra::Map::Map(const global_size_t numGlobalElements, const Teuchos::ArrayView< const GlobalOrdinal > &indexList, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=defaultArgNode< Node >());

%include "Tpetra_Map_decl.hpp"

%teuchos_rcp(Tpetra::Map<LO,GO,NO>)
%template(TpetraMap) Tpetra::Map<LO,GO,NO>;
