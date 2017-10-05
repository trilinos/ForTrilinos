// Dependencies
%include "Teuchos_RCP.i"
%import <std_string.i>
%import <Teuchos_ArrayView.i>
%import <Teuchos_Comm.i>

%{
#include "Teuchos_RCP.hpp"
#include "Tpetra_Map.hpp"
%}

// ignore Details namespace
%ignore Details;
%ignore Tpetra::Map::getNode;
// Ignore constructor that takes in Kokkos::View
%ignore Tpetra::Map::Map(const global_size_t numGlobalElements,
         const Kokkos::View<const GlobalOrdinal*, device_type>& indexList,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm);
// Ignore the Node versions
%ignore Tpetra::Map::Map(global_size_t numGlobalElements,
         GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
         LocalGlobal lg,
         const Teuchos::RCP<Node> &node);
%ignore Tpetra::Map::Map(global_size_t numGlobalElements,
         size_t numLocalElements,
         GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
         const Teuchos::RCP<Node> &node);
%ignore Tpetra::Map::Map(const global_size_t numGlobalElements,
         const Teuchos::ArrayView<const GlobalOrdinal>& indexList,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<Node>& node);

// POSTPONE
%ignore Tpetra::Map::getMyGlobalIndices; // Map does not expose global_indices_array_type externally, requires using `auto`
%ignore Tpetra::Map::describe;
%ignore Tpetra::Map::getLocalMap;
%ignore Tpetra::Map::removeEmptyProcesses;
%ignore Tpetra::Map::replaceCommWithSubset;

%teuchos_rcp(Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_Map_decl.hpp"

%template(TpetraMap) Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>;
