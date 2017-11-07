// Dependencies
%include "Teuchos_RCP.i"
%import <Teuchos_ArrayView.i>
%import <Teuchos_Comm.i>

%{
#include "Tpetra_CrsGraph.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================
%ignore Tpetra::CrsGraph::getNode;
%ignore Tpetra::CrsGraph::checkSizes;
%ignore Tpetra::CrsGraph::copyAndPermute;
%ignore Tpetra::CrsGraph::packAndPrepare;
%ignore Tpetra::CrsGraph::pack;
%ignore Tpetra::CrsGraph::pack_functor;
%ignore Tpetra::CrsGraph::unpackAndCombine;
%ignore Tpetra::CrsGraph::SLocalGlobalViews;
%ignore Tpetra::CrsGraph::SLocalGlobalNCViews;

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow, \
        const ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::DualView
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow, \
        ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::DualView
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::ArrayRCP< const size_t > &numEntPerRow, \
        const ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const Teuchos::ArrayRCP< const size_t > &numEntPerRow, \
        ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsGraph::CrsGraph(const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const typename local_graph_type::row_map_type &rowPointers, \
        const typename local_graph_type::entries_type::non_const_type &columnIndices, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::View; ±1 issue
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const Teuchos::ArrayRCP< size_t > &rowPointers, \
        const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Teuchos::ArrayRCP; ±1 issue
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const local_graph_type &lclGraph, \
        const Teuchos::RCP< Teuchos::ParameterList > &params);                  // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::CrsGraph (const local_graph_type &lclGraph, \
        const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap=Teuchos::null, \
        const Teuchos::RCP< const map_type > &domainMap=Teuchos::null, \
        const Teuchos::RCP< const map_type > &rangeMap=Teuchos::null, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::insertLocalIndices;                 // ±1 issue
%ignore Tpetra::CrsGraph::removeLocalIndices;                 // ±1 issue
%ignore Tpetra::CrsGraph::getNumEntriesInLocalRow;            // ±1 issue
%ignore Tpetra::CrsGraph::getNumAllocatedEntriesInLocalRow;   // ±1 issue
%ignore Tpetra::CrsGraph::getLocalRowCopy;                    // ±1 issue
%ignore Tpetra::CrsGraph::getLocalRowView;                    // ±1 issue
%ignore Tpetra::CrsGraph::describe;                           // needs Teuchos::FancyOStream
%ignore Tpetra::CrsGraph::getLocalDiagOffsets (const Kokkos::View< size_t *, device_type, Kokkos::MemoryUnmanaged > &offsets) const;    // needs Kokkos::View
%ignore Tpetra::CrsGraph::getLocalDiagOffsets (Teuchos::ArrayRCP< size_t > &offsets) const;     // needs Teuchos:ArrayRCP
%ignore Tpetra::CrsGraph::getNumEntriesPerLocalRowUpperBound; // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsGraph::setAllIndices;                      // needs Teuchos::ArrayRCP, or Kokkos::View
%ignore Tpetra::CrsGraph::getNodeRowPtrs;                     // needs Teuchos::RCP
%ignore Tpetra::CrsGraph::getNodePackedIndices;               // needs Teuchos::RCP
%ignore Tpetra::CrsGraph::getLocalGraph;                      // needs Kokkos::StaticCrsGraph



%teuchos_rcp(Tpetra::CrsGraph<LO,GO,NO,false>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_CrsGraph_decl.hpp"

%template(TpetraCrsGraph) Tpetra::CrsGraph<LO,GO,NO,false>;