/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%{
#include <Tpetra_CrsGraph.hpp>
%}

%ignore Tpetra::Details::HasDeprecatedMethods2630_WarningThisClassIsNotForUsers;

// =======================================================================
// Ignore permanently
// =======================================================================

%ignore Tpetra::CrsGraph::checkSizes;
%ignore Tpetra::CrsGraph::copyAndPermute;
%ignore Tpetra::CrsGraph::packAndPrepare;
%ignore Tpetra::CrsGraph::pack;
%ignore Tpetra::CrsGraph::pack_functor;
%ignore Tpetra::CrsGraph::unpackAndCombine;
%ignore Tpetra::CrsGraph::SLocalGlobalViews;
%ignore Tpetra::CrsGraph::SLocalGlobalNCViews;

// Depreceated functions
%ignore Tpetra::CrsGraph::getGlobalNumDiags;
%ignore Tpetra::CrsGraph::getNodeNumDiags;
%ignore Tpetra::CrsGraph::isLowerTriangular;
%ignore Tpetra::CrsGraph::isUpperTriangular;

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap,
        const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow,
        const ProfileType pftype = TPETRA_DEFAULT_PROFILE_TYPE,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::DualView
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
        const Teuchos::RCP<const map_type>& colMap,
        const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
        const ProfileType pftype = TPETRA_DEFAULT_PROFILE_TYPE,
        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);    // needs Kokkos::DualView
%ignore Tpetra::CrsGraph::CrsGraph(const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::RCP< const map_type > &colMap,
        const typename local_graph_type::row_map_type &rowPointers,
        const typename local_graph_type::entries_type::non_const_type &columnIndices,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::View; ±1 issue
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::RCP< const map_type > &colMap,
        const local_graph_type &lclGraph,
        const Teuchos::RCP< Teuchos::ParameterList > &params);                  // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::CrsGraph (const local_graph_type &lclGraph,
        const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::RCP< const map_type > &colMap=Teuchos::null,
        const Teuchos::RCP< const map_type > &domainMap=Teuchos::null,
        const Teuchos::RCP< const map_type > &rangeMap=Teuchos::null,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::CrsGraph (CrsGraph<local_ordinal_type, global_ordinal_type, node_type>&&);    // move constructor
%ignore Tpetra::CrsGraph::operator= (CrsGraph<local_ordinal_type, global_ordinal_type, node_type>&&);   // move assignment
%ignore Tpetra::CrsGraph::insertLocalIndices(const local_ordinal_type localRow, const local_ordinal_type numEnt, const local_ordinal_type inds[]);    // prefer ArrayView variant
%ignore Tpetra::CrsGraph::getLocalDiagOffsets (const Kokkos::View< size_t *, device_type, Kokkos::MemoryUnmanaged > &offsets) const;    // needs Kokkos::View
%ignore Tpetra::CrsGraph::getNumEntriesPerLocalRowUpperBound; // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsGraph::getLocalGraph;                      // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::getLocalRowView;                    // ±1 issue

// Returns ArrayView by reference
%ignore Tpetra::CrsGraph::getGlobalRowView;
// "should never be called by user code"
%ignore Tpetra::CrsGraph::setAllIndices;

// =======================================================================
// Fix ±1 issues
// =======================================================================
%apply int INDEX { int localRow, int LocalRow };

%apply const Teuchos::ArrayView<const int>& INDEX { const Teuchos::ArrayView<const LO>& indices }

%apply const Teuchos::ArrayView<int>& INDEX {
    const Teuchos::ArrayView<LO>& indices }

%apply const Teuchos::ArrayRCP<const int>& INDEX {
    const Teuchos::ArrayRCP<size_t>& rowPointers,
    const Teuchos::ArrayRCP<int>& columnIndices
}

%apply const Teuchos::ArrayView<int>& INDEX {
    const Teuchos::ArrayView<LO>& lclColInds
}

%apply int { size_t getNumEntriesInLocalRow,
             size_t getNumAllocatedEntriesInLocalRow}

// =======================================================================
// Make interface more Fortran friendly
// =======================================================================
%extend Tpetra::CrsGraph<LO,GO,NO> {
    void
    getLocalRowCopy (LO lclRow,
                     const Teuchos::ArrayView<LO>& lclColInds,
                     size_t& numColInds) const {
        Tpetra::CrsGraph<LO,GO,NO>::nonconst_local_inds_host_view_type temp_i(lclColInds.data(), lclColInds.size());
        return $self->getLocalRowCopy(lclRow, temp_i, numColInds);
    }
    void
    getGlobalRowCopy (GO gblRow,
                     const Teuchos::ArrayView<GO>& gblColInds,
                     size_t& numColInds) const {
        Tpetra::CrsGraph<LO,GO,NO>::nonconst_global_inds_host_view_type temp_i(gblColInds.data(), gblColInds.size());
        return $self->getGlobalRowCopy(gblRow, temp_i, numColInds);
    }
}

// Add doImport and doExport
%tpetra_extend_with_import_export(Tpetra::CrsGraph<LO,GO,NO>)

%ignore Tpetra::CrsGraph::getNodePackedIndices() const;
%ignore Tpetra::CrsGraph::getLocalDiagOffsets;
%ignore Tpetra::CrsGraph::setAllIndices (const typename local_graph_type::row_map_type &rowPointers, const typename local_graph_type::entries_type::non_const_type &columnIndices);


%include "Tpetra_CrsGraph_decl.hpp"

%teuchos_rcp(Tpetra::CrsGraph<LO,GO,NO>)
%template(TpetraCrsGraph) Tpetra::CrsGraph<LO,GO,NO>;
