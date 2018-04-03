# ForTrilinos roadmap

## Key

| Emoji            | Definition | Scope |
-------------------| ---------- | ----- |
:white_check_mark: | done                                            |
:star:             | done with modifications                         |
:star2:            | done with modifications that affect performance (e.g., array copy) |
:running:          | work in progress        | current release
:arrow_down:       | low priority            | current release
:x:                | ignored                 | future release/never

## Package status

Status | Package
-------|--------
:running: | [Teuchos](#teuchos)
:running: | [Tpetra](#tpetra)
:running: | [Belos](#belos)

## Teuchos

## Tpetra

Status | Class
-------|--------
:running: | [Map](#tpetramap)
:running: | [Export](#tpetraexport)
:running: | [Import](#tpetraimport)
:running: | [MultiVector](#tpetramultivector)
:running: | [CrsGraph](#tpetracrsgraph)
:running: | [CrsMatrix](#tpetracrsmatrix)
:running: | [Reader](#tpetrareader)
:running: | [Writer](#tpetrawriter)
:running: | [MatrixMatrix](#tpetramatrixmatrix)
:arrow_down: | Vector
:arrow_down: | RowMatrix
:arrow_down: | Operator
:arrow_down: | Distributor

### Tpetra::Map

**Constructor/destructor methods**

:grey_question: | Command | Comment
-------|--------|---------
:star: | `Map(numGlobalElements, indexBase, comm)` | no `indexBase`
:star: | `Map (global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, LocalGlobal lg=GloballyDistributed, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())` | no `indexBase`; no `node`
:star: | `Map (global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())` | no `indexBase`; no `node`
:x: | `Map (const global_size_t numGlobalElements, const Kokkos::View< const GlobalOrdinal *, device_type > &indexList, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)`
:x: | `Map (const global_size_t numGlobalElements, const GlobalOrdinal indexList[], const LocalOrdinal indexListSize, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)` | prefer `ArrayView` variant
:star: | `Map (const global_size_t numGlobalElements, const Teuchos::ArrayView< const GlobalOrdinal > &indexList, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())` | no `node`; no `indexBase`; Fortran arrays
:white_check_mark: | `Map ()`
:white_check_mark: | `~Map ()`

**Boolean tests**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `bool isNodeLocalElement (LocalOrdinal localIndex) const`
:white_check_mark: | `bool isNodeGlobalElement (GlobalOrdinal globalIndex) const`
:white_check_mark: | `bool isUniform () const`
:white_check_mark: | `bool isContiguous () const`
:white_check_mark: | `bool isDistributed () const`
:white_check_mark: | `bool isCompatible (const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const`
:white_check_mark: | `bool isSameAs (const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const`
:white_check_mark: | `bool locallySameAs (const Map< LocalOrdinal, GlobalOrdinal, node_type > &map) const`


:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `Teuchos::RCP< const Teuchos::Comm< int > > 	getComm () const`
:x: | `Teuchos::RCP< Node > 	getNode () const`
:arrow_down: | `std::string 	description () const`
:x: | `void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`
:x: | `template<class NodeOut > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, NodeOut > > 	clone (const Teuchos::RCP< NodeOut > &nodeOut) const`
:white_check_mark: | `Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	removeEmptyProcesses () const `
:white_check_mark: | `Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	replaceCommWithSubset (const Teuchos::RCP< const Teuchos::Comm< int > > &newComm) const`

**Related functions (non-members)**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `template<class LocalOrdinal , class GlobalOrdinal > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal > > 	createLocalMap (const size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)`
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createLocalMapWithNode (const size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::null)`
:x: | `template<class LocalOrdinal , class GlobalOrdinal > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal > > 	createUniformContigMap (const global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)`
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createUniformContigMapWithNode (const global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::null)`
:x: | `template<class LocalOrdinal , class GlobalOrdinal > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal > > 	createContigMap (const global_size_t numElements, const size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)`
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createContigMapWithNode (const global_size_t numElements, const size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())`
:x: | `template<class LocalOrdinal , class GlobalOrdinal > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal > > 	createNonContigMap (const Teuchos::ArrayView< const GlobalOrdinal > &elementList, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)`
:x: |  `emplate<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createNonContigMapWithNode (const Teuchos::ArrayView< const GlobalOrdinal > &elementList, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::null)`
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createWeightedContigMapWithNode (const int thisNodeWeight, const global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::null)`
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createOneToOne (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &M)`
:x: | ` template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createOneToOne (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &M, const Tpetra::Details::TieBreak< LocalOrdinal, GlobalOrdinal > &tie_break)`
:x: | ` template<class LocalOrdinal , class GlobalOrdinal , class Node >bool 	operator== (const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map1, const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map2)`
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > bool 	operator!= (const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map1, const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map2)`

**Attributes**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `bool 	isOneToOne () const`
:white_check_mark: | `global_size_t 	getGlobalNumElements () const`
:white_check_mark: | `size_t 	getNodeNumElements () const`
:x: | `GlobalOrdinal 	getIndexBase () const`
:star: | `LocalOrdinal 	getMinLocalIndex () const` | 1-based indices
:star: | `LocalOrdinal 	getMaxLocalIndex () const` | 1-based indices
:white_check_mark: | `GlobalOrdinal 	getMinGlobalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMaxGlobalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMinAllGlobalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMaxAllGlobalIndex () const`
:star: | `LocalOrdinal 	getLocalElement (GlobalOrdinal globalIndex) const` | 1-based indices
:white_check_mark: | `GlobalOrdinal 	getGlobalElement (LocalOrdinal localIndex) const`
:x: | `local_map_type 	getLocalMap () const`
:star2: | `LookupStatus 	getRemoteIndexList (const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< LocalOrdinal > &LIDList) const` | Fortran arrays; 1-based indices
:star2: | `LookupStatus 	getRemoteIndexList (const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const` | Fortran arrays; 1-based indices
:x: | `global_indices_array_type 	getMyGlobalIndices () const`
:star: | `Teuchos::ArrayView< const GlobalOrdinal > 	getNodeElementList () const` | Fortran arrays


### Tpetra::Export

**Constructor/Destructor Methods**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `Export (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target)`
:x:                | `Export (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::FancyOStream > &out)`
:white_check_mark: | `Export (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::ParameterList > &plist)`
:x:                | `Export (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::FancyOStream > &out, const Teuchos::RCP< Teuchos::ParameterList > &plist)`
:white_check_mark: | `Export (const Export< LocalOrdinal, GlobalOrdinal, Node > &rhs)`
:white_check_mark: | `Export (const Import< LocalOrdinal, GlobalOrdinal, Node > &importer)`
:white_check_mark: | `virtual 	~Export ()`
:white_check_mark: | `void 	setParameterList (const Teuchos::RCP< Teuchos::ParameterList > &plist)`

**Export Attribute Methods**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `size_t 	getNumSameIDs () const`
:white_check_mark: | `size_t 	getNumPermuteIDs () const`
:arrow_down: | `Teuchos::ArrayView< const LocalOrdinal > 	getPermuteFromLIDs () const`
:arrow_down: | `Teuchos::ArrayView< const LocalOrdinal > 	getPermuteToLIDs () const`
:white_check_mark: | `size_t 	getNumRemoteIDs () const`
:arrow_down: | `Teuchos::ArrayView< const LocalOrdinal > 	getRemoteLIDs () const`
:white_check_mark: | `size_t 	getNumExportIDs () const`
:arrow_down: | `Teuchos::ArrayView< const LocalOrdinal > 	getExportLIDs () const`
:arrow_down: | `Teuchos::ArrayView< const int > 	getExportPIDs () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getSourceMap () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getTargetMap () const`
:x: | `Distributor & 	getDistributor () const`
:white_check_mark: | `bool 	isLocallyComplete () const`
:arrow_down: | `Export< LocalOrdinal, GlobalOrdinal, Node > & 	operator= (const Export< LocalOrdinal, GlobalOrdinal, Node > &rhs)`

**I/O Methods**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `virtual void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`
:x: | `virtual void 	print (std::ostream &os) const`

**Related Functions (non-member functions)**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Export < LocalOrdinal, GlobalOrdinal, Node > > 	createExport (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &src, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &tgt)`

### Tpetra::Import

**Constructor/Destructor Methods**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target)`
:x: | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::FancyOStream > &out)`
:white_check_mark: | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::ParameterList > &plist)`
:x: | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::FancyOStream > &out, const Teuchos::RCP< Teuchos::ParameterList > &plist)`
:arrow_down: | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, Teuchos::Array< int > &remotePIDs)`
:white_check_mark: | `Import (const Import< LocalOrdinal, GlobalOrdinal, Node > &importer)`
:white_check_mark: | `Import (const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter)`
:arrow_down: | `Import (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &source, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &target, Teuchos::Array< int > &userRemotePIDs, Teuchos::Array< GlobalOrdinal > &remoteGIDs, const Teuchos::ArrayView< const LocalOrdinal > &userExportLIDs, const Teuchos::ArrayView< const int > &userExportPIDs, const bool useRemotePIDs, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &out=Teuchos::null)`
:white_check_mark: | `virtual 	~Import ()`
:white_check_mark: | `void 	setParameterList (const Teuchos::RCP< Teuchos::ParameterList > &plist)`

**Import Attribute Methods**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `size_t 	getNumSameIDs () const`
:white_check_mark: | `size_t 	getNumPermuteIDs () const`
:arrow_down: | `Teuchos::ArrayView< const LocalOrdinal > 	getPermuteFromLIDs () const`
:arrow_down: | `Teuchos::ArrayView< const LocalOrdinal > 	getPermuteToLIDs () const`
:white_check_mark: | `size_t 	getNumRemoteIDs () const`
:arrow_down: | `Teuchos::ArrayView< const LocalOrdinal > 	getRemoteLIDs () const`
:white_check_mark: | `size_t 	getNumExportIDs () const`
:arrow_down: | `Teuchos::ArrayView< const LocalOrdinal > 	getExportLIDs () const`
:arrow_down: | `Teuchos::ArrayView< const int > 	getExportPIDs () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getSourceMap () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getTargetMap () const`
:x: | `Distributor & 	getDistributor () const`
:white_check_mark: | `bool 	isLocallyComplete () const`
:x: | `Import< LocalOrdinal, GlobalOrdinal, Node > & 	operator= (const Import< LocalOrdinal, GlobalOrdinal, Node > &Source)`
:white_check_mark: | `Teuchos::RCP< const Import < LocalOrdinal, GlobalOrdinal, Node > > 	setUnion (const Import< LocalOrdinal, GlobalOrdinal, Node > &rhs) const`
:white_check_mark: | `Teuchos::RCP< const Import < LocalOrdinal, GlobalOrdinal, Node > > 	setUnion () const`
:white_check_mark: | `Teuchos::RCP< const Import < LocalOrdinal, GlobalOrdinal, Node > > 	createRemoteOnlyImport (const Teuchos::RCP< const map_type > &remoteTarget) const`

**I/O Methods**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `virtual void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`
:x: | `virtual void 	print (std::ostream &os) const`

**Related Functions (non-member functions)**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Import < LocalOrdinal, GlobalOrdinal, Node > > 	createImport (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &src, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &tgt)`
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Import < LocalOrdinal, GlobalOrdinal, Node > > 	createImport (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &src, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &tgt, const Teuchos::RCP< Teuchos::ParameterList > &plist)`

### Tpetra::MultiVector

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`
:white_check_mark: | `void 	setCopyOrView (const Teuchos::DataAccess copyOrView)`
:white_check_mark: | `Teuchos::DataAccess 	getCopyOrView () const`
:white_check_mark: | `void 	assign (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &src)`

**Constructors and destructor**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `MultiVector ()`
:white_check_mark: | `MultiVector (const Teuchos::RCP< const map_type > &map, const size_t numVecs, const bool zeroOut=true)`
:white_check_mark: | `MultiVector (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &source)`
:white_check_mark: | `MultiVector (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &source, const Teuchos::DataAccess copyOrView)`
:star:| `MultiVector (const Teuchos::RCP< const map_type > &map, const Teuchos::ArrayView< const Scalar > &A, const size_t LDA, const size_t NumVectors)` | Fortran arrays
:arrow_down: | `MultiVector (const Teuchos::RCP< const map_type > &map, const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs, const size_t NumVectors)`
:x: | `MultiVector (const Teuchos::RCP< const map_type > &map, const dual_view_type &view)`
:x: | `MultiVector (const Teuchos::RCP< const map_type > &map, const typename dual_view_type::t_dev &d_view)`
:x: | `MultiVector (const Teuchos::RCP< const map_type > &map, const dual_view_type &view, const dual_view_type &origView)`
:x: | `MultiVector (const Teuchos::RCP< const map_type > &map, const dual_view_type &view, const Teuchos::ArrayView< const size_t > &whichVectors)`
:x: | `MultiVector (const Teuchos::RCP< const map_type > &map, const dual_view_type &view, const dual_view_type &origView, const Teuchos::ArrayView< const size_t > &whichVectors)`
:white_check_mark: | `MultiVector (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, const map_type &subMap, const size_t offset=0)`
:x: | `template<class Node2 > Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node2 > > 	clone (const Teuchos::RCP< Node2 > &node2) const`
:white_check_mark: | `virtual 	~MultiVector () Destructor (virtual for memory safety of derived classes)`

**The following methods get either a (deep) copy or a view (shallow copy) of a subset of rows and/or columns of the MultiVector**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subCopy (const Teuchos::Range1D &colRng) const` | prefer `ArrayView` variant
:star2: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subCopy (const Teuchos::ArrayView< const size_t > &cols) const` | Fortran arrays; 1-based indices
:x: | `Teuchos::RCP< const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subView (const Teuchos::Range1D &colRng) const` | prefer `ArrayView` variant
:star2: | `Teuchos::RCP< const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subView (const Teuchos::ArrayView< const size_t > &cols) const` | Fortran arrays; 1-based indices
:x: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subViewNonConst (const Teuchos::Range1D &colRng)` | prefer `ArrayView` variant
:star2: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subViewNonConst (const Teuchos::ArrayView< const size_t > &cols)` | Fortran arrays; 1-based indices
:white_check_mark: | `Teuchos::RCP< const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	offsetView (const Teuchos::RCP< const map_type > &subMap, const size_t offset) const`|
:white_check_mark: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	offsetViewNonConst (const Teuchos::RCP< const map_type > &subMap, const size_t offset)`
:x: | `Teuchos::RCP< const Vector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	getVector (const size_t j) const`
:x: | `Teuchos::RCP< Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	getVectorNonConst (const size_t j)`
:star: | `Teuchos::ArrayRCP< const Scalar > 	getData (size_t j) const` | Fortran arrays; 1-based indices
:star: | `Teuchos::ArrayRCP< Scalar > 	getDataNonConst (size_t j)` | Fortran arrays; 1-based indices
:star: | `void 	get1dCopy (const Teuchos::ArrayView< Scalar > &A, const size_t LDA) const` | Fortran arrays
:arrow_down: | `void 	get2dCopy (const Teuchos::ArrayView< const Teuchos::ArrayView< Scalar > > &ArrayOfPtrs) const`
:star: | `Teuchos::ArrayRCP< const Scalar > 	get1dView () const` | Fortran arrays
:arrow_down: | `Teuchos::ArrayRCP < Teuchos::ArrayRCP< const Scalar > > 	get2dView () const`
:star: | `Teuchos::ArrayRCP< Scalar > 	get1dViewNonConst ()` | Fortran arrays
:arrow_down: | `Teuchos::ArrayRCP < Teuchos::ArrayRCP< Scalar > > 	get2dViewNonConst ()`
:x: | `dual_view_type 	getDualView () const`
:x: | `template<class TargetDeviceType > void 	sync ()`
:x: | `template<class TargetDeviceType > bool 	need_sync () const`
:x: | `template<class TargetDeviceType > void 	modify ()`
:x: | `template<class TargetDeviceType > Kokkos::Impl::if_c < std::is_same< typename device_type::memory_space, typename TargetDeviceType::memory_space > ::value, typename dual_view_type::t_dev, typename dual_view_type::t_host >::type 	getLocalView () const`

**Mathematical methods**

:grey_question: | Command | Comment
-------|--------|---------
:star: | `void 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Teuchos::ArrayView< dot_type > &dots) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < dot_type, T >::value), void > ::type 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Teuchos::ArrayView< T > &dots) const` | Fortran arrays
:x: | `template<typename T > std::enable_if< !(std::is_same < dot_type, T >::value), void > ::type 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, std::vector< T > &dots) const`
:x: | `void 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Kokkos::View< dot_type *, device_type > &dots) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < dot_type, T >::value), void > ::type 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Kokkos::View< T *, device_type > &dots) const`
:white_check_mark: | `void 	abs (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A)`
:white_check_mark: | `void 	reciprocal (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A)`
:white_check_mark: | `void 	scale (const Scalar &alpha)`
:star: | `void 	scale (const Teuchos::ArrayView< const Scalar > &alpha)` | Fortran array
:x: | `void 	scale (const Kokkos::View< const impl_scalar_type *, device_type > &alpha)`
:white_check_mark: | `void 	scale (const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A)`
:white_check_mark: | `void 	update (const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Scalar &beta)`
:white_check_mark: | `void 	update (const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Scalar &beta, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, const Scalar &gamma)`
:x: | `void 	norm1 (const Kokkos::View< mag_type *, device_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	norm1 (const Kokkos::View< T *, device_type > &norms) const`
:star: | `void 	norm1 (const Teuchos::ArrayView< mag_type > &norms) const` | Fortran arrays
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	norm1 (const Teuchos::ArrayView< T > &norms) const`
:x: | `void 	norm2 (const Kokkos::View< mag_type *, device_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	norm2 (const Kokkos::View< T *, device_type > &norms) const`
:star: | `void 	norm2 (const Teuchos::ArrayView< mag_type > &norms) const` | Fortran arrays
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	norm2 (const Teuchos::ArrayView< T > &norms) const`
:x: | `void 	normInf (const Kokkos::View< mag_type *, device_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	normInf (const Kokkos::View< T *, device_type > &norms) const`
:star: | `void 	normInf (const Teuchos::ArrayView< mag_type > &norms) const` | Fortran arrays
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	normInf (const Teuchos::ArrayView< T > &norms) const`
:x: | `void TPETRA_DEPRECATED 	normWeighted (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &weights, const Teuchos::ArrayView< mag_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type TPETRA_DEPRECATED 	normWeighted (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &weights, const Teuchos::ArrayView< T > &norms) const`
:star: | `void 	meanValue (const Teuchos::ArrayView< impl_scalar_type > &means) const` | Fortran arrays
:x: | `template<typename T > std::enable_if<!std::is_same < impl_scalar_type, T >::value, void >::type 	meanValue (const Teuchos::ArrayView< T > &means) const`
:white_check_mark: | `void 	multiply (Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, const Scalar &beta)`
:x: | `void 	elementWiseMultiply (Scalar scalarAB, const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, Scalar scalarThis)`

**Attribute access functions**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `size_t 	getNumVectors () const`
:white_check_mark: | `size_t 	getLocalLength () const`
:white_check_mark: | `global_size_t 	getGlobalLength () const`
:white_check_mark: | `size_t 	getStride () const`
:white_check_mark: | `bool 	isConstantStride () const`

**Overridden from Teuchos::Describable**

:grey_question: | Command | Comment
-------|--------|---------
:arrow_down: | `virtual std::string 	description () const`
:x: | `virtual void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`

**Public methods for redistributing data**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `void 	doImport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`
:white_check_mark: | `void 	doImport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:white_check_mark: | `void 	doExport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:white_check_mark: | `void 	doExport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`

**Attribute accessor methods**

:grey_question: | Command | Comment
-------|--------|---------
:arrow_down: | `bool 	isDistributed () const`
:white_check_mark: | `virtual Teuchos::RCP< const map_type > 	getMap () const`

**I/O methods**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `void 	print (std::ostream &os) const`

**Methods for use only by experts**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`

**Friends**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `template<class DS , class DL , class DG , class DN , const bool dstClassic, class SS , class SL , class SG , class SN , const bool srcClassic> void 	deep_copy (MultiVector< DS, DL, DG, DN, dstClassic > &dst, const MultiVector< SS, SL, SG, SN, srcClassic > &src)`

**Related Functions (non-member functions)**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `template<class DS , class DL , class DG , class DN , const bool dstClassic, class SS , class SL , class SG , class SN , const bool srcClassic> void 	deep_copy (MultiVector< DS, DL, DG, DN, dstClassic > &dst, const MultiVector< SS, SL, SG, SN, srcClassic > &src)`
:x: | `template<class ST , class LO , class GO , class NT , const bool classic = NT::classic> MultiVector< ST, LO, GO, NT, classic > 	createCopy (const MultiVector< ST, LO, GO, NT, classic > &src)`

**Post-construction modification routines**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `static const bool 	useAtomicUpdatesByDefault`
:star: | `void 	replaceGlobalValue (const GlobalOrdinal gblRow, const size_t col, const impl_scalar_type &value) const` | 1-based indices
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	replaceGlobalValue (GlobalOrdinal globalRow, size_t col, const T &value) const`
:star:| `void 	sumIntoGlobalValue (const GlobalOrdinal gblRow, const size_t col, const impl_scalar_type &value, const bool atomic=useAtomicUpdatesByDefault) const` | 1-based indices
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	sumIntoGlobalValue (const GlobalOrdinal gblRow, const size_t col, const T &val, const bool atomic=useAtomicUpdatesByDefault) const`
:star:| `void 	replaceLocalValue (const LocalOrdinal lclRow, const size_t col, const impl_scalar_type &value) const` | 1-based indices
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	replaceLocalValue (const LocalOrdinal lclRow, const size_t col, const T &val) const`
:star: | `void 	sumIntoLocalValue (const LocalOrdinal lclRow, const size_t col, const impl_scalar_type &val, const bool atomic=useAtomicUpdatesByDefault) const` | 1-based indices
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	sumIntoLocalValue (const LocalOrdinal lclRow, const size_t col, const T &val, const bool atomic=useAtomicUpdatesByDefault) const`
:white_check_mark: | `void 	putScalar (const Scalar &value)`
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	putScalar (const T &value)`
:white_check_mark: | `void 	randomize ()`
:white_check_mark: | `void 	randomize (const Scalar &minVal, const Scalar &maxVal)`
:white_check_mark: | `void 	replaceMap (const Teuchos::RCP< const map_type > &map)`
:white_check_mark: | `void 	reduce ()`
:x: | `MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > & 	operator= (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &source)`


### Tpetra::CrsGraph

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `bool 	haveGlobalConstants () const`
:white_check_mark: | `void 	computeGlobalConstants ()`
:x: | `local_graph_type 	getLocalGraph () const`

**Constructor/Destructor Methods**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:x: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow, const ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:star: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::ArrayRCP< const size_t > &numEntPerRow, const ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Fortran arrays
:white_check_mark: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const size_t maxNumEntriesPerRow, const ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:x: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:star: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const Teuchos::ArrayRCP< const size_t > &numEntPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Fortran arrays
:x: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const typename local_graph_type::row_map_type &rowPointers, const typename local_graph_type::entries_type::non_const_type &columnIndices, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:star2: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const Teuchos::ArrayRCP< size_t > &rowPointers, const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Fortran arrays; 1-based indices
:x: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const local_graph_type &lclGraph, const Teuchos::RCP< Teuchos::ParameterList > &params)`
:x: | `template<class Node2 > Teuchos::RCP< CrsGraph < LocalOrdinal, GlobalOrdinal, Node2, Node2::classic > > 	clone (const Teuchos::RCP< Node2 > &node2, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) const`
:white_check_mark: | `virtual 	~CrsGraph ()`
:white_check_mark: | `void 	setParameterList (const Teuchos::RCP< Teuchos::ParameterList > &params)`
:white_check_mark: | `Teuchos::RCP< const Teuchos::ParameterList > 	getValidParameters () const`

**Insertion/Removal Methods**

:grey_question: | Command | Comment
-------|--------|---------
:star: | `void 	insertGlobalIndices (const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &indices)` | Fortran arrays
:x: | `void 	insertGlobalIndices (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const GlobalOrdinal inds[])` | prefer `ArrayView` variant
:star2: | `void 	insertLocalIndices (const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &indices)` | Fortran arrays; 1-based indices
:x: | `void 	insertLocalIndices (const LocalOrdinal localRow, const LocalOrdinal numEnt, const LocalOrdinal inds[])` | prefer `ArrayView` variant
:white_check_mark: | `void 	removeLocalIndices (LocalOrdinal localRow)`

**Transformational Methods**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `void 	globalAssemble ()`
:white_check_mark: | `void 	resumeFill (const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	fillComplete (const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	fillComplete (const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	expertStaticFillComplete (const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< const import_type > &importer=Teuchos::null, const Teuchos::RCP< const export_type > &exporter=Teuchos::null, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`

**Methods implementing RowGraph.**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `Teuchos::RCP< const Teuchos::Comm< int > > 	getComm () const`
:x: | `Teuchos::RCP< node_type > 	getNode () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getRowMap () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getColMap () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getDomainMap () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getRangeMap () const`
:white_check_mark: | `Teuchos::RCP< const import_type > 	getImporter () const`
:white_check_mark: | `Teuchos::RCP< const export_type > 	getExporter () const`
:white_check_mark: | `global_size_t 	getGlobalNumRows () const`
:white_check_mark: | `global_size_t 	getGlobalNumCols () const`
:white_check_mark: | `size_t 	getNodeNumRows () const`
:white_check_mark: | `size_t 	getNodeNumCols () const`
:x: | `GlobalOrdinal 	getIndexBase () const`
:white_check_mark: | `global_size_t 	getGlobalNumEntries () const`
:white_check_mark: | `size_t 	getNodeNumEntries () const`
:white_check_mark: | `size_t 	getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const`
:star: | `size_t 	getNumEntriesInLocalRow (LocalOrdinal localRow) const` | 1-based indices
:white_check_mark: | `size_t 	getNodeAllocationSize () const`
:white_check_mark: | `size_t 	getNumAllocatedEntriesInGlobalRow (GlobalOrdinal globalRow) const`
:star: | `size_t 	getNumAllocatedEntriesInLocalRow (LocalOrdinal localRow) const`     | 1-based indices
:white_check_mark: | `global_size_t 	getGlobalNumDiags () const`
:white_check_mark: | `size_t 	getNodeNumDiags () const`
:white_check_mark: | `size_t 	getGlobalMaxNumRowEntries () const`
:white_check_mark: | `size_t 	getNodeMaxNumRowEntries () const`
:white_check_mark: | `bool 	hasColMap () const`
:white_check_mark: | `bool 	isLowerTriangular () const`
:white_check_mark: | `bool 	isUpperTriangular () const`
:white_check_mark: | `bool 	isLocallyIndexed () const`
:white_check_mark: | `bool 	isGloballyIndexed () const`
:white_check_mark: | `bool 	isFillComplete () const`
:white_check_mark: | `bool 	isFillActive () const`
:white_check_mark: | `bool 	isSorted () const`
:white_check_mark: | `bool 	isStorageOptimized () const`
:white_check_mark: | `ProfileType 	getProfileType () const`
:star: | `void 	getGlobalRowCopy (GlobalOrdinal GlobalRow, const Teuchos::ArrayView< GlobalOrdinal > &Indices, size_t &NumIndices) const` | Fortran arrays
:star2: | `void 	getLocalRowCopy (LocalOrdinal LocalRow, const Teuchos::ArrayView< LocalOrdinal > &indices, size_t &NumIndices) const` | Fortran arrays; 1-based indices
:star: | `void 	getGlobalRowView (const GlobalOrdinal gblRow, Teuchos::ArrayView< const GlobalOrdinal > &gblColInds) const` | Fortran arrays
:x: | `void 	getLocalRowView (const LocalOrdinal lclRow, Teuchos::ArrayView< const LocalOrdinal > &lclColInds) const` | cannot maintain the semantics
:white_check_mark: | `bool 	supportsRowViews () const`

**Overridden from Teuchos::Describable**

:grey_question: | Command | Comment
-------|--------|---------
:arrow_down: | `std::string 	description () const`
:x: | `void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`

**Implementation of DistObject**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `virtual bool 	checkSizes (const SrcDistObject &source)`
:x: | `virtual void 	copyAndPermute (const SrcDistObject &source, size_t numSameIDs, const Teuchos::ArrayView< const LocalOrdinal > &permuteToLIDs, const Teuchos::ArrayView< const LocalOrdinal > &permuteFromLIDs)`
:x: | `virtual void 	packAndPrepare (const SrcDistObject &source, const Teuchos::ArrayView< const LocalOrdinal > &exportLIDs, Teuchos::Array< GlobalOrdinal > &exports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t &constantNumPackets, Distributor &distor)`
:x: | `virtual void 	pack (const Teuchos::ArrayView< const LocalOrdinal > &exportLIDs, Teuchos::Array< GlobalOrdinal > &exports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t &constantNumPackets, Distributor &distor) const`
:x: | `virtual void 	unpackAndCombine (const Teuchos::ArrayView< const LocalOrdinal > &importLIDs, const Teuchos::ArrayView< const GlobalOrdinal > &imports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t constantNumPackets, Distributor &distor, CombineMode CM)`

**Advanced methods, at increased risk of deprecation.**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `void 	getLocalDiagOffsets (const Kokkos::View< size_t *, device_type, Kokkos::MemoryUnmanaged > &offsets) const`
:star: | `void 	getLocalDiagOffsets (Teuchos::ArrayRCP< size_t > &offsets) const` | Fortran arrays
:arrow_down: | `void 	getNumEntriesPerLocalRowUpperBound (Teuchos::ArrayRCP< const size_t > &boundPerLocalRow, size_t &boundForAllLocalRows, bool &boundSameForAllLocalRows) const`
:x: | `void 	setAllIndices (const typename local_graph_type::row_map_type &rowPointers, const typename local_graph_type::entries_type::non_const_type &columnIndices)`
:star2: | `void 	setAllIndices (const Teuchos::ArrayRCP< size_t > &rowPointers, const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices)` | Fortran arrays; 1-based indices
:star2: | `Teuchos::ArrayRCP< const size_t > 	getNodeRowPtrs () const` | Fortran arrays, different semantics
:star2: | `Teuchos::ArrayRCP< const LocalOrdinal > 	getNodePackedIndices () const` | Fortran arrays, different semantics
:white_check_mark: | `void 	replaceColMap (const Teuchos::RCP< const map_type > &newColMap)`
:white_check_mark: | `void 	reindexColumns (const Teuchos::RCP< const map_type > &newColMap, const Teuchos::RCP< const import_type > &newImport=Teuchos::null, const bool sortIndicesInEachRow=true)`
:white_check_mark: | `void 	replaceDomainMapAndImporter (const Teuchos::RCP< const map_type > &newDomainMap, const Teuchos::RCP< const import_type > &newImporter)`
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`
:white_check_mark: | `void 	doImport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`
:white_check_mark: | `void 	doImport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:white_check_mark: | `void 	doExport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:white_check_mark: | `void 	doExport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`

**Attribute accessor methods**

:grey_question: | Command | Comment
-------|--------|---------
:arrow_down: | `bool 	isDistributed () const`
:x: | `virtual Teuchos::RCP< const map_type > 	getMap () const`

**I/O methods**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `void 	print (std::ostream &os) const`

**Methods for use only by experts**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`

### Tpetra::CrsMatrix

**Public Member Functions**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `local_matrix_type::values_type 	getLocalValuesView () const`
:white_check_mark: | `void 	importAndFillComplete (Teuchos::RCP< CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > &destMatrix, const import_type &importer, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) const`
:white_check_mark: | `void 	importAndFillComplete (Teuchos::RCP< CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > &destMatrix, const import_type &rowImporter, const import_type &domainImporter, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params) const`
:white_check_mark: | `void 	exportAndFillComplete (Teuchos::RCP< CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > &destMatrix, const export_type &exporter, const Teuchos::RCP< const map_type > &domainMap=Teuchos::null, const Teuchos::RCP< const map_type > &rangeMap=Teuchos::null, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) const`
:white_check_mark: | `void 	exportAndFillComplete (Teuchos::RCP< CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > &destMatrix, const export_type &rowExporter, const export_type &domainExporter, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params) const`
:white_check_mark: | `bool 	haveGlobalConstants () const`

**Transformational methods**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `void 	globalAssemble ()`
:white_check_mark: | `void 	resumeFill (const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	fillComplete (const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	fillComplete (const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	expertStaticFillComplete (const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< const import_type > &importer=Teuchos::null, const Teuchos::RCP< const export_type > &exporter=Teuchos::null, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	replaceColMap (const Teuchos::RCP< const map_type > &newColMap)`
:white_check_mark: | `void 	reindexColumns (crs_graph_type *const graph, const Teuchos::RCP< const map_type > &newColMap, const Teuchos::RCP< const import_type > &newImport=Teuchos::null, const bool sortEachRow=true)`
:white_check_mark: | `void 	replaceDomainMapAndImporter (const Teuchos::RCP< const map_type > &newDomainMap, Teuchos::RCP< const import_type > &newImporter)`
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`

**Methods implementing Operator**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `void 	apply (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &Y, Teuchos::ETransp mode=Teuchos::NO_TRANS, Scalar alpha=Teuchos::ScalarTraits< Scalar >::one(), Scalar beta=Teuchos::ScalarTraits< Scalar >::zero()) const`
:white_check_mark: | `bool 	hasTransposeApply () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getDomainMap () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getRangeMap () const`

**Other "apply"-like methods**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `void 	gaussSeidel (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &D, const Scalar &dampingFactor, const ESweepDirection direction, const int numSweeps) const`
:arrow_down: | `void 	reorderedGaussSeidel (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &D, const Teuchos::ArrayView< LocalOrdinal > &rowIndices, const Scalar &dampingFactor, const ESweepDirection direction, const int numSweeps) const`
:white_check_mark: | `void 	gaussSeidelCopy (MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &D, const Scalar &dampingFactor, const ESweepDirection direction, const int numSweeps, const bool zeroInitialGuess) const`
:arrow_down: | `void 	reorderedGaussSeidelCopy (MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &D, const Teuchos::ArrayView< LocalOrdinal > &rowIndices, const Scalar &dampingFactor, const ESweepDirection direction, const int numSweeps, const bool zeroInitialGuess) const`
:arrow_down: | `virtual Teuchos::RCP < RowMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > 	add (const Scalar &alpha, const RowMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const Scalar &beta, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params) const`

**Implementation of Teuchos::Describable interface**

:grey_question: | Command | Comment
-------|--------|---------
:arrow_down: | `std::string 	description () const`
:x: | `void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`

**Extraction Methods**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `virtual void 	getLocalDiagCopy (Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, Node::classic > &diag) const =0`

**Mathematical methods**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `virtual void 	leftScale (const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, Node::classic > &x)=0`
:x: | `virtual void 	rightScale (const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, Node::classic > &x)=0`

**Pure virtual functions to be overridden by subclasses.**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `virtual void 	apply (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y, Teuchos::ETransp mode=Teuchos::NO_TRANS, Scalar alpha=Teuchos::ScalarTraits< Scalar >::one(), Scalar beta=Teuchos::ScalarTraits< Scalar >::zero()) const =0`

**Public methods for redistributing data**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `void 	doImport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`
:white_check_mark: | `void 	doImport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:white_check_mark: | `void 	doExport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:white_check_mark: | `void 	doExport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`

**Attribute accessor methods**

:grey_question: | Command | Comment
-------|--------|---------
:arrow_down: | `bool 	isDistributed () const`
:x: | `virtual Teuchos::RCP< const map_type > 	getMap () const`

**I/O methods**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `void 	print (std::ostream &os) const`

**Methods for use only by experts**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`

**Methods implemented by subclasses and used by doTransfer().**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `virtual bool 	useNewInterface ()`
:x: | `virtual void 	copyAndPermuteNew (const SrcDistObject &source, const size_t numSameIDs, const Kokkos::DualView< const local_ordinal_type *, device_type > &permuteToLIDs, const Kokkos::DualView< const local_ordinal_type *, device_type > &permuteFromLIDs)`
:x: | `virtual void 	packAndPrepare (const SrcDistObject &source, const Teuchos::ArrayView< const local_ordinal_type > &exportLIDs, Teuchos::Array< packet_type > &exports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t &constantNumPackets, Distributor &distor)`
:x: | `virtual void 	packAndPrepareNew (const SrcDistObject &source, const Kokkos::DualView< const local_ordinal_type *, device_type > &exportLIDs, Kokkos::DualView< packet_type *, device_type > &exports, const Kokkos::DualView< size_t *, device_type > &numPacketsPerLID, size_t &constantNumPackets, Distributor &distor)`
:x: | `virtual void 	unpackAndCombine (const Teuchos::ArrayView< const local_ordinal_type > &importLIDs, const Teuchos::ArrayView< const packet_type > &imports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t constantNumPackets, Distributor &distor, CombineMode CM)`
:x: | `virtual void 	unpackAndCombineNew (const Kokkos::DualView< const local_ordinal_type *, device_type > &importLIDs, const Kokkos::DualView< const packet_type *, device_type > &imports, const Kokkos::DualView< const size_t *, device_type > &numPacketsPerLID, const size_t constantNumPackets, Distributor &distor, const CombineMode CM)`

**Related Functions (non-member functions)**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `template<class Scalar , class LocalOrdinal , class GlobalOrdinal , class Node , const bool classic = Node::classic> Teuchos::RCP< CrsMatrix < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	createCrsMatrix (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map, size_t maxNumEntriesPerRow=0, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:x: | `template<class CrsMatrixType > Teuchos::RCP< CrsMatrixType > 	importAndFillCompleteCrsMatrix (const Teuchos::RCP< const CrsMatrixType > &sourceMatrix, const Import< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > &importer, const Teuchos::RCP< const Map< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > > &domainMap=Teuchos::null, const Teuchos::RCP< const Map< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > > &rangeMap=Teuchos::null, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:x: | `template<class CrsMatrixType > Teuchos::RCP< CrsMatrixType > 	importAndFillCompleteCrsMatrix (const Teuchos::RCP< const CrsMatrixType > &sourceMatrix, const Import< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > &rowImporter, const Import< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > &domainImporter, const Teuchos::RCP< const Map< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > > &domainMap, const Teuchos::RCP< const Map< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params)`
:x: | `template<class CrsMatrixType > Teuchos::RCP< CrsMatrixType > 	exportAndFillCompleteCrsMatrix (const Teuchos::RCP< const CrsMatrixType > &sourceMatrix, const Export< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > &exporter, const Teuchos::RCP< const Map< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > > &domainMap=Teuchos::null, const Teuchos::RCP< const Map< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > > &rangeMap=Teuchos::null, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:x: | `template<class CrsMatrixType > Teuchos::RCP< CrsMatrixType > 	exportAndFillCompleteCrsMatrix (const Teuchos::RCP< const CrsMatrixType > &sourceMatrix, const Export< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > &rowExporter, const Export< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > &domainExporter, const Teuchos::RCP< const Map< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > > &domainMap, const Teuchos::RCP< const Map< typename CrsMatrixType::local_ordinal_type, typename CrsMatrixType::global_ordinal_type, typename CrsMatrixType::node_type > > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params)`

**Constructors and destructor**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `CrsMatrix (const Teuchos::RCP< const map_type > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:star: | `CrsMatrix (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Fortran arrays
:star: | `CrsMatrix (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Fortran arrays
:star: | `CrsMatrix (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const Teuchos::ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Fortran arrays
:white_check_mark: | `CrsMatrix (const Teuchos::RCP< const crs_graph_type > &graph, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:x: | `CrsMatrix (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const typename local_matrix_type::row_map_type &rowPointers, const typename local_graph_type::entries_type::non_const_type &columnIndices, const typename local_matrix_type::values_type &values, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:star: | `CrsMatrix (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const Teuchos::ArrayRCP< size_t > &rowPointers, const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices, const Teuchos::ArrayRCP< Scalar > &values, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Fortran arrays
:x: | `CrsMatrix (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const local_matrix_type &lclMatrix, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:x: | `template<class Node2 > Teuchos::RCP< CrsMatrix < Scalar, LocalOrdinal, GlobalOrdinal, Node2, Node2::classic > > 	clone (const Teuchos::RCP< Node2 > &node2, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) const`
:white_check_mark: | `virtual 	~CrsMatrix ()`

**Methods for inserting, modifying, or removing entries**

:grey_question: | Command | Comment
-------|--------|---------
:star: | `void 	insertGlobalValues (const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals)` | Fortran arrays
:x: | `void 	insertGlobalValues (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const Scalar vals[], const GlobalOrdinal inds[])` | prefer `ArrayView` variant
:star2: | `void 	insertLocalValues (const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals)` | Fortran arrays; 1-based indices
:x: | `void 	insertLocalValues (const LocalOrdinal localRow, const LocalOrdinal numEnt, const Scalar vals[], const LocalOrdinal cols[])` | prefer `ArrayView` variant
:x: | `template<class GlobalIndicesViewType , class ImplScalarViewType > LocalOrdinal 	replaceGlobalValues (const GlobalOrdinal globalRow, const typename UnmanagedView< GlobalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals) const`
:star: | `LocalOrdinal 	replaceGlobalValues (const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals) const` | Fortran arrays
:x: | `LocalOrdinal 	replaceGlobalValues (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const Scalar vals[], const GlobalOrdinal cols[]) const` | prefer `ArrayView` variant
:x: | `template<class LocalIndicesViewType , class ImplScalarViewType > LocalOrdinal 	replaceLocalValues (const LocalOrdinal localRow, const typename UnmanagedView< LocalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals) const`
:star2: | `LocalOrdinal 	replaceLocalValues (const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals) const` | Fortran arrays; 1-based indices
:x: | `LocalOrdinal 	replaceLocalValues (const LocalOrdinal localRow, const LocalOrdinal numEnt, const Scalar inputVals[], const LocalOrdinal inputCols[]) const` | prefer `ArrayView` variant
:star: | `LocalOrdinal 	sumIntoGlobalValues (const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals, const bool atomic=useAtomicUpdatesByDefault)` | Fortran arrays; no `atomic`
:x: | `LocalOrdinal 	sumIntoGlobalValues (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const Scalar vals[], const GlobalOrdinal cols[], const bool atomic=useAtomicUpdatesByDefault)` | prefer `ArrayView` variant
:x: | `template<class LocalIndicesViewType , class ImplScalarViewType > LocalOrdinal 	sumIntoLocalValues (const LocalOrdinal localRow, const typename UnmanagedView< LocalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals, const bool atomic=useAtomicUpdatesByDefault) const`
:star2: | `LocalOrdinal 	sumIntoLocalValues (const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals, const bool atomic=useAtomicUpdatesByDefault) const` | Fortran arrays; 1-based indices; no `atomic`
:x: | `LocalOrdinal 	sumIntoLocalValues (const LocalOrdinal localRow, const LocalOrdinal numEnt, const Scalar vals[], const LocalOrdinal cols[], const bool atomic=useAtomicUpdatesByDefault) const` | prefer `ArrayView` variant
:x: | `template<class LocalIndicesViewType , class ImplScalarViewType , class BinaryFunction > LocalOrdinal 	transformLocalValues (const LocalOrdinal localRow, const typename UnmanagedView< LocalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals, BinaryFunction f, const bool atomic=useAtomicUpdatesByDefault) const`
:x: | `template<class BinaryFunction , class InputMemorySpace > LocalOrdinal 	transformGlobalValues (const GlobalOrdinal globalRow, const Kokkos::View< const GlobalOrdinal *, InputMemorySpace, Kokkos::MemoryUnmanaged > &inputInds, const Kokkos::View< const impl_scalar_type *, InputMemorySpace, Kokkos::MemoryUnmanaged > &inputVals, BinaryFunction f, const bool atomic=useAtomicUpdatesByDefault) const`
:white_check_mark: | `void 	setAllToScalar (const Scalar &alpha)`
:white_check_mark: | `void 	scale (const Scalar &alpha)`
:x: | `void 	setAllValues (const typename local_matrix_type::row_map_type &ptr, const typename local_graph_type::entries_type::non_const_type &ind, const typename local_matrix_type::values_type &val)`
:star2: | `void 	setAllValues (const Teuchos::ArrayRCP< size_t > &ptr, const Teuchos::ArrayRCP< LocalOrdinal > &ind, const Teuchos::ArrayRCP< Scalar > &val)` | Fortran arrays; 1-based indices
:star2: | `void 	getAllValues (Teuchos::ArrayRCP< const size_t > &rowPointers, Teuchos::ArrayRCP< const LocalOrdinal > &columnIndices, Teuchos::ArrayRCP< const Scalar > &values) const` | Fortran arrays, different semantics

**Methods implementing RowMatrix**

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `Teuchos::RCP< const Teuchos::Comm< int > > 	getComm () const`
:x: | `Teuchos::RCP< node_type > 	getNode () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getRowMap () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getColMap () const`
:x: | `Teuchos::RCP< const RowGraph < LocalOrdinal, GlobalOrdinal, Node > > 	getGraph () const`
:white_check_mark: | `Teuchos::RCP< const crs_graph_type > 	getCrsGraph () const`
:x: | `local_matrix_type 	getLocalMatrix () const`
:white_check_mark: | `global_size_t 	getGlobalNumRows () const`
:white_check_mark: | `global_size_t 	getGlobalNumCols () const`
:white_check_mark: | `size_t 	getNodeNumRows () const`
:white_check_mark: | `size_t 	getNodeNumCols () const`
:x: | `GlobalOrdinal 	getIndexBase () const`
:white_check_mark: | `global_size_t 	getGlobalNumEntries () const`
:white_check_mark: | `size_t 	getNodeNumEntries () const`
:white_check_mark: | `size_t 	getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const`
:star: | `size_t 	getNumEntriesInLocalRow (LocalOrdinal localRow) const` | 1-based indices
:white_check_mark: | `global_size_t 	getGlobalNumDiags () const`
:white_check_mark: | `size_t 	getNodeNumDiags () const`
:white_check_mark: | `size_t 	getGlobalMaxNumRowEntries () const`
:white_check_mark: | `size_t 	getNodeMaxNumRowEntries () const`
:white_check_mark: | `bool 	hasColMap () const`
:white_check_mark: | `bool 	isLowerTriangular () const`
:white_check_mark: | `bool 	isUpperTriangular () const`
:white_check_mark: | `bool 	isLocallyIndexed () const`
:white_check_mark: | `bool 	isGloballyIndexed () const`
:white_check_mark: | `bool 	isFillComplete () const`
:white_check_mark: | `bool 	isFillActive () const`
:white_check_mark: | `bool 	isStorageOptimized () const`
:white_check_mark: | `ProfileType 	getProfileType () const`
:white_check_mark: | `bool 	isStaticGraph () const`
:white_check_mark: | `mag_type 	getFrobeniusNorm () const`
:white_check_mark: | `virtual bool 	supportsRowViews () const`
:star: | `void 	getGlobalRowCopy (GlobalOrdinal GlobalRow, const Teuchos::ArrayView< GlobalOrdinal > &Indices, const Teuchos::ArrayView< Scalar > &Values, size_t &NumEntries) const` | Fortran arrays
:star2: | `void 	getLocalRowCopy (LocalOrdinal localRow, const Teuchos::ArrayView< LocalOrdinal > &colInds, const Teuchos::ArrayView< Scalar > &vals, size_t &numEntries) const` | Fortran arrays; 1-based indices
:star: | `void 	getGlobalRowView (GlobalOrdinal GlobalRow, Teuchos::ArrayView< const GlobalOrdinal > &indices, Teuchos::ArrayView< const Scalar > &values) const` | Fortran arrays
:x: | `void 	getLocalRowView (LocalOrdinal LocalRow, Teuchos::ArrayView< const LocalOrdinal > &indices, Teuchos::ArrayView< const Scalar > &values) const` | use `getLocalRowCopy` as we cannot maintain the semantics
:x: | `LocalOrdinal 	getLocalRowViewRaw (const LocalOrdinal lclRow, LocalOrdinal &numEnt, const LocalOrdinal *&lclColInds, const Scalar *&vals) const` | prefer `ArrayView` variant
:x: | `LocalOrdinal 	getLocalRowView (const LocalOrdinal lclRow, LocalOrdinal &numEnt, const impl_scalar_type *&val, const LocalOrdinal *&ind) const` | prefer `ArrayView` variant
:x: | `template<class OutputScalarType > std::enable_if<!std::is_same < OutputScalarType, impl_scalar_type >::value &&std::is_convertible < impl_scalar_type, OutputScalarType >::value, LocalOrdinal >::type 	getLocalRowView (const LocalOrdinal lclRow, LocalOrdinal &numEnt, const OutputScalarType *&val, const LocalOrdinal *&ind) const` | prefer `ArrayView` variant
:arrow_down: | `void 	getLocalDiagCopy (Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &diag) const`
:star: | `void 	getLocalDiagOffsets (Teuchos::ArrayRCP< size_t > &offsets) const` | Fortran arrays
:x: | `void 	getLocalDiagCopy (Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &diag, const Kokkos::View< const size_t *, device_type, Kokkos::MemoryUnmanaged > &offsets) const`
:arrow_down: | `void 	getLocalDiagCopy (Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &diag, const Teuchos::ArrayView< const size_t > &offsets) const`
:arrow_down: | `void 	leftScale (const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &x)`
:arrow_down: | `void 	rightScale (const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &x)`

**Advanced templated methods**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `template<class DomainScalar , class RangeScalar > void 	localMultiply (const MultiVector< DomainScalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, MultiVector< RangeScalar, LocalOrdinal, GlobalOrdinal, Node, classic > &Y, Teuchos::ETransp mode, RangeScalar alpha, RangeScalar beta) const`
:x: | `template<class DomainScalar , class RangeScalar > void 	localGaussSeidel (const MultiVector< DomainScalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, MultiVector< RangeScalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &D, const RangeScalar &dampingFactor, const KokkosClassic::ESweepDirection direction) const`
:x: | `template<class DomainScalar , class RangeScalar > void 	reorderedLocalGaussSeidel (const MultiVector< DomainScalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, MultiVector< RangeScalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &D, const Teuchos::ArrayView< LocalOrdinal > &rowIndices, const RangeScalar &dampingFactor, const KokkosClassic::ESweepDirection direction) const`
:x: | `template<class DomainScalar , class RangeScalar > void 	localSolve (const MultiVector< RangeScalar, LocalOrdinal, GlobalOrdinal, Node, classic > &Y, MultiVector< DomainScalar, LocalOrdinal, GlobalOrdinal, Node, classic > &X, Teuchos::ETransp mode) const`
:x: | `template<class T > Teuchos::RCP< CrsMatrix< T, LocalOrdinal, GlobalOrdinal, Node, classic > > 	convert () const`

**Implementation of DistObject interface**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `virtual bool 	checkSizes (const SrcDistObject &source)`
:x: | `virtual void 	copyAndPermute (const SrcDistObject &source, size_t numSameIDs, const Teuchos::ArrayView< const LocalOrdinal > &permuteToLIDs, const Teuchos::ArrayView< const LocalOrdinal > &permuteFromLIDs)`
:x: | `virtual void 	packAndPrepare (const SrcDistObject &source, const Teuchos::ArrayView< const LocalOrdinal > &exportLIDs, Teuchos::Array< char > &exports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t &constantNumPackets, Distributor &distor)`
:x: | `void 	unpackAndCombine (const Teuchos::ArrayView< const LocalOrdinal > &importLIDs, const Teuchos::ArrayView< const char > &imports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t constantNumPackets, Distributor &distor, CombineMode combineMode)`

**Implementation of Packable interface**

:grey_question: | Command | Comment
-------|--------|---------
:x: | `virtual void 	pack (const Teuchos::ArrayView< const LocalOrdinal > &exportLIDs, Teuchos::Array< char > &exports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t &constantNumPackets, Distributor &distor) const`
:x: | `void 	packNonStatic (const Teuchos::ArrayView< const LocalOrdinal > &exportLIDs, Teuchos::Array< char > &exports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t &constantNumPackets, Distributor &distor) const`

### Tpetra::MatrixMatrix

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `template<class Scalar , class LocalOrdinal , class GlobalOrdinal , class Node > void 	Multiply (const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, bool transposeA, const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, bool transposeB, CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &C, bool call_FillComplete_on_result=true, const std::string &label=std::string(), const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `template<class Scalar , class LocalOrdinal , class GlobalOrdinal , class Node > void 	Add (const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, bool transposeA, Scalar scalarA, CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, Scalar scalarB)`
:x: | `template<class Scalar , class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< CrsMatrix < Scalar, LocalOrdinal, GlobalOrdinal, Node > > 	add (const Scalar &alpha, const bool transposeA, const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const Scalar &beta, const bool transposeB, const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap=Teuchos::null, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap=Teuchos::null, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `template<class Scalar , class LocalOrdinal , class GlobalOrdinal , class Node > void 	Add (const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, bool transposeA, Scalar scalarA, const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, bool transposeB, Scalar scalarB, Teuchos::RCP< CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > > C)`
:arrow_down: | `template<class Scalar , class LocalOrdinal , class GlobalOrdinal , class Node > void 	Jacobi (Scalar omega, const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Dinv, const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > &C, bool call_FillComplete_on_result=true, const std::string &label=std::string(), const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`

### Tpetra::Reader

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraphFile (const std::string &filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraphFile (const std::string &filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< node_type > &pNode, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:white_check_mark: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraphFile (const std::string &filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraphFile (const std::string &filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< node_type > &pNode, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false)`
:white_check_mark: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraphFile (const std::string &filename, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraph (std::istream &in, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraph (std::istream &in, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< node_type > &pNode, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraph (std::istream &in, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraph (std::istream &in, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< node_type > &pNode, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_graph_type > 	readSparseGraph (std::istream &in, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:white_check_mark: | `static Teuchos::RCP < sparse_matrix_type > 	readSparseFile (const std::string &filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_matrix_type > 	readSparseFile (const std::string &filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< node_type > &pNode, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:white_check_mark: | `static Teuchos::RCP < sparse_matrix_type > 	readSparseFile (const std::string &filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_matrix_type > 	readSparseFile (const std::string &filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< node_type > &pNode, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false)`
:white_check_mark: | `static Teuchos::RCP < sparse_matrix_type > 	readSparseFile (const std::string &filename, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_matrix_type > 	readSparse (std::istream &in, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_matrix_type > 	readSparse (std::istream &in, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< node_type > &pNode, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_matrix_type > 	readSparse (std::istream &in, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_matrix_type > 	readSparse (std::istream &in, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< node_type > &pNode, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < sparse_matrix_type > 	readSparse (std::istream &in, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false)`
:white_check_mark: | `static Teuchos::RCP < multivector_type > 	readDenseFile (const std::string &filename, const Teuchos::RCP< const comm_type > &comm, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < multivector_type > 	readDenseFile (const std::string &filename, const Teuchos::RCP< const comm_type > &comm, const Teuchos::RCP< node_type > &node, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false)`
:white_check_mark: | `static Teuchos::RCP< vector_type > 	readVectorFile (const std::string &filename, const Teuchos::RCP< const comm_type > &comm, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP< vector_type > 	readVectorFile (const std::string &filename, const Teuchos::RCP< const comm_type > &comm, const Teuchos::RCP< node_type > &node, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < multivector_type > 	readDense (std::istream &in, const Teuchos::RCP< const comm_type > &comm, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP < multivector_type > 	readDense (std::istream &in, const Teuchos::RCP< const comm_type > &comm, const Teuchos::RCP< node_type > &node, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP< vector_type > 	readVector (std::istream &in, const Teuchos::RCP< const comm_type > &comm, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP< vector_type > 	readVector (std::istream &in, const Teuchos::RCP< const comm_type > &comm, const Teuchos::RCP< node_type > &node, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false)`
:white_check_mark: | `static Teuchos::RCP< const map_type > 	readMapFile (const std::string &filename, const Teuchos::RCP< const comm_type > &comm, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP< const map_type > 	readMapFile (const std::string &filename, const Teuchos::RCP< const comm_type > &comm, const Teuchos::RCP< node_type > &node, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP< const map_type > 	readMap (std::istream &in, const Teuchos::RCP< const comm_type > &comm, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP< const map_type > 	readMap (std::istream &in, const Teuchos::RCP< const comm_type > &comm, const Teuchos::RCP< node_type > &node, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP< const map_type > 	readMap (std::istream &in, const Teuchos::RCP< const comm_type > &comm, const Teuchos::RCP< Teuchos::FancyOStream > &err, const bool tolerant=false, const bool debug=false)`
:x: | `static Teuchos::RCP< const map_type > 	readMap (std::istream &in, const Teuchos::RCP< const comm_type > &comm, const Teuchos::RCP< node_type > &node, const Teuchos::RCP< Teuchos::FancyOStream > &err, const bool tolerant=false, const bool debug=false)`

### Tpetra::Writer

:grey_question: | Command | Comment
-------|--------|---------
:white_check_mark: | `static void 	writeSparseFile (const std::string &filename, const Teuchos::RCP< const sparse_matrix_type > &pMatrix, const std::string &matrixName, const std::string &matrixDescription, const bool debug=false)`
:white_check_mark: | `static void 	writeSparseFile (const std::string &filename, const Teuchos::RCP< const sparse_matrix_type > &pMatrix, const bool debug=false)`
:x: | `static void 	writeSparse (std::ostream &out, const Teuchos::RCP< const sparse_matrix_type > &pMatrix, const std::string &matrixName, const std::string &matrixDescription, const bool debug=false)`
:x: | `static void 	writeSparseGraph (std::ostream &out, const crs_graph_type &graph, const std::string &graphName, const std::string &graphDescription, const bool debug=false)`
:x: | `static void 	writeSparseGraph (std::ostream &out, const crs_graph_type &graph, const bool debug=false)`
:white_check_mark: | `static void 	writeSparseGraphFile (const std::string &filename, const crs_graph_type &graph, const std::string &graphName, const std::string &graphDescription, const bool debug=false)`
:white_check_mark: | `static void 	writeSparseGraphFile (const std::string &filename, const crs_graph_type &graph, const bool debug=false)`
:white_check_mark: | `static void 	writeSparseGraphFile (const std::string &filename, const Teuchos::RCP< const crs_graph_type > &pGraph, const std::string &graphName, const std::string &graphDescription, const bool debug=false)`
:white_check_mark: | `static void 	writeSparseGraphFile (const std::string &filename, const Teuchos::RCP< const crs_graph_type > &pGraph, const bool debug=false)`
:x: | `static void 	writeSparse (std::ostream &out, const Teuchos::RCP< const sparse_matrix_type > &pMatrix, const bool debug=false)`
:x: | `static void 	writeDenseFile (const std::string &filename, const multivector_type &X, const std::string &matrixName, const std::string &matrixDescription, const Teuchos::RCP< Teuchos::FancyOStream > &err=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &dbg=Teuchos::null)`
:x: | `static void 	writeDenseFile (const std::string &filename, const Teuchos::RCP< const multivector_type > &X, const std::string &matrixName, const std::string &matrixDescription, const Teuchos::RCP< Teuchos::FancyOStream > &err=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &dbg=Teuchos::null)`
:x: | `static void 	writeDenseFile (const std::string &filename, const multivector_type &X, const Teuchos::RCP< Teuchos::FancyOStream > &err=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &dbg=Teuchos::null)`
:x: | `static void 	writeDenseFile (const std::string &filename, const Teuchos::RCP< const multivector_type > &X, const Teuchos::RCP< Teuchos::FancyOStream > &err=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &dbg=Teuchos::null)`
:x: | `static void 	writeDense (std::ostream &out, const multivector_type &X, const std::string &matrixName, const std::string &matrixDescription, const Teuchos::RCP< Teuchos::FancyOStream > &err=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &dbg=Teuchos::null)`
:x: | `static void 	writeDense (std::ostream &out, const Teuchos::RCP< const multivector_type > &X, const std::string &matrixName, const std::string &matrixDescription, const Teuchos::RCP< Teuchos::FancyOStream > &err=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &dbg=Teuchos::null)`
:x: | `static void 	writeDense (std::ostream &out, const multivector_type &X, const Teuchos::RCP< Teuchos::FancyOStream > &err=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &dbg=Teuchos::null)`
:x: | `static void 	writeDense (std::ostream &out, const Teuchos::RCP< const multivector_type > &X, const Teuchos::RCP< Teuchos::FancyOStream > &err=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &dbg=Teuchos::null)`
:x: | `static void 	writeMap (std::ostream &out, const map_type &map, const bool debug=false)`
:x: | `static void 	writeMap (std::ostream &out, const map_type &map, const Teuchos::RCP< Teuchos::FancyOStream > &err, const bool debug=false)`
:white_check_mark: | `static void 	writeMapFile (const std::string &filename, const map_type &map)`
:arrow_down: | `static void 	writeOperator (const std::string &fileName, operator_type const &A)`
:x: | `static void 	writeOperator (std::ostream &out, const operator_type &A)`
:arrow_down: | `static void 	writeOperator (const std::string &fileName, const operator_type &A, const Teuchos::ParameterList &params)`
:x: | `static void 	writeOperator (std::ostream &out, const operator_type &A, const Teuchos::ParameterList &params)`

## Belos
