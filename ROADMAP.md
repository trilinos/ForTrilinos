# ForTrilinos roadmap

## Key

:white_check_mark: done

:white_check_mark: :star: done with some modifications

:running: - work in progress

:warning: - something is missing

:arrow_down: low priority in the current cycle

:x: ignored

## Package status

Status | Command
---|--------
:running: | Teuchos
:running: | Tpetra
:running: | Belos

Now is the extensive list.

## Teuchos

## Tpetra

### Tpetra::Map

**Constructor/destructor methods**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `Map(numGlobalElements, indexBase, comm)` | No `indexBase`
:white_check_mark: :star: | `Map (global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, LocalGlobal lg=GloballyDistributed, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())` | No `Node` version, no `indexBase`
:white_check_mark: :star: | `Map (global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())` | No `Node` version, no `indexBase`
:arrow_down: | `Map (const global_size_t numGlobalElements, const Kokkos::View< const GlobalOrdinal *, device_type > &indexList, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)` | `Kokkos::View` is planned in the future
:x: | `Map (const global_size_t numGlobalElements, const GlobalOrdinal indexList[], const LocalOrdinal indexListSize, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)` | The `ArrayView` version is used
:white_check_mark: :star: | `Map (const global_size_t numGlobalElements, const Teuchos::ArrayView< const GlobalOrdinal > &indexList, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())` | No `Node` version. Takes in Fortran array instead of `ArrayView`. No `indexBase`.
:white_check_mark: | `Map ()`
:white_check_mark: | `~Map ()`

**Boolean tests**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `bool isNodeLocalElement (LocalOrdinal localIndex) const`
:white_check_mark: | `bool isNodeGlobalElement (GlobalOrdinal globalIndex) const`
:white_check_mark: | `bool isUniform () const`
:white_check_mark: | `bool isContiguous () const`
:white_check_mark: | `bool isDistributed () const`
:white_check_mark: | `bool isCompatible (const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const`
:white_check_mark: | `bool isSameAs (const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const`
:white_check_mark: | `bool locallySameAs (const Map< LocalOrdinal, GlobalOrdinal, node_type > &map) const`


Status | Command| Comment
-------|--------|---------
:white_check_mark: | `Teuchos::RCP< const Teuchos::Comm< int > > 	getComm () const`
:x: | `Teuchos::RCP< Node > 	getNode () const`
:running: | `std::string 	description () const`
:x: | `void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`
:x: | `template<class NodeOut > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, NodeOut > > 	clone (const Teuchos::RCP< NodeOut > &nodeOut) const`
:white_check_mark: | `Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	removeEmptyProcesses () const `
:white_check_mark: | `Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	replaceCommWithSubset (const Teuchos::RCP< const Teuchos::Comm< int > > &newComm) const`

**Related functions (non-members)**

Status | Command| Comment
-------|--------|---------
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal > > 	createLocalMap (const size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createLocalMapWithNode (const size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::null)`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal > > 	createUniformContigMap (const global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createUniformContigMapWithNode (const global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::null)`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal > > 	createContigMap (const global_size_t numElements, const size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createContigMapWithNode (const global_size_t numElements, const size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal > > 	createNonContigMap (const Teuchos::ArrayView< const GlobalOrdinal > &elementList, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)`
:arrow_down: |  `emplate<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createNonContigMapWithNode (const Teuchos::ArrayView< const GlobalOrdinal > &elementList, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::null)`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createWeightedContigMapWithNode (const int thisNodeWeight, const global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=Teuchos::null)`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createOneToOne (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &M)`
:arrow_down: | ` template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Map < LocalOrdinal, GlobalOrdinal, Node > > 	createOneToOne (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &M, const Tpetra::Details::TieBreak< LocalOrdinal, GlobalOrdinal > &tie_break)`
:arrow_down: | ` template<class LocalOrdinal , class GlobalOrdinal , class Node >bool 	operator== (const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map1, const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map2)`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > bool 	operator!= (const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map1, const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > &map2)`

**Attributes**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `bool 	isOneToOne () const`
:white_check_mark: | `global_size_t 	getGlobalNumElements () const`
:white_check_mark: | `size_t 	getNodeNumElements () const`
:x: | `GlobalOrdinal 	getIndexBase () const`
:white_check_mark: :star: | `LocalOrdinal 	getMinLocalIndex () const` | Returns 1-based index
:white_check_mark: :star: | `LocalOrdinal 	getMaxLocalIndex () const` | Returns 1-based index
:white_check_mark: | `GlobalOrdinal 	getMinGlobalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMaxGlobalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMinAllGlobalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMaxAllGlobalIndex () const`
:white_check_mark: :star: | `LocalOrdinal 	getLocalElement (GlobalOrdinal globalIndex) const` | Returns 1-based index
:white_check_mark: | `GlobalOrdinal 	getGlobalElement (LocalOrdinal localIndex) const`
:white_check_mark: | `local_map_type 	getLocalMap () const`
:white_check_mark: :star: | `LookupStatus 	getRemoteIndexList (const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< LocalOrdinal > &LIDList) const` | Using Fortran arrays instead of `ArrayView`, 1-based arrays
:white_check_mark: :star: | `LookupStatus 	getRemoteIndexList (const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const` | Using Fortran arrays instead of `ArrayView`, 1-based arrays
:arrow_down: | `global_indices_array_type 	getMyGlobalIndices () const` | Return type is not exposed, needs to be implemented.
:white_check_mark: :star: | `Teuchos::ArrayView< const GlobalOrdinal > 	getNodeElementList () const` | converted to subroutine


### Tpetra::Export

**Constructor/Destructor Methods**

Status | Command| Comment
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

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `size_t 	getNumSameIDs () const`
:white_check_mark: | `size_t 	getNumPermuteIDs () const`
:running: | `Teuchos::ArrayView< const LocalOrdinal > 	getPermuteFromLIDs () const`
:running: | `Teuchos::ArrayView< const LocalOrdinal > 	getPermuteToLIDs () const`
:white_check_mark: | `size_t 	getNumRemoteIDs () const`
:running: | `Teuchos::ArrayView< const LocalOrdinal > 	getRemoteLIDs () const`
:white_check_mark: | `size_t 	getNumExportIDs () const`
:running: | `Teuchos::ArrayView< const LocalOrdinal > 	getExportLIDs () const`
:running: | `Teuchos::ArrayView< const int > 	getExportPIDs () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getSourceMap () const`
:white_check_mark: | `Teuchos::RCP< const map_type > 	getTargetMap () const`
:x: | `Distributor & 	getDistributor () const`
:white_check_mark: | `bool 	isLocallyComplete () const`
:arrow_down: | `Export< LocalOrdinal, GlobalOrdinal, Node > & 	operator= (const Export< LocalOrdinal, GlobalOrdinal, Node > &rhs)`

**I/O Methods**

Status | Command| Comment
-------|--------|---------
:x: | `virtual void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`
:x: | `virtual void 	print (std::ostream &os) const`

**Related Functions (non-member functions)**

Status | Command| Comment
-------|--------|---------
:x: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Export < LocalOrdinal, GlobalOrdinal, Node > > 	createExport (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &src, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &tgt)`

### Tpetra::Import

**Constructor/Destructor Methods**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target)`
:x: | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::FancyOStream > &out)`
:white_check_mark: | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::ParameterList > &plist)`
:x: | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::FancyOStream > &out, const Teuchos::RCP< Teuchos::ParameterList > &plist)`
? | `Import (const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, Teuchos::Array< int > &remotePIDs)`
:white_check_mark: | `Import (const Import< LocalOrdinal, GlobalOrdinal, Node > &importer)`
:white_check_mark: | `Import (const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter)`
? | `Import (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &source, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &target, Teuchos::Array< int > &userRemotePIDs, Teuchos::Array< GlobalOrdinal > &remoteGIDs, const Teuchos::ArrayView< const LocalOrdinal > &userExportLIDs, const Teuchos::ArrayView< const int > &userExportPIDs, const bool useRemotePIDs, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null, const Teuchos::RCP< Teuchos::FancyOStream > &out=Teuchos::null)`
:white_check_mark: | `virtual 	~Import ()`
:white_check_mark: | `void 	setParameterList (const Teuchos::RCP< Teuchos::ParameterList > &plist)`

**Import Attribute Methods**

Status | Command| Comment
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

Status | Command| Comment
-------|--------|---------
:x: | `virtual void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`
:x: | `virtual void 	print (std::ostream &os) const`

**Related Functions (non-member functions)**

Status | Command| Comment
-------|--------|---------
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Import < LocalOrdinal, GlobalOrdinal, Node > > 	createImport (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &src, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &tgt)`
:arrow_down: | `template<class LocalOrdinal , class GlobalOrdinal , class Node > Teuchos::RCP< const Import < LocalOrdinal, GlobalOrdinal, Node > > 	createImport (const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &src, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &tgt, const Teuchos::RCP< Teuchos::ParameterList > &plist)`

### Tpetra::MultiVector

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`
:white_check_mark: | `void 	setCopyOrView (const Teuchos::DataAccess copyOrView)`
:white_check_mark: | `Teuchos::DataAccess 	getCopyOrView () const`
:white_check_mark: | `void 	assign (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &src)`

**Constructors and destructor**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `MultiVector ()`
:white_check_mark: | `MultiVector (const Teuchos::RCP< const map_type > &map, const size_t numVecs, const bool zeroOut=true)`
:white_check_mark: | `MultiVector (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &source)`
:white_check_mark: | `MultiVector (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &source, const Teuchos::DataAccess copyOrView)`
:white_check_mark: :star:| `MultiVector (const Teuchos::RCP< const map_type > &map, const Teuchos::ArrayView< const Scalar > &A, const size_t LDA, const size_t NumVectors)` | Fortran array version
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

Status | Command| Comment
-------|--------|---------
:x: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subCopy (const Teuchos::Range1D &colRng) const` | Prefer `ArrayView` version
:white_check_mark: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subCopy (const Teuchos::ArrayView< const size_t > &cols) const` | 1-based
:arrow_down: | `Teuchos::RCP< const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subView (const Teuchos::Range1D &colRng) const` | Prefer `ArrayView` version
:white_check_mark: | `Teuchos::RCP< const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subView (const Teuchos::ArrayView< const size_t > &cols) const` | 1-based
:arrow_down: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subViewNonConst (const Teuchos::Range1D &colRng)` | Prefer `ArrayView` version
:white_check_mark: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	subViewNonConst (const Teuchos::ArrayView< const size_t > &cols)` | 1-based
:white_check_mark: | `Teuchos::RCP< const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	offsetView (const Teuchos::RCP< const map_type > &subMap, const size_t offset) const`
:white_check_mark: | `Teuchos::RCP< MultiVector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	offsetViewNonConst (const Teuchos::RCP< const map_type > &subMap, const size_t offset)`
:white_check_mark: :star: | `Teuchos::RCP< const Vector < Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	getVector (const size_t j) const` | 1-based
:white_check_mark: :star:  | `Teuchos::RCP< Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > > 	getVectorNonConst (const size_t j)` | 1-based
:running: | `Teuchos::ArrayRCP< const Scalar > 	getData (size_t j) const`
:running: | `Teuchos::ArrayRCP< Scalar > 	getDataNonConst (size_t j)`
:arrow_down: | `void 	get1dCopy (const Teuchos::ArrayView< Scalar > &A, const size_t LDA) const`
:arrow_down: | `void 	get2dCopy (const Teuchos::ArrayView< const Teuchos::ArrayView< Scalar > > &ArrayOfPtrs) const`
:arrow_down: | `Teuchos::ArrayRCP< const Scalar > 	get1dView () const`
:arrow_down: | `Teuchos::ArrayRCP < Teuchos::ArrayRCP< const Scalar > > 	get2dView () const`
:arrow_down: | `Teuchos::ArrayRCP< Scalar > 	get1dViewNonConst ()`
:arrow_down: | `Teuchos::ArrayRCP < Teuchos::ArrayRCP< Scalar > > 	get2dViewNonConst ()`
:x: | `dual_view_type 	getDualView () const`
:x: | `template<class TargetDeviceType > void 	sync ()`
:x: | `template<class TargetDeviceType > bool 	need_sync () const`
:x: | `template<class TargetDeviceType > void 	modify ()`
:x: | `template<class TargetDeviceType > Kokkos::Impl::if_c < std::is_same< typename device_type::memory_space, typename TargetDeviceType::memory_space > ::value, typename dual_view_type::t_dev, typename dual_view_type::t_host >::type 	getLocalView () const`

**Mathematical methods**

Status | Command| Comment
-------|--------|---------
:white_check_mark: :star: | `void 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Teuchos::ArrayView< dot_type > &dots) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < dot_type, T >::value), void > ::type 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Teuchos::ArrayView< T > &dots) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < dot_type, T >::value), void > ::type 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, std::vector< T > &dots) const`
:x: | `void 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Kokkos::View< dot_type *, device_type > &dots) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < dot_type, T >::value), void > ::type 	dot (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Kokkos::View< T *, device_type > &dots) const`
:white_check_mark: | `void 	abs (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A)`
:white_check_mark: | `void 	reciprocal (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A)`
:white_check_mark: | `void 	scale (const Scalar &alpha)`
:white_check_mark: :star: | `void 	scale (const Teuchos::ArrayView< const Scalar > &alpha)` | Fortran array
:x: | `void 	scale (const Kokkos::View< const impl_scalar_type *, device_type > &alpha)`
:white_check_mark: | `void 	scale (const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A)`
:white_check_mark: | `void 	update (const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Scalar &beta)`
:white_check_mark: | `void 	update (const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Scalar &beta, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, const Scalar &gamma)`
:x: | `void 	norm1 (const Kokkos::View< mag_type *, device_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	norm1 (const Kokkos::View< T *, device_type > &norms) const`
:white_check_mark: :star: | `void 	norm1 (const Teuchos::ArrayView< mag_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	norm1 (const Teuchos::ArrayView< T > &norms) const`
:x: | `void 	norm2 (const Kokkos::View< mag_type *, device_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	norm2 (const Kokkos::View< T *, device_type > &norms) const`
:white_check_mark: :star: | `void 	norm2 (const Teuchos::ArrayView< mag_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	norm2 (const Teuchos::ArrayView< T > &norms) const`
:x: | `void 	normInf (const Kokkos::View< mag_type *, device_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	normInf (const Kokkos::View< T *, device_type > &norms) const`
:white_check_mark: :star: | `void 	normInf (const Teuchos::ArrayView< mag_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type 	normInf (const Teuchos::ArrayView< T > &norms) const`
:x: | `void TPETRA_DEPRECATED 	normWeighted (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &weights, const Teuchos::ArrayView< mag_type > &norms) const`
:x: | `template<typename T > std::enable_if< !(std::is_same < mag_type, T >::value), void > ::type TPETRA_DEPRECATED 	normWeighted (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &weights, const Teuchos::ArrayView< T > &norms) const`
:white_check_mark: :star: | `void 	meanValue (const Teuchos::ArrayView< impl_scalar_type > &means) const` | Fortran array
:x: | `template<typename T > std::enable_if<!std::is_same < impl_scalar_type, T >::value, void >::type 	meanValue (const Teuchos::ArrayView< T > &means) const`
:white_check_mark: | `void 	multiply (Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, const Scalar &beta)`
:arrow_down: | `void 	elementWiseMultiply (Scalar scalarAB, const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &B, Scalar scalarThis)`

**Attribute access functions**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `size_t 	getNumVectors () const`
:white_check_mark: | `size_t 	getLocalLength () const`
:white_check_mark: | `global_size_t 	getGlobalLength () const`
:white_check_mark: | `size_t 	getStride () const`
:white_check_mark: | `bool 	isConstantStride () const`

**Overridden from Teuchos::Describable**

Status | Command| Comment
-------|--------|---------
:running: | `virtual std::string 	description () const`
:x: | `virtual void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`

**Public methods for redistributing data**

Status | Command| Comment
-------|--------|---------
:x: | `void 	doImport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`
:x: | `void 	doImport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:x: | `void 	doExport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:x: | `void 	doExport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`

**Attribute accessor methods**

Status | Command| Comment
-------|--------|---------
:arrow_down: | `bool 	isDistributed () const`
:arrow_down: | `virtual Teuchos::RCP< const map_type > 	getMap () const`

**I/O methods**

Status | Command| Comment
-------|--------|---------
:x: | `void 	print (std::ostream &os) const`

**Methods for use only by experts**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`

**Friends**

Status | Command| Comment
-------|--------|---------
:x: | `template<class DS , class DL , class DG , class DN , const bool dstClassic, class SS , class SL , class SG , class SN , const bool srcClassic> void 	deep_copy (MultiVector< DS, DL, DG, DN, dstClassic > &dst, const MultiVector< SS, SL, SG, SN, srcClassic > &src)`

**Related Functions (non-member functions)**

Status | Command| Comment
-------|--------|---------
:x: | `template<class DS , class DL , class DG , class DN , const bool dstClassic, class SS , class SL , class SG , class SN , const bool srcClassic> void 	deep_copy (MultiVector< DS, DL, DG, DN, dstClassic > &dst, const MultiVector< SS, SL, SG, SN, srcClassic > &src)`
:x: | `template<class ST , class LO , class GO , class NT , const bool classic = NT::classic> MultiVector< ST, LO, GO, NT, classic > 	createCopy (const MultiVector< ST, LO, GO, NT, classic > &src)`

**Post-construction modification routines**

Status | Command| Comment
-------|--------|---------
:x: | `static const bool 	useAtomicUpdatesByDefault`
:white_check_mark: :star: | `void 	replaceGlobalValue (const GlobalOrdinal gblRow, const size_t col, const impl_scalar_type &value) const` | 1-based
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	replaceGlobalValue (GlobalOrdinal globalRow, size_t col, const T &value) const`
:white_check_mark: | `void 	sumIntoGlobalValue (const GlobalOrdinal gblRow, const size_t col, const impl_scalar_type &value, const bool atomic=useAtomicUpdatesByDefault) const` | 1-based
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	sumIntoGlobalValue (const GlobalOrdinal gblRow, const size_t col, const T &val, const bool atomic=useAtomicUpdatesByDefault) const`
:white_check_mark: | `void 	replaceLocalValue (const LocalOrdinal lclRow, const size_t col, const impl_scalar_type &value) const` | 1-based
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	replaceLocalValue (const LocalOrdinal lclRow, const size_t col, const T &val) const`
:white_check_mark: | `void 	sumIntoLocalValue (const LocalOrdinal lclRow, const size_t col, const impl_scalar_type &val, const bool atomic=useAtomicUpdatesByDefault) const`
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	sumIntoLocalValue (const LocalOrdinal lclRow, const size_t col, const T &val, const bool atomic=useAtomicUpdatesByDefault) const`
:white_check_mark: | `void 	putScalar (const Scalar &value)`
:x: | `template<typename T > std::enable_if<!std::is_same < T, impl_scalar_type >::value &&std::is_convertible< T, impl_scalar_type >::value, void >::type 	putScalar (const T &value)`
:white_check_mark: | `void 	randomize ()`
:white_check_mark: | `void 	randomize (const Scalar &minVal, const Scalar &maxVal)`
:white_check_mark: | `void 	replaceMap (const Teuchos::RCP< const map_type > &map)`
:white_check_mark: | `void 	reduce ()`
:x: | `MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > & 	operator= (const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &source)`

### Tpetra::CrsGraph

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `bool 	haveGlobalConstants () const`
:white_check_mark: | `void 	computeGlobalConstants ()`
:x: | `local_graph_type 	getLocalGraph () const`
 
**Constructor/Destructor Methods**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:x: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow, const ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: :star: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::ArrayRCP< const size_t > &numEntPerRow, const ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Use Fortran arrays
:white_check_mark: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const size_t maxNumEntriesPerRow, const ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:x: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: :star: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const Teuchos::ArrayRCP< const size_t > &numEntPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Use Fortran arrays
:x: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const typename local_graph_type::row_map_type &rowPointers, const typename local_graph_type::entries_type::non_const_type &columnIndices, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: :x: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const Teuchos::ArrayRCP< size_t > &rowPointers, const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)` | Use Fortran arrays; 1-based; may have ownership problems
:x: | `CrsGraph (const Teuchos::RCP< const map_type > &rowMap, const Teuchos::RCP< const map_type > &colMap, const local_graph_type &lclGraph, const Teuchos::RCP< Teuchos::ParameterList > &params)`
:x: | `template<class Node2 > Teuchos::RCP< CrsGraph < LocalOrdinal, GlobalOrdinal, Node2, Node2::classic > > 	clone (const Teuchos::RCP< Node2 > &node2, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) const`
:white_check_mark: | `virtual 	~CrsGraph ()`
:white_check_mark: | `void 	setParameterList (const Teuchos::RCP< Teuchos::ParameterList > &params)`
:white_check_mark: | `Teuchos::RCP< const Teuchos::ParameterList > 	getValidParameters () const`

**Insertion/Removal Methods**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `void 	insertGlobalIndices (const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &indices)`
:white_check_mark: | `void 	insertGlobalIndices (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const GlobalOrdinal inds[])`
:white_check_mark: | `void 	insertLocalIndices (const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &indices)` | Use Fortran arrays; 1-based
:x: | `void 	insertLocalIndices (const LocalOrdinal localRow, const LocalOrdinal numEnt, const LocalOrdinal inds[])` | Prefer `ArrayView` variant
:white_check_mark: | `void 	removeLocalIndices (LocalOrdinal localRow)`

**Transformational Methods**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `void 	globalAssemble ()`
:white_check_mark: | `void 	resumeFill (const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	fillComplete (const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	fillComplete (const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`
:white_check_mark: | `void 	expertStaticFillComplete (const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const Teuchos::RCP< const import_type > &importer=Teuchos::null, const Teuchos::RCP< const export_type > &exporter=Teuchos::null, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)`

**Methods implementing RowGraph.**

Status | Command| Comment
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
:white_check_mark: | `size_t 	getNumEntriesInLocalRow (LocalOrdinal localRow) const`
:white_check_mark: | `size_t 	getNodeAllocationSize () const`
:white_check_mark: | `size_t 	getNumAllocatedEntriesInGlobalRow (GlobalOrdinal globalRow) const`
:white_check_mark: | `size_t 	getNumAllocatedEntriesInLocalRow (LocalOrdinal localRow) const`
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
:white_check_mark: | `void 	getGlobalRowCopy (GlobalOrdinal GlobalRow, const Teuchos::ArrayView< GlobalOrdinal > &Indices, size_t &NumIndices) const`
:running: | `void 	getLocalRowCopy (LocalOrdinal LocalRow, const Teuchos::ArrayView< LocalOrdinal > &indices, size_t &NumIndices) const`
:white_check_mark: | `void 	getGlobalRowView (const GlobalOrdinal gblRow, Teuchos::ArrayView< const GlobalOrdinal > &gblColInds) const`
:white_check_mark: | `bool 	supportsRowViews () const`
:running: | `void 	getLocalRowView (const LocalOrdinal lclRow, Teuchos::ArrayView< const LocalOrdinal > &lclColInds) const`
 
**Overridden from Teuchos::Describable**

Status | Command| Comment
-------|--------|---------
:running: | `std::string 	description () const`
:x: | `void 	describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const`
 
**Implementation of DistObject**

Status | Command| Comment
-------|--------|---------
:x: | `virtual bool 	checkSizes (const SrcDistObject &source)`
:x: | `virtual void 	copyAndPermute (const SrcDistObject &source, size_t numSameIDs, const Teuchos::ArrayView< const LocalOrdinal > &permuteToLIDs, const Teuchos::ArrayView< const LocalOrdinal > &permuteFromLIDs)`
:x: | `virtual void 	packAndPrepare (const SrcDistObject &source, const Teuchos::ArrayView< const LocalOrdinal > &exportLIDs, Teuchos::Array< GlobalOrdinal > &exports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t &constantNumPackets, Distributor &distor)`
:x: | `virtual void 	pack (const Teuchos::ArrayView< const LocalOrdinal > &exportLIDs, Teuchos::Array< GlobalOrdinal > &exports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t &constantNumPackets, Distributor &distor) const`
:x: | `virtual void 	unpackAndCombine (const Teuchos::ArrayView< const LocalOrdinal > &importLIDs, const Teuchos::ArrayView< const GlobalOrdinal > &imports, const Teuchos::ArrayView< size_t > &numPacketsPerLID, size_t constantNumPackets, Distributor &distor, CombineMode CM)`

**Advanced methods, at increased risk of deprecation.**

Status | Command| Comment
-------|--------|---------
:x: | `void 	getLocalDiagOffsets (const Kokkos::View< size_t *, device_type, Kokkos::MemoryUnmanaged > &offsets) const`
:running: | `void 	getLocalDiagOffsets (Teuchos::ArrayRCP< size_t > &offsets) const`
:running: | `void 	getNumEntriesPerLocalRowUpperBound (Teuchos::ArrayRCP< const size_t > &boundPerLocalRow, size_t &boundForAllLocalRows, bool &boundSameForAllLocalRows) const`
:x: | `void 	setAllIndices (const typename local_graph_type::row_map_type &rowPointers, const typename local_graph_type::entries_type::non_const_type &columnIndices)`
:running: | `void 	setAllIndices (const Teuchos::ArrayRCP< size_t > &rowPointers, const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices)`
:running: | `Teuchos::ArrayRCP< const size_t > 	getNodeRowPtrs () const`
:running: | `Teuchos::ArrayRCP< const LocalOrdinal > 	getNodePackedIndices () const`
:white_check_mark: | `void 	replaceColMap (const Teuchos::RCP< const map_type > &newColMap)`
:white_check_mark: | `void 	reindexColumns (const Teuchos::RCP< const map_type > &newColMap, const Teuchos::RCP< const import_type > &newImport=Teuchos::null, const bool sortIndicesInEachRow=true)`
:white_check_mark: | `void 	replaceDomainMapAndImporter (const Teuchos::RCP< const map_type > &newDomainMap, const Teuchos::RCP< const import_type > &newImporter)`
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`
:x: | `void 	doImport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`
:x: | `void 	doImport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:x: | `void 	doExport (const SrcDistObject &source, const Export< LocalOrdinal, GlobalOrdinal, Node > &exporter, CombineMode CM)`
:x: | `void 	doExport (const SrcDistObject &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)`

**Attribute accessor methods**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `bool 	isDistributed () const`
:arrow_down: | `virtual Teuchos::RCP< const map_type > 	getMap () const`

**I/O methods**

Status | Command| Comment
-------|--------|---------
:x: | `void 	print (std::ostream &os) const`

**Methods for use only by experts**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `virtual void 	removeEmptyProcessesInPlace (const Teuchos::RCP< const map_type > &newMap)`

### Tpetra::CrsMatrix

## Belos
