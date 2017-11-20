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

**Public Member Functions**

Status | Command| Comment
-------|--------|---------
:white_check_mark: | `Map(numGlobalElements, indexBase, comm)`
:white_check_mark: :star: | `Map (global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, LocalGlobal lg=GloballyDistributed, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())` | No `Node` version
:white_check_mark: :star: | `Map (global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())` | No `Node` version
:arrow_down: | `Map (const global_size_t numGlobalElements, const Kokkos::View< const GlobalOrdinal *, device_type > &indexList, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)` | `Kokkos::View` is planned in the future
:x: | `Map (const global_size_t numGlobalElements, const GlobalOrdinal indexList[], const LocalOrdinal indexListSize, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm)` | The `ArrayView` version is used
:white_check_mark: :star: | `Map (const global_size_t numGlobalElements, const Teuchos::ArrayView< const GlobalOrdinal > &indexList, const GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=defaultArgNode< Node >())` | No `Node` version. Takes in Fortran array instead of `ArrayView`
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
:white_check_mark: | `GlobalOrdinal 	getIndexBase () const`
:white_check_mark: | `LocalOrdinal 	getMinLocalIndex () const`
:white_check_mark: | `LocalOrdinal 	getMaxLocalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMinGlobalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMaxGlobalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMinAllGlobalIndex () const`
:white_check_mark: | `GlobalOrdinal 	getMaxAllGlobalIndex () const`
:white_check_mark: | `LocalOrdinal 	getLocalElement (GlobalOrdinal globalIndex) const`
:white_check_mark: | `GlobalOrdinal 	getGlobalElement (LocalOrdinal localIndex) const`
:white_check_mark: | `local_map_type 	getLocalMap () const`
:white_check_mark: :star: | `LookupStatus 	getRemoteIndexList (const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< LocalOrdinal > &LIDList) const` | Using Fortran arrays instead of `ArrayView`
:white_check_mark: :star: | `LookupStatus 	getRemoteIndexList (const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const` | Using Fortran arrays insted of `ArrayView`
:x: | `global_indices_array_type 	getMyGlobalIndices () const` | Return type is not exposed, needs to be implemented.
:white_check_mark: | `Teuchos::ArrayView< const GlobalOrdinal > 	getNodeElementList () const`


Belos
-----
