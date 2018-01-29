/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"
%import <std_string.i>

%{
#include "Teuchos_RCP.hpp"
#include "MatrixMarket_Tpetra.hpp"
%}

%inline %{
  typedef Tpetra::CrsMatrix<SC,LO,GO,NO> CMT;
%}

// Declare these instantiations as RCPs
%teuchos_rcp(Tpetra::MatrixMarket::Reader<CMT>)
%teuchos_rcp(Tpetra::MatrixMarket::Writer<CMT>)

// =======================================================================
// Declare selected typedefs and functions
//
// (Easiest way to ignore all the default arguments)
// =======================================================================

namespace Tpetra {
namespace MatrixMarket {
template<class SparseMatrixType>
class Reader {
private:
    Reader();
public:
    typedef Teuchos::Comm<int> comm_type;
    typedef SparseMatrixType sparse_matrix_type;
    typedef Teuchos::RCP<sparse_matrix_type> sparse_matrix_ptr;

    typedef typename SparseMatrixType::scalar_type scalar_type;
    typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename SparseMatrixType::node_type node_type;

    typedef CrsGraph<local_ordinal_type, global_ordinal_type, node_type> sparse_graph_type;
    typedef MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> multivector_type;
    typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

%define TPETRA_READ_METHOD(RTYPE, METHOD)
    static Teuchos::RCP < RTYPE > METHOD(const std::string& filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false);
    static Teuchos::RCP < RTYPE > METHOD(const std::string& filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false);
    static Teuchos::RCP < RTYPE > METHOD(const std::string& filename, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false);
%enddef
    TPETRA_READ_METHOD(sparse_graph_type, readSparseGraphFile)
    TPETRA_READ_METHOD(sparse_matrix_type, readSparseFile)
    static Teuchos::RCP < multivector_type > readDenseFile (const std::string& filename, const Teuchos::RCP< const comm_type > &comm, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false);
    static Teuchos::RCP< const map_type > readMapFile (const std::string& filename, const Teuchos::RCP< const comm_type > &comm, const bool tolerant=false, const bool debug=false);
#undef TPETRA_READ_METHOD
};

template<class SparseMatrixType>
class Writer {
public:
    typedef Teuchos::Comm<int> comm_type;
    typedef SparseMatrixType sparse_matrix_type;
    typedef typename SparseMatrixType::scalar_type scalar_type;
    typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename SparseMatrixType::node_type node_type;

    typedef CrsGraph<local_ordinal_type, global_ordinal_type, node_type> crs_graph_type;
    typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

%define TPETRA_WRITE_METHOD(INTYPE, METHOD)
    static void METHOD (const std::string& filename, const Teuchos::RCP< const INTYPE > &input, const std::string &name, const std::string &description, const bool debug=false);
    static void METHOD (const std::string& filename, const Teuchos::RCP< const INTYPE > &input, const bool debug=false);
%enddef

    TPETRA_WRITE_METHOD(sparse_matrix_type, writeSparseFile)
    TPETRA_WRITE_METHOD(crs_graph_type, writeSparseGraphFile)
    static void writeMapFile (const std::string& filename, const map_type &map);

 #undef TPETRA_WRITE_METHOD
};
} // end namespace Tpetra::MatrixMarket
} // end namespace Tpetra

%template(TpetraReader) Tpetra::MatrixMarket::Reader<CMT>;
%template(TpetraWriter) Tpetra::MatrixMarket::Writer<CMT>;
