/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%{
#include "MatrixMarket_Tpetra.hpp"
%}

// Function signatures for local quantities are incorrectly declared as size_t
%apply LO { size_t getNumEntriesInLocalRow };

// =======================================================================
// Just declare a few methods
// (There are many overloads that we don't want)
// =======================================================================

namespace Tpetra {
namespace MatrixMarket {
template<class SparseMatrixType>
class Reader {
  public:
    typedef SparseMatrixType sparse_matrix_type;

    typedef typename SparseMatrixType::scalar_type scalar_type;
    typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename SparseMatrixType::node_type node_type;

    typedef CrsGraph<local_ordinal_type,
                     global_ordinal_type,
                     node_type> sparse_graph_type;
    typedef MultiVector<scalar_type,
                        local_ordinal_type,
                        global_ordinal_type,
                        node_type> multivector_type;
    typedef Teuchos::Comm<int> comm_type;
    typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    static Teuchos::RCP < sparse_graph_type > readSparseGraphFile (const std::string& filename, const Teuchos::RCP< const comm_type > &pComm, const bool callFillComplete=true);
    static Teuchos::RCP < sparse_graph_type > readSparseGraphFile (const std::string& filename, const Teuchos::RCP< const comm_type > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams);
    static Teuchos::RCP < sparse_graph_type > readSparseGraphFile (const std::string& filename, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true);
    static Teuchos::RCP < sparse_matrix_type > readSparseFile (const std::string& filename, const Teuchos::RCP< const comm_type > &pComm, const bool callFillComplete=true);
    static Teuchos::RCP < sparse_matrix_type > readSparseFile (const std::string& filename, const Teuchos::RCP< const comm_type > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams);
    static Teuchos::RCP < sparse_matrix_type > readSparseFile (const std::string& filename, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true);
    static Teuchos::RCP < multivector_type > readDenseFile (const std::string& filename, const Teuchos::RCP< const comm_type > &comm, Teuchos::RCP< const map_type > &map);
    static Teuchos::RCP< const map_type > readMapFile (const std::string& filename, const Teuchos::RCP< const comm_type > &comm);
  private:
    Reader();
};

template<class SparseMatrixType>
class Writer {
  public:
    typedef SparseMatrixType sparse_matrix_type;

    typedef typename SparseMatrixType::scalar_type scalar_type;
    typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename SparseMatrixType::node_type node_type;

    typedef CrsGraph<local_ordinal_type,
                     global_ordinal_type,
                     node_type> crs_graph_type;
    typedef MultiVector<scalar_type,
                        local_ordinal_type,
                        global_ordinal_type,
                        node_type> multivector_type;
    typedef Teuchos::Comm<int> comm_type;
    typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    static void writeSparseFile (const std::string& filename, const Teuchos::RCP< const sparse_matrix_type > &pMatrix, const std::string &matrixName, const std::string &matrixDescription);
    static void writeSparseFile (const std::string& filename, const Teuchos::RCP< const sparse_matrix_type > &pMatrix);
    static void writeSparseGraphFile (const std::string& filename, const Teuchos::RCP< const crs_graph_type > &pGraph, const std::string &graphName, const std::string &graphDescription);
    static void writeSparseGraphFile (const std::string& filename, const Teuchos::RCP< const crs_graph_type > &pGraph);
    static void writeDenseFile (const std::string &filename, const Teuchos::RCP< const multivector_type > &X, const std::string &matrixName, const std::string &matrixDescription);
    static void writeDenseFile (const std::string &filename, const Teuchos::RCP< const multivector_type > &X);
  private:
    Writer();
};
} // namespace MatrixMarket
} // namespace Tpetra

%teuchos_rcp(Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<SC,LO,GO,NO> >)
%teuchos_rcp(Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<SC,LO,GO,NO> >)

%template(TpetraReader) Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<SC,LO,GO,NO> >;
%template(TpetraWriter) Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<SC,LO,GO,NO> >;
