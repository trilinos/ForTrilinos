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

// =======================================================================
// Ignore permanently
// =======================================================================
// Ignore everything as all require std::string, and carefully extend
%ignore Tpetra::MatrixMarket::Reader::readSparseGraph;
%ignore Tpetra::MatrixMarket::Reader::readSparse;
%ignore Tpetra::MatrixMarket::Reader::readDense;
%ignore Tpetra::MatrixMarket::Reader::readVector;
%ignore Tpetra::MatrixMarket::Reader::readMap;
%ignore Tpetra::MatrixMarket::Writer::writeSparse;
%ignore Tpetra::MatrixMarket::Writer::writeSparseGraph;
%ignore Tpetra::MatrixMarket::Writer::writeDense;
%ignore Tpetra::MatrixMarket::Writer::writeMap;
%ignore Tpetra::MatrixMarket::Writer::writeOperator;

%inline %{
  typedef Tpetra::CrsMatrix<SC,LO,GO,NO> CMT;
%}

// =======================================================================
// Make interface more Fortran friendly
// =======================================================================
%extend Tpetra::MatrixMarket::Reader<CMT> {
    static Teuchos::RCP < sparse_graph_type > readSparseGraphFile (const std::string& filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false) {
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseGraphFile(filename, pComm, callFillComplete, tolerant, debug);
    }
    static Teuchos::RCP < sparse_graph_type > readSparseGraphFile (const std::string& filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false) {
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseGraphFile(filename, pComm, constructorParams, fillCompleteParams, tolerant, debug);
    }
    static Teuchos::RCP < sparse_graph_type > readSparseGraphFile (const std::string& filename, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false) {
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseGraphFile(filename, rowMap, colMap, domainMap, rangeMap, callFillComplete, tolerant, debug);
    }
    static Teuchos::RCP < sparse_matrix_type > readSparseFile (const std::string& filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false) {
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseFile(filename, pComm, callFillComplete, tolerant, debug);
    }
    static Teuchos::RCP < sparse_matrix_type > readSparseFile (const std::string& filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false) {
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseFile(filename, pComm, constructorParams, fillCompleteParams, tolerant, debug);
    }
    static Teuchos::RCP < sparse_matrix_type > readSparseFile (const std::string& filename, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false) {
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseFile(filename, rowMap, colMap, domainMap, rangeMap, callFillComplete, tolerant, debug);
    }
    static Teuchos::RCP < multivector_type > readDenseFile (const std::string& filename, const Teuchos::RCP< const comm_type > &comm, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false) {
      return Tpetra::MatrixMarket::Reader<CMT>::readDenseFile(filename, comm, map, tolerant, debug);
    }
    static Teuchos::RCP< const map_type > readMapFile (const std::string& filename, const Teuchos::RCP< const comm_type > &comm, const bool tolerant=false, const bool debug=false) {
      return Tpetra::MatrixMarket::Reader<CMT>::readMapFile(filename, comm, tolerant, debug);
    }
}
%extend Tpetra::MatrixMarket::Writer<CMT> {
    static void writeSparseFile (const std::string& filename, const Teuchos::RCP< const sparse_matrix_type > &pMatrix, const std::string &matrixName, const std::string &matrixDescription, const bool debug=false) {
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseFile(filename, pMatrix, matrixName, matrixDescription, debug);
    }
    static void writeSparseFile (const std::string& filename, const Teuchos::RCP< const sparse_matrix_type > &pMatrix, const bool debug=false) {
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseFile(filename, pMatrix, debug);
    }
    static void writeSparseGraphFile (const std::string& filename, const Teuchos::RCP< const crs_graph_type > &pGraph, const std::string &graphName, const std::string &graphDescription, const bool debug=false) {
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseGraphFile(filename, pGraph, graphName, graphDescription, debug);
    }
    static void writeSparseGraphFile (const std::string& filename, const Teuchos::RCP< const crs_graph_type > &pGraph, const bool debug=false) {
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseGraphFile(filename, pGraph, debug);
    }
    static void writeMapFile (const std::string& filename, const map_type &map) {
      Tpetra::MatrixMarket::Writer<CMT>::writeMapFile(filename, map);
    }
}
%ignore Tpetra::MatrixMarket::Reader::readSparseGraphFile;
%ignore Tpetra::MatrixMarket::Reader::readSparseFile;
%ignore Tpetra::MatrixMarket::Reader::readDenseFile;
%ignore Tpetra::MatrixMarket::Reader::readVectorFile;
%ignore Tpetra::MatrixMarket::Reader::readMapFile;
%ignore Tpetra::MatrixMarket::Reader::readVectorVectorFile;
%ignore Tpetra::MatrixMarket::Writer::writeSparseFile;
%ignore Tpetra::MatrixMarket::Writer::writeSparseGraphFile;
%ignore Tpetra::MatrixMarket::Writer::writeDenseFile;
%ignore Tpetra::MatrixMarket::Writer::writeOperator;

// The introduction of typedef seems to work. No clue why.
%teuchos_rcp(Tpetra::MatrixMarket::Reader<CMT>)
%teuchos_rcp(Tpetra::MatrixMarket::Writer<CMT>)

%include "MatrixMarket_Tpetra.hpp"

%template(TpetraReader) Tpetra::MatrixMarket::Reader<CMT>;
%template(TpetraWriter) Tpetra::MatrixMarket::Writer<CMT>;
