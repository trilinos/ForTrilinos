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
    static Teuchos::RCP < sparse_graph_type > readSparseGraphFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseGraphFile(filenameStr, pComm, callFillComplete, tolerant, debug);
    }
    static Teuchos::RCP < sparse_graph_type > readSparseGraphFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseGraphFile(filenameStr, pComm, constructorParams, fillCompleteParams, tolerant, debug);
    }
    static Teuchos::RCP < sparse_graph_type > readSparseGraphFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseGraphFile(filenameStr, rowMap, colMap, domainMap, rangeMap, callFillComplete, tolerant, debug);
    }
    static Teuchos::RCP < sparse_matrix_type > readSparseFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseFile(filenameStr, pComm, callFillComplete, tolerant, debug);
    }
    static Teuchos::RCP < sparse_matrix_type > readSparseFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const Teuchos::Comm< int > > &pComm, const Teuchos::RCP< Teuchos::ParameterList > &constructorParams, const Teuchos::RCP< Teuchos::ParameterList > &fillCompleteParams, const bool tolerant=false, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseFile(filenameStr, pComm, constructorParams, fillCompleteParams, tolerant, debug);
    }
    static Teuchos::RCP < sparse_matrix_type > readSparseFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const map_type > &rowMap, Teuchos::RCP< const map_type > &colMap, const Teuchos::RCP< const map_type > &domainMap, const Teuchos::RCP< const map_type > &rangeMap, const bool callFillComplete=true, const bool tolerant=false, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      return Tpetra::MatrixMarket::Reader<CMT>::readSparseFile(filenameStr, rowMap, colMap, domainMap, rangeMap, callFillComplete, tolerant, debug);
    }
    static Teuchos::RCP < multivector_type > readDenseFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const comm_type > &comm, Teuchos::RCP< const map_type > &map, const bool tolerant=false, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      return Tpetra::MatrixMarket::Reader<CMT>::readDenseFile(filenameStr, comm, map, tolerant, debug);
    }
    static Teuchos::RCP< const map_type > readMapFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const comm_type > &comm, const bool tolerant=false, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      return Tpetra::MatrixMarket::Reader<CMT>::readMapFile(filenameStr, comm, tolerant, debug);
    }
}
%extend Tpetra::MatrixMarket::Writer<CMT> {
    static void writeSparseFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const sparse_matrix_type > &pMatrix, const std::string &matrixName, const std::string &matrixDescription, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseFile(filenameStr, pMatrix, matrixName, matrixDescription, debug);
    }
    static void writeSparseFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const sparse_matrix_type > &pMatrix, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseFile(filenameStr, pMatrix, debug);
    }
    static void writeSparseGraphFile (std::pair<const char*, size_t> filename, const crs_graph_type &graph, const std::string &graphName, const std::string &graphDescription, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseGraphFile(filenameStr, graph, graphName, graphDescription, debug);
    }
    static void writeSparseGraphFile (std::pair<const char*, size_t> filename, const crs_graph_type &graph, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseGraphFile(filenameStr, graph, debug);
    }
    static void writeSparseGraphFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const crs_graph_type > &pGraph, const std::string &graphName, const std::string &graphDescription, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseGraphFile(filenameStr, pGraph, graphName, graphDescription, debug);
    }
    static void writeSparseGraphFile (std::pair<const char*, size_t> filename, const Teuchos::RCP< const crs_graph_type > &pGraph, const bool debug=false) {
      std::string filenameStr(filename.first, filename.second);
      Tpetra::MatrixMarket::Writer<CMT>::writeSparseGraphFile(filenameStr, pGraph, debug);
    }
    static void writeMapFile (std::pair<const char*, size_t> filename, const map_type &map) {
      std::string filenameStr(filename.first, filename.second);
      Tpetra::MatrixMarket::Writer<CMT>::writeMapFile(filenameStr, map);
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
