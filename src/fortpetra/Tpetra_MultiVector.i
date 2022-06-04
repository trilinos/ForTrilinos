/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%{
#include <Tpetra_MultiVector.hpp>
%}

// Treat array RCP return values as array views
// (they will not reference count, though, of course)
%apply Teuchos::ArrayView<double> { Teuchos::ArrayRCP<double> };
%apply Teuchos::ArrayView<const double> { Teuchos::ArrayRCP<const double> };

// Function signatures for local quantities are incorrectly declared as size_t
%apply LO { size_t getLocalLength };

// Add doImport and doExport
%tpetra_extend_with_import_export(Tpetra::MultiVector<SC,LO,GO,NO>)

// Fix ±1 issues
%apply int INDEX { size_t j, size_t col, int lclRow }
%apply const Teuchos::ArrayView<const int>& INDEX { const Teuchos::ArrayView<const size_t>& cols }

/*
 * Note: the MultiVector here are derived directly from Trilinos 13.2:
 * - Protected and private methods removed
 * - Inline definitions removed
 * - Previously %ignored functions removed
 * - Unused templated methods removed
 * - Documentation removed
 *
 * using the following license.
 *
 * ---
 *
 *          Tpetra: Templated Linear Algebra Services Package
 *                 Copyright (2008) Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Michael A. Heroux (maherou@sandia.gov)
 */
namespace Tpetra {
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class MultiVector :
    public DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
  public:
    using scalar_type = Scalar;
    using impl_scalar_type =
      typename Kokkos::Details::ArithTraits<Scalar>::val_type;
    using map_type = Map<LocalOrdinal, GlobalOrdinal, Node>;
    using local_ordinal_type = typename map_type::local_ordinal_type;
    using global_ordinal_type = typename map_type::global_ordinal_type;
    using device_type = typename map_type::device_type;
    using node_type = typename map_type::node_type;

    using dot_type =
      typename Kokkos::Details::InnerProductSpaceTraits<impl_scalar_type>::dot_type;

    using mag_type = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;

    using execution_space = typename device_type::execution_space;

    MultiVector();
    MultiVector(const Teuchos::RCP<const map_type>& map,
                const size_t numVecs,
                const bool zeroOut = true);
    MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
                 const Teuchos::DataAccess copyOrView);
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const Teuchos::ArrayView<const Scalar>& A,
                 const size_t LDA,
                 const size_t NumVectors);
    MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                 const Teuchos::RCP<const map_type>& subMap,
                 const local_ordinal_type rowOffset = 0);
    MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&) = default;

    virtual ~MultiVector () = default;

    //! Swap contents of \c mv with contents of \c *this.
    void swap (MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& mv);

    void
    replaceGlobalValue (const GlobalOrdinal gblRow,
                        const size_t col,
                        const impl_scalar_type& value);

    void
    sumIntoGlobalValue (const GlobalOrdinal gblRow,
                        const size_t col,
                        const impl_scalar_type& value,
                        const bool atomic = useAtomicUpdatesByDefault);

    void
    replaceLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const impl_scalar_type& value);

    void
    sumIntoLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const impl_scalar_type& val,
                       const bool atomic = useAtomicUpdatesByDefault);

    void putScalar (const Scalar& value);

    void randomize();

    void randomize (const Scalar& minVal, const Scalar& maxVal);

    void replaceMap (const Teuchos::RCP<const map_type>& map);

    void reduce();
    //! Return a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subCopy (const Teuchos::ArrayView<const size_t>& cols) const;

    //! Return a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subView (const Teuchos::ArrayView<const size_t>& cols) const;

    //! Return a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subViewNonConst (const Teuchos::ArrayView<const size_t>& cols);

    Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    offsetView (const Teuchos::RCP<const map_type>& subMap,
                const size_t offset) const;

    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                        const size_t offset);

    //! Const view of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP<const Scalar> getData (size_t j) const;

    //! View of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP<Scalar> getDataNonConst (size_t j);

    void
    get1dCopy (const Teuchos::ArrayView<Scalar>& A,
               const size_t LDA) const;

    Teuchos::ArrayRCP<const Scalar> get1dView () const;

    Teuchos::ArrayRCP<Scalar> get1dViewNonConst ();

    void
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Teuchos::ArrayView<dot_type>& dots) const;

    void
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Kokkos::View<dot_type*, Kokkos::HostSpace>& norms) const;
    //! Put element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

    //! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

    void scale (const Scalar& alpha);

    void scale (const Teuchos::ArrayView<const Scalar>& alpha);

    void scale (const Kokkos::View<const impl_scalar_type*, device_type>& alpha);

    void
    scale (const Scalar& alpha,
           const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

    void
    update (const Scalar& alpha,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            const Scalar& beta);

    void
    update (const Scalar& alpha,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            const Scalar& beta,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
            const Scalar& gamma);

    void
    norm1 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;

    void norm1 (const Teuchos::ArrayView<mag_type>& norms) const;

    void
    norm2 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;

    void norm2 (const Teuchos::ArrayView<mag_type>& norms) const;

    void normInf (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;

    void normInf (const Teuchos::ArrayView<mag_type>& norms) const;

    void meanValue (const Teuchos::ArrayView<impl_scalar_type>& means) const;

    void
    multiply (Teuchos::ETransp transA,
              Teuchos::ETransp transB,
              const Scalar& alpha,
              const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
              const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
              const Scalar& beta);

    size_t getNumVectors() const;
    size_t getLocalLength() const;
    global_size_t getGlobalLength() const;
    size_t getStride() const;
    bool isConstantStride() const;
    bool aliases(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& other) const;
    virtual std::string description() const;
    virtual void removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap);
    void setCopyOrView (const Teuchos::DataAccess copyOrView);
    bool isSameSize(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & vec) const;
    bool importsAreAliased();

  }; // class MultiVector
}
/* End Sandia copyright */

%teuchos_rcp(Tpetra::MultiVector<SC,LO,GO,NO>)
%template(TpetraMultiVector) Tpetra::MultiVector<SC,LO,GO,NO>;

%clear Teuchos::ArrayRCP<double>;
%clear Teuchos::ArrayRCP<const double>;
