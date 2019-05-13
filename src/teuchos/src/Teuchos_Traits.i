/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Declare and ignore some traits classes
namespace Teuchos {
  template<typename O, typename T> class SerializationTraits;
  template<typename T> class TypeNameTraits;
}
%ignore Teuchos::SerializationTraits;
%ignore Teuchos::TypeNameTraits;

