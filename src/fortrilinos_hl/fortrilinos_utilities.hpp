/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#ifndef FORTRILINOS_UTILITIES_HPP
#define FORTRILINOS_UTILITIES_HPP

#include <Teuchos_Comm.hpp>

namespace ForTrilinos {

  void SetRandomSeed(const Teuchos::Comm<int> &comm);

}

#endif // FORTRILINOS_UTILITIES_HPP
