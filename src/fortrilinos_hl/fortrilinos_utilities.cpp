#include "fortrilinos_utilities.hpp"

#include <Teuchos_ScalarTraits.hpp>

#include <climits>

namespace ForTrilinos {

  // Taken from MueLu
  void SetRandomSeed(const Teuchos::Comm<int> &comm) {
    // Distribute the seeds evenly in [1,maxint-1].  This guarantees nothing
    // about where in random number stream we are, but avoids overflow situations
    // in parallel when multiplying by a PID.  It would be better to use
    // a good parallel random number generator.
    double one = 1.0;
    int maxint = INT_MAX; //= 2^31-1 = 2147483647 for 32-bit integers
    int mySeed = Teuchos::as<int>((maxint-1) * (one -(comm.getRank()+1)/(comm.getSize()+one)) );
    if (mySeed < 1 || mySeed == maxint) {
      std::ostringstream errStr;
      errStr << "Error detected with random seed = " << mySeed << ". It should be in the interval [1, 2^31-2].";
      throw std::runtime_error(errStr.str());
    }

    std::srand(mySeed);
    // For Tpetra, we could use Kokkos' random number generator here.
    Teuchos::ScalarTraits<double>::seedrandom(mySeed);
  }

}
