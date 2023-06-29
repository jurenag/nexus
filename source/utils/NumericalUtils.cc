#include "NumericalUtils.h"

namespace nexus {

  G4double FindClosestNumber(G4double x, std::vector<G4double>* numbers)
  {

    auto CompareDistance = [x](G4double a, G4double b) -> G4double {  // The comparison function depends on the input value, x.
        return std::abs(x-a) < std::abs(x-b);                   // This lambda function let us capture such x from the context.
    };

    auto closest = std::min_element(numbers->begin(), numbers->end(), CompareDistance); // Returns an iterator to the minimum element
                                                                                      // within numbers, up to the defined comparison
    return *closest;
  }

}
