#ifndef LAR_SPHERE_H
#define LAR_SPHERE_H

#include "GeometryBase.h"

class G4Material;
class G4GenericMessenger;
namespace nexus { class SpherePointSampler; }

namespace nexus {

  /// Spherical chamber filled with liquid argon

  class LArSphere: public GeometryBase
  {
  public:
    /// Default constructor
    LArSphere();
    /// Alternative constructor
    LArSphere(G4double);
    /// Destructor
    ~LArSphere();

    /// Return vertex within region <region> of the chamber
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

  private:
    G4double radius_;   ///< Radius of the sphere

    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    /// Vertexes random generator
    SpherePointSampler* sphere_vertex_gen_;

  };

} // end namespace nexus

#endif
