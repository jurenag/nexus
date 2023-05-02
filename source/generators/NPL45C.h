
#ifndef NPL45C_H
#define NPL45C_H

#include <G4VPrimaryGenerator.hh>

#include <random>
#include <initializer_list>
#include <vector>

class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;


namespace nexus {

  class GeometryBase;
  // Thorlabs pulsed laser: thorlabs.com/thorproduct.cfm?partnumber=NPL45C
  class NPL45C: public G4VPrimaryGenerator
  {
  public:
    /// Constructor
    NPL45C();
    /// Destructor
    ~NPL45C();

    /// This method is invoked at the beginning of the event. It sets
    /// a primary vertex (that is, a photon in a given position and time)
    /// in the event.
    void GeneratePrimaryVertex(G4Event*);

  private:    
    G4GenericMessenger* msg_;

    const GeometryBase* geom_;      ///< Pointer to the detector geometry
    
    // For further implementation:
    G4double ed_x_, ed_y_, ed_z_;   ///< Emission direction

    G4String region_;               ///< Region label that is given to the geometry to generate a vertex

    std::random_device rd_;     ///< Weak random generator (used to seed mt19937 random gen.)
    std::mt19937 gen_;          ///< Random number generator

    std::vector<G4double> wl_bin_edges_;                          ///< Bin edges for the laser emission spectrum
    std::vector<G4double> emission_spectrum_;                     ///< Laser emission spectrum
    std::piecewise_constant_distribution<G4double>* wl_sampler_;  ///< Laser emission spectrum histogram sampler

    void LoadNCheckSpectrumData();  ///< Loads LED emission spectrum data into wl_bin_edges_ and emission_spectrum_
    G4bool SpectrumDataIsIllFormed();   ///< Checks whether the loaded data is ill-formed
    G4double RandomEnergy();    ///< Generate a random kinetic energy according to LED emission spectrum.
  };

} // end namespace nexus

#endif
