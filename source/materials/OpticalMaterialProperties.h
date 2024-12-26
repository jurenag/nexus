// ----------------------------------------------------------------------------
// nexus | OpticalMaterialProperties.h
//
// Optical properties of relevant materials.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef OPTICAL_MATERIAL_PROPERTIES_H
#define OPTICAL_MATERIAL_PROPERTIES_H

#include <globals.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

class G4MaterialPropertiesTable;


namespace opticalprops {

  using namespace CLHEP;


  G4MaterialPropertiesTable* Vacuum();

  G4MaterialPropertiesTable* Air();

  G4MaterialPropertiesTable* TunableRIMat(G4double rindex=1.5);

  G4MaterialPropertiesTable* GlassEpoxy();

  G4MaterialPropertiesTable* HamamatsuEpoxy();

  G4MaterialPropertiesTable* FusedSilica();

  G4MaterialPropertiesTable* FakeFusedSilica(G4double transparency = .9,
                                             G4double thickness    = 1. * mm);

  G4MaterialPropertiesTable* Epoxy();

  G4MaterialPropertiesTable* ImperfectDielectricDielectricSurface(G4double escaping_probability=.01);

  G4MaterialPropertiesTable* ITO();

  G4MaterialPropertiesTable* PEDOT();

  G4MaterialPropertiesTable* Sapphire();

  G4MaterialPropertiesTable* OptCoupler();

  G4MaterialPropertiesTable* GAr(G4double sc_yield,
                                 G4double e_lifetime=1000.*ms);

  G4MaterialPropertiesTable* GXe(G4double pressure=1.*bar,
                                 G4double temperature=STP_Temperature,
                                 G4int sc_yield=25510/MeV,
                                 G4double e_lifetime=1000.*ms);

  G4MaterialPropertiesTable* LXe();

  G4MaterialPropertiesTable* LAr();

  G4MaterialPropertiesTable* PTP(G4double refractive_index = 1.65);

  G4MaterialPropertiesTable* LArPTPArtifact();
  
  G4MaterialPropertiesTable* PerfectPhotonCollector();

  G4MaterialPropertiesTable* PerfectPhotonReflector();

  G4MaterialPropertiesTable* PerfectPolishedSurfaceTransmitter();

  G4MaterialPropertiesTable* FakeGrid(G4double pressure=1.*bar,
                                      G4double temperature=STP_Temperature,
                                      G4double transparency=.9,
                                      G4double thickness=1.*mm,
                                      G4int sc_yield=25510/MeV,
                                      G4double e_lifetime=1000.*ms,
                                      G4double photoe_p=0);

  G4MaterialPropertiesTable* TPB();

  G4MaterialPropertiesTable* DegradedTPB(G4double wls_eff);

  G4MaterialPropertiesTable* TPH();

  G4MaterialPropertiesTable* PTFE();

  G4MaterialPropertiesTable* PolishedAl();

  G4MaterialPropertiesTable* EJ280();

  G4MaterialPropertiesTable* EJ282();

  G4MaterialPropertiesTable* EJ286(G4double attenuation_length=1.*m);

  G4MaterialPropertiesTable* G2P_FB118(G4double cromophore_concentration, G4double rindex=1.502, G4bool verbosity = true);

  G4MaterialPropertiesTable* FakeG2P_FB118(G4double cromophore_concentration, G4double rindex=1.502, G4bool verbosity = true);

  G4MaterialPropertiesTable* SCHOTT_B270();

  G4MaterialPropertiesTable* SCHOTT_BOROFLOAT_33();

  G4MaterialPropertiesTable* Y11();

  G4MaterialPropertiesTable* B2();

  G4MaterialPropertiesTable* Pethylene();

  G4MaterialPropertiesTable* FPethylene();

  G4MaterialPropertiesTable* PMMA();

  G4MaterialPropertiesTable* Copper();

  G4MaterialPropertiesTable* Steel();

  G4MaterialPropertiesTable* XXX();

  G4MaterialPropertiesTable* Vikuiti(G4int reflectivity_type=0, G4double reflectivity_factor_scale=1.);


  constexpr G4double optPhotMinE_ =  0.2  * eV;
  constexpr G4double optPhotMaxE_ = 11.5  * eV;
  constexpr G4double optPhotFusedSilicaMaxE_ = 10.7  * eV; // formulas for fused silica are valid up to this energy
  constexpr G4double optPhotSapphireMaxE_ = 10.3  * eV; // formulas for sapphire are valid up to this energy
  constexpr G4double noAbsLength_ = 1.e8  * m;

  constexpr G4double hc_ = h_Planck * c_light;


} // end namespace opticalprops

#endif
