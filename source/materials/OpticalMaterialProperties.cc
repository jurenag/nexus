// ----------------------------------------------------------------------------
// nexus | OpticalMaterialProperties.cc
//
// Optical properties of relevant materials.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "OpticalMaterialProperties.h"
#include "XenonProperties.h"
#include "SellmeierEquation.h"

#include <G4MaterialPropertiesTable.hh>

#include <assert.h>

using namespace nexus;
using namespace CLHEP;

namespace opticalprops {
  /// Vacuum ///
  G4MaterialPropertiesTable* Vacuum()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> photEnergy = {optPhotMinE_, optPhotMaxE_};

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex = {1., 1.};
    mpt->AddProperty("RINDEX", photEnergy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> absLength = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", photEnergy, absLength);

    return mpt;
  }



  /// Fused Silica ///
  G4MaterialPropertiesTable* FusedSilica()
  {
    // Optical properties of Suprasil 311/312(c) synthetic fused silica.
    // Obtained from http://heraeus-quarzglas.com

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    // The range is chosen to be up to ~10.7 eV because Sellmeier's equation
    // for fused silica is valid only in that range
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    // The following values for the refractive index have been calculated
    // using Sellmeier's equation:
    //    n^2 - 1 = B_1 * \lambda^2 / (\lambda^2 - C_1) +
    //            + B_2 * \lambda^2 / (\lambda^2 - C_2) +
    //            + B_3 * \lambda^2 / (\lambda^2 - C_3),
    // with wavelength \lambda in micrometers and
    //    B_1 = 4.73E-1, B_2 = 6.31E-1, B_3 = 9.06E-1
    //    C_1 = 1.30E-2, C_2 = 4.13E-3, C_3 = 9.88E+1.

    G4double B_1 = 4.73e-1;
    G4double B_2 = 6.31e-1;
    G4double B_3 = 9.06e-1;
    G4double C_1 = 1.30e-2;
    G4double C_2 = 4.13e-3;
    G4double C_3 = 9.88e+1;

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double lambda = h_Planck*c_light/ri_energy[i]*1000; // in micron
      G4double n2 = 1 + B_1*pow(lambda,2)/(pow(lambda,2)-C_1)
        + B_2*pow(lambda,2)/(pow(lambda,2)-C_2)
        + B_3*pow(lambda,2)/(pow(lambda,2)-C_3);
      rIndex.push_back(sqrt(n2));
      // G4cout << "* FusedSilica rIndex:  " << std::setw(5) << ri_energy[i]/eV
      //       << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,  6.46499 * eV,
      6.54000 * eV,  6.59490 * eV,  6.64000 * eV,  6.72714 * eV,
      6.73828 * eV,  6.75000 * eV,  6.82104 * eV,  6.86000 * eV,
      6.88000 * eV,  6.89000 * eV,  7.00000 * eV,  7.01000 * eV,
      7.01797 * eV,  7.05000 * eV,  7.08000 * eV,  7.08482 * eV,
      7.30000 * eV,  7.36000 * eV,  7.40000 * eV,  7.48000 * eV,
      7.52000 * eV,  7.58000 * eV,  7.67440 * eV,  7.76000 * eV,
      7.89000 * eV,  7.93000 * eV,  8.00000 * eV,
      optPhotMaxE_
    };

    std::vector<G4double> absLength = {
      noAbsLength_, noAbsLength_,
      200.0 * cm,   200.0 * cm,  90.0 * cm,  45.0 * cm,
      45.0 * cm,    30.0 * cm,  24.0 * cm,  21.0 * cm,
      20.0 * cm,    19.0 * cm,  16.0 * cm,  14.0 * cm,
      13.0 * cm,     8.5 * cm,   8.0 * cm,   6.0 * cm,
       1.5 * cm,     1.2 * cm,   1.0 * cm,   .65 * cm,
        .4 * cm,     .37 * cm,   .32 * cm,   .28 * cm,
        .22 * cm,    .215 * cm,  .00005*cm,
      .00005* cm
    };

    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }



  /// Fake Fused Silica ///
  G4MaterialPropertiesTable* FakeFusedSilica(G4double transparency,
                                            G4double thickness)
  {
    // Optical properties of Suprasil 311/312(c) synthetic fused silica.
    // Obtained from http://heraeus-quarzglas.com

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    // The range is chosen to be up to ~10.7 eV because Sellmeier's equation
    // for fused silica is valid only in that range
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    // The following values for the refractive index have been calculated
    // using Sellmeier's equation:
    //    n^2 - 1 = B_1 * \lambda^2 / (\lambda^2 - C_1) +
    //            + B_2 * \lambda^2 / (\lambda^2 - C_2) +
    //            + B_3 * \lambda^2 / (\lambda^2 - C_3),
    // with wavelength \lambda in micrometers and
    //    B_1 = 4.73E-1, B_2 = 6.31E-1, B_3 = 9.06E-1
    //    C_1 = 1.30E-2, C_2 = 4.13E-3, C_3 = 9.88E+1.

    G4double B_1 = 4.73e-1;
    G4double B_2 = 6.31e-1;
    G4double B_3 = 9.06e-1;
    G4double C_1 = 1.30e-2;
    G4double C_2 = 4.13e-3;
    G4double C_3 = 9.88e+1;

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double lambda = h_Planck*c_light/ri_energy[i]*1000; // in micron
      G4double n2 = 1 + B_1*pow(lambda,2)/(pow(lambda,2)-C_1)
        + B_2*pow(lambda,2)/(pow(lambda,2)-C_2)
        + B_3*pow(lambda,2)/(pow(lambda,2)-C_3);
      rIndex.push_back(sqrt(n2));
      //G4cout << "* FakeFusedSilica rIndex:  " << std::setw(5)
      //       << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH (Set to match the transparency)
    G4double abs_length     = -thickness / log(transparency);
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> abs_l      = {abs_length, abs_length};
    mpt->AddProperty("ABSLENGTH", abs_energy, abs_l);

    return mpt;
  }



  /// ITO ///
  G4MaterialPropertiesTable* ITO()
  {
    // Input data: complex refraction index obtained from:
    // https://refractiveindex.info/?shelf=other&book=In2O3-SnO2&page=Moerland
    // Only valid in [1000 - 400] nm

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> energies = {
      optPhotMinE_,
      h_Planck * c_light / (1000. * nm),  h_Planck * c_light / (800. * nm),
      h_Planck * c_light / ( 700. * nm),  h_Planck * c_light / (600. * nm),
      h_Planck * c_light / ( 580. * nm),  h_Planck * c_light / (560. * nm),
      h_Planck * c_light / ( 540. * nm),  h_Planck * c_light / (520. * nm),
      h_Planck * c_light / ( 500. * nm),  h_Planck * c_light / (480. * nm),
      h_Planck * c_light / ( 460. * nm),  h_Planck * c_light / (440. * nm),
      h_Planck * c_light / ( 420. * nm),  h_Planck * c_light / (400. * nm),
      optPhotMaxE_ };

    std::vector<G4double> rIndex = {
      1.635,
      1.635, 1.775,
      1.835, 1.894,
      1.906, 1.919,
      1.931, 1.945,
      1.960, 1.975,
      1.993, 2.012,
      2.036, 2.064,
      2.064 };
    mpt->AddProperty("RINDEX", energies, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_length = {
      (1000. * nm) / (4*pi * 0.0103),
      (1000. * nm) / (4*pi * 0.0103),  (800. * nm) / (4*pi * 0.0049),
      ( 700. * nm) / (4*pi * 0.0033),  (600. * nm) / (4*pi * 0.0023),
      ( 580. * nm) / (4*pi * 0.0022),  (560. * nm) / (4*pi * 0.0022),
      ( 540. * nm) / (4*pi * 0.0022),  (520. * nm) / (4*pi * 0.0023),
      ( 500. * nm) / (4*pi * 0.0026),  (480. * nm) / (4*pi * 0.0031),
      ( 460. * nm) / (4*pi * 0.0039),  (440. * nm) / (4*pi * 0.0053),
      ( 420. * nm) / (4*pi * 0.0080),  (400. * nm) / (4*pi * 0.0125),
      (400. * nm) / (4*pi * 0.0125) };
    mpt->AddProperty("ABSLENGTH", energies, abs_length);

    //G4cout << "*** ITO properties ...  " << G4endl;
    //for (int i=0; i<num_energies; i++) {
    //  G4cout << "* Energy: " << std::setw(5) << energies[i]/eV << " eV  ->  " << std::setw(5)
    //         << (1240./ (energies[i] / eV)) << " nm" << G4endl;
    //  G4cout << "  rIndex    : " << std::setw(5) << rIndex[i]
    //         << "  Abs Length: " << std::setw(5) << abs_length[i] / nm << " nm" << G4endl;
    //}

    return mpt;
  }



  /// PEDOT ///
  G4MaterialPropertiesTable* PEDOT()
  {
    // Input data: complex refraction index obtained from:
    // https://refractiveindex.info/?shelf=other&book=PEDOT-PSS&page=Chen
    // Only valid in [1097 - 302] nm

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> energies = {
      optPhotMinE_,
      h_Planck * c_light / (1097. * nm),  h_Planck * c_light / (1000. * nm),
      h_Planck * c_light / ( 950. * nm),  h_Planck * c_light / ( 900. * nm),
      h_Planck * c_light / ( 800. * nm),  h_Planck * c_light / ( 700. * nm),
      h_Planck * c_light / ( 600. * nm),  h_Planck * c_light / ( 550. * nm),
      h_Planck * c_light / ( 500. * nm),  h_Planck * c_light / ( 450. * nm),
      h_Planck * c_light / ( 420. * nm),  h_Planck * c_light / ( 400. * nm),
      h_Planck * c_light / ( 370. * nm),  h_Planck * c_light / ( 350. * nm),
      h_Planck * c_light / ( 302. * nm),  optPhotMaxE_ };

    std::vector<G4double> rIndex = {
      1.4760,
      1.4760, 1.4662,
      1.4665, 1.4693,
      1.4802, 1.4935,
      1.5080, 1.5155,
      1.5235, 1.5328,
      1.5391, 1.5439,
      1.5522, 1.5587,
      1.5805, 1.5805 };
    mpt->AddProperty("RINDEX", energies, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_length = {
      (1097. * nm) / (4*pi * 0.1191),
      (1097. * nm) / (4*pi * 0.1191),   (1000. * nm) / (4*pi * 0.0859),
      ( 950. * nm) / (4*pi * 0.0701),   ( 900. * nm) / (4*pi * 0.0561),
      ( 800. * nm) / (4*pi * 0.0340),   ( 700. * nm) / (4*pi * 0.0197),
      ( 600. * nm) / (4*pi * 0.0107),   ( 550. * nm) / (4*pi * 0.0076),
      ( 500. * nm) / (4*pi * 0.0051),   ( 450. * nm) / (4*pi * 0.0035),
      ( 420. * nm) / (4*pi * 0.0025),   ( 400. * nm) / (4*pi * 0.00194),
      ( 370. * nm) / (4*pi * 0.00135),  ( 350. * nm) / (4*pi * 0.00103),
      ( 302. * nm) / (4*pi * 0.0004),   ( 302. * nm) / (4*pi * 0.0004) };
    mpt->AddProperty("ABSLENGTH", energies, abs_length);

    //G4cout << "*** PEDOT properties ...  " << G4endl;
    //for (int i=0; i<num_energies; i++) {
    //  G4cout << "* Energy: " << std::setw(5) << energies[i]/eV << " eV  ->  " << std::setw(5)
    //         << (1240./ (energies[i] / eV)) << " nm" << G4endl;
    //  G4cout << "  rIndex    : " << std::setw(5) << rIndex[i]
    //         << "  Abs Length: " << std::setw(5) << abs_length[i] / nm << " nm" << G4endl;
    //}

    return mpt;
  }



  /// Glass Epoxy ///
  G4MaterialPropertiesTable* GlassEpoxy()
  {
    // Optical properties of Optorez 1330 glass epoxy.
    // Obtained from http://refractiveindex.info and
    // https://www.zeonex.com/Optics.aspx.html#glass-like

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    // The range is chosen to be up to ~10.7 eV because Sellmeier's equation
    // for fused silica is valid only in that range
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double lambda = h_Planck*c_light/ri_energy[i]*1000; // in micron
      G4double n2 = 2.291142 - 3.311944E-2*pow(lambda,2) - 1.630099E-2*pow(lambda,-2) +
                    7.265983E-3*pow(lambda,-4) - 6.806145E-4*pow(lambda,-6) +
                    1.960732E-5*pow(lambda,-8);
      rIndex.push_back(sqrt(n2));
      // G4cout << "* GlassEpoxy rIndex:  " << std::setw(5)
      //        << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_, 2.000 * eV,
      2.132 * eV,   2.735 * eV,  2.908 * eV,  3.119 * eV,
      3.320 * eV,   3.476 * eV,  3.588 * eV,  3.749 * eV,
      3.869 * eV,   3.973 * eV,  4.120 * eV,  4.224 * eV,
      4.320 * eV,   4.420 * eV,  5.018 * eV
    };

    std::vector<G4double> absLength = {
      noAbsLength_, noAbsLength_,
      326.00 * mm,  117.68 * mm,  85.89 * mm,  50.93 * mm,
      31.25 * mm,   17.19 * mm,  10.46 * mm,   5.26 * mm,
        3.77 * mm,    2.69 * mm,   1.94 * mm,   1.33 * mm,
        0.73 * mm,    0.32 * mm,   0.10 * mm
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }

  /// Epoxy resin of Hamamatsu S13360 MPPC's ///
  G4MaterialPropertiesTable* HamamatsuEpoxy()
  {
    // See hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s13360-2050ve_etc_kapd1053e.pdf
    // There's very little data regarding this MPPC window. I am gonna stick to the available data.
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    G4double energy_span[] =  {optPhotMinE_,  optPhotMaxE_};
    G4double rindex[] =       {1.55        ,  1.55        };
    G4double abslength[] =    {noAbsLength_,  noAbsLength_};

    mpt->AddProperty("RINDEX",    energy_span, rindex,    2);
    mpt->AddProperty("ABSLENGTH", energy_span, abslength, 2);

    return mpt;
  }


  /// Sapphire ///
  G4MaterialPropertiesTable* Sapphire()
  {
    // Input data: Sellmeier equation coeficients extracted from:
    // https://refractiveindex.info/?shelf=3d&book=crystals&page=sapphire
    //C[i] coeficients at line 362 are squared.

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    G4double um2 = micrometer*micrometer;
    G4double B[3] = {1.4313493, 0.65054713, 5.3414021};
    G4double C[3] = {0.0052799261 * um2, 0.0142382647 * um2, 325.017834 * um2};
    SellmeierEquation seq(B, C);

    const G4int ri_entries = 100;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      rIndex.push_back(seq.RefractiveIndex(h_Planck*c_light/ri_energy[i]));
      //G4cout << "* Sapphire rIndex:  " << std::setw(5)
      //       << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_, 0.900 * eV,
      1.000 * eV,   1.296 * eV,  1.683 * eV,  2.075 * eV,
      2.585 * eV,   3.088 * eV,  3.709 * eV,  4.385 * eV,
      4.972 * eV,   5.608 * eV,  6.066 * eV,  6.426 * eV,
      6.806 * eV,   7.135 * eV,  7.401 * eV,  7.637 * eV,
      7.880 * eV,   8.217 * eV
    };

    std::vector<G4double> absLength = {
      noAbsLength_, noAbsLength_,
      3455.0  * mm,  3455.0  * mm,  3455.0  * mm,  3455.0  * mm,
      3455.0  * mm,  3140.98 * mm,  2283.30 * mm,  1742.11 * mm,
      437.06 * mm,   219.24 * mm,  117.773 * mm,   80.560 * mm,
      48.071 * mm,   28.805 * mm,   17.880 * mm,   11.567 * mm,
        7.718 * mm,    4.995 * mm
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }



  /// Optical Coupler ///
  G4MaterialPropertiesTable* OptCoupler()
  {
    // gel NyoGel OCK-451
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    G4double um2 = micrometer*micrometer;

    G4double constTerm  = 1.4954;
    G4double squareTerm = 0.008022 * um2;
    G4double quadTerm   = 0.;

    const G4int ri_entries = 100;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double wl = h_Planck * c_light / ri_energy[i];
      rIndex.push_back(constTerm + squareTerm/(wl*wl) + quadTerm/pow(wl,4));
      //G4cout << "* OptCoupler rIndex:  " << std::setw(5)
      //       << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    // Values estimated from printed plot (to be improved).
    std::vector<G4double> abs_energy = {
      optPhotMinE_,  1.70 * eV,
      1.77 * eV,     2.07 * eV,  2.48 * eV,  2.76 * eV,
      2.92 * eV,     3.10 * eV,  3.31 * eV,  3.54 * eV,
      3.81 * eV,     4.13 * eV
    };

    std::vector<G4double> absLength = {
      noAbsLength_, noAbsLength_,
      1332.8 * mm,  1332.8 * mm,  1332.8 * mm,  666.17 * mm,
      499.5 * mm,   399.5 * mm,   199.5 * mm,  132.83 * mm,
        99.5 * mm,     4.5 * mm
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }



  /// Gaseous Argon ///
  G4MaterialPropertiesTable* GAr(G4double sc_yield,
                                G4double e_lifetime)
  {
    // An argon gas proportional scintillation counter with UV avalanche photodiode scintillation
    // readout C.M.B. Monteiro, J.A.M. Lopes, P.C.P.S. Simoes, J.M.F. dos Santos, C.A.N. Conde

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double wl = h_Planck * c_light / ri_energy[i] * 1000; // in micron
      // From refractiveindex.info
      rIndex.push_back(1 + 0.012055*(0.2075*pow(wl,2)/(91.012*pow(wl,2)-1) +
                                     0.0415*pow(wl,2)/(87.892*pow(wl,2)-1) +
                                     4.3330*pow(wl,2)/(214.02*pow(wl,2)-1)));
      //G4cout << "* GAr rIndex:  " << std::setw(5) << ri_energy[i]/eV
      //       << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // EMISSION SPECTRUM
    G4double Wavelength_peak  = 128.000 * nm;
    G4double Wavelength_sigma =   2.929 * nm;
    G4double Energy_peak  = (h_Planck*c_light / Wavelength_peak);
    G4double Energy_sigma = (h_Planck*c_light * Wavelength_sigma / pow(Wavelength_peak,2));
    //G4cout << "*** GAr Energy_peak: " << Energy_peak/eV << " eV   Energy_sigma: "
    //       << Energy_sigma/eV << " eV" << G4endl;

    // Sampling from ~110 nm to 150 nm <----> from ~11.236 eV to 8.240 eV
    const G4int sc_entries = 380;
    std::vector<G4double> sc_energy;
    std::vector<G4double> intensity;
    for (int i=0; i<sc_entries; i++){
      sc_energy.push_back(8.240*eV + 0.008*i*eV);
      intensity.push_back(exp(-pow(Energy_peak/eV-sc_energy[i]/eV,2) /
                              (2*pow(Energy_sigma/eV, 2)))/(Energy_sigma/eV*sqrt(pi*2.)));
      //G4cout << "* GAr energy: " << std::setw(6) << sc_energy[i]/eV << " eV  ->  "
      //       << std::setw(6) << intensity[i] << G4endl;
    }
    mpt->AddProperty("SCINTILLATIONCOMPONENT1", sc_energy, intensity);
    mpt->AddProperty("SCINTILLATIONCOMPONENT2", sc_energy, intensity);
    mpt->AddProperty("ELSPECTRUM"             , sc_energy, intensity, 1);

    // CONST PROPERTIES
    mpt->AddConstProperty("SCINTILLATIONYIELD", sc_yield);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",   6.*ns);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",   37.*ns);
    mpt->AddConstProperty("SCINTILLATIONYIELD1", .342);
    mpt->AddConstProperty("SCINTILLATIONYIELD2", .658);
    mpt->AddConstProperty("RESOLUTIONSCALE",    1.0);
    mpt->AddConstProperty("ATTACHMENT",         e_lifetime, 1);

    return mpt;
  }



  /// Gaseous Xenon ///
  G4MaterialPropertiesTable* GXe(G4double pressure,
                                G4double temperature,
                                G4int    sc_yield,
                                G4double e_lifetime)
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    G4double density = GXeDensity(pressure);
    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      rIndex.push_back(XenonRefractiveIndex(ri_energy[i], density));
      // G4cout << "* GXe rIndex:  " << std::setw(7)
      //        << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex, ri_entries);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // EMISSION SPECTRUM
    // Sampling from ~150 nm to 200 nm <----> from 6.20625 eV to 8.20625 eV
    const G4int sc_entries = 200;
    std::vector<G4double> sc_energy;
    for (int i=0; i<sc_entries; i++){
      sc_energy.push_back(6.20625 * eV + 0.01 * i * eV);
    }
    std::vector<G4double> intensity;
    for (G4int i=0; i<sc_entries; i++) {
      intensity.push_back(GXeScintillation(sc_energy[i], pressure));
    }
    //for (int i=0; i<sc_entries; i++) {
    //  G4cout << "* GXe Scint:  " << std::setw(7) << sc_energy[i]/eV
    //         << " eV -> " << intensity[i] << G4endl;
    //}
    mpt->AddProperty("SCINTILLATIONCOMPONENT1", sc_energy, intensity);
    mpt->AddProperty("SCINTILLATIONCOMPONENT2", sc_energy, intensity);
    mpt->AddProperty("ELSPECTRUM"             , sc_energy, intensity, 1);

    // CONST PROPERTIES
    mpt->AddConstProperty("SCINTILLATIONYIELD", sc_yield);
    mpt->AddConstProperty("RESOLUTIONSCALE",    1.0);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",   4.5  * ns);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",   100. * ns);
    mpt->AddConstProperty("SCINTILLATIONYIELD1", .1);
    mpt->AddConstProperty("SCINTILLATIONYIELD2", .9);
    mpt->AddConstProperty("ATTACHMENT",         e_lifetime, 1);

    return mpt;
  }

  G4MaterialPropertiesTable* LXe()
  {
    /// The time constants are taken from E. Hogenbirk et al 2018 JINST 13 P10031
    G4MaterialPropertiesTable* LXe_mpt = new G4MaterialPropertiesTable();

    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    G4double density = LXeDensity();
    std::vector<G4double> ri_index;
    for (G4int i=0; i<ri_entries; i++) {
      ri_index.push_back(XenonRefractiveIndex(ri_energy[i], density));
    }
    LXe_mpt->AddProperty("RINDEX", ri_energy, ri_index);

    // for (G4int i=ri_entries-1; i>=0; i--) {
    //   G4cout << h_Planck*c_light/ri_energy[i]/nanometer << " nm, " << rindex[i] << G4endl;
    // }

    // Sampling from ~151 nm to 200 nm <----> from 6.20625 eV to 8.21 eV
    const G4int sc_entries = 500;
    const G4double minE = 6.20625*eV;
    eWidth = (optPhotMaxE_ - minE) / sc_entries;

    std::vector<G4double> sc_energy;
    for (G4int j=0; j<sc_entries; j++){
      sc_energy.push_back(minE + j * eWidth);
    }
    std::vector<G4double> intensity;
    for (G4int i=0; i<sc_entries; i++) {
      intensity.push_back(LXeScintillation(sc_energy[i]));
    }

    LXe_mpt->AddProperty("SCINTILLATIONCOMPONENT1", sc_energy, intensity);
    LXe_mpt->AddProperty("SCINTILLATIONCOMPONENT2", sc_energy, intensity);

    LXe_mpt->AddConstProperty("SCINTILLATIONYIELD", 58708./MeV);
    LXe_mpt->AddConstProperty("RESOLUTIONSCALE", 1);
    LXe_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.*ns);
    LXe_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 43.5*ns);
    LXe_mpt->AddConstProperty("SCINTILLATIONYIELD1", .03);
    LXe_mpt->AddConstProperty("SCINTILLATIONYIELD2", .97);
    LXe_mpt->AddConstProperty("ATTACHMENT", 1000.*ms, 1);

    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> abs_length = {noAbsLength_, noAbsLength_};

    LXe_mpt->AddProperty("ABSLENGTH", abs_energy, abs_length);

    std::vector<G4double> rayleigh_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> rayleigh_length = {36.*cm, 36.*cm};

    LXe_mpt->AddProperty("RAYLEIGH", rayleigh_energy, rayleigh_length);

    return LXe_mpt;
  }

  G4MaterialPropertiesTable* LAr()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    G4double energies_rindex[] = {2.043*eV, 2.2109*eV, 2.3789*eV, 2.5468*eV, 2.7148*eV, 2.8828*eV, 3.0507*eV, 
                          3.2187*eV, 3.3867*eV, 3.5546*eV, 3.7226*eV, 3.8905*eV, 4.0585*eV, 4.2265*eV, 
                          4.3944*eV, 4.5624*eV, 4.7303*eV, 4.8983*eV, 5.0663*eV, 5.2342*eV, 5.4022*eV, 
                          5.5702*eV, 5.7381*eV, 5.9061*eV, 6.074*eV, 6.242*eV, 6.41*eV, 6.5779*eV, 
                          6.7459*eV, 6.9139*eV, 7.0818*eV, 7.2498*eV, 7.4177*eV, 7.5857*eV, 7.7537*eV, 
                          7.9216*eV, 8.0896*eV, 8.2575*eV, 8.4255*eV, 8.5935*eV, 8.7614*eV, 8.9294*eV, 
                          9.0974*eV, 9.2577*eV, 9.3951*eV, 9.5173*eV, 9.6471*eV, 9.7616*eV, 9.8532*eV, 
                          9.9448*eV, 10.029*eV, 10.097*eV, 10.159*eV, 10.235*eV, 10.288*eV, 10.328*eV, 
                          10.383*eV, 10.461*eV, 10.475*eV, 10.502*eV};

    // Before considering any change on this LAr rindex, please take into account that LArPTPArtifact
    // is meant to have the exact same rindex as LAr. So, please, introduce the same changes in both
    // mpts.
    G4double rindex[] = {1.2288, 1.233, 1.2332, 1.2344, 1.235, 1.2362, 1.2373, 1.2379, 1.2393, 1.2405, 
                        1.2414, 1.2432, 1.2441, 1.2452, 1.2471, 1.2495, 1.2511, 1.2527, 1.2548, 1.2571, 
                        1.2592, 1.2614, 1.2641, 1.2667, 1.2694, 1.2725, 1.2758, 1.2794, 1.2835, 1.2877, 
                        1.2928, 1.2979, 1.3031, 1.3085, 1.3146, 1.3218, 1.3291, 1.337, 1.3454, 1.3556, 
                        1.3673, 1.3796, 1.3921, 1.4071, 1.4229, 1.4381, 1.4535, 1.4691, 1.485, 1.4999, 
                        1.5149, 1.5332, 1.5506, 1.5672, 1.5859, 1.5997, 1.6171, 1.6301, 1.6464, 1.6623};

    mpt->AddProperty("RINDEX", energies_rindex, rindex, 60);

    // Absorption length (ABSLENGTH)
    G4double energies_abslength[]  = {optPhotMinE_, optPhotMaxE_};
      
    G4double abslength[] = {noAbsLength_, noAbsLength_};

    mpt->AddProperty("ABSLENGTH", energies_abslength, abslength, 2);

    return mpt;
  }

  G4MaterialPropertiesTable* PTP(){
    // This material is meant to model p-Terphenyl, which is a wavelength shifter commonly used to WLS the 
    // LAr VUV scintillation light. The emission spectrum was taken from
    // sciencedirect.com/science/article/abs/pii/016890029390701I ,
    // whereas the PLQY (photo luminiscence quantum yield = #photons emitted/#photons absorbed), the rindex
    // and the WLS delay time was taken from 
    // mdpi-res.com/d_attachment/instruments/instruments-05-00004/article_deploy/instruments-05-00004-v2.pdf?version=1609810328
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // Refractive index (RINDEX)
    G4double energies_rindex[] =    {optPhotMinE_, optPhotMaxE_};
    G4double rindex[] =             {1.65, 1.65};
    mpt->AddProperty("RINDEX", energies_rindex, rindex, 2);

    // Absorption length (ABSLENGTH)
    G4double energies_abslength[]  = {optPhotMinE_, optPhotMaxE_};
    G4double abslength[] = {noAbsLength_, noAbsLength_};

    mpt->AddProperty("ABSLENGTH", energies_abslength, abslength, 2);

    // WLS ABSORPTION LENGTH
    // Getting the LAr scintillation spectrum from researchgate.net/figure/Scintillation-light-spectrum-of-LAr-recorded-with-sulfur-beam-excitation-blue-in-color_fig7_258169905
    // , it is clear that there's almost no LAr scintillation for wavelengths bigger than 145nm. The following WLS absorption length is set so as to absorb every LAr scintillation
    // photon in the range [118, 145]nm, and absorb nothing outside that range.
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_, 
      h_Planck*c_light/(146.25*nm), h_Planck*c_light/(145.75*nm), // These points are meant to fix wlsabslength=infinite at 146nm
      h_Planck*c_light/(145.25*nm), h_Planck*c_light/(144.75*nm), // These points are meant to fix wlsabslength=very small at 145nm
      h_Planck*c_light/(118.25*nm), h_Planck*c_light/(117.75*nm), // These points are meant to fix wlsabslength=very small at 118nm
      h_Planck*c_light/(117.25*nm), h_Planck*c_light/(116.75*nm), // These points are meant to fix wlsabslength=infinite at 117nm   
      optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_  ,
      noAbsLength_  , noAbsLength_  ,
      1.*nm         , 1.*nm         ,
      1.*nm         , 1.*nm         ,
      noAbsLength_  , noAbsLength_  ,
      noAbsLength_  
    };
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_, 
      h_Planck*c_light/(404.*nm), h_Planck*c_light/(403.*nm), h_Planck*c_light/(402.*nm),
      h_Planck*c_light/(401.7764*nm), h_Planck*c_light/(399.4336*nm), h_Planck*c_light/(396.6098*nm), h_Planck*c_light/(395.2065*nm), h_Planck*c_light/(393.0391*nm), 
      h_Planck*c_light/(391.6609*nm), h_Planck*c_light/(390.1205*nm), h_Planck*c_light/(389.0432*nm), h_Planck*c_light/(388.1055*nm), h_Planck*c_light/(387.257*nm), 
      h_Planck*c_light/(386.2195*nm), h_Planck*c_light/(385.8109*nm), h_Planck*c_light/(384.8051*nm), h_Planck*c_light/(383.7688*nm), h_Planck*c_light/(382.7972*nm), 
      h_Planck*c_light/(382.0305*nm), h_Planck*c_light/(380.9974*nm), h_Planck*c_light/(379.9932*nm), h_Planck*c_light/(379.087*nm), h_Planck*c_light/(378.0467*nm), 
      h_Planck*c_light/(377.0923*nm), h_Planck*c_light/(376.0287*nm), h_Planck*c_light/(375.1186*nm), h_Planck*c_light/(374.9371*nm), h_Planck*c_light/(374.2467*nm), 
      h_Planck*c_light/(373.9871*nm), h_Planck*c_light/(373.5364*nm), h_Planck*c_light/(372.9633*nm), h_Planck*c_light/(372.7951*nm), h_Planck*c_light/(372.5039*nm), 
      h_Planck*c_light/(372.1126*nm), h_Planck*c_light/(371.8224*nm), h_Planck*c_light/(371.4771*nm), h_Planck*c_light/(371.1768*nm), h_Planck*c_light/(370.9658*nm), 
      h_Planck*c_light/(370.1462*nm), h_Planck*c_light/(369.9143*nm), h_Planck*c_light/(369.4954*nm), h_Planck*c_light/(368.9566*nm), h_Planck*c_light/(368.3537*nm), 
      h_Planck*c_light/(367.4803*nm), h_Planck*c_light/(366.0482*nm), h_Planck*c_light/(363.8461*nm), h_Planck*c_light/(362.5375*nm), h_Planck*c_light/(361.607*nm), 
      h_Planck*c_light/(360.7652*nm), h_Planck*c_light/(360.2621*nm), h_Planck*c_light/(359.6768*nm), h_Planck*c_light/(359.364*nm), h_Planck*c_light/(359.1766*nm), 
      h_Planck*c_light/(358.9998*nm), h_Planck*c_light/(358.8232*nm), h_Planck*c_light/(358.6468*nm), h_Planck*c_light/(358.4601*nm), h_Planck*c_light/(358.284*nm), 
      h_Planck*c_light/(358.1081*nm), h_Planck*c_light/(357.6845*nm), h_Planck*c_light/(357.2105*nm), h_Planck*c_light/(357.0665*nm), h_Planck*c_light/(356.1536*nm), 
      h_Planck*c_light/(354.5748*nm), h_Planck*c_light/(353.2414*nm), h_Planck*c_light/(353.0*nm), h_Planck*c_light/(352.5283*nm), h_Planck*c_light/(352.398*nm), 
      h_Planck*c_light/(351.9379*nm), h_Planck*c_light/(351.8779*nm), h_Planck*c_light/(351.6284*nm), h_Planck*c_light/(351.2499*nm), h_Planck*c_light/(351.1007*nm), 
      h_Planck*c_light/(350.6142*nm), h_Planck*c_light/(350.5052*nm), h_Planck*c_light/(350.0203*nm), h_Planck*c_light/(349.2217*nm), h_Planck*c_light/(348.2604*nm), 
      h_Planck*c_light/(347.1974*nm), h_Planck*c_light/(346.6925*nm), h_Planck*c_light/(346.5278*nm), h_Planck*c_light/(345.6487*nm), h_Planck*c_light/(344.7644*nm), 
      h_Planck*c_light/(344.0374*nm), h_Planck*c_light/(343.5226*nm), h_Planck*c_light/(343.0569*nm), h_Planck*c_light/(342.924*nm), h_Planck*c_light/(342.5829*nm), 
      h_Planck*c_light/(342.2519*nm), h_Planck*c_light/(341.9122*nm), h_Planck*c_light/(341.8744*nm), h_Planck*c_light/(341.8085*nm), h_Planck*c_light/(341.7519*nm), 
      h_Planck*c_light/(341.6672*nm), h_Planck*c_light/(341.338*nm), h_Planck*c_light/(341.1971*nm), h_Planck*c_light/(341.0188*nm), h_Planck*c_light/(340.9437*nm), 
      h_Planck*c_light/(340.5878*nm), h_Planck*c_light/(340.5317*nm), h_Planck*c_light/(340.298*nm), h_Planck*c_light/(340.214*nm), h_Planck*c_light/(340.046*nm), 
      h_Planck*c_light/(339.9714*nm), h_Planck*c_light/(339.5618*nm), h_Planck*c_light/(339.4967*nm), h_Planck*c_light/(339.1623*nm), h_Planck*c_light/(338.9213*nm), 
      h_Planck*c_light/(338.8472*nm), h_Planck*c_light/(338.6066*nm), h_Planck*c_light/(338.4495*nm), h_Planck*c_light/(338.2832*nm), h_Planck*c_light/(338.1264*nm), 
      h_Planck*c_light/(337.7856*nm), h_Planck*c_light/(337.6476*nm), h_Planck*c_light/(337.4913*nm), h_Planck*c_light/(337.2435*nm), h_Planck*c_light/(336.8586*nm), 
      h_Planck*c_light/(336.6483*nm), h_Planck*c_light/(336.3743*nm), h_Planck*c_light/(336.0278*nm), h_Planck*c_light/(335.7094*nm), h_Planck*c_light/(335.328*nm), 
      h_Planck*c_light/(334.9023*nm), h_Planck*c_light/(334.0541*nm), h_Planck*c_light/(332.7988*nm), 
      h_Planck*c_light/(331.*nm), h_Planck*c_light/(330.*nm), h_Planck*c_light/(329.*nm),
      optPhotMaxE_
    };
    std::vector<G4double> WLS_emiSpectrum = {
      0.0, 
      0.0, 0.0, 0.0,
      0.036, 0.051, 0.066, 0.076, 0.092, 0.108, 0.123, 0.138, 0.153, 0.168, 0.186, 0.199, 0.214, 0.232, 0.25, 0.269, 0.285, 0.3, 0.314, 0.328, 0.341, 0.358, 0.381, 
      0.399, 0.408, 0.426, 0.448, 0.472, 0.489, 0.507, 0.542, 0.529, 0.554, 0.575, 0.591, 0.603, 0.614, 0.631, 0.648, 0.663, 0.683, 0.701, 0.709, 0.723, 0.744, 0.765, 
      0.783, 0.798, 0.824, 0.813, 0.854, 0.839, 0.887, 0.87, 0.917, 0.9, 0.931, 0.951, 0.968, 0.988, 1.0, 0.978, 0.961, 0.943, 0.923, 0.894, 0.908, 0.879, 0.862, 
      0.844, 0.824, 0.803, 0.784, 0.76, 0.742, 0.757, 0.776, 0.79, 0.807, 0.813, 0.796, 0.779, 0.759, 0.742, 0.72, 0.698, 0.66, 0.679, 0.645, 0.61, 0.627, 0.571, 
      0.592, 0.538, 0.555, 0.511, 0.526, 0.466, 0.487, 0.435, 0.451, 0.404, 0.404, 0.368, 0.349, 0.328, 0.29, 0.307, 0.258, 0.274, 0.222, 0.239, 0.194, 0.16, 0.143, 
      0.125, 0.105, 0.082, 0.061, 0.04, 0.022, 0.007, 0.005, 
      0.0, 0.0, 0.0,
      0.0
    };
    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay 
    mpt->AddConstProperty("WLSTIMECONSTANT", 1. * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.82);

    return mpt;
  };

  G4MaterialPropertiesTable* LArPTPArtifact(){
    // This material is model p-Terphenyl, but with the same refractive index as opticalprops::LAr(). 
    // That's why this material must have the same refractive indexs as opticalprops::LAr(), and, on top of 
    // that, it must have also a defined a Optical WLS process (OpWLS) with a WLS absorption length, such 
    // that every LAr scintillation photon (127-128nm) is quickly absorbed and reemitted according to the 
    // PTP emission spectrum. The parameter WLSMEANNUMBERPHOTONS of OpWLS can be used to simulate the 
    // overall efficiency of PTP, both due to an actual-too-long-wlsabslength or due to a conversion efficiency.
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    G4double energies_rindex[] = {2.043*eV, 2.2109*eV, 2.3789*eV, 2.5468*eV, 2.7148*eV, 2.8828*eV, 3.0507*eV, 
                          3.2187*eV, 3.3867*eV, 3.5546*eV, 3.7226*eV, 3.8905*eV, 4.0585*eV, 4.2265*eV, 
                          4.3944*eV, 4.5624*eV, 4.7303*eV, 4.8983*eV, 5.0663*eV, 5.2342*eV, 5.4022*eV, 
                          5.5702*eV, 5.7381*eV, 5.9061*eV, 6.074*eV, 6.242*eV, 6.41*eV, 6.5779*eV, 
                          6.7459*eV, 6.9139*eV, 7.0818*eV, 7.2498*eV, 7.4177*eV, 7.5857*eV, 7.7537*eV, 
                          7.9216*eV, 8.0896*eV, 8.2575*eV, 8.4255*eV, 8.5935*eV, 8.7614*eV, 8.9294*eV, 
                          9.0974*eV, 9.2577*eV, 9.3951*eV, 9.5173*eV, 9.6471*eV, 9.7616*eV, 9.8532*eV, 
                          9.9448*eV, 10.029*eV, 10.097*eV, 10.159*eV, 10.235*eV, 10.288*eV, 10.328*eV, 
                          10.383*eV, 10.461*eV, 10.475*eV, 10.502*eV};

            
    // Before considering any change on this LArPTPArtifact rindex, please take into account that
    // this mpt is meant to have the exact same rindex as LAr's mpt.
    G4double rindex[] = {1.2288, 1.233, 1.2332, 1.2344, 1.235, 1.2362, 1.2373, 1.2379, 1.2393, 1.2405, 
                        1.2414, 1.2432, 1.2441, 1.2452, 1.2471, 1.2495, 1.2511, 1.2527, 1.2548, 1.2571, 
                        1.2592, 1.2614, 1.2641, 1.2667, 1.2694, 1.2725, 1.2758, 1.2794, 1.2835, 1.2877, 
                        1.2928, 1.2979, 1.3031, 1.3085, 1.3146, 1.3218, 1.3291, 1.337, 1.3454, 1.3556, 
                        1.3673, 1.3796, 1.3921, 1.4071, 1.4229, 1.4381, 1.4535, 1.4691, 1.485, 1.4999, 
                        1.5149, 1.5332, 1.5506, 1.5672, 1.5859, 1.5997, 1.6171, 1.6301, 1.6464, 1.6623};

    mpt->AddProperty("RINDEX", energies_rindex, rindex, 60);

    // Absorption length (ABSLENGTH)
    G4double energies_abslength[]  = {optPhotMinE_, optPhotMaxE_};
    G4double abslength[] = {noAbsLength_, noAbsLength_};

    mpt->AddProperty("ABSLENGTH", energies_abslength, abslength, 2);

    // WLS ABSORPTION LENGTH
    // Getting the LAr scintillation spectrum from researchgate.net/figure/Scintillation-light-spectrum-of-LAr-recorded-with-sulfur-beam-excitation-blue-in-color_fig7_258169905
    // , it is clear that there's almost no LAr scintillation for wavelengths bigger than 145nm. The following WLS absorption length is set so as to absorb every LAr scintillation
    // photon in the range [118, 145]nm, and absorb nothing outside that range.
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_, 
      h_Planck*c_light/(146.25*nm), h_Planck*c_light/(145.75*nm), // These points are meant to fix wlsabslength=infinite at 146nm
      h_Planck*c_light/(145.25*nm), h_Planck*c_light/(144.75*nm), // These points are meant to fix wlsabslength=very small at 145nm
      h_Planck*c_light/(118.25*nm), h_Planck*c_light/(117.75*nm), // These points are meant to fix wlsabslength=very small at 118nm
      h_Planck*c_light/(117.25*nm), h_Planck*c_light/(116.75*nm), // These points are meant to fix wlsabslength=infinite at 117nm   
      optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_  ,
      noAbsLength_  , noAbsLength_  ,
      1.*nm         , 1.*nm         ,
      1.*nm         , 1.*nm         ,
      noAbsLength_  , noAbsLength_  ,
      noAbsLength_  
    };
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_, 
      h_Planck*c_light/(404.*nm), h_Planck*c_light/(403.*nm), h_Planck*c_light/(402.*nm),
      h_Planck*c_light/(401.7764*nm), h_Planck*c_light/(399.4336*nm), h_Planck*c_light/(396.6098*nm), h_Planck*c_light/(395.2065*nm), h_Planck*c_light/(393.0391*nm), 
      h_Planck*c_light/(391.6609*nm), h_Planck*c_light/(390.1205*nm), h_Planck*c_light/(389.0432*nm), h_Planck*c_light/(388.1055*nm), h_Planck*c_light/(387.257*nm), 
      h_Planck*c_light/(386.2195*nm), h_Planck*c_light/(385.8109*nm), h_Planck*c_light/(384.8051*nm), h_Planck*c_light/(383.7688*nm), h_Planck*c_light/(382.7972*nm), 
      h_Planck*c_light/(382.0305*nm), h_Planck*c_light/(380.9974*nm), h_Planck*c_light/(379.9932*nm), h_Planck*c_light/(379.087*nm), h_Planck*c_light/(378.0467*nm), 
      h_Planck*c_light/(377.0923*nm), h_Planck*c_light/(376.0287*nm), h_Planck*c_light/(375.1186*nm), h_Planck*c_light/(374.9371*nm), h_Planck*c_light/(374.2467*nm), 
      h_Planck*c_light/(373.9871*nm), h_Planck*c_light/(373.5364*nm), h_Planck*c_light/(372.9633*nm), h_Planck*c_light/(372.7951*nm), h_Planck*c_light/(372.5039*nm), 
      h_Planck*c_light/(372.1126*nm), h_Planck*c_light/(371.8224*nm), h_Planck*c_light/(371.4771*nm), h_Planck*c_light/(371.1768*nm), h_Planck*c_light/(370.9658*nm), 
      h_Planck*c_light/(370.1462*nm), h_Planck*c_light/(369.9143*nm), h_Planck*c_light/(369.4954*nm), h_Planck*c_light/(368.9566*nm), h_Planck*c_light/(368.3537*nm), 
      h_Planck*c_light/(367.4803*nm), h_Planck*c_light/(366.0482*nm), h_Planck*c_light/(363.8461*nm), h_Planck*c_light/(362.5375*nm), h_Planck*c_light/(361.607*nm), 
      h_Planck*c_light/(360.7652*nm), h_Planck*c_light/(360.2621*nm), h_Planck*c_light/(359.6768*nm), h_Planck*c_light/(359.364*nm), h_Planck*c_light/(359.1766*nm), 
      h_Planck*c_light/(358.9998*nm), h_Planck*c_light/(358.8232*nm), h_Planck*c_light/(358.6468*nm), h_Planck*c_light/(358.4601*nm), h_Planck*c_light/(358.284*nm), 
      h_Planck*c_light/(358.1081*nm), h_Planck*c_light/(357.6845*nm), h_Planck*c_light/(357.2105*nm), h_Planck*c_light/(357.0665*nm), h_Planck*c_light/(356.1536*nm), 
      h_Planck*c_light/(354.5748*nm), h_Planck*c_light/(353.2414*nm), h_Planck*c_light/(353.0*nm), h_Planck*c_light/(352.5283*nm), h_Planck*c_light/(352.398*nm), 
      h_Planck*c_light/(351.9379*nm), h_Planck*c_light/(351.8779*nm), h_Planck*c_light/(351.6284*nm), h_Planck*c_light/(351.2499*nm), h_Planck*c_light/(351.1007*nm), 
      h_Planck*c_light/(350.6142*nm), h_Planck*c_light/(350.5052*nm), h_Planck*c_light/(350.0203*nm), h_Planck*c_light/(349.2217*nm), h_Planck*c_light/(348.2604*nm), 
      h_Planck*c_light/(347.1974*nm), h_Planck*c_light/(346.6925*nm), h_Planck*c_light/(346.5278*nm), h_Planck*c_light/(345.6487*nm), h_Planck*c_light/(344.7644*nm), 
      h_Planck*c_light/(344.0374*nm), h_Planck*c_light/(343.5226*nm), h_Planck*c_light/(343.0569*nm), h_Planck*c_light/(342.924*nm), h_Planck*c_light/(342.5829*nm), 
      h_Planck*c_light/(342.2519*nm), h_Planck*c_light/(341.9122*nm), h_Planck*c_light/(341.8744*nm), h_Planck*c_light/(341.8085*nm), h_Planck*c_light/(341.7519*nm), 
      h_Planck*c_light/(341.6672*nm), h_Planck*c_light/(341.338*nm), h_Planck*c_light/(341.1971*nm), h_Planck*c_light/(341.0188*nm), h_Planck*c_light/(340.9437*nm), 
      h_Planck*c_light/(340.5878*nm), h_Planck*c_light/(340.5317*nm), h_Planck*c_light/(340.298*nm), h_Planck*c_light/(340.214*nm), h_Planck*c_light/(340.046*nm), 
      h_Planck*c_light/(339.9714*nm), h_Planck*c_light/(339.5618*nm), h_Planck*c_light/(339.4967*nm), h_Planck*c_light/(339.1623*nm), h_Planck*c_light/(338.9213*nm), 
      h_Planck*c_light/(338.8472*nm), h_Planck*c_light/(338.6066*nm), h_Planck*c_light/(338.4495*nm), h_Planck*c_light/(338.2832*nm), h_Planck*c_light/(338.1264*nm), 
      h_Planck*c_light/(337.7856*nm), h_Planck*c_light/(337.6476*nm), h_Planck*c_light/(337.4913*nm), h_Planck*c_light/(337.2435*nm), h_Planck*c_light/(336.8586*nm), 
      h_Planck*c_light/(336.6483*nm), h_Planck*c_light/(336.3743*nm), h_Planck*c_light/(336.0278*nm), h_Planck*c_light/(335.7094*nm), h_Planck*c_light/(335.328*nm), 
      h_Planck*c_light/(334.9023*nm), h_Planck*c_light/(334.0541*nm), h_Planck*c_light/(332.7988*nm), 
      h_Planck*c_light/(331.*nm), h_Planck*c_light/(330.*nm), h_Planck*c_light/(329.*nm),
      optPhotMaxE_
    };
    std::vector<G4double> WLS_emiSpectrum = {
      0.0, 
      0.0, 0.0, 0.0,
      0.036, 0.051, 0.066, 0.076, 0.092, 0.108, 0.123, 0.138, 0.153, 0.168, 0.186, 0.199, 0.214, 0.232, 0.25, 0.269, 0.285, 0.3, 0.314, 0.328, 0.341, 0.358, 0.381, 
      0.399, 0.408, 0.426, 0.448, 0.472, 0.489, 0.507, 0.542, 0.529, 0.554, 0.575, 0.591, 0.603, 0.614, 0.631, 0.648, 0.663, 0.683, 0.701, 0.709, 0.723, 0.744, 0.765, 
      0.783, 0.798, 0.824, 0.813, 0.854, 0.839, 0.887, 0.87, 0.917, 0.9, 0.931, 0.951, 0.968, 0.988, 1.0, 0.978, 0.961, 0.943, 0.923, 0.894, 0.908, 0.879, 0.862, 
      0.844, 0.824, 0.803, 0.784, 0.76, 0.742, 0.757, 0.776, 0.79, 0.807, 0.813, 0.796, 0.779, 0.759, 0.742, 0.72, 0.698, 0.66, 0.679, 0.645, 0.61, 0.627, 0.571, 
      0.592, 0.538, 0.555, 0.511, 0.526, 0.466, 0.487, 0.435, 0.451, 0.404, 0.404, 0.368, 0.349, 0.328, 0.29, 0.307, 0.258, 0.274, 0.222, 0.239, 0.194, 0.16, 0.143, 
      0.125, 0.105, 0.082, 0.061, 0.04, 0.022, 0.007, 0.005, 
      0.0, 0.0, 0.0,
      0.0
    };
    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // The following data was taken from mdpi-res.com/d_attachment/instruments/instruments-05-00004/article_deploy/instruments-05-00004-v2.pdf?version=1609810328
    // WLS Delay 
    mpt->AddConstProperty("WLSTIMECONSTANT", 1. * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.82);

    return mpt;
  };

  /// 100% Efficient photon collector ///
  G4MaterialPropertiesTable* PerfectPhotonCollector()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> energy =          {optPhotMinE_   ,   optPhotMaxE_};
    std::vector<G4double> efficiency =      {1.             ,   1.          };
    std::vector<G4double> reflectivity =    {0.             ,   0.          };
    std::vector<G4double> rindex =          {1.             ,   1.          };

    mpt->AddProperty("EFFICIENCY",      energy.data(), efficiency.data(),   energy.size());
    mpt->AddProperty("REFLECTIVITY",    energy.data(), reflectivity.data(), energy.size());
    mpt->AddProperty("RINDEX",          energy.data(), rindex.data(),       energy.size());

    return mpt;
  }

  G4MaterialPropertiesTable* PerfectPhotonReflector()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> energy =          {optPhotMinE_   ,   optPhotMaxE_};
    std::vector<G4double> efficiency =      {0.             ,   0.          };
    std::vector<G4double> reflectivity =    {0.85             ,   0.85          };
    std::vector<G4double> rindex =          {1.             ,   1.          };

    mpt->AddProperty("EFFICIENCY",      energy.data(), efficiency.data(),   energy.size());
    mpt->AddProperty("REFLECTIVITY",    energy.data(), reflectivity.data(), energy.size());
    mpt->AddProperty("RINDEX",          energy.data(), rindex.data(),       energy.size());

    return mpt;
  }

  /// Fake Grid ///
  G4MaterialPropertiesTable* FakeGrid(G4double pressure,
                                      G4double temperature,
                                      G4double transparency,
                                      G4double thickness,
                                      G4int    sc_yield,
                                      G4double e_lifetime,
                                      G4double photoe_p)
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // PROPERTIES FROM XENON
    G4MaterialPropertiesTable* xenon_pt = opticalprops::GXe(pressure, temperature, sc_yield, e_lifetime);

    mpt->AddProperty("RINDEX",        xenon_pt->GetProperty("RINDEX"));
    mpt->AddProperty("SCINTILLATIONCOMPONENT1", xenon_pt->GetProperty("SCINTILLATIONCOMPONENT1"));
    mpt->AddProperty("SCINTILLATIONCOMPONENT2", xenon_pt->GetProperty("SCINTILLATIONCOMPONENT2"));

    mpt->AddConstProperty("SCINTILLATIONYIELD", xenon_pt->GetConstProperty("SCINTILLATIONYIELD"));
    mpt->AddConstProperty("RESOLUTIONSCALE",    xenon_pt->GetConstProperty("RESOLUTIONSCALE"));
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",   xenon_pt->GetConstProperty("SCINTILLATIONTIMECONSTANT1"));
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",   xenon_pt->GetConstProperty("SCINTILLATIONTIMECONSTANT2"));
    mpt->AddConstProperty("SCINTILLATIONYIELD1", xenon_pt->GetConstProperty("SCINTILLATIONYIELD1"));
    mpt->AddConstProperty("SCINTILLATIONYIELD2", xenon_pt->GetConstProperty("SCINTILLATIONYIELD2"));
    mpt->AddConstProperty("ATTACHMENT",         xenon_pt->GetConstProperty("ATTACHMENT"), 1);

    // ABSORPTION LENGTH
    G4double abs_length   = -thickness/log(transparency);
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {abs_length, abs_length};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // PHOTOELECTRIC REEMISSION
    // https://aip.scitation.org/doi/10.1063/1.1708797
    G4double stainless_wf = 4.3 * eV; // work function
    mpt->AddConstProperty("WORK_FUNCTION", stainless_wf, 1);
    mpt->AddConstProperty("OP_PHOTOELECTRIC_PROBABILITY", photoe_p, 1);

    return mpt;
  }



  /// PTFE (== TEFLON) ///
  G4MaterialPropertiesTable* PTFE()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFLECTIVITY
    std::vector<G4double> ENERGIES = {
      optPhotMinE_,  2.8 * eV,  3.5 * eV,  4. * eV,
      6. * eV,       7.2 * eV,  optPhotMaxE_
    };
    std::vector<G4double> REFLECTIVITY = {
      .98,  .98,  .98,  .98,
      .72,  .72,  .72
    };
    mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY);

    // REFLEXION BEHAVIOR
    std::vector<G4double> ENERGIES_2    = {optPhotMinE_, optPhotMaxE_};
    // Specular reflection about the normal to a microfacet.
    // Such a vector is chosen according to a gaussian distribution with
    // sigma = SigmaAlhpa (in rad) and centered in the average normal.
    std::vector<G4double> specularlobe  = {0., 0.};
    // specular reflection about the average normal
    std::vector<G4double> specularspike = {0., 0.};
    // 180 degrees reflection.
    std::vector<G4double> backscatter   = {0., 0.};
    // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

    mpt->AddProperty("SPECULARLOBECONSTANT", ENERGIES_2, specularlobe);
    mpt->AddProperty("SPECULARSPIKECONSTANT",ENERGIES_2, specularspike);
    mpt->AddProperty("BACKSCATTERCONSTANT",  ENERGIES_2, backscatter);

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex = {1.41, 1.41};
    mpt->AddProperty("RINDEX", ENERGIES_2, rIndex);

    return mpt;
  }

  /// PolishedAl ///
  G4MaterialPropertiesTable* PolishedAl()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> ENERGIES = {
       h_Planck * c_light / (2456.42541 * nm), h_Planck * c_light / (2396.60266 * nm),
       h_Planck * c_light / (2276.95716 * nm), h_Planck * c_light / (2159.52733 * nm),
       h_Planck * c_light / (2037.66617 * nm), h_Planck * c_light / (1918.02068 * nm),
       h_Planck * c_light / (1798.37518 * nm), h_Planck * c_light / (1676.51403 * nm),
       h_Planck * c_light / (1559.08419 * nm), h_Planck * c_light / (1437.22304 * nm),
       h_Planck * c_light / (1319.79321 * nm), h_Planck * c_light / (1197.93205 * nm),
       h_Planck * c_light / (1078.28656 * nm), h_Planck * c_light / (956.42541 * nm),
       h_Planck * c_light / (838.99557 * nm), h_Planck * c_light / (717.13442 * nm),
       h_Planck * c_light / (597.48892 * nm), h_Planck * c_light / (477.84343 * nm),
       h_Planck * c_light / (418.02068 * nm), h_Planck * c_light / (358.19793 * nm),
       h_Planck * c_light / (293.94387 * nm)
    };
    std::vector<G4double> REFLECTIVITY = {
      .99088, .99082, .98925, .98623, .98611,
      .98163, .98006, .97849, .97401, .97098,
      .96941, .96784, .96481, .96033, .96167,
      .96301, .96289, .96278, .96126, .95830,
      .94224
    };
    // DOI:10.4236/ampc.2015.511046
    mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY);

    // REFLEXION BEHAVIOR
    std::vector<G4double> ENERGIES_2    = {optPhotMinE_, optPhotMaxE_};
    // Specular reflection about the normal to a microfacet.
    // Such a vector is chosen according to a gaussian distribution with
    // sigma = SigmaAlhpa (in rad) and centered in the average normal.
    std::vector<G4double> specularlobe  = {0., 0.};
    // specular reflection about the average normal
    std::vector<G4double> specularspike = {0., 0.};
    // 180 degrees reflection.
    std::vector<G4double> backscatter   = {0., 0.};
    // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

    mpt->AddProperty("SPECULARLOBECONSTANT", ENERGIES_2, specularlobe);
    mpt->AddProperty("SPECULARSPIKECONSTANT",ENERGIES_2, specularspike);
    mpt->AddProperty("BACKSCATTERCONSTANT",  ENERGIES_2, backscatter);

    // REFRACTIVE INDEX
    std::vector<G4double> ENERGIES_3    = {
      0.005 * eV, 0.19581 * eV, 0.43227 * eV,
      0.84211 * eV, 1.2254 * eV, 1.4477 * eV,
      1.7831 * eV, 2.8203 * eV, 3.6216 * eV,
      5.0548 * eV, 7.0554 * eV, 9.4450 * eV,
      12.645 * eV, 14.939 * eV, 16.238 * eV,
      18.4 * eV, 20. * eV
    };
    std::vector<G4double> rIndex = {
      473.49, 12.843, 3.8841, 1.437, 1.4821, 2.4465, 1.6203, 0.58336, 0.32634, 0.1686,
      0.089866, 0.051461, 0.039232, 0.11588, 0.39013, 0.58276, 0.66415
    };
    // from https://refractiveindex.info/?shelf=3d&book=metals&page=aluminium
    mpt->AddProperty("RINDEX", ENERGIES_3, rIndex);

    return mpt;
  }

  /// TPB (tetraphenyl butadiene) ///
  G4MaterialPropertiesTable* TPB()
  {
    // Data from https://doi.org/10.1140/epjc/s10052-018-5807-z
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex_energies = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> TPB_rIndex      = {1.67    , 1.67};
    mpt->AddProperty("RINDEX", rIndex_energies, TPB_rIndex);

    // ABSORPTION LENGTH
    // Assuming no absorption except WLS
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    // This is a combination of figure 11 (for wavelength > 270 nm) and
    // figure 20 (for 50 nm < wavelength < 250 nm).
    // XXX There is a mismatch in the border of the figures that anyway, we implement.
    // Fig 20 -> WLS_absLength = 400 nm for wavelength = 250 nm
    // Fig 11 -> WLS_absLength = 100 nm for wavelength = 270 nm
    // Values for wavelength shorter than 100 nm NOT included as they fit outside
    // the simulation energy limits set in the header.

    //G4double WLS_abs_energy[] = {
    //  optPhotMinE_,                      h_Planck * c_light / (450. * nm),
    //  h_Planck * c_light / (440. * nm),  h_Planck * c_light / (430. * nm),
    //  h_Planck * c_light / (420. * nm),  h_Planck * c_light / (410. * nm),
    //  h_Planck * c_light / (400. * nm),  h_Planck * c_light / (390. * nm),
    //  h_Planck * c_light / (380. * nm),  h_Planck * c_light / (370. * nm),
    //  h_Planck * c_light / (360. * nm),  h_Planck * c_light / (330. * nm),
    //  h_Planck * c_light / (320. * nm),  h_Planck * c_light / (310. * nm),
    //  h_Planck * c_light / (300. * nm),  h_Planck * c_light / (270. * nm),
    //  h_Planck * c_light / (250. * nm),  h_Planck * c_light / (230. * nm),
    //  h_Planck * c_light / (210. * nm),  h_Planck * c_light / (190. * nm),
    //  h_Planck * c_light / (170. * nm),  h_Planck * c_light / (150. * nm),
    //  h_Planck * c_light / (100. * nm),  optPhotMaxE_
    //};
    //const G4int WLS_abs_entries = sizeof(WLS_abs_energy) / sizeof(G4double);
  //
    //G4double WLS_absLength[] = {
    //  noAbsLength_,  noAbsLength_,  //       450 nm
    //  1.e6 * nm,     1.e5 * nm,     // 440 , 430 nm
    //  2.2e4 * nm,     7.e3 * nm,     // 420 , 410 nm
    //  2.2e3 * nm,     700. * nm,     // 400 , 390 nm
    //  200. * nm,      50. * nm,     // 380 , 370 nm
    //   30. * nm,      30. * nm,     // 360 , 330 nm
    //   50. * nm,      80. * nm,     // 320 , 310 nm
    //  100. * nm,     100. * nm,     // 300 , 270 nm
    //  400. * nm,     400. * nm,     // 250 , 230 nm
    //  350. * nm,     250. * nm,     // 210 , 190 nm
    //  350. * nm,     400. * nm,     // 170 , 150 nm
    //  400. * nm,     noAbsLength_   // 100 nm
    //};

    // WLS ABSORPTION LENGTH (Version NoSecWLS)
    // The NoSecWLS is forced by setting the WLS_absLength to noAbsLength_
    // for wavelengths higher than 380 nm where the WLS emission spectrum starts.
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,
      h_Planck * c_light / (380. * nm),  h_Planck * c_light / (370. * nm),
      h_Planck * c_light / (360. * nm),  h_Planck * c_light / (330. * nm),
      h_Planck * c_light / (320. * nm),  h_Planck * c_light / (310. * nm),
      h_Planck * c_light / (300. * nm),  h_Planck * c_light / (270. * nm),
      h_Planck * c_light / (250. * nm),  h_Planck * c_light / (230. * nm),
      h_Planck * c_light / (210. * nm),  h_Planck * c_light / (190. * nm),
      h_Planck * c_light / (170. * nm),  h_Planck * c_light / (150. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_,                 // ~6200 nm
      noAbsLength_,   50. * nm,     // 380 , 370 nm
      30. * nm,      30. * nm,      // 360 , 330 nm
      50. * nm,      80. * nm,      // 320 , 310 nm
      100. * nm,     100. * nm,     // 300 , 270 nm
      400. * nm,     400. * nm,     // 250 , 230 nm
      350. * nm,     250. * nm,     // 210 , 190 nm
      350. * nm,     400. * nm,     // 170 , 150 nm
      400. * nm                     // ~108 nm
    };

    //for (int i=0; i<WLS_abs_energy.size(); i++)
    //  G4cout << "* TPB WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (h_Planck * c_light / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / nm << " nm" << G4endl;
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    // Implemented with formula (7), with parameter values in table (3)
    // Sampling from ~380 nm to 600 nm <--> from 2.06 to 3.26 eV
    const G4int WLS_emi_entries = 120;
    std::vector<G4double> WLS_emi_energy;
    for (int i=0; i<WLS_emi_entries; i++)
      WLS_emi_energy.push_back(2.06 * eV + 0.01 * i * eV);

    std::vector<G4double> WLS_emiSpectrum;
    G4double A      = 0.782;
    G4double alpha  = 3.7e-2;
    G4double sigma1 = 15.43;
    G4double mu1    = 418.10;
    G4double sigma2 = 9.72;
    G4double mu2    = 411.2;

    for (int i=0; i<WLS_emi_entries; i++) {
      G4double wl = (h_Planck * c_light / WLS_emi_energy[i]) / nm;
      WLS_emiSpectrum.push_back(A * (alpha/2.) * exp((alpha/2.) *
                          (2*mu1 + alpha*pow(sigma1,2) - 2*wl)) *
                          erfc((mu1 + alpha*pow(sigma1,2) - wl) / (sqrt(2)*sigma1)) +
                          (1-A) * (1 / sqrt(2*pow(sigma2,2)*3.1416)) *
                                exp((-pow(wl-mu2,2)) / (2*pow(sigma2,2))));
      // G4cout << "* TPB WLSemi:  " << std::setw(4)
      //        << wl << " nm -> " << WLS_emiSpectrum[i] << G4endl;
    };
    mpt->AddProperty("WLSCOMPONENT", WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 1.2 * ns);

    // WLS Quantum Efficiency
    // According to the paper, the QE of TPB depends on the incident wavelength.
    // As Geant4 doesn't allow this possibility, it is set to the value corresponding
    // to Xe scintillation spectrum peak.
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.65);

    return mpt;
  }



  /// Degraded TPB ///
  G4MaterialPropertiesTable* DegradedTPB(G4double wls_eff)
  {
    // It has all the same properties of TPB except the WaveLengthShifting robability
    // that is set by parameter, trying to model a degraded behaviour of the TPB coating

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // All Optical Material Properties from normal TPB ...
    mpt->AddProperty("RINDEX",       opticalprops::TPB()->GetProperty("RINDEX"));
    mpt->AddProperty("ABSLENGTH",    opticalprops::TPB()->GetProperty("ABSLENGTH"));
    mpt->AddProperty("WLSABSLENGTH", opticalprops::TPB()->GetProperty("WLSABSLENGTH"));
    mpt->AddProperty("WLSCOMPONENT", opticalprops::TPB()->GetProperty("WLSCOMPONENT"));
    mpt->AddConstProperty("WLSTIMECONSTANT", opticalprops::TPB()->GetConstProperty("WLSTIMECONSTANT"));

    // Except WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", wls_eff);

    return mpt;
  }


  /// TPH (p-terphenyl) ///
  G4MaterialPropertiesTable* TPH()
  {
    // Data from https://doi.org/10.1016/j.nima.2011.12.036
    // and https://iopscience.iop.org/article/10.1088/1748-0221/5/04/P04007/
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex_energies = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> PTP_rIndex      = {1.65    , 1.65};
    mpt->AddProperty("RINDEX", rIndex_energies, PTP_rIndex);

    // ABSORPTION LENGTH
    // Assuming no absorption except WLS
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);


    // WLS ABSORPTION LENGTH
    // There are no tabulated values in the literature for the PTP absorption
    // length as a function of the wavelength.
    // However in https://iopscience.iop.org/article/10.1088/1748-0221/5/04/P04007/
    // it says that "the thickness >=150 nm of a polycrystalline p-terphenyl layer
    // is enough for absorption of >=99.9% of the xenon light". Thus, using the
    // Beer-Lambert law 1 - P = e^(-x/lambda), where lambda is the absorption length
    // and P is the absorption probability, we can place an upper limit of 21 nm
    // on the absorption length at 175 nm. This is reasonably close to the TPB value.
    // Then, we scale accordingly to the absorption spectrum (which is in a.u.).
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,
      h_Planck * c_light / (337. * nm),
      h_Planck * c_light / (318. * nm), h_Planck * c_light / (292. * nm),
      h_Planck * c_light / (276. * nm), h_Planck * c_light / (253. * nm),
      h_Planck * c_light / (238. * nm), h_Planck * c_light / (222. * nm),
      h_Planck * c_light / (204. * nm), h_Planck * c_light / (190. * nm),
      h_Planck * c_light / (175. * nm), h_Planck * c_light / (169. * nm),
      optPhotMaxE_
    };

    float XePeakAbsValue = 1.879;
    float XePeakAbsLength = 21 * nm;

    std::vector<float> PTP_absorption = {
      0.002, // 337
      0.174, 0.414, // 318, 292
      0.949, 0.540, // 276, 253
      1.218, 1.858, // 238, 222
      1.429, 1.716, // 204, 190
      1.879, 1.803, // 175, 169
    };

    std::vector<G4double> WLS_absLength = {noAbsLength_};
    for (auto& abs_value : PTP_absorption)
      WLS_absLength.push_back(XePeakAbsLength / (abs_value / XePeakAbsValue));

    WLS_absLength.push_back(noAbsLength_);

    //for (int i=0; i<WLS_abs_entries; i++)
    //  G4cout << "* TPB WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (h_Planck * c_light / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / nm << " nm" << G4endl;
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_,
      h_Planck * c_light / (452. * nm), h_Planck * c_light / (430. * nm),
      h_Planck * c_light / (412. * nm), h_Planck * c_light / (398. * nm),
      h_Planck * c_light / (385. * nm), h_Planck * c_light / (371. * nm),
      h_Planck * c_light / (361. * nm), h_Planck * c_light / (354. * nm),
      h_Planck * c_light / (336. * nm), h_Planck * c_light / (317. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.,
      0.044, 0.179,
      0.351, 0.514,
      0.849, 0.993,
      0.745, 0.421,
      0.173, 0.022,
      0.
    };

    mpt->AddProperty("WLSCOMPONENT", WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    // https://www.sciencedirect.com/science/article/pii/S0168900219301652
    mpt->AddConstProperty("WLSTIMECONSTANT", 2.4 * ns);

    // WLS Quantum Efficiency
    // This is set to QE at the Xenon peak, which the paper claims to be >90%
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.9);

    return mpt;
  }



  /// EJ-280 ///
  G4MaterialPropertiesTable* EJ280()
  {
    // https://eljentechnology.com/products/wavelength-shifting-plastics/ej-280-ej-282-ej-284-ej-286
    // and data sheets from the provider.
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> ri_energy = {
      optPhotMinE_,
      h_Planck * c_light / (609. * nm),  h_Planck * c_light / (589.26 * nm),
      h_Planck * c_light / (550. * nm),  h_Planck * c_light / (530.   * nm),
      h_Planck * c_light / (500. * nm),  h_Planck * c_light / (490.   * nm),
      h_Planck * c_light / (481. * nm),  h_Planck * c_light / (460.   * nm),
      h_Planck * c_light / (435. * nm),  h_Planck * c_light / (425.   * nm),
      optPhotMaxE_
    };

    std::vector<G4double> rIndex = {
      1.5780,
      1.5780,  1.5800,   // 609 , 589.26 nm
      1.5845,  1.5870,   // 550 , 530 nm
      1.5913,  1.5929,   // 500 , 490 nm
      1.5945,  1.5986,   // 481 , 460 nm
      1.6050,  1.6080,   // 435 , 425 nm
      1.608
    };
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,                      h_Planck * c_light / (750. * nm),
      h_Planck * c_light / (740. * nm),  h_Planck * c_light / (380. * nm),
      h_Planck * c_light / (370. * nm),  optPhotMaxE_
    };
    std::vector<G4double> absLength = {
      noAbsLength_,  noAbsLength_,
      3.0 * m,       3.0 * m,
      noAbsLength_,  noAbsLength_
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,                      h_Planck * c_light / (500. * nm),
      h_Planck * c_light / (495. * nm),  h_Planck * c_light / (490. * nm),
      h_Planck * c_light / (485. * nm),  h_Planck * c_light / (480. * nm),
      h_Planck * c_light / (475. * nm),  h_Planck * c_light / (470. * nm),
      h_Planck * c_light / (465. * nm),  h_Planck * c_light / (460. * nm),
      h_Planck * c_light / (455. * nm),  h_Planck * c_light / (450. * nm),
      h_Planck * c_light / (445. * nm),  h_Planck * c_light / (440. * nm),
      h_Planck * c_light / (435. * nm),  h_Planck * c_light / (430. * nm),
      h_Planck * c_light / (425. * nm),  h_Planck * c_light / (420. * nm),
      h_Planck * c_light / (415. * nm),  h_Planck * c_light / (410. * nm),
      h_Planck * c_light / (405. * nm),  h_Planck * c_light / (400. * nm),
      h_Planck * c_light / (395. * nm),  h_Planck * c_light / (390. * nm),
      h_Planck * c_light / (385. * nm),  h_Planck * c_light / (380. * nm),
      h_Planck * c_light / (375. * nm),  h_Planck * c_light / (370. * nm),
      h_Planck * c_light / (365. * nm),  h_Planck * c_light / (360. * nm),
      h_Planck * c_light / (355. * nm),  h_Planck * c_light / (350. * nm),
      h_Planck * c_light / (345. * nm),  optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_,        noAbsLength_,
      (1. / 0.0080) * cm,  (1. / 0.0165) * cm,    // 495 , 490 nm
      (1. / 0.0325) * cm,  (1. / 0.0815) * cm,    // 485 , 480 nm
      (1. / 0.2940) * cm,  (1. / 0.9640) * cm,    // 475 , 470 nm
      (1. / 2.8600) * cm,  (1. / 6.3900) * cm,    // 465 , 460 nm
      (1. / 9.9700) * cm,  (1. / 11.0645)* cm,    // 455 , 450 nm
      (1. / 10.198) * cm,  (1. / 9.4465) * cm,    // 445 , 440 nm
      (1. / 10.2145)* cm,  (1. / 12.240) * cm,    // 435 , 430 nm
      (1. / 12.519) * cm,  (1. / 10.867) * cm,    // 425 , 420 nm
      (1. / 9.0710) * cm,  (1. / 8.0895) * cm,    // 415 , 410 nm
      (1. / 7.6650) * cm,  (1. / 6.7170) * cm,    // 405 , 400 nm
      (1. / 5.2460) * cm,  (1. / 4.1185) * cm,    // 395 , 390 nm
      (1. / 3.3175) * cm,  (1. / 2.6800) * cm,    // 385 , 380 nm
      (1. / 1.9610) * cm,  (1. / 1.4220) * cm,    // 375 , 370 nm
      (1. / 1.0295) * cm,  (1. / 0.7680) * cm,    // 365 , 360 nm
      (1. / 0.6865) * cm,  (1. / 0.5885) * cm,    // 355 , 350 nm
      noAbsLength_,        noAbsLength_
    };
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);
    //for (int i=0; i<WLS_abs_entries; i++)
    //  G4cout << "* EJ280 WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (h_Planck * c_light / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / mm << " mm" << G4endl;

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_,                      h_Planck * c_light / (610. * nm),
      h_Planck * c_light / (605. * nm),  h_Planck * c_light / (600. * nm),
      h_Planck * c_light / (595. * nm),  h_Planck * c_light / (590. * nm),
      h_Planck * c_light / (585. * nm),  h_Planck * c_light / (580. * nm),
      h_Planck * c_light / (575. * nm),  h_Planck * c_light / (570. * nm),
      h_Planck * c_light / (565. * nm),  h_Planck * c_light / (560. * nm),
      h_Planck * c_light / (555. * nm),  h_Planck * c_light / (550. * nm),
      h_Planck * c_light / (545. * nm),  h_Planck * c_light / (540. * nm),
      h_Planck * c_light / (535. * nm),  h_Planck * c_light / (530. * nm),
      h_Planck * c_light / (525. * nm),  h_Planck * c_light / (520. * nm),
      h_Planck * c_light / (515. * nm),  h_Planck * c_light / (510. * nm),
      h_Planck * c_light / (505. * nm),  h_Planck * c_light / (500. * nm),
      h_Planck * c_light / (495. * nm),  h_Planck * c_light / (490. * nm),
      h_Planck * c_light / (485. * nm),  h_Planck * c_light / (480. * nm),
      h_Planck * c_light / (475. * nm),  h_Planck * c_light / (470. * nm),
      h_Planck * c_light / (465. * nm),  h_Planck * c_light / (460. * nm),
      h_Planck * c_light / (455. * nm),  h_Planck * c_light / (450. * nm),
      h_Planck * c_light / (445. * nm),  h_Planck * c_light / (440. * nm),
      h_Planck * c_light / (435. * nm),  optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.000,    0.000,   //     , 610 nm
      0.003,    0.006,   // 605 , 600 nm
      0.007,    0.009,   // 595 , 590 nm
      0.014,    0.017,   // 585 , 580 nm
      0.024,    0.033,   // 575 , 570 nm
      0.042,    0.051,   // 565 , 560 nm
      0.063,    0.081,   // 555 , 550 nm
      0.112,    0.157,   // 545 , 540 nm
      0.211,    0.274,   // 535 , 530 nm
      0.329,    0.341,   // 525 , 520 nm
      0.325,    0.346,   // 515 , 510 nm
      0.433,    0.578,   // 505 , 500 nm
      0.792,    1.000,   // 495 , 490 nm
      0.966,    0.718,   // 485 , 480 nm
      0.604,    0.681,   // 475 , 470 nm
      0.708,    0.525,   // 465 , 460 nm
      0.242,    0.046,   // 455 , 450 nm
      0.012,    0.003,   // 445 , 440 nm
      0.000,    0.000    // 435 ,     nm
    };
    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 8.5 * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.86);

    return mpt;
  }

  /// EJ-282 ///
  G4MaterialPropertiesTable* EJ282()
  {
    // https://eljentechnology.com/products/wavelength-shifting-plastics/ej-280-ej-282-ej-284-ej-286
    // and data sheets from the provider.
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX (assuming the same as for EJ286)
    std::vector<G4double> ri_energy = {
      optPhotMinE_,
      h_Planck * c_light / (609. * nm),    h_Planck * c_light / (589.26 * nm),
      h_Planck * c_light / (550. * nm),    h_Planck * c_light / (530.   * nm),
      h_Planck * c_light / (500. * nm),    h_Planck * c_light / (490.   * nm),
      h_Planck * c_light / (481. * nm),    h_Planck * c_light / (460.   * nm),
      h_Planck * c_light / (435. * nm),    h_Planck * c_light / (425.   * nm),
      optPhotMaxE_
    };

    std::vector<G4double> rIndex = {
      1.5780,
      1.5780,  1.5800,   // 609 , 589.26 nm
      1.5845,  1.5870,   // 550 , 530 nm
      1.5913,  1.5929,   // 500 , 490 nm
      1.5945,  1.5986,   // 481 , 460 nm
      1.6050,  1.6080,   // 435 , 425 nm
      1.6080
    };

    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,                      h_Planck * c_light / (750. * nm),
      h_Planck * c_light / (740. * nm),  h_Planck * c_light / (380. * nm),
      h_Planck * c_light / (370. * nm),  optPhotMaxE_
    };

    std::vector<G4double> absLength = {
      noAbsLength_,  noAbsLength_,
      3.0 * m,       3.0 * m,//1.0 * m,       1.0 * m,
      noAbsLength_,  noAbsLength_
    };

    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_, h_Planck * c_light / (490.84 * nm), h_Planck * c_light / (486.15 * nm), 
                    h_Planck * c_light / (479.11 * nm), h_Planck * c_light / (473.0  * nm), 
                    h_Planck * c_light / (466.9 *  nm), h_Planck * c_light / (462.66 * nm), 
                    h_Planck * c_light / (459.36 * nm), h_Planck * c_light / (456.04 * nm), 
                    h_Planck * c_light / (452.72 * nm), h_Planck * c_light / (448.93 * nm), 
                    h_Planck * c_light / (444.63 * nm), h_Planck * c_light / (439.89 * nm), 
                    h_Planck * c_light / (434.69 * nm), h_Planck * c_light / (429.5  * nm), 
                    h_Planck * c_light / (426.19 * nm), h_Planck * c_light / (422.41 * nm), 
                    h_Planck * c_light / (417.7 *  nm), h_Planck * c_light / (412.52 * nm), 
                    h_Planck * c_light / (405.93 * nm), h_Planck * c_light / (400.75 * nm), 
                    h_Planck * c_light / (396.97 * nm), h_Planck * c_light / (393.2  * nm), 
                    h_Planck * c_light / (387.09 * nm), h_Planck * c_light / (381.94 * nm), 
                    h_Planck * c_light / (377.26 * nm), h_Planck * c_light / (372.12 * nm), 
                    h_Planck * c_light / (367.91 * nm), h_Planck * c_light / (364.65 * nm), 
                    h_Planck * c_light / (360.92 * nm), h_Planck * c_light / (357.66 * nm), 
                    h_Planck * c_light / (352.99 * nm), h_Planck * c_light / (349.26 * nm), 
                    h_Planck * c_light / (344.6 *  nm), h_Planck * c_light / (340.41 * nm), 
                    h_Planck * c_light / (336.69 * nm), h_Planck * c_light / (332.49 * nm), 
                    h_Planck * c_light / (327.35 * nm),  
      optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_, noAbsLength_, noAbsLength_, noAbsLength_, noAbsLength_, 
                    504.1636*mm, 149.4848*mm, 81.9183*mm, 36.417*mm, 
                    21.9209*mm, 14.8313*mm, 9.0541*mm, 6.754*mm, 
                    5.4275*mm, 4.6198*mm, 4.0332*mm, 3.5726*mm, 
                    3.113*mm, 2.7394*mm, 2.5869*mm, 2.1279*mm, 
                    1.6981*mm, 1.2173*mm, 1.0827*mm, 1.3746*mm, 
                    1.8709*mm, 2.2824*mm, 2.7394*mm, 3.2736*mm, 
                    3.7805*mm, 4.4167*mm, 5.0977*mm, 5.8328*mm, 
                    7.1454*mm, 9.3416*mm, 11.5256*mm, 13.872*mm, 
                    17.3803*mm, 
      noAbsLength_
    };

    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);


    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_, h_Planck * c_light / (598.72 * nm), h_Planck * c_light / (586.83 * nm), 
                    h_Planck * c_light / (575.91 * nm), h_Planck * c_light / (563.06 * nm), 
                    h_Planck * c_light / (552.29 * nm), h_Planck * c_light / (543.68 * nm), 
                    h_Planck * c_light / (536.24 * nm), h_Planck * c_light / (529.98 * nm), 
                    h_Planck * c_light / (523.32 * nm), h_Planck * c_light / (516.67 * nm), 
                    h_Planck * c_light / (511.19 * nm), h_Planck * c_light / (506.89 * nm), 
                    h_Planck * c_light / (503.36 * nm), h_Planck * c_light / (500.43 * nm), 
                    h_Planck * c_light / (498.28 * nm), h_Planck * c_light / (495.93 * nm), 
                    h_Planck * c_light / (493.58 * nm), h_Planck * c_light / (490.84 * nm), 
                    h_Planck * c_light / (487.58 * nm), h_Planck * c_light / (480.27 * nm), 
                    h_Planck * c_light / (472.37 * nm), h_Planck * c_light / (466.97 * nm), 
                    h_Planck * c_light / (463.84 * nm), h_Planck * c_light / (461.96 * nm), 
                    h_Planck * c_light / (460.23 * nm), h_Planck * c_light / (458.36 * nm), 
                    h_Planck * c_light / (456.79 * nm), h_Planck * c_light / (455.23 * nm), 
                    h_Planck * c_light / (452.57 * nm), h_Planck * c_light / (450.53 * nm),  
      optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.0000,   0.0064, 0.0136, 0.0297, 0.0504, 0.0809, 0.1286, 0.1888, 0.2484, 0.3114, 0.3709, 
                0.4305, 0.4964, 0.5704, 0.6392, 0.7046, 0.7734, 0.834, 0.8994, 0.9624, 0.9949, 
                0.9606, 0.8981, 0.8358, 0.7695, 0.698, 0.6205, 0.5449, 0.4618, 0.3111, 0.2238,
      0.0000
    };

    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 1.9 * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.93);

    return mpt;
  }

  /// EJ-286 ///
  G4MaterialPropertiesTable* EJ286(G4double attenuation_length)
  {
    // https://eljentechnology.com/products/wavelength-shifting-plastics/ej-280-ej-282-ej-284-ej-286
    // and data sheets from the provider.
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> ri_energy = {
      optPhotMinE_,
      h_Planck * c_light / (609. * nm),    h_Planck * c_light / (589.26 * nm),
      h_Planck * c_light / (550. * nm),    h_Planck * c_light / (530.   * nm),
      h_Planck * c_light / (500. * nm),    h_Planck * c_light / (490.   * nm),
      h_Planck * c_light / (481. * nm),    h_Planck * c_light / (460.   * nm),
      h_Planck * c_light / (435. * nm),    h_Planck * c_light / (425.   * nm),
      optPhotMaxE_
    };

    std::vector<G4double> rIndex = {
      1.5780,
      1.5780,  1.5800,   // 609 , 589.26 nm
      1.5845,  1.5870,   // 550 , 530 nm
      1.5913,  1.5929,   // 500 , 490 nm
      1.5945,  1.5986,   // 481 , 460 nm
      1.6050,  1.6080,   // 435 , 425 nm
      1.6080
    };
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,                      h_Planck * c_light / (750. * nm),
      h_Planck * c_light / (740. * nm),  h_Planck * c_light / (380. * nm),
      h_Planck * c_light / (370. * nm),  optPhotMaxE_
    };
    std::vector<G4double> absLength = {
      attenuation_length,  attenuation_length,
      attenuation_length,  attenuation_length,
      attenuation_length,  attenuation_length
    };

    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,                      h_Planck * c_light / (445. * nm),
                                        h_Planck * c_light / (440. * nm),
      h_Planck * c_light / (435. * nm),  h_Planck * c_light / (430. * nm),
      h_Planck * c_light / (425. * nm),  h_Planck * c_light / (420. * nm),
      h_Planck * c_light / (415. * nm),  h_Planck * c_light / (410. * nm),
      h_Planck * c_light / (405. * nm),  h_Planck * c_light / (400. * nm),
      h_Planck * c_light / (395. * nm),  h_Planck * c_light / (390. * nm),
      h_Planck * c_light / (385. * nm),  h_Planck * c_light / (380. * nm),
      h_Planck * c_light / (375. * nm),  h_Planck * c_light / (370. * nm),
      h_Planck * c_light / (365. * nm),  h_Planck * c_light / (360. * nm),
      h_Planck * c_light / (355. * nm),  h_Planck * c_light / (350. * nm),
      h_Planck * c_light / (345. * nm),  h_Planck * c_light / (340. * nm),
      h_Planck * c_light / (335. * nm),  h_Planck * c_light / (330. * nm),
      h_Planck * c_light / (325. * nm),  h_Planck * c_light / (320. * nm),
      h_Planck * c_light / (315. * nm),  h_Planck * c_light / (310. * nm),
      h_Planck * c_light / (305. * nm),  h_Planck * c_light / (300. * nm),
      h_Planck * c_light / (295. * nm),  h_Planck * c_light / (290. * nm),
      h_Planck * c_light / (285. * nm),  h_Planck * c_light / (280. * nm),
      h_Planck * c_light / (275. * nm),  optPhotMaxE_
    };


    std::vector<G4double> WLS_absLength = {
      noAbsLength_,         noAbsLength_,
                            (1. / 0.00007) * cm,    //       440 nm
      (1. /  0.0003) * cm,  (1. / 0.00104) * cm,    // 435 , 430 nm
      (1. / 0.00223) * cm,  (1. / 0.00408) * cm,    // 425 , 420 nm
      (1. /  0.0104) * cm,  (1. / 0.18544) * cm,    // 415 , 410 nm
      (1. /  1.4094) * cm,  (1. /  3.7088) * cm,    // 405 , 400 nm
      (1. /  7.4176) * cm,  (1. / 11.8682) * cm,    // 395 , 390 nm
      (1. / 16.6155) * cm,  (1. / 22.2529) * cm,    // 385 , 380 nm
      (1. / 27.8162) * cm,  (1. / 33.3794) * cm,    // 375 , 370 nm
      (1. / 37.8671) * cm,  (1. / 40.4262) * cm,    // 365 , 360 nm
      (1. / 41.5388) * cm,  (1. / 41.1679) * cm,    // 355 , 350 nm
      (1. / 38.9426) * cm,  (1. / 35.0113) * cm,    // 345 , 340 nm
      (1. / 31.1541) * cm,  (1. / 27.4453) * cm,    // 335 , 330 nm
      (1. / 23.4398) * cm,  (1. / 20.0276) * cm,    // 325 , 320 nm
      (1. / 16.3188) * cm,  (1. / 13.3518) * cm,    // 315 , 310 nm
      (1. / 10.5331) * cm,  (1. /  8.1594) * cm,    // 305 , 300 nm
      (1. /  6.1196) * cm,  (1. /  4.6731) * cm,    // 295 , 290 nm
      (1. /  3.6346) * cm,  (1. /  3.0412) * cm,    // 285 , 280 nm
      noAbsLength_,         noAbsLength_
    };
    // XXX We are assuming that EJ286 doesn't absorb wave lengths shorter than 280 nm
    // although the spectrum continues ...
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);
    //for (int i=0; i<WLS_abs_entries; i++)
    //  G4cout << "* EJ286 WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (h_Planck * c_light / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / mm << " mm" << G4endl;

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_,
      h_Planck * c_light / (533. * nm), h_Planck * c_light / (532. * nm), h_Planck * c_light / (531. * nm),
      h_Planck * c_light / (530. * nm), h_Planck * c_light / (525. * nm),  
      h_Planck * c_light / (520. * nm), h_Planck * c_light / (515. * nm),  
      h_Planck * c_light / (510. * nm), h_Planck * c_light / (505. * nm),  
      h_Planck * c_light / (500. * nm), h_Planck * c_light / (495. * nm),  
      h_Planck * c_light / (490. * nm), h_Planck * c_light / (485. * nm),  
      h_Planck * c_light / (480. * nm), h_Planck * c_light / (475. * nm),  
      h_Planck * c_light / (470. * nm), h_Planck * c_light / (465. * nm),  
      h_Planck * c_light / (460. * nm), h_Planck * c_light / (455. * nm),  
      h_Planck * c_light / (450. * nm), h_Planck * c_light / (445. * nm),  
      h_Planck * c_light / (440. * nm), h_Planck * c_light / (435. * nm),  
      h_Planck * c_light / (430. * nm), h_Planck * c_light / (425. * nm),  
      h_Planck * c_light / (420. * nm), h_Planck * c_light / (415. * nm),  
      h_Planck * c_light / (410. * nm), h_Planck * c_light / (405. * nm),  
      h_Planck * c_light / (400. * nm), h_Planck * c_light / (395. * nm),  
      h_Planck * c_light / (390. * nm), h_Planck * c_light / (385. * nm),  
      h_Planck * c_light / (380. * nm), h_Planck * c_light / (375. * nm),  
      h_Planck * c_light / (370. * nm), h_Planck * c_light / (365. * nm),  
      h_Planck * c_light / (364. * nm), h_Planck * c_light / (363. * nm), h_Planck * c_light / (362. * nm),  
      optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.0000,
      0.0000, 0.0000, 0.0000,
      0.0089, 0.0100,  
      0.0181, 0.0210,  
      0.0270, 0.0380,
      0.0496, 0.0600,  
      0.0721, 0.0900,
      0.1125, 0.1500,  
      0.1848, 0.2100,  
      0.2388, 0.2800,  
      0.3289, 0.4000,  
      0.4956, 0.5700,  
      0.6230, 0.6450,  
      0.6667, 0.8000,  
      0.9800, 0.9900,  
      0.8559, 0.7118,  
      0.7400, 0.8000,  
      0.6702, 0.3800,  
      0.1082, 0.0400,  
      0.0000, 0.0000,  0.0000,
      0.0000
    };
    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 1.2 * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.92);

    return mpt;
  }

  G4MaterialPropertiesTable* G2P_FB118(G4double attenuation_length)
  {
    // iopscience.iop.org/article/10.1088/1748-0221/16/09/P09027
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> ri_energy = {
      optPhotMinE_,
      h_Planck * c_light / (609. * nm),    h_Planck * c_light / (589.26 * nm),
      h_Planck * c_light / (550. * nm),    h_Planck * c_light / (530.   * nm),
      h_Planck * c_light / (500. * nm),    h_Planck * c_light / (490.   * nm),
      h_Planck * c_light / (481. * nm),    h_Planck * c_light / (460.   * nm),
      h_Planck * c_light / (435. * nm),    h_Planck * c_light / (425.   * nm),
      optPhotMaxE_
    };

    std::vector<G4double> rIndex = {
      1.502,
      1.502,  1.502,   // 609 , 589.26 nm
      1.502,  1.502,   // 550 , 530 nm
      1.502,  1.502,   // 500 , 490 nm
      1.502,  1.502,   // 481 , 460 nm
      1.502,  1.502,   // 435 , 425 nm
      1.502
    };
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,                      h_Planck * c_light / (750. * nm),
      h_Planck * c_light / (740. * nm),  h_Planck * c_light / (380. * nm),
      h_Planck * c_light / (370. * nm),  optPhotMaxE_
    };
    std::vector<G4double> absLength = {
      attenuation_length,  attenuation_length,
      attenuation_length,  attenuation_length,
      attenuation_length,  attenuation_length
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_, h_Planck*c_light/(404.1469*nm), h_Planck*c_light/(403.7884*nm), h_Planck*c_light/(403.4298*nm), h_Planck*c_light/(403.0713*nm), 
                    h_Planck*c_light/(402.7128*nm), h_Planck*c_light/(402.3542*nm), h_Planck*c_light/(401.9957*nm), h_Planck*c_light/(401.6372*nm), 
                    h_Planck*c_light/(401.2786*nm), h_Planck*c_light/(400.9201*nm), h_Planck*c_light/(400.5616*nm), h_Planck*c_light/(400.203*nm), 
                    h_Planck*c_light/(399.8445*nm), h_Planck*c_light/(399.486*nm), h_Planck*c_light/(399.1274*nm), h_Planck*c_light/(398.7689*nm), 
                    h_Planck*c_light/(398.4104*nm), h_Planck*c_light/(398.0518*nm), h_Planck*c_light/(397.6933*nm), h_Planck*c_light/(397.3348*nm), 
                    h_Planck*c_light/(396.9762*nm), h_Planck*c_light/(396.6177*nm), h_Planck*c_light/(396.2592*nm), h_Planck*c_light/(395.9006*nm), 
                    h_Planck*c_light/(395.5421*nm), h_Planck*c_light/(395.1836*nm), h_Planck*c_light/(394.825*nm), h_Planck*c_light/(394.4665*nm), 
                    h_Planck*c_light/(394.108*nm), h_Planck*c_light/(393.7494*nm), h_Planck*c_light/(393.3909*nm), h_Planck*c_light/(393.0324*nm), 
                    h_Planck*c_light/(392.6738*nm), h_Planck*c_light/(392.3153*nm), h_Planck*c_light/(391.9568*nm), h_Planck*c_light/(391.5982*nm), 
                    h_Planck*c_light/(391.2397*nm), h_Planck*c_light/(390.8812*nm), h_Planck*c_light/(390.5226*nm), h_Planck*c_light/(390.1641*nm), 
                    h_Planck*c_light/(389.8056*nm), h_Planck*c_light/(389.447*nm), h_Planck*c_light/(389.0885*nm), h_Planck*c_light/(388.73*nm), 
                    h_Planck*c_light/(388.3714*nm), h_Planck*c_light/(388.0129*nm), h_Planck*c_light/(387.6544*nm), h_Planck*c_light/(387.2958*nm), 
                    h_Planck*c_light/(386.9373*nm), h_Planck*c_light/(386.5788*nm), h_Planck*c_light/(386.2202*nm), h_Planck*c_light/(385.8617*nm), 
                    h_Planck*c_light/(385.5032*nm), h_Planck*c_light/(385.1446*nm), h_Planck*c_light/(384.7861*nm), h_Planck*c_light/(384.4276*nm), 
                    h_Planck*c_light/(384.069*nm), h_Planck*c_light/(383.7105*nm), h_Planck*c_light/(383.352*nm), h_Planck*c_light/(382.9934*nm), 
                    h_Planck*c_light/(382.6349*nm), h_Planck*c_light/(382.2764*nm), h_Planck*c_light/(381.9178*nm), h_Planck*c_light/(381.5593*nm), 
                    h_Planck*c_light/(381.2008*nm), h_Planck*c_light/(380.8422*nm), h_Planck*c_light/(380.4837*nm), h_Planck*c_light/(380.1252*nm), 
                    h_Planck*c_light/(379.7666*nm), h_Planck*c_light/(379.4081*nm), h_Planck*c_light/(379.0496*nm), h_Planck*c_light/(378.691*nm), 
                    h_Planck*c_light/(378.3325*nm), h_Planck*c_light/(377.974*nm), h_Planck*c_light/(377.6154*nm), h_Planck*c_light/(377.2569*nm), 
                    h_Planck*c_light/(376.8984*nm), h_Planck*c_light/(376.5398*nm), h_Planck*c_light/(376.1813*nm), h_Planck*c_light/(375.8228*nm), 
                    h_Planck*c_light/(375.4642*nm), h_Planck*c_light/(375.1057*nm), h_Planck*c_light/(374.7472*nm), h_Planck*c_light/(374.3886*nm), 
                    h_Planck*c_light/(374.0301*nm), h_Planck*c_light/(373.6716*nm), h_Planck*c_light/(373.313*nm), h_Planck*c_light/(372.9545*nm), 
                    h_Planck*c_light/(372.596*nm), h_Planck*c_light/(372.2374*nm), h_Planck*c_light/(371.8789*nm), h_Planck*c_light/(371.5204*nm), 
                    h_Planck*c_light/(371.1618*nm), h_Planck*c_light/(370.8033*nm), h_Planck*c_light/(370.4448*nm), h_Planck*c_light/(370.0862*nm), 
                    h_Planck*c_light/(369.7277*nm), h_Planck*c_light/(369.3692*nm), h_Planck*c_light/(369.0106*nm), h_Planck*c_light/(368.6521*nm), 
                    h_Planck*c_light/(368.2936*nm), h_Planck*c_light/(367.935*nm), h_Planck*c_light/(367.5765*nm), h_Planck*c_light/(367.218*nm), 
                    h_Planck*c_light/(366.8594*nm), h_Planck*c_light/(366.5009*nm), h_Planck*c_light/(366.1424*nm), h_Planck*c_light/(365.7838*nm), 
                    h_Planck*c_light/(365.4253*nm), h_Planck*c_light/(365.0668*nm), h_Planck*c_light/(364.7082*nm), h_Planck*c_light/(364.3497*nm), 
                    h_Planck*c_light/(363.9912*nm), h_Planck*c_light/(363.6326*nm), h_Planck*c_light/(363.2741*nm), h_Planck*c_light/(362.9156*nm), 
                    h_Planck*c_light/(362.557*nm), h_Planck*c_light/(362.1985*nm), h_Planck*c_light/(361.84*nm), h_Planck*c_light/(361.4814*nm), 
                    h_Planck*c_light/(361.1229*nm), h_Planck*c_light/(360.7644*nm), h_Planck*c_light/(360.4058*nm), h_Planck*c_light/(360.0473*nm), 
                    h_Planck*c_light/(359.6888*nm), h_Planck*c_light/(359.3302*nm), h_Planck*c_light/(358.9717*nm), h_Planck*c_light/(358.6132*nm), 
                    h_Planck*c_light/(358.2546*nm), h_Planck*c_light/(357.8961*nm), h_Planck*c_light/(357.5376*nm), h_Planck*c_light/(357.179*nm), 
                    h_Planck*c_light/(356.8205*nm), h_Planck*c_light/(356.462*nm), h_Planck*c_light/(356.1034*nm), h_Planck*c_light/(355.7449*nm), 
                    h_Planck*c_light/(355.3864*nm), h_Planck*c_light/(355.0278*nm), h_Planck*c_light/(354.6693*nm), h_Planck*c_light/(354.3108*nm), 
                    h_Planck*c_light/(353.9522*nm), h_Planck*c_light/(353.5937*nm), h_Planck*c_light/(353.2352*nm), h_Planck*c_light/(352.8766*nm), 
                    h_Planck*c_light/(352.5181*nm), h_Planck*c_light/(352.1596*nm), h_Planck*c_light/(351.801*nm), h_Planck*c_light/(351.4425*nm), 
                    h_Planck*c_light/(351.084*nm), h_Planck*c_light/(350.7254*nm), h_Planck*c_light/(350.3669*nm), h_Planck*c_light/(350.0084*nm), 
                    h_Planck*c_light/(349.6498*nm), h_Planck*c_light/(349.2913*nm), h_Planck*c_light/(348.9328*nm), h_Planck*c_light/(348.5742*nm), 
                    h_Planck*c_light/(348.2157*nm), h_Planck*c_light/(347.8572*nm), h_Planck*c_light/(347.4986*nm), h_Planck*c_light/(347.1401*nm), 
                    h_Planck*c_light/(346.7816*nm), h_Planck*c_light/(346.423*nm), h_Planck*c_light/(346.0645*nm), h_Planck*c_light/(345.706*nm), 
                    h_Planck*c_light/(345.3474*nm), h_Planck*c_light/(344.9889*nm), h_Planck*c_light/(344.6304*nm), h_Planck*c_light/(344.2718*nm), 
                    h_Planck*c_light/(343.9133*nm), h_Planck*c_light/(343.5548*nm), h_Planck*c_light/(343.1963*nm), h_Planck*c_light/(342.8377*nm), 
                    h_Planck*c_light/(342.4792*nm), h_Planck*c_light/(342.1207*nm), h_Planck*c_light/(341.7621*nm), h_Planck*c_light/(341.4036*nm), 
                    h_Planck*c_light/(341.0451*nm), h_Planck*c_light/(340.6865*nm), h_Planck*c_light/(340.328*nm), h_Planck*c_light/(339.9695*nm), 
                    h_Planck*c_light/(339.6109*nm), h_Planck*c_light/(339.2524*nm), h_Planck*c_light/(338.8939*nm), h_Planck*c_light/(338.5353*nm), 
                    h_Planck*c_light/(338.1768*nm), h_Planck*c_light/(337.8183*nm), h_Planck*c_light/(337.4597*nm), h_Planck*c_light/(337.1012*nm), 
                    h_Planck*c_light/(336.7427*nm), h_Planck*c_light/(336.3841*nm), h_Planck*c_light/(336.0256*nm), h_Planck*c_light/(335.6671*nm), 
                    h_Planck*c_light/(335.3085*nm), h_Planck*c_light/(334.95*nm), h_Planck*c_light/(334.5915*nm), h_Planck*c_light/(334.2329*nm), 
                    h_Planck*c_light/(333.8744*nm), h_Planck*c_light/(333.5159*nm), h_Planck*c_light/(333.1573*nm), 
      optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_, 4.511 *mm, 4.132 *mm, 3.823 *mm, 3.564 *mm, 3.343 *mm, 3.148 *mm, 2.97 *mm, 2.803 *mm, 2.64 *mm, 2.476 *mm, 2.31 *mm, 2.147 *mm, 1.99 *mm, 
                    1.842 *mm, 1.706 *mm, 1.586 *mm, 1.485 *mm, 1.406 *mm, 1.345 *mm, 1.3 *mm, 1.269 *mm, 1.247 *mm, 1.23 *mm, 1.216 *mm, 1.203 *mm, 1.191 *mm, 
                    1.18 *mm, 1.171 *mm, 1.163 *mm, 1.155 *mm, 1.148 *mm, 1.142 *mm, 1.136 *mm, 1.131 *mm, 1.124 *mm, 1.117 *mm, 1.107 *mm, 1.094 *mm, 1.077 *mm, 
                    1.054 *mm, 1.025 *mm, 0.986 *mm, 0.935 *mm, 0.863 *mm, 0.751 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 
                    0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 
                    0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 
                    0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 
                    0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 
                    0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 
                    0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 
                    0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 
                    0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 
                    0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.248 *mm, 0.48 *mm, 
                    0.842 *mm, 0.98 *mm, 1.086 *mm, 1.175 *mm, 1.256 *mm, 1.33 *mm, 1.401 *mm, 1.469 *mm, 1.535 *mm, 1.6 *mm, 1.665 *mm, 1.729 *mm, 1.794 *mm, 
                    1.859 *mm, 1.926 *mm, 1.993 *mm, 2.063 *mm, 2.134 *mm, 2.207 *mm, 2.283 *mm, 2.36 *mm, 2.44 *mm, 2.521 *mm, 2.604 *mm, 2.689 *mm, 2.775 *mm, 
                    2.864 *mm, 2.954 *mm, 3.045 *mm, 3.138 *mm,   
      noAbsLength_
    };
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_, h_Planck * c_light / (602.00 * nm), h_Planck * c_light / (601.00 * nm),
                    h_Planck * c_light / (599.44 * nm), h_Planck * c_light / (579.10 * nm), 
                    h_Planck * c_light / (547.19 * nm), h_Planck * c_light / (522.45 * nm), 
                    h_Planck * c_light / (502.87 * nm), h_Planck * c_light / (490.87 * nm), 
                    h_Planck * c_light / (481.27 * nm), h_Planck * c_light / (472.86 * nm), 
                    h_Planck * c_light / (467.23 * nm), h_Planck * c_light / (464.41 * nm), 
                    h_Planck * c_light / (462.00 * nm), h_Planck * c_light / (457.20 * nm), 
                    h_Planck * c_light / (451.61 * nm), h_Planck * c_light / (446.42 * nm), 
                    h_Planck * c_light / (443.61 * nm), h_Planck * c_light / (440.39 * nm), 
                    h_Planck * c_light / (437.96 * nm), h_Planck * c_light / (436.34 * nm), 
                    h_Planck * c_light / (435.13 * nm), h_Planck * c_light / (433.52 * nm), 
                    h_Planck * c_light / (430.73 * nm), h_Planck * c_light / (429.14 * nm), 
                    h_Planck * c_light / (427.57 * nm), h_Planck * c_light / (425.21 * nm), 
                    h_Planck * c_light / (422.08 * nm), h_Planck * c_light / (418.54 * nm), 
                    h_Planck * c_light / (414.61 * nm), h_Planck * c_light / (411.48 * nm), 
                    h_Planck * c_light / (409.54 * nm), h_Planck * c_light / (406.81 * nm), 
                    h_Planck * c_light / (404.45 * nm), h_Planck * c_light / (401.67 * nm), 
                    h_Planck * c_light / (394.50 * nm), h_Planck * c_light / (380.54 * nm),
                    h_Planck * c_light / (379.00 * nm), h_Planck * c_light / (378.00 * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.0000,   
                0.0000, 0.0000,
                0.0044, 0.0081, 
                0.0328, 0.0786, 
                0.1782, 0.2586, 
                0.3237, 0.4214, 
                0.5210, 0.5766, 
                0.6188, 0.6533, 
                0.6666, 0.6915, 
                0.7298, 0.8046, 
                0.8755, 0.9407, 
                0.9732, 0.9943, 
                1.0000, 0.9789, 
                0.9329, 0.8562, 
                0.7373, 0.6242, 
                0.5034, 0.3672, 
                0.2484, 0.0988, 
                0.0336, 0.0106, 
                0.0029, 0.0009,
                0.0000, 0.0000,
      0.0000
    };
    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 1.2 * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.90);

    return mpt;

  }

  G4MaterialPropertiesTable* SCHOTT_B270()
  {
      // refractiveindex.info/?shelf=glass&book=SCHOTT-multipurpose&page=B270   <-- They got this refractive index from:
      // refractiveindex.info/download/data/1993/Schott_B270.pdf                <-- Which is a SCHOTT datasheet from 1993. The data in these one matches the data in this newer datasheet from SCHOTT
      // media.schott.com/api/public/content/bd52af0b4b4e4fe7b73817a5c92eb71c?v=5f2e41a4&download=true&_ga=2.119782792.1059616897.1654520734-1060119037.1654183676
      // schott.com/es-es/products/b-270-p1000313/technical-details             <-- Here, they also provide 250-450nm transmitance curves


    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> ri_energy =   {optPhotMinE_, 2.7675*eV, 2.7799*eV, 2.7924*eV, 2.8051*eV, 2.8178*eV, 2.8307*eV, 
                                        2.8437*eV, 2.8568*eV, 2.87*eV, 2.8834*eV, 2.8968*eV, 2.9104*eV, 2.9242*eV, 2.938*eV, 
                                        2.952*eV, 2.9661*eV, 2.9804*eV, 2.9948*eV, 3.0093*eV, 3.024*eV, 3.0388*eV, 3.0538*eV, 
                                        3.0689*eV, 3.0842*eV, 3.0996*eV, 3.1152*eV, 3.1309*eV, 3.1468*eV, 3.1629*eV, 3.1791*eV, 
                                        3.1955*eV, 3.212*eV, 3.2288*eV, 3.2457*eV, 3.2627*eV, 3.28*eV, 3.2975*eV, 3.3151*eV, 
                                        3.3329*eV, 3.3509*eV, 3.3691*eV, 3.3875*eV, 3.4062*eV, 3.425*eV, 3.444*eV, 3.4632*eV, 
                                        3.4827*eV, 3.5024*eV, 3.5223*eV, 3.5424*eV, 3.5628*eV, 3.5834*eV, 3.6042*eV, 3.6253*eV, 
                                        3.6466*eV, 3.6682*eV, 3.69*eV, 3.7121*eV, 3.7345*eV, 3.7571*eV, 3.78*eV, 3.8032*eV, 
                                        .8267*eV, 3.8504*eV, 3.8745*eV, 3.8989*eV, 3.9236*eV, 3.9485*eV, 3.9739*eV, 3.9995*eV, 
                                        4.0255*eV, 4.0518*eV, 4.0784*eV, optPhotMaxE_};
    std::vector<G4double> rIndex =  {1.5328, 1.5328, 1.533, 1.5332, 1.5334, 1.5336, 1.5338, 1.534, 1.5342, 1.5344, 1.5346, 1.5348, 
                                    1.535, 1.5352, 1.5354, 1.5356, 1.5358, 1.5361, 1.5363, 1.5365, 1.5367, 1.5369, 1.5371, 1.5374, 
                                    1.5376, 1.5378, 1.538, 1.5383, 1.5385, 1.5387, 1.539, 1.5392, 1.5394, 1.5397, 1.5399, 1.5401, 
                                    1.5404, 1.5406, 1.5409, 1.5411, 1.5413, 1.5416, 1.5418, 1.5421, 1.5423, 1.5426, 1.5429, 1.5431, 
                                    1.5434, 1.5436, 1.5439, 1.5441, 1.5444, 1.5447, 1.5449, 1.5452, 1.5455, 1.5457, 1.546, 1.5463, 
                                    1.5466, 1.5468, 1.5471, 1.5474, 1.5477, 1.5479, 1.5482, 1.5485, 1.5488, 1.5491, 1.5494, 1.5497, 
                                    1.5499, 1.5502, 1.5502};
    mpt->AddProperty("RINDEX", ri_energy, rIndex);


    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy =  {optPhotMinE_, 2.7588*eV, 2.7733*eV, 2.7879*eV, 2.8027*eV, 2.8227*eV, 2.8407*eV, 2.8532*eV, 
                                        2.8687*eV, 2.8818*eV, 2.8949*eV, 2.9149*eV, 2.9257*eV, 2.9366*eV, 2.9488*eV, 2.9641*eV, 
                                        2.9822*eV, 2.9949*eV, 3.0068*eV, 3.0175*eV, 3.0282*eV, 3.0423*eV, 3.0614*eV, 3.0763*eV, 
                                        3.0928*eV, 3.111*eV, 3.1282*eV, 3.1398*eV, 3.1514*eV, 3.1631*eV, 3.1748*eV, 3.1974*eV, 
                                        3.2199*eV, 3.2449*eV, 3.2671*eV, 3.2787*eV, 3.3024*eV, 3.3232*eV, 3.3451*eV, 3.3709*eV, 
                                        3.3895*eV, 3.4163*eV, 3.442*eV, 3.4641*eV, 3.4818*eV, 3.5056*eV, 3.5165*eV, 3.5393*eV, 
                                        3.5568*eV, 3.5747*eV, 3.5944*eV, 3.6118*eV, 3.6314*eV, 3.6495*eV, 3.6612*eV, 3.6782*eV, 
                                        3.6997*eV, 3.7171*eV, 3.7345*eV, 3.7485*eV, 3.7637*eV, 3.7804*eV, 3.8014*eV, 3.8142*eV, 
                                        3.8328*eV, 3.8472*eV, 3.8574*eV, 3.866*eV, 3.8852*eV, 3.8888*eV, 3.9015*eV, 3.9104*eV, 
                                        3.9169*eV, 3.9267*eV, 3.9389*eV, 3.9497*eV, 3.9634*eV, 3.9742*eV, 3.9851*eV, 3.9945*eV, 
                                        4.0007*eV, 4.0085*eV, 4.018*eV, 4.029*eV, 4.0402*eV, 4.0563*eV, 4.0679*eV, 4.0756*eV, 
                                        4.0886*eV, 4.0993*eV, 4.1133*eV, 4.1299*eV, 4.1433*eV, 4.1534*eV, 4.1635*eV, 4.1805*eV, 
                                        4.2011*eV, 4.222*eV, 4.2421*eV, 4.2655*eV, 4.3001*eV, 4.338*eV, 4.38*eV, 4.4166*eV, 
                                        4.4665*eV, 4.511*eV, 4.5629*eV, 4.6161*eV, 4.6774*eV, 4.7404*eV, 4.7943*eV, 4.8604*eV, optPhotMaxE_};
    std::vector<G4double> absLength =   {23.3403*mm, 23.3403*mm, 22.8518*mm, 22.8518*mm, 22.8518*mm, 23.3403*mm, 23.0327*mm, 22.8518*mm, 22.8518*mm, 
                                        22.8518*mm, 23.3403*mm, 22.8518*mm, 22.3851*mm, 23.3403*mm, 23.6499*mm, 23.6499*mm, 23.3403*mm, 23.3403*mm, 
                                        23.6499*mm, 23.6499*mm, 23.6499*mm, 22.8518*mm, 23.3403*mm, 23.3403*mm, 23.3403*mm, 23.3403*mm, 23.0327*mm, 
                                        23.0327*mm, 23.0327*mm, 23.3374*mm, 23.6499*mm, 23.3403*mm, 22.6792*mm, 22.8518*mm, 21.889*mm, 22.2598*mm, 
                                        21.9337*mm, 21.9337*mm, 21.5324*mm, 21.0836*mm, 20.9742*mm, 20.6806*mm, 20.4822*mm, 19.8994*mm, 19.4537*mm, 
                                        18.6188*mm, 18.6188*mm, 17.8475*mm, 17.136*mm, 16.4761*mm, 16.0121*mm, 15.1974*mm, 14.1325*mm, 13.3397*mm, 
                                        12.8225*mm, 12.0863*mm, 11.1549*mm, 10.4281*mm, 9.8144*mm, 9.205*mm, 8.3988*mm, 7.7092*mm, 7.1125*mm, 
                                        6.7153*mm, 6.1858*mm, 5.8705*mm, 5.4229*mm, 5.1275*mm, 4.6954*mm, 4.571*mm, 4.257*mm, 4.0262*mm, 
                                        3.8241*mm, 3.6253*mm, 3.354*mm, 3.1466*mm, 2.864*mm, 2.6756*mm, 2.5024*mm, 2.3898*mm, 2.2897*mm, 
                                        2.188*mm, 2.0444*mm, 1.9268*mm, 1.8156*mm, 1.6477*mm, 1.5535*mm, 1.4743*mm, 1.3886*mm, 1.3037*mm, 
                                        1.1916*mm, 1.1031*mm, 1.0222*mm, 0.9416*mm, 0.8991*mm, 0.8175*mm, 0.734*mm, 0.6535*mm, 0.5682*mm, 
                                        0.5059*mm, 0.4079*mm, 0.3696*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 
                                        0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);
    return mpt;
  }

  /// Y-11 ///
  G4MaterialPropertiesTable* Y11()
  {
    // http://kuraraypsf.jp/psf/index.html
    // http://kuraraypsf.jp/psf/ws.html
    // Excel provided by kuraray with Tabulated WLS absorption lengths
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> ri_energy = {
      optPhotMinE_,  optPhotMaxE_
    };
    std::vector<G4double> rIndex = {
      1.59,  1.59
    };
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,                      h_Planck * c_light / (750. * nm),
      h_Planck * c_light / (740. * nm),  h_Planck * c_light / (380. * nm),
      h_Planck * c_light / (370. * nm),  optPhotMaxE_
    };
    std::vector<G4double> absLength = {
      noAbsLength_,  noAbsLength_,
      3.5 * m,       3.5 * m,
      noAbsLength_,  noAbsLength_
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,                      h_Planck * c_light / (490. * nm),
      h_Planck * c_light / (485. * nm),  h_Planck * c_light / (475. * nm),
      h_Planck * c_light / (454. * nm),  h_Planck * c_light / (443. * nm),
      h_Planck * c_light / (430. * nm),  h_Planck * c_light / (410. * nm),
      h_Planck * c_light / (405. * nm),  h_Planck * c_light / (359. * nm),
      h_Planck * c_light / (350. * nm),  h_Planck * c_light / (345. * nm),
      optPhotMaxE_
    };
    std::vector<G4double> WLS_absLength = {
      noAbsLength_,  noAbsLength_,    //     , 490 nm
      44.2  * mm,    5.39 * mm,       // 485 , 475 nm
      0.395 * mm,    0.462 * mm,      // 454 , 443 nm
      0.354 * mm,    0.571 * mm,      // 430 , 410 nm
      0.612 * mm,    4.51 * mm,       // 405 , 359 nm
      4.81  * mm,    noAbsLength_,    // 350 , 345 nm
      noAbsLength_
    };
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);
    //for (int i=0; i<WLS_abs_entries; i++)
    //  G4cout << "* Y11 WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (h_Planck * c_light / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / mm << " mm" << G4endl;

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_,                      h_Planck * c_light / (580. * nm),
      h_Planck * c_light / (550. * nm),  h_Planck * c_light / (530. * nm),
      h_Planck * c_light / (525. * nm),  h_Planck * c_light / (520. * nm),
      h_Planck * c_light / (515. * nm),  h_Planck * c_light / (510. * nm),
      h_Planck * c_light / (505. * nm),  h_Planck * c_light / (500. * nm),
      h_Planck * c_light / (495. * nm),  h_Planck * c_light / (490. * nm),
      h_Planck * c_light / (485. * nm),  h_Planck * c_light / (480. * nm),
      h_Planck * c_light / (475. * nm),  h_Planck * c_light / (470. * nm),
      h_Planck * c_light / (465. * nm),  h_Planck * c_light / (460. * nm),
      h_Planck * c_light / (455. * nm),  h_Planck * c_light / (450. * nm),
      h_Planck * c_light / (445. * nm),  optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.000,    0.000,   //     , 580 nm
      0.200,    0.300,   // 550 , 530 nm
      0.400,    0.600,   // 525 , 520 nm
      0.750,    0.750,   // 515 , 510 nm
      0.720,    0.700,   // 505 , 500 nm
      0.680,    0.650,   // 495 , 490 nm
      0.700,    0.900,   // 485 , 480 nm
      1.000,    0.950,   // 475 , 470 nm
      0.500,    0.300,   // 465 , 460 nm
      0.100,    0.050,   // 455 , 450 nm
      0.000,    0.000    // 445 ,     nm
    };
    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 8.5 * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.87);

    return mpt;
  }


/// B-2 ///
  G4MaterialPropertiesTable* B2()
  {
    // http://kuraraypsf.jp/psf/index.html
    // http://kuraraypsf.jp/psf/ws.html
    // Excel provided by kuraray with Tabulated WLS absorption lengths
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> ri_energy = {
      optPhotMinE_,  optPhotMaxE_
    };
    std::vector<G4double> rIndex = {
      1.59,  1.59
    };
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,                      h_Planck * c_light / (750. * nm),
      h_Planck * c_light / (740. * nm),  h_Planck * c_light / (380. * nm),
      h_Planck * c_light / (370. * nm),  optPhotMaxE_
    };
    std::vector<G4double> absLength = {
      noAbsLength_,  noAbsLength_,
      3.5 * m,       3.5 * m,
      noAbsLength_,  noAbsLength_
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    // For B2 fibers Kuraray provides absorption spectrum and not
    // absorption length. We assume that the absorption length at the
    // absorption maximum is the same as with the Y11 fiber and
    // scale according to the absorption spectrum. This is not perfect
    // but it was verified to be a good approximation with the Y11 fiber,
    // for which Kuraray did provide the absorption length.

    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,
      h_Planck * c_light / (418. * nm), h_Planck * c_light / (412. * nm),
      h_Planck * c_light / (405. * nm), h_Planck * c_light / (400. * nm),
      h_Planck * c_light / (394. * nm), h_Planck * c_light / (387. * nm),
      h_Planck * c_light / (384. * nm), h_Planck * c_light / (382. * nm),
      h_Planck * c_light / (378. * nm), h_Planck * c_light / (370. * nm),
      h_Planck * c_light / (361. * nm), h_Planck * c_light / (353. * nm),
      h_Planck * c_light / (345. * nm), h_Planck * c_light / (341. * nm),
      h_Planck * c_light / (336. * nm), h_Planck * c_light / (331. * nm),
      h_Planck * c_light / (316. * nm), h_Planck * c_light / (301. * nm),
      h_Planck * c_light / (280. * nm),
      optPhotMaxE_
    };

    float minAbsLength = 0.395 * mm;

    std::vector<float> B2_absorption {
      -0.01,        // 418
      -0.06, -0.26, // 405, 412
      -0.44, -0.59, // 394, 400
      -0.59, -0.64, // 384, 387
      -0.77, -0.92, // 378, 382
      -1.00, -0.93, // 361, 370
      -0.85, -0.87, // 345, 353
      -0.87, -0.77, // 336, 341
      -0.56, -0.35, // 316, 331
      -0.22, -0.12  // 280, 301
    };

    std::vector<G4double> WLS_absLength {noAbsLength_};

    for (auto &abs_value : B2_absorption)
      WLS_absLength.push_back(- minAbsLength / abs_value);

    WLS_absLength.push_back(noAbsLength_);

    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_,
      h_Planck * c_light / (542 * nm), h_Planck * c_light / (525 * nm),
      h_Planck * c_light / (508 * nm), h_Planck * c_light / (497 * nm),
      h_Planck * c_light / (488 * nm), h_Planck * c_light / (479 * nm),
      h_Planck * c_light / (473 * nm), h_Planck * c_light / (467 * nm),
      h_Planck * c_light / (463 * nm), h_Planck * c_light / (458 * nm),
      h_Planck * c_light / (454 * nm), h_Planck * c_light / (449 * nm),
      h_Planck * c_light / (445 * nm), h_Planck * c_light / (442 * nm),
      h_Planck * c_light / (440 * nm), h_Planck * c_light / (438 * nm),
      h_Planck * c_light / (433 * nm), h_Planck * c_light / (429 * nm),
      h_Planck * c_light / (424 * nm), h_Planck * c_light / (420 * nm),
      h_Planck * c_light / (418 * nm), h_Planck * c_light / (416 * nm),
      h_Planck * c_light / (411 * nm), h_Planck * c_light / (404 * nm),
      h_Planck * c_light / (402 * nm), h_Planck * c_light / (399 * nm),
      h_Planck * c_light / (398 * nm), h_Planck * c_light / (396 * nm),
      h_Planck * c_light / (395 * nm), h_Planck * c_light / (394 * nm),
      h_Planck * c_light / (392 * nm), h_Planck * c_light / (391 * nm),
      h_Planck * c_light / (386 * nm), h_Planck * c_light / (380 * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.000,
      0.053, 0.070, // 542, 525
      0.109, 0.143, // 508, 497
      0.199, 0.270, // 488, 479
      0.337, 0.423, // 473, 467
      0.497, 0.582, // 463, 458
      0.615, 0.645, // 454, 449
      0.679, 0.750, // 445, 442
      0.801, 0.857, // 440, 438
      0.957, 0.999, // 433, 429
      0.949, 0.906, // 424, 420
      0.855, 0.809, // 418, 416
      0.750, 0.750, // 411, 404
      0.719, 0.671, // 402, 399
      0.590, 0.500, // 398, 396
      0.421, 0.327, // 395, 394
      0.217, 0.138, // 392, 391
      0.065, 0.023, // 386, 380
      0.000
    };

    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 8.5 * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.87);

    return mpt;
  }



  /// Pethylene ///
  G4MaterialPropertiesTable* Pethylene()
  {
    // Fiber cladding material.
    // Properties from geant4/examples/extended/optical/wls
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex_energies = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> rIndex          = {1.49, 1.49};
    mpt->AddProperty("RINDEX", rIndex_energies, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }



  /// FPethylene ///
  G4MaterialPropertiesTable* FPethylene()
  {
    // Fiber cladding material.
    // Properties from geant4/examples/extended/optical/wls
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex_energies = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> rIndex          = {1.42, 1.42};
    mpt->AddProperty("RINDEX", rIndex_energies, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }



  /// PMMA == PolyMethylmethacrylate ///
  G4MaterialPropertiesTable* PMMA()
  {
    // Fiber cladding material.
    // Properties from geant4/examples/extended/optical/wls
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex_energies = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> rIndex          = {1.49, 1.49};
    mpt->AddProperty("RINDEX", rIndex_energies, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,
      2.722 * eV,  3.047 * eV,  3.097 * eV,  3.136 * eV,  3.168 * eV,  3.229 * eV,  3.291 * eV,
      3.323 * eV,  3.345 * eV,  3.363 * eV,  3.397 * eV,  3.451 * eV,  3.511 * eV,  3.590 * eV,
      optPhotMaxE_
    };
    std::vector<G4double> abslength = {
      noAbsLength_,
      noAbsLength_,  4537. * mm,  329.7 * mm,  98.60 * mm,  36.94 * mm,  10.36 * mm,  4.356 * mm,
      2.563 * mm,    1.765 * mm,  1.474 * mm,  1.153 * mm,  0.922 * mm,  0.765 * mm,  0.671 * mm,
      0.671 * mm
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, abslength);

    return mpt;
  }



  /// XXX ///
  G4MaterialPropertiesTable* XXX()
  {
    // Playing material properties
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex_energies = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> rIndex          = {1.0    , 1.0};
    mpt->AddProperty("RINDEX", rIndex_energies, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {10.*cm, 10.*cm};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    std::vector<G4double> WLS_abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> WLS_absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy  = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> WLS_emiSpectrum = {1.0, 1.0};
    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 1. * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 1.);

    return mpt;
  }

  /// VIKUITI ///
  G4MaterialPropertiesTable* specularspikeVIKUITI()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    //Info from David Warner 11/2019 talk
    //It is said that the reflective foils are coated with TPB,
    //shall we care about this?

    //Reflectivity data taken from https://indico.fnal.gov/event/24273/contributions/188657/attachments/130083/158244/DUNE_60Review1.pdf for 45º curve
    const G4int entries = 120;
    G4double energy[]       = {1.0334*eV, 1.0378*eV, 1.0454*eV, 1.0532*eV, 1.0645*eV, 1.0731*eV, 
                              1.0786*eV, 1.0875*eV, 1.0984*eV, 1.1081*eV, 1.1209*eV, 1.1289*eV, 
                              1.1415*eV, 1.156*eV, 1.1658*eV, 1.1746*eV, 1.1836*eV, 1.1986*eV, 
                              1.2126*eV, 1.2222*eV, 1.2304*eV, 1.2424*eV, 1.2531*eV, 1.2721*eV, 
                              1.2832*eV, 1.2944*eV, 1.3053*eV, 1.3175*eV, 1.3379*eV, 1.3521*eV, 
                              1.3586*eV, 1.364*eV, 1.3684*eV, 1.3729*eV, 1.3774*eV, 1.3819*eV, 
                              1.3864*eV, 1.3899*eV, 1.3922*eV, 1.3945*eV, 1.3968*eV, 1.3991*eV, 
                              1.4015*eV, 1.4038*eV, 1.4062*eV, 1.4086*eV, 1.4109*eV, 1.4133*eV, 
                              1.4157*eV, 1.4181*eV, 1.4205*eV, 1.4229*eV, 1.4253*eV, 1.4277*eV, 
                              1.4301*eV, 1.4338*eV, 1.4387*eV, 1.4436*eV, 1.4486*eV, 1.4535*eV, 
                              1.4611*eV, 1.4725*eV, 1.4905*eV, 1.505*eV, 1.5253*eV, 1.5461*eV, 
                              1.5703*eV, 1.5864*eV, 1.6014*eV, 1.629*eV, 1.6415*eV, 1.6608*eV, 
                              1.6805*eV, 1.7006*eV, 1.7213*eV, 1.7424*eV, 1.777*eV, 1.8015*eV, 
                              1.8321*eV, 1.8553*eV, 1.8811*eV, 1.9064*eV, 1.9324*eV, 1.9625*eV, 
                              2.0148*eV, 2.0748*eV, 2.1314*eV, 2.1634*eV, 2.1969*eV, 2.2286*eV, 
                              2.3104*eV, 2.3982*eV, 2.436*eV, 2.4786*eV, 2.5227*eV, 2.6192*eV, 
                              2.716*eV, 2.8056*eV, 2.9177*eV, 3.0176*eV, 3.1041*eV, 3.1423*eV, 
                              3.1852*eV, 3.2051*eV, 3.219*eV, 3.2309*eV, 3.2369*eV, 3.2488*eV, 
                              3.2489*eV, 3.2608*eV, 3.2609*eV, 3.2728*eV, 3.2847*eV, 3.2968*eV, 
                              3.3083*eV, 3.3088*eV, 3.3216*eV, 3.373*eV, 3.4437*eV, 3.5064*eV};

    // The following array contains the first two (and last two) entries of energy[])
    // This one is for the sake of comfortably adding constant properties through the 
    // whole energy range given by energy[]
    G4int energy_span_entries = 4;
    G4double energy_span[]  = {1.0334*eV, 1.0378*eV, 3.4437*eV, 3.5064*eV};
                              
    G4double reflectivity[] = {0.1848, 0.1848, 0.1848, 0.1848, 0.1823, 0.1849, 0.1849, 0.187, 
                              0.1892, 0.1865, 0.1789, 0.1772, 0.1866, 0.1927, 0.1851, 0.1825, 
                              0.1835, 0.1963, 0.199, 0.1903, 0.1826, 0.1852, 0.1927, 0.1948, 
                              0.1895, 0.1965, 0.2045, 0.2077, 0.212, 0.2274, 0.2502, 0.2718, 
                              0.2925, 0.3165, 0.3421, 0.3701, 0.4013, 0.4269, 0.4445, 0.4653, 
                              0.4917, 0.5261, 0.5517, 0.5788, 0.6028, 0.6252, 0.646, 0.67, 
                              0.694, 0.7132, 0.7372, 0.7596, 0.7788, 0.798, 0.8171, 0.8411, 
                              0.8755, 0.9019, 0.9243, 0.9435, 0.9631, 0.9811, 0.9894, 0.9857, 
                              0.9852, 0.9852, 0.9852, 0.9788, 0.9797, 0.9853, 0.9869, 0.9885, 
                              0.9901, 0.9886, 0.9896, 0.9918, 0.9854, 0.9854, 0.9804, 0.9784, 
                              0.9796, 0.9791, 0.9799, 0.9832, 0.984, 0.9823, 0.9806, 0.9828, 
                              0.9809, 0.9825, 0.9859, 0.9806, 0.981, 0.9818, 0.9816, 0.986, 
                              0.9819, 0.9796, 0.9741, 0.9724, 0.9706, 0.9608, 0.9325, 0.9072, 
                              0.875, 0.839, 0.8166, 0.751, 0.7702, 0.6919, 0.7094, 0.6183, 
                              0.512, 0.4376, 0.2596, 0.3233, 0.1908, 0.145, 0.1463, 0.1393};
                            
    G4double test_energy[] = {opticalprops::optPhotMinE_,
                            opticalprops::optPhotMinE_ +(1*eV), 
                            opticalprops::optPhotMaxE_ -(1*eV),
                            opticalprops::optPhotMaxE_};
    G4double test_reflectivity[] = {1., 1., 1., 1.};

    mpt->AddProperty("REFLECTIVITY", energy, reflectivity, entries);

    // From geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html#optical-photon-processes
    // The specular lobe constant (material property name SPECULARLOBECONSTANT) represents the reflection probability 
    // about the normal of a micro facet. The specular spike constant (material property name SPECULARSPIKECONSTANT), 
    // in turn, illustrates the probability of reflection about the average surface normal. The diffuse lobe constant 
    // is for the probability of internal Lambertian reflection, and finally the back-scatter spike constant (material 
    // property name BACKSCATTERCONSTANT) is for the case of several reflections within a deep groove with the ultimate 
    // result of exact back-scattering. The four probabilities add up to one, with the diffuse lobe constant being calculated 
    // by the code from other other three values that the user entered.
    
    G4double specularlobe[]     = {0., 0., 0., 0.};
    G4double specularspike[]    = {1., 1., 1., 1.};
    G4double backscatter[]      = {0., 0., 0., 0.};
    G4double efficiency[]       = {0., 0., 0., 0.};

    // By default, transmission equals to 0, so no need to explictly set it.
    // "The material properties REFLECTIVITY and TRANSMISSION are used. 
    // By default, REFLECTIVITY equals 1 and TRANSMISSION equals 0."

    mpt->AddProperty("SPECULARLOBECONSTANT", energy_span, specularlobe, energy_span_entries);
    mpt->AddProperty("SPECULARSPIKECONSTANT", energy_span, specularspike, energy_span_entries);
    mpt->AddProperty("BACKSCATTERCONSTANT", energy_span, backscatter, energy_span_entries);
    mpt->AddProperty("EFFICIENCY", energy_span, efficiency, energy_span_entries);
    
    return mpt;
  }

  G4MaterialPropertiesTable* specularlobeVIKUITI()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    //Info from David Warner 11/2019 talk
    //It is said that the reflective foils are coated with TPB,
    //shall we care about this?

    //Reflectivity data taken from https://indico.fnal.gov/event/24273/contributions/188657/attachments/130083/158244/DUNE_60Review1.pdf for 45º curve
    const G4int entries = 120;
    G4double energy[]       = {1.0334*eV, 1.0378*eV, 1.0454*eV, 1.0532*eV, 1.0645*eV, 1.0731*eV, 
                              1.0786*eV, 1.0875*eV, 1.0984*eV, 1.1081*eV, 1.1209*eV, 1.1289*eV, 
                              1.1415*eV, 1.156*eV, 1.1658*eV, 1.1746*eV, 1.1836*eV, 1.1986*eV, 
                              1.2126*eV, 1.2222*eV, 1.2304*eV, 1.2424*eV, 1.2531*eV, 1.2721*eV, 
                              1.2832*eV, 1.2944*eV, 1.3053*eV, 1.3175*eV, 1.3379*eV, 1.3521*eV, 
                              1.3586*eV, 1.364*eV, 1.3684*eV, 1.3729*eV, 1.3774*eV, 1.3819*eV, 
                              1.3864*eV, 1.3899*eV, 1.3922*eV, 1.3945*eV, 1.3968*eV, 1.3991*eV, 
                              1.4015*eV, 1.4038*eV, 1.4062*eV, 1.4086*eV, 1.4109*eV, 1.4133*eV, 
                              1.4157*eV, 1.4181*eV, 1.4205*eV, 1.4229*eV, 1.4253*eV, 1.4277*eV, 
                              1.4301*eV, 1.4338*eV, 1.4387*eV, 1.4436*eV, 1.4486*eV, 1.4535*eV, 
                              1.4611*eV, 1.4725*eV, 1.4905*eV, 1.505*eV, 1.5253*eV, 1.5461*eV, 
                              1.5703*eV, 1.5864*eV, 1.6014*eV, 1.629*eV, 1.6415*eV, 1.6608*eV, 
                              1.6805*eV, 1.7006*eV, 1.7213*eV, 1.7424*eV, 1.777*eV, 1.8015*eV, 
                              1.8321*eV, 1.8553*eV, 1.8811*eV, 1.9064*eV, 1.9324*eV, 1.9625*eV, 
                              2.0148*eV, 2.0748*eV, 2.1314*eV, 2.1634*eV, 2.1969*eV, 2.2286*eV, 
                              2.3104*eV, 2.3982*eV, 2.436*eV, 2.4786*eV, 2.5227*eV, 2.6192*eV, 
                              2.716*eV, 2.8056*eV, 2.9177*eV, 3.0176*eV, 3.1041*eV, 3.1423*eV, 
                              3.1852*eV, 3.2051*eV, 3.219*eV, 3.2309*eV, 3.2369*eV, 3.2488*eV, 
                              3.2489*eV, 3.2608*eV, 3.2609*eV, 3.2728*eV, 3.2847*eV, 3.2968*eV, 
                              3.3083*eV, 3.3088*eV, 3.3216*eV, 3.373*eV, 3.4437*eV, 3.5064*eV};

    // The following array contains the first two (and last two) entries of energy[])
    // This one is for the sake of comfortably adding constant properties through the 
    // whole energy range given by energy[]
    G4int energy_span_entries = 4;
    G4double energy_span[]  = {1.0334*eV, 1.0378*eV, 3.4437*eV, 3.5064*eV};
                              
    G4double reflectivity[] = {0.1848, 0.1848, 0.1848, 0.1848, 0.1823, 0.1849, 0.1849, 0.187, 
                              0.1892, 0.1865, 0.1789, 0.1772, 0.1866, 0.1927, 0.1851, 0.1825, 
                              0.1835, 0.1963, 0.199, 0.1903, 0.1826, 0.1852, 0.1927, 0.1948, 
                              0.1895, 0.1965, 0.2045, 0.2077, 0.212, 0.2274, 0.2502, 0.2718, 
                              0.2925, 0.3165, 0.3421, 0.3701, 0.4013, 0.4269, 0.4445, 0.4653, 
                              0.4917, 0.5261, 0.5517, 0.5788, 0.6028, 0.6252, 0.646, 0.67, 
                              0.694, 0.7132, 0.7372, 0.7596, 0.7788, 0.798, 0.8171, 0.8411, 
                              0.8755, 0.9019, 0.9243, 0.9435, 0.9631, 0.9811, 0.9894, 0.9857, 
                              0.9852, 0.9852, 0.9852, 0.9788, 0.9797, 0.9853, 0.9869, 0.9885, 
                              0.9901, 0.9886, 0.9896, 0.9918, 0.9854, 0.9854, 0.9804, 0.9784, 
                              0.9796, 0.9791, 0.9799, 0.9832, 0.984, 0.9823, 0.9806, 0.9828, 
                              0.9809, 0.9825, 0.9859, 0.9806, 0.981, 0.9818, 0.9816, 0.986, 
                              0.9819, 0.9796, 0.9741, 0.9724, 0.9706, 0.9608, 0.9325, 0.9072, 
                              0.875, 0.839, 0.8166, 0.751, 0.7702, 0.6919, 0.7094, 0.6183, 
                              0.512, 0.4376, 0.2596, 0.3233, 0.1908, 0.145, 0.1463, 0.1393};
                            
    G4double test_energy[] = {opticalprops::optPhotMinE_,
                            opticalprops::optPhotMinE_ +(1*eV), 
                            opticalprops::optPhotMaxE_ -(1*eV),
                            opticalprops::optPhotMaxE_};
    G4double test_reflectivity[] = {1., 1., 1., 1.};

    mpt->AddProperty("REFLECTIVITY", energy, reflectivity, entries);

    // From geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html#optical-photon-processes
    // The specular lobe constant (material property name SPECULARLOBECONSTANT) represents the reflection probability 
    // about the normal of a micro facet. The specular spike constant (material property name SPECULARSPIKECONSTANT), 
    // in turn, illustrates the probability of reflection about the average surface normal. The diffuse lobe constant 
    // is for the probability of internal Lambertian reflection, and finally the back-scatter spike constant (material 
    // property name BACKSCATTERCONSTANT) is for the case of several reflections within a deep groove with the ultimate 
    // result of exact back-scattering. The four probabilities add up to one, with the diffuse lobe constant being calculated 
    // by the code from other other three values that the user entered.
    
    G4double specularlobe[]     = {1., 1., 1., 1.};
    G4double specularspike[]    = {0., 0., 0., 0.};
    G4double backscatter[]      = {0., 0., 0., 0.};
    G4double efficiency[]       = {0., 0., 0., 0.};

    // By default, transmission equals to 0, so no need to explictly set it.
    // "The material properties REFLECTIVITY and TRANSMISSION are used. 
    // By default, REFLECTIVITY equals 1 and TRANSMISSION equals 0."

    mpt->AddProperty("SPECULARLOBECONSTANT", energy_span, specularlobe, energy_span_entries);
    mpt->AddProperty("SPECULARSPIKECONSTANT", energy_span, specularspike, energy_span_entries);
    mpt->AddProperty("BACKSCATTERCONSTANT", energy_span, backscatter, energy_span_entries);
    mpt->AddProperty("EFFICIENCY", energy_span, efficiency, energy_span_entries);
    
    return mpt;
  }

  G4MaterialPropertiesTable* diffusiveVIKUITI()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    //Info from David Warner 11/2019 talk
    //It is said that the reflective foils are coated with TPB,
    //shall we care about this?

    //Reflectivity data taken from https://indico.fnal.gov/event/24273/contributions/188657/attachments/130083/158244/DUNE_60Review1.pdf for 45º curve
    const G4int entries = 120;
    G4double energy[]       = {1.0334*eV, 1.0378*eV, 1.0454*eV, 1.0532*eV, 1.0645*eV, 1.0731*eV, 
                              1.0786*eV, 1.0875*eV, 1.0984*eV, 1.1081*eV, 1.1209*eV, 1.1289*eV, 
                              1.1415*eV, 1.156*eV, 1.1658*eV, 1.1746*eV, 1.1836*eV, 1.1986*eV, 
                              1.2126*eV, 1.2222*eV, 1.2304*eV, 1.2424*eV, 1.2531*eV, 1.2721*eV, 
                              1.2832*eV, 1.2944*eV, 1.3053*eV, 1.3175*eV, 1.3379*eV, 1.3521*eV, 
                              1.3586*eV, 1.364*eV, 1.3684*eV, 1.3729*eV, 1.3774*eV, 1.3819*eV, 
                              1.3864*eV, 1.3899*eV, 1.3922*eV, 1.3945*eV, 1.3968*eV, 1.3991*eV, 
                              1.4015*eV, 1.4038*eV, 1.4062*eV, 1.4086*eV, 1.4109*eV, 1.4133*eV, 
                              1.4157*eV, 1.4181*eV, 1.4205*eV, 1.4229*eV, 1.4253*eV, 1.4277*eV, 
                              1.4301*eV, 1.4338*eV, 1.4387*eV, 1.4436*eV, 1.4486*eV, 1.4535*eV, 
                              1.4611*eV, 1.4725*eV, 1.4905*eV, 1.505*eV, 1.5253*eV, 1.5461*eV, 
                              1.5703*eV, 1.5864*eV, 1.6014*eV, 1.629*eV, 1.6415*eV, 1.6608*eV, 
                              1.6805*eV, 1.7006*eV, 1.7213*eV, 1.7424*eV, 1.777*eV, 1.8015*eV, 
                              1.8321*eV, 1.8553*eV, 1.8811*eV, 1.9064*eV, 1.9324*eV, 1.9625*eV, 
                              2.0148*eV, 2.0748*eV, 2.1314*eV, 2.1634*eV, 2.1969*eV, 2.2286*eV, 
                              2.3104*eV, 2.3982*eV, 2.436*eV, 2.4786*eV, 2.5227*eV, 2.6192*eV, 
                              2.716*eV, 2.8056*eV, 2.9177*eV, 3.0176*eV, 3.1041*eV, 3.1423*eV, 
                              3.1852*eV, 3.2051*eV, 3.219*eV, 3.2309*eV, 3.2369*eV, 3.2488*eV, 
                              3.2489*eV, 3.2608*eV, 3.2609*eV, 3.2728*eV, 3.2847*eV, 3.2968*eV, 
                              3.3083*eV, 3.3088*eV, 3.3216*eV, 3.373*eV, 3.4437*eV, 3.5064*eV};

    // The following array contains the first two (and last two) entries of energy[])
    // This one is for the sake of comfortably adding constant properties through the 
    // whole energy range given by energy[]
    G4int energy_span_entries = 4;
    G4double energy_span[]  = {1.0334*eV, 1.0378*eV, 3.4437*eV, 3.5064*eV};
                              
    G4double reflectivity[] = {0.1848, 0.1848, 0.1848, 0.1848, 0.1823, 0.1849, 0.1849, 0.187, 
                              0.1892, 0.1865, 0.1789, 0.1772, 0.1866, 0.1927, 0.1851, 0.1825, 
                              0.1835, 0.1963, 0.199, 0.1903, 0.1826, 0.1852, 0.1927, 0.1948, 
                              0.1895, 0.1965, 0.2045, 0.2077, 0.212, 0.2274, 0.2502, 0.2718, 
                              0.2925, 0.3165, 0.3421, 0.3701, 0.4013, 0.4269, 0.4445, 0.4653, 
                              0.4917, 0.5261, 0.5517, 0.5788, 0.6028, 0.6252, 0.646, 0.67, 
                              0.694, 0.7132, 0.7372, 0.7596, 0.7788, 0.798, 0.8171, 0.8411, 
                              0.8755, 0.9019, 0.9243, 0.9435, 0.9631, 0.9811, 0.9894, 0.9857, 
                              0.9852, 0.9852, 0.9852, 0.9788, 0.9797, 0.9853, 0.9869, 0.9885, 
                              0.9901, 0.9886, 0.9896, 0.9918, 0.9854, 0.9854, 0.9804, 0.9784, 
                              0.9796, 0.9791, 0.9799, 0.9832, 0.984, 0.9823, 0.9806, 0.9828, 
                              0.9809, 0.9825, 0.9859, 0.9806, 0.981, 0.9818, 0.9816, 0.986, 
                              0.9819, 0.9796, 0.9741, 0.9724, 0.9706, 0.9608, 0.9325, 0.9072, 
                              0.875, 0.839, 0.8166, 0.751, 0.7702, 0.6919, 0.7094, 0.6183, 
                              0.512, 0.4376, 0.2596, 0.3233, 0.1908, 0.145, 0.1463, 0.1393};
                            
    G4double test_energy[] = {opticalprops::optPhotMinE_,
                            opticalprops::optPhotMinE_ +(1*eV), 
                            opticalprops::optPhotMaxE_ -(1*eV),
                            opticalprops::optPhotMaxE_};
    G4double test_reflectivity[] = {1., 1., 1., 1.};

    mpt->AddProperty("REFLECTIVITY", energy, reflectivity, entries);

    // From geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html#optical-photon-processes
    // The specular lobe constant (material property name SPECULARLOBECONSTANT) represents the reflection probability 
    // about the normal of a micro facet. The specular spike constant (material property name SPECULARSPIKECONSTANT), 
    // in turn, illustrates the probability of reflection about the average surface normal. The diffuse lobe constant 
    // is for the probability of internal Lambertian reflection, and finally the back-scatter spike constant (material 
    // property name BACKSCATTERCONSTANT) is for the case of several reflections within a deep groove with the ultimate 
    // result of exact back-scattering. The four probabilities add up to one, with the diffuse lobe constant being calculated 
    // by the code from other other three values that the user entered.
    
    G4double specularlobe[]     = {0., 0., 0., 0.};
    G4double specularspike[]    = {0., 0., 0., 0.};
    G4double backscatter[]      = {0., 0., 0., 0.};
    G4double efficiency[]       = {0., 0., 0., 0.};

    // By default, transmission equals to 0, so no need to explictly set it.
    // "The material properties REFLECTIVITY and TRANSMISSION are used. 
    // By default, REFLECTIVITY equals 1 and TRANSMISSION equals 0."

    mpt->AddProperty("SPECULARLOBECONSTANT", energy_span, specularlobe, energy_span_entries);
    mpt->AddProperty("SPECULARSPIKECONSTANT", energy_span, specularspike, energy_span_entries);
    mpt->AddProperty("BACKSCATTERCONSTANT", energy_span, backscatter, energy_span_entries);
    mpt->AddProperty("EFFICIENCY", energy_span, efficiency, energy_span_entries);
    
    return mpt;
  }

}
