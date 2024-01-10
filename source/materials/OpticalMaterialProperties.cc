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
#include "NumericalUtils.h"

#include <G4MaterialPropertiesTable.hh>


#include <map>

using namespace nexus;
using namespace CLHEP;

namespace opticalprops {

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




  /// Air /// < For the time being, the code for this one is copied from
  // Vacuum(). We should look for further information on the air optical
  // properties and tune this. 
  G4MaterialPropertiesTable* Air()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();


    // REFRACTIVE INDEX
    std::vector<G4double> photEnergy_rindex = { h_Planck*c_light/(1681.7055*nm), h_Planck*c_light/(1649.4682*nm), h_Planck*c_light/(1610.0732*nm), 
                                                h_Planck*c_light/(1570.6809*nm), h_Planck*c_light/(1531.2798*nm), h_Planck*c_light/(1491.8881*nm), 
                                                h_Planck*c_light/(1452.4924*nm), h_Planck*c_light/(1413.0967*nm), h_Planck*c_light/(1373.705*nm), 
                                                h_Planck*c_light/(1334.31*nm), h_Planck*c_light/(1294.9176*nm), h_Planck*c_light/(1255.5279*nm), 
                                                h_Planck*c_light/(1216.1356*nm), h_Planck*c_light/(1176.7459*nm), h_Planck*c_light/(1137.359*nm), 
                                                h_Planck*c_light/(1097.972*nm), h_Planck*c_light/(1058.5863*nm), h_Planck*c_light/(1019.2041*nm), 
                                                h_Planck*c_light/(979.8211*nm), h_Planck*c_light/(940.4435*nm), h_Planck*c_light/(901.068*nm), 
                                                h_Planck*c_light/(861.6971*nm), h_Planck*c_light/(822.3296*nm), h_Planck*c_light/(782.9674*nm), 
                                                h_Planck*c_light/(743.612*nm), h_Planck*c_light/(704.2647*nm), h_Planck*c_light/(664.9287*nm), 
                                                h_Planck*c_light/(625.6028*nm), h_Planck*c_light/(586.2964*nm), h_Planck*c_light/(547.0115*nm), 
                                                h_Planck*c_light/(507.7601*nm), h_Planck*c_light/(472.1026*nm), h_Planck*c_light/(443.6127*nm), 
                                                h_Planck*c_light/(420.4943*nm), h_Planck*c_light/(400.9555*nm), h_Planck*c_light/(383.2187*nm), 
                                                h_Planck*c_light/(367.2825*nm), h_Planck*c_light/(353.1358*nm), h_Planck*c_light/(340.7813*nm), 
                                                h_Planck*c_light/(330.2152*nm), h_Planck*c_light/(319.664*nm), h_Planck*c_light/(310.895*nm), 
                                                h_Planck*c_light/(303.89*nm), h_Planck*c_light/(296.896*nm), h_Planck*c_light/(289.9242*nm), 
                                                h_Planck*c_light/(282.9709*nm), h_Planck*c_light/(276.0397*nm), h_Planck*c_light/(269.127*nm), 
                                                h_Planck*c_light/(262.7743*nm), h_Planck*c_light/(258.7874*nm), h_Planck*c_light/(255.368*nm), 
                                                h_Planck*c_light/(251.9522*nm), h_Planck*c_light/(248.5697*nm)};

    
    std::vector<G4double> rIndex = {  1.000273071, 1.000273075, 1.000273116, 1.000273175, 1.000273179, 1.000273241, 1.000273279, 1.000273316, 
                                      1.000273379, 1.00027342, 1.000273478, 1.000273553, 1.000273612, 1.000273687, 1.000273778, 1.00027387, 
                                      1.00027397, 1.000274091, 1.000274207, 1.000274358, 1.00027452, 1.000274712, 1.000274924, 1.00027517, 
                                      1.000275458, 1.000275796, 1.000276205, 1.000276676, 1.000277268, 1.000277994, 1.000278928, 1.000279936, 
                                      1.000280966, 1.000281988, 1.000282995, 1.00028407, 1.000285205, 1.000286331, 1.000287464, 1.000288581, 
                                      1.000289791, 1.000290946, 1.000291934, 1.00029299, 1.000294183, 1.000295492, 1.000296938, 1.000298498, 
                                      1.00029983, 1.000301023, 1.000302033, 1.000303065, 1.000304305};

    mpt->AddProperty("RINDEX", photEnergy_rindex, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> photEnergy_abslength = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", photEnergy_abslength, absLength);

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
    G4double um2 = micrometer*micrometer;
    G4double B[3] = {4.73e-1, 6.31e-1, 9.06e-1};
    G4double C[3] = {1.30e-2 * um2, 4.13e-3 * um2, 9.88e+1 * um2};
    SellmeierEquation seq(B, C);


    const G4int ri_entries = 200;
    G4double eWidth = (optPhotFusedSilicaMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      rIndex.push_back(seq.RefractiveIndex(hc_/ri_energy[i]));
    }
    ri_energy.push_back(optPhotMaxE_);          // This sets the refractive index between optPhotFusedSilicaMaxE_ and
    rIndex.push_back(rIndex[rIndex.size()-1]);  // optPhotMaxE_ to the value obtained at optPhotFusedSilicaMaxE_

    // for (unsigned int i=0; i<ri_energy.size(); i++) {
    // G4cout << "* FusedSilica rIndex:  " << std::setw(5) << ri_energy[i]/eV
    //       << " eV -> " << rIndex[i] << G4endl;
    // }
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



  G4MaterialPropertiesTable* FakeFusedSilica(G4double transparency,
                                             G4double thickness)
  {
    // Optical properties of Suprasil 311/312(c) synthetic fused silica.
    // Obtained from http://heraeus-quarzglas.com

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    G4MaterialPropertiesTable* fused_sil_pt = opticalprops::FusedSilica();
    mpt->AddProperty("RINDEX", fused_sil_pt->GetProperty("RINDEX"));

    // ABSORPTION LENGTH (Set to match the transparency)
    G4double abs_length     = -thickness / log(transparency);
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> abs_l      = {abs_length, abs_length};
    mpt->AddProperty("ABSLENGTH", abs_energy, abs_l);

    return mpt;
  }


  G4MaterialPropertiesTable* Epoxy()
  {
    // This is the material used as a window for NEXT-100 SiPMs.

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    // The only information we have is that n = 1.55 at a
    // not specified wavelength. It seems that it's common to measure
    // the refractive index at the D line of sodium (590 nm),
    // therefore we assume that.
    // The dependence of n on the wavelength is taken from
    // https://www.epotek.com/docs/en/Related/Tech%20Tip%2018%20Understanding%20Optical%20Properties%20of%20Epoxy%20Applications.pdf,
    // shifting the whole graphs in y down to match 1.55 at 590 nm.
    // We fill the values outside the range of the plot with the
    // minimum and maximum values.
     std::vector<G4double> ri_energy = {
      optPhotMinE_,
      hc_ / (2451.63 * nm), hc_ / (2430.88 * nm), hc_ / (2405.51 * nm),
      hc_ / (2380.14 * nm), hc_ / (2354.77 * nm), hc_ / (2329.41 * nm),
      hc_ / (2304.04 * nm), hc_ / (2278.67 * nm), hc_ / (2253.3 * nm),
      hc_ / (2227.94 * nm), hc_ / (2202.57 * nm), hc_ / (2177.2 * nm),
      hc_ / (2151.83 * nm), hc_ / (2126.47 * nm), hc_ / (2101.1 * nm),
      hc_ / (2075.73 * nm), hc_ / (2050.36 * nm), hc_ / (2025.0 * nm),
      hc_ / (1968.5 * nm),  hc_ / (1951.2 * nm),  hc_ / (1925.83 * nm),
      hc_ / (1900.47 * nm), hc_ / (1875.1 * nm),  hc_ / (1849.73 * nm),
      hc_ / (1824.36 * nm), hc_ / (1799.0 * nm),  hc_ / (1773.63 * nm),
      hc_ / (1747.6 * nm),  hc_ / (1722.89 * nm), hc_ / (1697.53 * nm),
      hc_ / (1672.16 * nm), hc_ / (1646.79 * nm), hc_ / (1621.42 * nm),
      hc_ / (1596.06 * nm), hc_ / (1570.69 * nm), hc_ / (1545.32 * nm),
      hc_ / (1519.95 * nm), hc_ / (1494.59 * nm), hc_ / (1469.22 * nm),
      hc_ / (1443.85 * nm), hc_ / (1418.49 * nm), hc_ / (1393.12 * nm),
      hc_ / (1367.75 * nm), hc_ / (1342.38 * nm), hc_ / (1317.02 * nm),
      hc_ / (1291.65 * nm), hc_ / (1266.28 * nm), hc_ / (1240.91 * nm),
      hc_ / (1215.55 * nm), hc_ / (1190.18 * nm), hc_ / (1164.81 * nm),
      hc_ / (1139.44 * nm), hc_ / (1114.08 * nm), hc_ / (1088.71 * nm),
      hc_ / (1063.34 * nm), hc_ / (1037.97 * nm), hc_ / (1012.61 * nm),
      hc_ / (987.24 * nm),  hc_ / (961.87 * nm),  hc_ / (936.5 * nm),
      hc_ / (911.14 * nm),  hc_ / (885.77 * nm),  hc_ / (860.4 * nm),
      hc_ / (835.03 * nm),  hc_ / (809.67 * nm),  hc_ / (784.3 * nm),
      hc_ / (758.93 * nm),  hc_ / (733.57 * nm),  hc_ / (708.2 * nm),
      hc_ / (682.83 * nm),  hc_ / (657.46 * nm),  hc_ / (624.02 * nm),
      hc_ / (606.73 * nm),  hc_ / (587.13 * nm),  hc_ / (569.83 * nm),
      hc_ / (554.84 * nm),  hc_ / (541.0 * nm),   hc_ / (519.1 * nm),
      hc_ / (509.87 * nm),  hc_ / (499.49 * nm),  hc_ / (490.27 * nm),
      hc_ / (481.04 * nm),  hc_ / (470.67 * nm),  hc_ / (456.83 * nm),
      hc_ / (451.06 * nm),  hc_ / (442.99 * nm),  hc_ / (434.92 * nm),
      hc_ / (426.85 * nm),  hc_ / (417.63 * nm),  hc_ / (401.48 * nm),
      hc_ / (395.95 * nm),
      optPhotMaxE_
     };

     std::vector<G4double> rIndex = {
       1.524,
       1.524, 1.525, 1.525,
       1.525, 1.525, 1.526,
       1.526, 1.526, 1.526,
       1.526, 1.527, 1.527,
       1.527, 1.527, 1.527,
       1.528, 1.528, 1.528,
       1.528, 1.528, 1.528,
       1.529, 1.529, 1.529,
       1.529, 1.529, 1.53,
       1.53,  1.53,  1.53,
       1.531, 1.531, 1.531,
       1.531, 1.532, 1.532,
       1.532, 1.532, 1.533,
       1.533, 1.533, 1.533,
       1.534, 1.534, 1.534,
       1.534, 1.535, 1.535,
       1.535, 1.535, 1.536,
       1.536, 1.536, 1.536,
       1.536, 1.537, 1.537,
       1.537, 1.538, 1.538,
       1.539, 1.539, 1.54,
       1.54,  1.541, 1.542,
       1.543, 1.544, 1.545,
       1.546, 1.547, 1.549,
       1.55,  1.551, 1.552,
       1.554, 1.555, 1.558,
       1.559, 1.56,  1.562,
       1.563, 1.565, 1.567,
       1.568, 1.569, 1.57,
       1.572, 1.573, 1.576,
       1.577,
       1.577
     };

    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    // We don't have information about the absorption length,
    // therefore we use that of GlassEpoxy since it is virtually infinite
    // for the TPB emission wavelengths.
    std::vector<G4double> abs_energy = {
      optPhotMinE_, 2.000 * eV,
      2.132 * eV,   2.735 * eV,  2.908 * eV,  3.119 * eV,
      3.320 * eV,   3.476 * eV,  3.588 * eV,  3.749 * eV,
      3.869 * eV,   3.973 * eV,  4.120 * eV,  4.224 * eV,
      4.320 * eV,   4.420 * eV,  5.018 * eV,
      optPhotMaxE_
    };

    std::vector<G4double> absLength = {
      noAbsLength_, noAbsLength_,
      326.00 * mm,  117.68 * mm,  85.89 * mm,  50.93 * mm,
      31.25 * mm,   17.19 * mm,  10.46 * mm,   5.26 * mm,
       3.77 * mm,    2.69 * mm,   1.94 * mm,   1.33 * mm,
       0.73 * mm,    0.32 * mm,   0.10 * mm,
       0.10 * mm
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }



  /// Imperfect Dielectric-dielectric surface ///
  G4MaterialPropertiesTable* ImperfectDielectricDielectricSurface(G4double tunneling_probability)
  {
    // This surface has a certain probability for photons, regarless their energy, to be 
    // straightforward transmitted without being fresnel-refracted or -reflected. Such 
    // probability is the one that is given to the tunneling_probability parameter. This 
    // is based on the following chunk of documentation of Geant4, which can be found in:
    // https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html#optical-photon-processes
    // 
    // "It is possible to specify that a given fraction of photons are absorbed at the surface, 
    // or transmitted without change in direction or polarization. This is applicable for 
    // dielectric_dielectric interfaces that are not backpainted. The material properties 
    // REFLECTIVITY and TRANSMITTANCE are used. By default, REFLECTIVITY equals 1 and 
    // TRANSMITTANCE equals 0. At a surface interaction, a random number is chosen. If the 
    // random number is greater than the sum of the values of REFLECTIVITY and TRANSMITTANCE 
    // at the photon energy, the photon is absorbed. Otherwise, if the random number is greater 
    // than the REFLECTIVITY value, the photon is transmitted. Otherwise, the usual calculation 
    // of scattering takes place."
    
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    G4double energy[] =         { optPhotMinE_,             optPhotMinE_+(0.1*eV),    optPhotMinE_+(0.2*eV),    optPhotMinE_+(0.3*eV), 
                                  optPhotMaxE_-(0.3*eV),    optPhotMaxE_-(0.2*eV),    optPhotMaxE_-(0.1*eV),    optPhotMaxE_            };
    G4double transmittance[] =  { tunneling_probability,    tunneling_probability,    tunneling_probability,    tunneling_probability,
                                  tunneling_probability,    tunneling_probability,    tunneling_probability,    tunneling_probability    };
    G4double reflectivity[] =   { 1.-tunneling_probability, 1.-tunneling_probability, 1.-tunneling_probability, 1.-tunneling_probability,
                                  1.-tunneling_probability, 1.-tunneling_probability, 1.-tunneling_probability, 1.-tunneling_probability };

    mpt->AddProperty("TRANSMITTANCE", energy, transmittance,  8);
    mpt->AddProperty("REFLECTIVITY",  energy, reflectivity,   8);
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
      hc_ / (1000. * nm), hc_ / (800. * nm), hc_ / (700. * nm),
      hc_ / ( 600. * nm), hc_ / (580. * nm), hc_ / (560. * nm),
      hc_ / ( 540. * nm), hc_ / (520. * nm), hc_ / (500. * nm),
      hc_ / ( 480. * nm), hc_ / (460. * nm), hc_ / (440. * nm),
      hc_ / ( 420. * nm), hc_ / (400. * nm),
      optPhotMaxE_ };

    std::vector<G4double> rIndex = {
      1.635,
      1.635, 1.775, 1.835,
      1.894, 1.906, 1.919,
      1.931, 1.945, 1.960,
      1.975, 1.993, 2.012,
      2.036, 2.064,
      2.064 };
    mpt->AddProperty("RINDEX", energies, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_length = {
      (1000. * nm) / (4*pi * 0.0103),
      (1000. * nm) / (4*pi * 0.0103), (800. * nm) / (4*pi * 0.0049),
      ( 700. * nm) / (4*pi * 0.0033), (600. * nm) / (4*pi * 0.0023),
      ( 580. * nm) / (4*pi * 0.0022), (560. * nm) / (4*pi * 0.0022),
      ( 540. * nm) / (4*pi * 0.0022), (520. * nm) / (4*pi * 0.0023),
      ( 500. * nm) / (4*pi * 0.0026), (480. * nm) / (4*pi * 0.0031),
      ( 460. * nm) / (4*pi * 0.0039), (440. * nm) / (4*pi * 0.0053),
      ( 420. * nm) / (4*pi * 0.0080), (400. * nm) / (4*pi * 0.0125),
      ( 400. * nm) / (4*pi * 0.0125) };
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



  G4MaterialPropertiesTable* PEDOT()
  {
    // Input data: complex refraction index obtained from:
    // https://refractiveindex.info/?shelf=other&book=PEDOT-PSS&page=Chen
    // Only valid in [1097 - 302] nm

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> energies = {
      optPhotMinE_,
      hc_ / (1097. * nm),  hc_ / (1000. * nm),  hc_ / ( 950. * nm),
      hc_ / ( 900. * nm),  hc_ / ( 800. * nm),  hc_ / ( 700. * nm),
      hc_ / ( 600. * nm),  hc_ / ( 550. * nm),  hc_ / ( 500. * nm),
      hc_ / ( 450. * nm),  hc_ / ( 420. * nm),  hc_ / ( 400. * nm),
      hc_ / ( 370. * nm),  hc_ / ( 350. * nm),  hc_ / ( 302. * nm),
      optPhotMaxE_ };

    std::vector<G4double> rIndex = {
      1.4760,
      1.4760, 1.4662, 1.4665,
      1.4693, 1.4802, 1.4935,
      1.5080, 1.5155, 1.5235,
      1.5328, 1.5391, 1.5439,
      1.5522, 1.5587, 1.5805,
      1.5805 };
    mpt->AddProperty("RINDEX", energies, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_length = {
      (1097. * nm) / (4*pi * 0.1191),
      (1097. * nm) / (4*pi * 0.1191),  (1000. * nm) / (4*pi * 0.0859),
      ( 950. * nm) / (4*pi * 0.0701),  ( 900. * nm) / (4*pi * 0.0561),
      ( 800. * nm) / (4*pi * 0.0340),  ( 700. * nm) / (4*pi * 0.0197),
      ( 600. * nm) / (4*pi * 0.0107),  ( 550. * nm) / (4*pi * 0.0076),
      ( 500. * nm) / (4*pi * 0.0051),  ( 450. * nm) / (4*pi * 0.0035),
      ( 420. * nm) / (4*pi * 0.0025),  ( 400. * nm) / (4*pi * 0.00194),
      ( 370. * nm) / (4*pi * 0.00135), ( 350. * nm) / (4*pi * 0.00103),
      ( 302. * nm) / (4*pi * 0.0004),  ( 302. * nm) / (4*pi * 0.0004) };
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



    /// Non-absorbent and r. index-tunable material ///
  G4MaterialPropertiesTable* TunableRIMat(G4double refractive_index)
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // Refractive index (RINDEX)
    G4double energies_rindex[] =    { optPhotMinE_,           optPhotMinE_+(0.1*eV),  optPhotMinE_+(0.2*eV),  optPhotMinE_+(0.3*eV), 
                                      optPhotMaxE_-(0.3*eV),  optPhotMaxE_-(0.2*eV),  optPhotMaxE_-(0.1*eV),  optPhotMaxE_            };
    G4double rindex[] =             { refractive_index,       refractive_index,       refractive_index,       refractive_index,
                                      refractive_index,       refractive_index,       refractive_index,       refractive_index        };
    mpt->AddProperty("RINDEX", energies_rindex, rindex, 8);

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
      G4double lambda = hc_/ri_energy[i]*1000; // in micron
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
    G4double rindex[] =       {1.55        ,  1.55        };  // The windows for Hamamatsu S13360-6050VE  
                                                              // and Hamamatsu S13360-5075HD-HQR
                                                              // are made out of epoxy resin. The r. index
                                                              // value for the epoxy resin of HS13360-6050VE 
                                                              // is 1.55. I am assuming this one for the 
                                                              // S13360-5075HD-HQR SiPM as well. Also, we know
                                                              // the FBK_NUV-HD-CRYO triple trench window
                                                              // is made out of epoxy resin. Since this is
                                                              // everything we know, we will also use this
                                                              // MPT for the window of the FBK SiPM.
    G4double abslength[] =    {noAbsLength_,  noAbsLength_};

    mpt->AddProperty("RINDEX",    energy_span, rindex,    2);
    mpt->AddProperty("ABSLENGTH", energy_span, abslength, 2);

    return mpt;
  }


  G4MaterialPropertiesTable* Sapphire()
  {
    // Input data: Sellmeier equation coeficients extracted from:
    // https://refractiveindex.info/?shelf=3d&book=crystals&page=sapphire
    // C[i] coeficients at line 362 are squared.

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    G4double um2 = micrometer*micrometer;
    G4double B[3] = {1.4313493, 0.65054713, 5.3414021};
    G4double C[3] = {0.0052799261 * um2, 0.0142382647 * um2, 325.017834 * um2};
    SellmeierEquation seq(B, C);

    const G4int ri_entries = 100;
    G4double eWidth = (optPhotSapphireMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      rIndex.push_back(seq.RefractiveIndex(hc_/ri_energy[i]));
    }
    // This sets the refractive index between optPhotSapphireMaxE_ and
    // optPhotMaxE_ to the value obtained at optPhotSapphireMaxE_
    ri_energy.push_back(optPhotMaxE_);
    rIndex.push_back(rIndex[rIndex.size()-1]);
    // for (unsigned int i=0; i<ri_energy.size(); i++) {
    //   G4cout << "* Sapphire rIndex:  " << std::setw(5)
    //          << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    // }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_, 0.900 * eV,
      1.000 * eV,   1.296 * eV,  1.683 * eV,  2.075 * eV,
      2.585 * eV,   3.088 * eV,  3.709 * eV,  4.385 * eV,
      4.972 * eV,   5.608 * eV,  6.066 * eV,  6.426 * eV,
      6.806 * eV,   7.135 * eV,  7.401 * eV,  7.637 * eV,
      7.880 * eV,   8.217 * eV,
      optPhotMaxE_
    };

    std::vector<G4double> absLength = {
      noAbsLength_, noAbsLength_,
      3455.0  * mm,  3455.0  * mm,  3455.0  * mm,  3455.0  * mm,
      3455.0  * mm,  3140.98 * mm,  2283.30 * mm,  1742.11 * mm,
      437.06 * mm,   219.24 * mm,  117.773 * mm,   80.560 * mm,
      48.071 * mm,   28.805 * mm,   17.880 * mm,   11.567 * mm,
      7.718 * mm,    4.995 * mm,
      4.995 * mm
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    return mpt;
  }



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
      G4double wl = hc_ / ri_energy[i];
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



  /// Gaseous argon ///
  G4MaterialPropertiesTable* GAr(G4double sc_yield,
                                G4double e_lifetime)
  {
    // An argon gas proportional scintillation counter with UV avalanche photodiode scintillation
    // readout C.M.B. Monteiro, J.A.M. Lopes, P.C.P.S. Simoes, J.M.F. dos Santos, C.A.N. Conde
    //
    // May 2023:
    // Updated scintillation decay and yields from:
    // Triplet Lifetime in Gaseous Argon. Michael Akashi-Ronquest et al.

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
      G4double wl = hc_ / ri_energy[i] * 1000; // in micron
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
    G4double Energy_peak  = (hc_ / Wavelength_peak);
    G4double Energy_sigma = (hc_ * Wavelength_sigma / pow(Wavelength_peak,2));
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
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",   3480.*ns);
    mpt->AddConstProperty("SCINTILLATIONYIELD1", .136);
    mpt->AddConstProperty("SCINTILLATIONYIELD2", .864);
    mpt->AddConstProperty("RESOLUTIONSCALE",    1.0);
    mpt->AddConstProperty("ATTACHMENT",         e_lifetime, 1);

    return mpt;
  }



  /// Gaseous xenon ///
  G4MaterialPropertiesTable* GXe(G4double pressure,
                                 G4double /*temperature*/,
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



  /// Liquid xenon ///
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
    //   G4cout << hc_/ri_energy[i]/nanometer << " nm, " << rindex[i] << G4endl;
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

  G4MaterialPropertiesTable* PTP(G4double refractive_index){
    // This material is meant to model p-Terphenyl, which is a wavelength shifter commonly used to WLS the 
    // LAr VUV scintillation light. The emission spectrum was taken from
    // sciencedirect.com/science/article/abs/pii/016890029390701I ,
    // whereas the PLQY (photo luminiscence quantum yield = #photons emitted/#photons absorbed), the rindex
    // and the WLS delay time was taken from 
    // mdpi-res.com/d_attachment/instruments/instruments-05-00004/article_deploy/instruments-05-00004-v2.pdf?version=1609810328
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // Refractive index (RINDEX)
    G4double energies_rindex[] =    {optPhotMinE_, optPhotMaxE_};
    G4double rindex[] =             {refractive_index, refractive_index};
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
    // This data was digitized from the 100 ºK curve of Fig. 1 of 
    // 'FLUORESCENCE AND PHASE TRANSITION OF p-TERPHENYL CRYSTAL' by Nobuko I. WAKAYAMA.
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_, 
      h_Planck*c_light/(437.0000*nm), h_Planck*c_light/(436.0000*nm), h_Planck*c_light/(435.0000*nm),
      h_Planck*c_light/(434.58*nm), h_Planck*c_light/(430.84*nm), h_Planck*c_light/(426.96*nm), h_Planck*c_light/(423.91*nm), 
      h_Planck*c_light/(420.87*nm), h_Planck*c_light/(418.1*nm), h_Planck*c_light/(415.05*nm), h_Planck*c_light/(412.42*nm), 
      h_Planck*c_light/(408.95*nm), h_Planck*c_light/(406.04*nm), h_Planck*c_light/(403.68*nm), h_Planck*c_light/(402.44*nm), 
      h_Planck*c_light/(401.21*nm), h_Planck*c_light/(399.7*nm), h_Planck*c_light/(398.18*nm), h_Planck*c_light/(397.22*nm), 
      h_Planck*c_light/(396.41*nm), h_Planck*c_light/(395.32*nm), h_Planck*c_light/(394.64*nm), h_Planck*c_light/(393.4*nm), 
      h_Planck*c_light/(391.74*nm), h_Planck*c_light/(390.64*nm), h_Planck*c_light/(389.82*nm), h_Planck*c_light/(388.57*nm), 
      h_Planck*c_light/(387.72*nm), h_Planck*c_light/(386.46*nm), h_Planck*c_light/(385.48*nm), h_Planck*c_light/(384.23*nm), 
      h_Planck*c_light/(382.83*nm), h_Planck*c_light/(381.72*nm), h_Planck*c_light/(380.9*nm), h_Planck*c_light/(379.95*nm), 
      h_Planck*c_light/(379.27*nm), h_Planck*c_light/(378.73*nm), h_Planck*c_light/(378.2*nm), h_Planck*c_light/(377.67*nm), 
      h_Planck*c_light/(377.14*nm), h_Planck*c_light/(376.49*nm), h_Planck*c_light/(375.83*nm), h_Planck*c_light/(374.74*nm), 
      h_Planck*c_light/(373.64*nm), h_Planck*c_light/(373.11*nm), h_Planck*c_light/(372.58*nm), h_Planck*c_light/(372.18*nm), 
      h_Planck*c_light/(371.8*nm), h_Planck*c_light/(371.54*nm), h_Planck*c_light/(370.44*nm), h_Planck*c_light/(369.86*nm), 
      h_Planck*c_light/(369.42*nm), h_Planck*c_light/(368.85*nm), h_Planck*c_light/(367.88*nm), h_Planck*c_light/(366.89*nm), 
      h_Planck*c_light/(366.47*nm), h_Planck*c_light/(366.17*nm), h_Planck*c_light/(365.72*nm), h_Planck*c_light/(365.28*nm), 
      h_Planck*c_light/(364.7*nm), h_Planck*c_light/(363.84*nm), h_Planck*c_light/(363.26*nm), h_Planck*c_light/(362.68*nm), 
      h_Planck*c_light/(362.12*nm), h_Planck*c_light/(361.41*nm), h_Planck*c_light/(360.3*nm), h_Planck*c_light/(359.9*nm), 
      h_Planck*c_light/(359.5*nm), h_Planck*c_light/(359.11*nm), h_Planck*c_light/(358.59*nm), h_Planck*c_light/(358.2*nm), 
      h_Planck*c_light/(357.95*nm), h_Planck*c_light/(357.29*nm), h_Planck*c_light/(356.91*nm), h_Planck*c_light/(356.54*nm), 
      h_Planck*c_light/(356.33*nm), h_Planck*c_light/(356.1*nm), h_Planck*c_light/(355.72*nm), h_Planck*c_light/(355.61*nm), 
      h_Planck*c_light/(355.23*nm), h_Planck*c_light/(354.56*nm), h_Planck*c_light/(354.26*nm), h_Planck*c_light/(353.83*nm), 
      h_Planck*c_light/(353.53*nm), h_Planck*c_light/(353.1*nm), h_Planck*c_light/(352.4*nm), h_Planck*c_light/(352.13*nm), 
      h_Planck*c_light/(351.72*nm), h_Planck*c_light/(351.03*nm), h_Planck*c_light/(350.74*nm), h_Planck*c_light/(350.43*nm), 
      h_Planck*c_light/(349.98*nm), h_Planck*c_light/(349.51*nm), h_Planck*c_light/(349.06*nm), h_Planck*c_light/(348.73*nm), 
      h_Planck*c_light/(348.28*nm), h_Planck*c_light/(347.95*nm), h_Planck*c_light/(347.49*nm), h_Planck*c_light/(347.04*nm), 
      h_Planck*c_light/(346.59*nm), h_Planck*c_light/(346.01*nm), h_Planck*c_light/(345.54*nm), h_Planck*c_light/(344.96*nm), 
      h_Planck*c_light/(344.51*nm), h_Planck*c_light/(343.79*nm), h_Planck*c_light/(343.49*nm), h_Planck*c_light/(342.92*nm), 
      h_Planck*c_light/(342.21*nm), h_Planck*c_light/(340.53*nm), h_Planck*c_light/(339.41*nm), h_Planck*c_light/(338.7*nm), 
      h_Planck*c_light/(338.26*nm), h_Planck*c_light/(337.96*nm), h_Planck*c_light/(337.52*nm), h_Planck*c_light/(336.96*nm), 
      h_Planck*c_light/(336.25*nm), h_Planck*c_light/(335.54*nm), h_Planck*c_light/(334.56*nm), h_Planck*c_light/(333.31*nm), 
      h_Planck*c_light/(331.5*nm), h_Planck*c_light/(329.84*nm), h_Planck*c_light/(328.31*nm),
      h_Planck*c_light/(327.0000*nm), h_Planck*c_light/(326.0000*nm), h_Planck*c_light/(325.0000*nm),
      optPhotMaxE_
    };
    std::vector<G4double> WLS_emiSpectrum = {
      0.0, 
      0.0, 0.0, 0.0,
      0.03085, 0.03413, 0.04332, 0.05612, 0.06892, 0.08055, 0.09335, 0.10026, 0.10238, 0.10808, 0.1162, 0.12913, 0.14561, 
      0.16442, 0.18207, 0.19856, 0.21743, 0.2422, 0.26463, 0.2811, 0.28572, 0.29865, 0.31042, 0.31034, 0.29607, 0.28415, 
      0.27581, 0.26743, 0.25905, 0.26253, 0.27785, 0.30026, 0.31559, 0.33448, 0.36993, 0.40539, 0.43611, 0.49048, 0.52948, 
      0.55187, 0.56954, 0.59671, 0.62034, 0.64871, 0.69246, 0.709, 0.71957, 0.69113, 0.66389, 0.65438, 0.6484, 0.63532, 
      0.62346, 0.59859, 0.56425, 0.53346, 0.50384, 0.46947, 0.43986, 0.41734, 0.40192, 0.38885, 0.39233, 0.40531, 0.42776, 
      0.45612, 0.49986, 0.54242, 0.57552, 0.61334, 0.65827, 0.7174, 0.79783, 0.84514, 0.89125, 0.9291, 0.97403, 1.0, 
      0.97869, 0.9621, 0.93723, 0.91827, 0.90521, 0.91584, 0.92646, 0.93113, 0.91928, 0.88377, 0.83524, 0.77723, 0.73224, 
      0.67425, 0.6269, 0.56654, 0.51919, 0.47656, 0.43512, 0.40432, 0.3475, 0.31552, 0.27763, 0.248, 0.22077, 0.19943, 
      0.17572, 0.16259, 0.1495, 0.13407, 0.11393, 0.07959, 0.05945, 0.04285, 0.02624, 0.01317, 0.00482, 0.0, 0.00105, 0.00094, 0.00083,
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
    std::vector<G4double> reflectivity =    {.98             ,   .98          };
    std::vector<G4double> rindex =          {1.             ,   1.          };

    mpt->AddProperty("EFFICIENCY",      energy.data(), efficiency.data(),   energy.size());
    mpt->AddProperty("REFLECTIVITY",    energy.data(), reflectivity.data(), energy.size());
    mpt->AddProperty("RINDEX",          energy.data(), rindex.data(),       energy.size());

    return mpt;
  }

  G4MaterialPropertiesTable* PerfectPolishedSurfaceTransmitter()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> energy =          {optPhotMinE_   ,   optPhotMaxE_};
    std::vector<G4double> reflectivity =    {0.             ,   0.          };
    std::vector<G4double> transmission =    {1.             ,   1.          };

    mpt->AddProperty("REFLECTIVITY",    energy.data(), reflectivity.data(), energy.size());
    mpt->AddProperty("TRANSMITTANCE",    energy.data(), transmission.data(), energy.size());

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



  G4MaterialPropertiesTable* PolishedAl()
  {
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    std::vector<G4double> ENERGIES = {
       hc_ / (2456.42541 * nm), hc_ / (2396.60266 * nm), hc_ / (2276.95716 * nm),
       hc_ / (2159.52733 * nm), hc_ / (2037.66617 * nm), hc_ / (1918.02068 * nm),
       hc_ / (1798.37518 * nm), hc_ / (1676.51403 * nm), hc_ / (1559.08419 * nm),
       hc_ / (1437.22304 * nm), hc_ / (1319.79321 * nm), hc_ / (1197.93205 * nm),
       hc_ / (1078.28656 * nm), hc_ / (956.42541 * nm),  hc_ / (838.99557 * nm),
       hc_ / (717.13442 * nm),  hc_ / (597.48892 * nm),  hc_ / (477.84343 * nm),
       hc_ / (418.02068 * nm),  hc_ / (358.19793 * nm),  hc_ / (293.94387 * nm)
    };
    std::vector<G4double> REFLECTIVITY = {
      .99088, .99082, .98925,
      .98623, .98611, .98163,
      .98006, .97849, .97401,
      .97098, .96941, .96784,
      .96481, .96033, .96167,
      .96301, .96289, .96278,
      .96126, .95830, .94224
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
    //  optPhotMinE_,                      hc_ / (450. * nm),
    //  hc_ / (440. * nm),  hc_ / (430. * nm),
    //  hc_ / (420. * nm),  hc_ / (410. * nm),
    //  hc_ / (400. * nm),  hc_ / (390. * nm),
    //  hc_ / (380. * nm),  hc_ / (370. * nm),
    //  hc_ / (360. * nm),  hc_ / (330. * nm),
    //  hc_ / (320. * nm),  hc_ / (310. * nm),
    //  hc_ / (300. * nm),  hc_ / (270. * nm),
    //  hc_ / (250. * nm),  hc_ / (230. * nm),
    //  hc_ / (210. * nm),  hc_ / (190. * nm),
    //  hc_ / (170. * nm),  hc_ / (150. * nm),
    //  hc_ / (100. * nm),  optPhotMaxE_
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
      hc_ / (380. * nm),  hc_ / (370. * nm), hc_ / (360. * nm),
      hc_ / (330. * nm),  hc_ / (320. * nm), hc_ / (310. * nm),
      hc_ / (300. * nm),  hc_ / (270. * nm), hc_ / (250. * nm),
      hc_ / (230. * nm),  hc_ / (210. * nm), hc_ / (190. * nm),
      hc_ / (170. * nm),  hc_ / (150. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_,                 // ~6200 nm
      noAbsLength_, 50. * nm, 30. * nm, // 380, 370, 360 nm
      30. * nm, 50. * nm, 80. * nm,     // 330, 320, 310 nm
      100. * nm, 100. * nm, 400. * nm,  // 300, 270, 250 nm
      400. * nm, 350. * nm, 250. * nm,  // 230, 210, 190 nm
      350. * nm, 400. * nm, 400. * nm   // 170, 150, ~108 nm
    };

    //for (int i=0; i<WLS_abs_energy.size(); i++)
    //  G4cout << "* TPB WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (hc_ / WLS_abs_energy[i]) / nm
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
      G4double wl = (hc_ / WLS_emi_energy[i]) / nm;
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



  G4MaterialPropertiesTable* DegradedTPB(G4double wls_eff)
  {
    // It has all the same properties of TPB except the WaveLengthShifting probability
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
      hc_ / (337. * nm), hc_ / (318. * nm), hc_ / (292. * nm),
      hc_ / (276. * nm), hc_ / (253. * nm), hc_ / (238. * nm),
      hc_ / (222. * nm), hc_ / (204. * nm), hc_ / (190. * nm),
      hc_ / (175. * nm), hc_ / (169. * nm),
      optPhotMaxE_
    };

    float XePeakAbsValue = 1.879;
    float XePeakAbsLength = 21 * nm;

    std::vector<float> PTP_absorption = {
      0.002, 0.174, 0.414, // 337, 318, 292
      0.949, 0.540, 1.218, // 276, 253, 238
      1.858, 1.429, 1.716, // 222, 204, 190
      1.879, 1.803 // 175, 169
    };

    std::vector<G4double> WLS_absLength = {noAbsLength_};
    for (auto& abs_value : PTP_absorption)
      WLS_absLength.push_back(XePeakAbsLength / (abs_value / XePeakAbsValue));

    WLS_absLength.push_back(noAbsLength_);

    //for (int i=0; i<WLS_abs_entries; i++)
    //  G4cout << "* TPB WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (hc_ / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / nm << " nm" << G4endl;
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_,
      hc_ / (452. * nm), hc_ / (430. * nm), hc_ / (412. * nm), hc_ / (398. * nm),
      hc_ / (385. * nm), hc_ / (371. * nm), hc_ / (361. * nm), hc_ / (354. * nm),
      hc_ / (336. * nm), hc_ / (317. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.,
      0.044, 0.179, 0.351, 0.514,
      0.849, 0.993, 0.745, 0.421,
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



  G4MaterialPropertiesTable* EJ280()
  {
    // https://eljentechnology.com/products/wavelength-shifting-plastics/ej-280-ej-282-ej-284-ej-286
    // and data sheets from the provider.
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> ri_energy = {
      optPhotMinE_,
      hc_ / (609. * nm),  hc_ / (589.26 * nm), hc_ / (550. * nm),
      hc_ / (530. * nm),  hc_ / (500. * nm),   hc_ / (490. * nm),
      hc_ / (481. * nm),  hc_ / (460. * nm),   hc_ / (435. * nm),
      hc_ / (425. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> rIndex = {
      1.5780,
      1.5780, 1.5800, 1.5845,  // 609 , 589.26, 550 nm
      1.5870, 1.5913, 1.5929,  // 530, 500, 490 nm
      1.5945, 1.5986, 1.6050,  // 481, 460, 435 nm
      1.6080,                  // 425 nm
      1.608
    };
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,
      hc_ / (750. * nm), hc_ / (740. * nm), hc_ / (380. * nm), hc_ / (370. * nm),
      optPhotMaxE_
    };
    std::vector<G4double> absLength = {
      noAbsLength_,
      noAbsLength_, 3.0 * m, 3.0 * m, noAbsLength_,
      noAbsLength_
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,
      hc_ / (500. * nm), hc_ / (495. * nm), hc_ / (490. * nm),
      hc_ / (485. * nm), hc_ / (480. * nm), hc_ / (475. * nm),
      hc_ / (470. * nm), hc_ / (465. * nm), hc_ / (460. * nm),
      hc_ / (455. * nm), hc_ / (450. * nm), hc_ / (445. * nm),
      hc_ / (440. * nm), hc_ / (435. * nm), hc_ / (430. * nm),
      hc_ / (425. * nm), hc_ / (420. * nm), hc_ / (415. * nm),
      hc_ / (410. * nm), hc_ / (405. * nm), hc_ / (400. * nm),
      hc_ / (395. * nm), hc_ / (390. * nm), hc_ / (385. * nm),
      hc_ / (380. * nm), hc_ / (375. * nm), hc_ / (370. * nm),
      hc_ / (365. * nm), hc_ / (360. * nm), hc_ / (355. * nm),
      hc_ / (350. * nm), hc_ / (345. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_,
      noAbsLength_, (1. / 0.0080) * cm, (1. / 0.0165) * cm,
      (1. / 0.0325) * cm, (1. / 0.0815) * cm, (1. / 0.2940) * cm,
      (1. / 0.9640) * cm, (1. / 2.8600) * cm, (1. / 6.3900) * cm,
      (1. / 9.9700) * cm, (1. / 11.0645)* cm, (1. / 10.198) * cm,
      (1. / 9.4465) * cm, (1. / 10.2145)* cm, (1. / 12.240) * cm,
      (1. / 12.519) * cm, (1. / 10.867) * cm, (1. / 9.0710) * cm,
      (1. / 8.0895) * cm, (1. / 7.6650) * cm, (1. / 6.7170) * cm,
      (1. / 5.2460) * cm, (1. / 4.1185) * cm, (1. / 3.3175) * cm,
      (1. / 2.6800) * cm, (1. / 1.9610) * cm, (1. / 1.4220) * cm,
      (1. / 1.0295) * cm, (1. / 0.7680) * cm, (1. / 0.6865) * cm,
      (1. / 0.5885) * cm, noAbsLength_,
      noAbsLength_
    };
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);
    //for (int i=0; i<WLS_abs_entries; i++)
    //  G4cout << "* EJ280 WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (hc_ / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / mm << " mm" << G4endl;

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_,
      hc_ / (610. * nm), hc_ / (605. * nm), hc_ / (600. * nm),
      hc_ / (595. * nm), hc_ / (590. * nm), hc_ / (585. * nm),
      hc_ / (580. * nm), hc_ / (575. * nm), hc_ / (570. * nm),
      hc_ / (565. * nm), hc_ / (560. * nm), hc_ / (555. * nm),
      hc_ / (550. * nm), hc_ / (545. * nm), hc_ / (540. * nm),
      hc_ / (535. * nm), hc_ / (530. * nm), hc_ / (525. * nm),
      hc_ / (520. * nm), hc_ / (515. * nm), hc_ / (510. * nm),
      hc_ / (505. * nm), hc_ / (500. * nm), hc_ / (495. * nm),
      hc_ / (490. * nm), hc_ / (485. * nm), hc_ / (480. * nm),
      hc_ / (475. * nm), hc_ / (470. * nm), hc_ / (465. * nm),
      hc_ / (460. * nm), hc_ / (455. * nm), hc_ / (450. * nm),
      hc_ / (445. * nm), hc_ / (440. * nm), hc_ / (435. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.000,
      0.000, 0.003, 0.006,  // 610, 605, 600 nm
      0.007, 0.009, 0.014,  // 595, 590, 585 nm
      0.017, 0.024, 0.033,  // 580, 575, 570 nm
      0.042, 0.051, 0.063,  // 565, 560, 555 nm
      0.081, 0.112, 0.157,  // 550, 545, 540 nm
      0.211, 0.274, 0.329,  // 535, 530, 525 nm
      0.341, 0.325, 0.346,  // 520, 515, 510 nm
      0.433, 0.578, 0.792,  // 505, 500, 495 nm
      1.000, 0.966, 0.718,  // 490, 485, 480 nm
      0.604, 0.681, 0.708,  // 475, 470, 465 nm
      0.525, 0.242, 0.046,  // 460, 455, 450 nm
      0.012, 0.003, 0.000,  // 445, 440, 435 nm
      0.000
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

  G4MaterialPropertiesTable* EJ286()
  {
    // https://eljentechnology.com/products/wavelength-shifting-plastics/ej-280-ej-282-ej-284-ej-286
    // and data sheets from the provider.
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> ri_energy = {
      optPhotMinE_,
      hc_ / (609. * nm), hc_ / (589.26 * nm), hc_ / (550. * nm),
      hc_ / (530. * nm), hc_ / (500. * nm),   hc_ / (490. * nm),
      hc_ / (481. * nm), hc_ / (460. * nm),   hc_ / (435. * nm),
      hc_ / (425. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> rIndex = {
      1.5780,
      1.5780, 1.5800, 1.5845,  // 609, 589.26, 550 nm
      1.5870, 1.5913, 1.5929,  // 530, 500, 490 nm
      1.5945, 1.5986, 1.6050,  // 481, 460, 435 nm
      1.6080,                  // 425 nm
      1.6080
    };
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {
      optPhotMinE_,
      hc_ / (750. * nm), hc_ / (740. * nm), hc_ / (380. * nm), hc_ / (370. * nm),
      optPhotMaxE_
    };
    std::vector<G4double> absLength = {
      attenuation_length,  attenuation_length,
      attenuation_length,  attenuation_length,
      attenuation_length,  attenuation_length
    };

    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,
      hc_ / (445. * nm), hc_ / (440. * nm), hc_ / (435. * nm),
      hc_ / (430. * nm), hc_ / (425. * nm), hc_ / (420. * nm),
      hc_ / (415. * nm), hc_ / (410. * nm), hc_ / (405. * nm),
      hc_ / (400. * nm), hc_ / (395. * nm), hc_ / (390. * nm),
      hc_ / (385. * nm), hc_ / (380. * nm), hc_ / (375. * nm),
      hc_ / (370. * nm), hc_ / (365. * nm), hc_ / (360. * nm),
      hc_ / (355. * nm), hc_ / (350. * nm), hc_ / (345. * nm),
      hc_ / (340. * nm), hc_ / (335. * nm), hc_ / (330. * nm),
      hc_ / (325. * nm), hc_ / (320. * nm), hc_ / (315. * nm),
      hc_ / (310. * nm), hc_ / (305. * nm), hc_ / (300. * nm),
      hc_ / (295. * nm), hc_ / (290. * nm), hc_ / (285. * nm),
      hc_ / (280. * nm), hc_ / (275. * nm),
      optPhotMaxE_
    };


    std::vector<G4double> WLS_absLength = {
      noAbsLength_,
      noAbsLength_,        (1. / 0.00007) * cm, (1. /  0.0003) * cm,
      (1. / 0.00104) * cm, (1. / 0.00223) * cm, (1. / 0.00408) * cm,
      (1. /  0.0104) * cm, (1. / 0.18544) * cm, (1. /  1.4094) * cm,
      (1. /  3.7088) * cm, (1. /  7.4176) * cm, (1. / 11.8682) * cm,
      (1. / 16.6155) * cm, (1. / 22.2529) * cm, (1. / 27.8162) * cm,
      (1. / 33.3794) * cm, (1. / 37.8671) * cm, (1. / 40.4262) * cm,
      (1. / 41.5388) * cm, (1. / 41.1679) * cm, (1. / 38.9426) * cm,
      (1. / 35.0113) * cm, (1. / 31.1541) * cm, (1. / 27.4453) * cm,
      (1. / 23.4398) * cm, (1. / 20.0276) * cm, (1. / 16.3188) * cm,
      (1. / 13.3518) * cm, (1. / 10.5331) * cm, (1. /  8.1594) * cm,
      (1. /  6.1196) * cm, (1. /  4.6731) * cm, (1. /  3.6346) * cm,
      (1. /  3.0412) * cm,  noAbsLength_,
      noAbsLength_
    };
    // XXX We are assuming that EJ286 doesn't absorb wave lengths shorter than 280 nm
    // although the spectrum continues ...
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);
    //for (int i=0; i<WLS_abs_entries; i++)
    //  G4cout << "* EJ286 WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (hc_ / WLS_abs_energy[i]) / nm
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

  G4MaterialPropertiesTable* G2P_FB118(G4double cromophore_concentration, G4double rindex, G4bool verbosity)
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

    G4cout << "rindex=" << rindex << G4endl;

    std::vector<G4double> rIndex = {
      rindex,
      rindex,  rindex,   // 609 , 589.26 nm
      rindex,  rindex,   // 550 , 530 nm
      rindex,  rindex,   // 500 , 490 nm
      rindex,  rindex,   // 481 , 460 nm
      rindex,  rindex,   // 435 , 425 nm
      rindex
    };
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // abs_energy got 324 entries
    std::vector<G4double> abs_energy = {optPhotMinE_, h_Planck*c_light/(601.073*nm), h_Planck*c_light/(599.073*nm), h_Planck*c_light/(598.034*nm), h_Planck*c_light/(596.994*nm), h_Planck*c_light/(595.954*nm), h_Planck*c_light/(595.063*nm), h_Planck*c_light/(594.022*nm), h_Planck*c_light/(592.982*nm), h_Planck*c_light/(591.941*nm), h_Planck*c_light/(591.048*nm), h_Planck*c_light/(590.007*nm), h_Planck*c_light/(588.965*nm), h_Planck*c_light/(588.072*nm), h_Planck*c_light/(587.029*nm), h_Planck*c_light/(585.987*nm), h_Planck*c_light/(584.944*nm), h_Planck*c_light/(584.05*nm), h_Planck*c_light/(583.006*nm), h_Planck*c_light/(581.963*nm), h_Planck*c_light/(581.068*nm), h_Planck*c_light/(580.024*nm), h_Planck*c_light/(578.979*nm), h_Planck*c_light/(577.935*nm), h_Planck*c_light/(577.039*nm), h_Planck*c_light/(575.994*nm), h_Planck*c_light/(574.948*nm), h_Planck*c_light/(574.052*nm), h_Planck*c_light/(573.006*nm), h_Planck*c_light/(571.96*nm), h_Planck*c_light/(571.063*nm), h_Planck*c_light/(570.016*nm), h_Planck*c_light/(568.969*nm), h_Planck*c_light/(568.071*nm), h_Planck*c_light/(567.024*nm), h_Planck*c_light/(565.976*nm), h_Planck*c_light/(564.928*nm), h_Planck*c_light/(564.029*nm), h_Planck*c_light/(562.98*nm), h_Planck*c_light/(561.932*nm), h_Planck*c_light/(561.032*nm), h_Planck*c_light/(559.983*nm), h_Planck*c_light/(558.934*nm), h_Planck*c_light/(558.034*nm), h_Planck*c_light/(556.984*nm), h_Planck*c_light/(555.933*nm), h_Planck*c_light/(555.033*nm), h_Planck*c_light/(553.982*nm), h_Planck*c_light/(552.931*nm), h_Planck*c_light/(552.03*nm), h_Planck*c_light/(550.978*nm), h_Planck*c_light/(549.926*nm), h_Planck*c_light/(549.024*nm), h_Planck*c_light/(547.972*nm), h_Planck*c_light/(547.07*nm), h_Planck*c_light/(546.017*nm), h_Planck*c_light/(544.964*nm), h_Planck*c_light/(544.061*nm), h_Planck*c_light/(543.008*nm), h_Planck*c_light/(541.954*nm), h_Planck*c_light/(541.05*nm), h_Planck*c_light/(539.996*nm), h_Planck*c_light/(538.942*nm), h_Planck*c_light/(538.038*nm), h_Planck*c_light/(536.982*nm), h_Planck*c_light/(535.927*nm), h_Planck*c_light/(535.022*nm), h_Planck*c_light/(533.967*nm), h_Planck*c_light/(533.062*nm), h_Planck*c_light/(532.006*nm), h_Planck*c_light/(530.949*nm), h_Planck*c_light/(530.043*nm), h_Planck*c_light/(528.986*nm), h_Planck*c_light/(527.929*nm), h_Planck*c_light/(527.023*nm), h_Planck*c_light/(525.965*nm), h_Planck*c_light/(525.059*nm), h_Planck*c_light/(524.0*nm), h_Planck*c_light/(522.942*nm), h_Planck*c_light/(522.035*nm), h_Planck*c_light/(520.976*nm), h_Planck*c_light/(520.068*nm), h_Planck*c_light/(519.009*nm), h_Planck*c_light/(517.949*nm), h_Planck*c_light/(517.041*nm), h_Planck*c_light/(515.981*nm), h_Planck*c_light/(515.073*nm), h_Planck*c_light/(514.012*nm), h_Planck*c_light/(512.952*nm), h_Planck*c_light/(512.042*nm), h_Planck*c_light/(510.981*nm), h_Planck*c_light/(510.071*nm), h_Planck*c_light/(509.01*nm), h_Planck*c_light/(507.948*nm), h_Planck*c_light/(507.038*nm), h_Planck*c_light/(505.976*nm), h_Planck*c_light/(505.065*nm), h_Planck*c_light/(504.002*nm), h_Planck*c_light/(502.939*nm), h_Planck*c_light/(502.028*nm), h_Planck*c_light/(500.965*nm), h_Planck*c_light/(500.053*nm), h_Planck*c_light/(498.989*nm), h_Planck*c_light/(497.925*nm), h_Planck*c_light/(497.013*nm), h_Planck*c_light/(495.949*nm), h_Planck*c_light/(495.036*nm), h_Planck*c_light/(493.971*nm), h_Planck*c_light/(493.058*nm), h_Planck*c_light/(491.993*nm), h_Planck*c_light/(490.927*nm), h_Planck*c_light/(490.014*nm), h_Planck*c_light/(488.948*nm), h_Planck*c_light/(488.034*nm), h_Planck*c_light/(486.968*nm), h_Planck*c_light/(486.053*nm), h_Planck*c_light/(484.986*nm), h_Planck*c_light/(484.072*nm), h_Planck*c_light/(483.004*nm), h_Planck*c_light/(481.937*nm), h_Planck*c_light/(481.022*nm), h_Planck*c_light/(479.954*nm), h_Planck*c_light/(479.038*nm), h_Planck*c_light/(477.97*nm), h_Planck*c_light/(477.054*nm), h_Planck*c_light/(475.985*nm), h_Planck*c_light/(475.069*nm), h_Planck*c_light/(474.0*nm), h_Planck*c_light/(472.93*nm), h_Planck*c_light/(472.014*nm), h_Planck*c_light/(470.944*nm), h_Planck*c_light/(470.027*nm), h_Planck*c_light/(468.956*nm), h_Planck*c_light/(468.039*nm), h_Planck*c_light/(466.968*nm), h_Planck*c_light/(466.05*nm), h_Planck*c_light/(464.979*nm), h_Planck*c_light/(464.061*nm), h_Planck*c_light/(462.99*nm), h_Planck*c_light/(462.071*nm), h_Planck*c_light/(460.999*nm), h_Planck*c_light/(459.927*nm), h_Planck*c_light/(459.008*nm), h_Planck*c_light/(457.936*nm), h_Planck*c_light/(457.016*nm), h_Planck*c_light/(455.943*nm), h_Planck*c_light/(455.023*nm), h_Planck*c_light/(453.95*nm), h_Planck*c_light/(453.03*nm), h_Planck*c_light/(451.956*nm), h_Planck*c_light/(451.036*nm), h_Planck*c_light/(449.962*nm), h_Planck*c_light/(449.041*nm), h_Planck*c_light/(447.966*nm), h_Planck*c_light/(447.045*nm), h_Planck*c_light/(445.97*nm), h_Planck*c_light/(445.049*nm), h_Planck*c_light/(443.974*nm), h_Planck*c_light/(443.052*nm), h_Planck*c_light/(441.976*nm), h_Planck*c_light/(441.054*nm), h_Planck*c_light/(439.978*nm), h_Planck*c_light/(439.056*nm), h_Planck*c_light/(437.979*nm), h_Planck*c_light/(437.056*nm), h_Planck*c_light/(435.98*nm), h_Planck*c_light/(435.056*nm), h_Planck*c_light/(433.979*nm), h_Planck*c_light/(433.056*nm), h_Planck*c_light/(431.978*nm), h_Planck*c_light/(431.054*nm), h_Planck*c_light/(429.976*nm), h_Planck*c_light/(429.044*nm), h_Planck*c_light/(427.967*nm), h_Planck*c_light/(427.043*nm), h_Planck*c_light/(425.965*nm), h_Planck*c_light/(425.041*nm), h_Planck*c_light/(423.963*nm), h_Planck*c_light/(423.039*nm), h_Planck*c_light/(421.96*nm), h_Planck*c_light/(421.036*nm), h_Planck*c_light/(419.957*nm), h_Planck*c_light/(419.032*nm), h_Planck*c_light/(417.952*nm), h_Planck*c_light/(417.027*nm), h_Planck*c_light/(415.947*nm), h_Planck*c_light/(415.022*nm), h_Planck*c_light/(413.942*nm), h_Planck*c_light/(413.016*nm), h_Planck*c_light/(411.936*nm), h_Planck*c_light/(411.009*nm), h_Planck*c_light/(409.928*nm), h_Planck*c_light/(409.002*nm), h_Planck*c_light/(408.075*nm), h_Planck*c_light/(406.994*nm), h_Planck*c_light/(406.067*nm), h_Planck*c_light/(404.985*nm), h_Planck*c_light/(404.058*nm), h_Planck*c_light/(402.976*nm), h_Planck*c_light/(402.049*nm), h_Planck*c_light/(400.966*nm), h_Planck*c_light/(400.038*nm), h_Planck*c_light/(398.956*nm), h_Planck*c_light/(398.027*nm), h_Planck*c_light/(396.944*nm), h_Planck*c_light/(396.016*nm), h_Planck*c_light/(394.932*nm), h_Planck*c_light/(394.004*nm), h_Planck*c_light/(393.075*nm), h_Planck*c_light/(391.991*nm), h_Planck*c_light/(391.062*nm), h_Planck*c_light/(389.977*nm), h_Planck*c_light/(389.048*nm), h_Planck*c_light/(387.963*nm), h_Planck*c_light/(387.034*nm), h_Planck*c_light/(385.948*nm), h_Planck*c_light/(385.018*nm), h_Planck*c_light/(383.933*nm), h_Planck*c_light/(383.003*nm), h_Planck*c_light/(382.072*nm), h_Planck*c_light/(380.986*nm), h_Planck*c_light/(380.056*nm), h_Planck*c_light/(378.97*nm), h_Planck*c_light/(378.038*nm), h_Planck*c_light/(376.952*nm), h_Planck*c_light/(376.021*nm), h_Planck*c_light/(374.934*nm), h_Planck*c_light/(374.002*nm), h_Planck*c_light/(373.071*nm), h_Planck*c_light/(371.983*nm), h_Planck*c_light/(371.051*nm), h_Planck*c_light/(369.964*nm), h_Planck*c_light/(369.031*nm), h_Planck*c_light/(367.944*nm), h_Planck*c_light/(367.011*nm), h_Planck*c_light/(365.923*nm), h_Planck*c_light/(364.99*nm), h_Planck*c_light/(364.057*nm), h_Planck*c_light/(362.968*nm), h_Planck*c_light/(362.035*nm), h_Planck*c_light/(360.946*nm), h_Planck*c_light/(360.012*nm), h_Planck*c_light/(358.923*nm), h_Planck*c_light/(357.989*nm), h_Planck*c_light/(357.055*nm), h_Planck*c_light/(355.965*nm), h_Planck*c_light/(355.031*nm), h_Planck*c_light/(353.941*nm), h_Planck*c_light/(353.006*nm), h_Planck*c_light/(352.072*nm), h_Planck*c_light/(350.981*nm), h_Planck*c_light/(350.046*nm), h_Planck*c_light/(348.956*nm), h_Planck*c_light/(348.02*nm), h_Planck*c_light/(346.929*nm), h_Planck*c_light/(345.994*nm), h_Planck*c_light/(345.058*nm), h_Planck*c_light/(343.967*nm), h_Planck*c_light/(343.031*nm), h_Planck*c_light/(341.939*nm), h_Planck*c_light/(341.003*nm), h_Planck*c_light/(340.067*nm), h_Planck*c_light/(338.974*nm), h_Planck*c_light/(338.038*nm), h_Planck*c_light/(336.945*nm), h_Planck*c_light/(336.009*nm), h_Planck*c_light/(335.072*nm), h_Planck*c_light/(333.979*nm), h_Planck*c_light/(333.042*nm), h_Planck*c_light/(331.948*nm), h_Planck*c_light/(331.011*nm), h_Planck*c_light/(330.074*nm), h_Planck*c_light/(328.98*nm), h_Planck*c_light/(328.042*nm), h_Planck*c_light/(326.948*nm), h_Planck*c_light/(326.01*nm), h_Planck*c_light/(325.072*nm), h_Planck*c_light/(323.978*nm), h_Planck*c_light/(323.04*nm), h_Planck*c_light/(321.945*nm), h_Planck*c_light/(321.006*nm), h_Planck*c_light/(320.068*nm), h_Planck*c_light/(318.972*nm), h_Planck*c_light/(318.034*nm), h_Planck*c_light/(316.938*nm), h_Planck*c_light/(315.999*nm), h_Planck*c_light/(315.06*nm), h_Planck*c_light/(313.964*nm), h_Planck*c_light/(313.025*nm), h_Planck*c_light/(311.928*nm), h_Planck*c_light/(310.989*nm), h_Planck*c_light/(310.049*nm), h_Planck*c_light/(308.952*nm), h_Planck*c_light/(308.012*nm), h_Planck*c_light/(307.072*nm), h_Planck*c_light/(305.975*nm), h_Planck*c_light/(305.035*nm), h_Planck*c_light/(303.938*nm), h_Planck*c_light/(302.997*nm), h_Planck*c_light/(302.056*nm), h_Planck*c_light/(300.959*nm), h_Planck*c_light/(300.018*nm), h_Planck*c_light/(299.077*nm), h_Planck*c_light/(297.979*nm), h_Planck*c_light/(297.038*nm), h_Planck*c_light/(295.94*nm), h_Planck*c_light/(294.998*nm), h_Planck*c_light/(294.056*nm), h_Planck*c_light/(292.958*nm), h_Planck*c_light/(292.016*nm), h_Planck*c_light/(291.074*nm), h_Planck*c_light/(289.975*nm), h_Planck*c_light/(289.033*nm), h_Planck*c_light/(287.934*nm), h_Planck*c_light/(286.992*nm), h_Planck*c_light/(286.049*nm), h_Planck*c_light/(284.95*nm), h_Planck*c_light/(284.007*nm), h_Planck*c_light/(283.064*nm), h_Planck*c_light/(281.964*nm), h_Planck*c_light/(281.021*nm), h_Planck*c_light/(280.078*nm), h_Planck*c_light/(278.078*nm), optPhotMaxE_};
    
    // WLS_abs_energy got 102 entries
    std::vector<G4double> WLS_abs_energy = {optPhotMinE_, h_Planck*c_light/(431.052*nm), h_Planck*c_light/(429.052*nm), h_Planck*c_light/(427.974*nm), h_Planck*c_light/(427.05*nm), h_Planck*c_light/(425.971*nm), h_Planck*c_light/(425.046*nm), h_Planck*c_light/(423.967*nm), h_Planck*c_light/(423.042*nm), h_Planck*c_light/(421.962*nm), h_Planck*c_light/(421.037*nm), h_Planck*c_light/(419.957*nm), h_Planck*c_light/(419.032*nm), h_Planck*c_light/(417.951*nm), h_Planck*c_light/(417.025*nm), h_Planck*c_light/(415.945*nm), h_Planck*c_light/(415.018*nm), h_Planck*c_light/(413.937*nm), h_Planck*c_light/(413.011*nm), h_Planck*c_light/(411.93*nm), h_Planck*c_light/(411.002*nm), h_Planck*c_light/(410.075*nm), h_Planck*c_light/(408.994*nm), h_Planck*c_light/(408.066*nm), h_Planck*c_light/(406.984*nm), h_Planck*c_light/(406.056*nm), h_Planck*c_light/(404.974*nm), h_Planck*c_light/(404.046*nm), h_Planck*c_light/(402.963*nm), h_Planck*c_light/(402.034*nm), h_Planck*c_light/(400.951*nm), h_Planck*c_light/(400.023*nm), h_Planck*c_light/(398.939*nm), h_Planck*c_light/(398.01*nm), h_Planck*c_light/(396.926*nm), h_Planck*c_light/(395.997*nm), h_Planck*c_light/(395.068*nm), h_Planck*c_light/(393.983*nm), h_Planck*c_light/(393.054*nm), h_Planck*c_light/(391.969*nm), h_Planck*c_light/(391.039*nm), h_Planck*c_light/(389.954*nm), h_Planck*c_light/(389.023*nm), h_Planck*c_light/(387.938*nm), h_Planck*c_light/(387.007*nm), h_Planck*c_light/(386.077*nm), h_Planck*c_light/(384.991*nm), h_Planck*c_light/(384.06*nm), h_Planck*c_light/(382.973*nm), h_Planck*c_light/(382.042*nm), h_Planck*c_light/(380.956*nm), h_Planck*c_light/(380.024*nm), h_Planck*c_light/(378.937*nm), h_Planck*c_light/(378.005*nm), h_Planck*c_light/(377.073*nm), h_Planck*c_light/(375.986*nm), h_Planck*c_light/(375.053*nm), h_Planck*c_light/(373.966*nm), h_Planck*c_light/(373.033*nm), h_Planck*c_light/(371.945*nm), h_Planck*c_light/(371.012*nm), h_Planck*c_light/(369.924*nm), h_Planck*c_light/(368.991*nm), h_Planck*c_light/(368.057*nm), h_Planck*c_light/(366.968*nm), h_Planck*c_light/(366.035*nm), h_Planck*c_light/(364.946*nm), h_Planck*c_light/(364.012*nm), h_Planck*c_light/(362.922*nm), h_Planck*c_light/(361.988*nm), h_Planck*c_light/(361.054*nm), h_Planck*c_light/(359.964*nm), h_Planck*c_light/(359.03*nm), h_Planck*c_light/(357.939*nm), h_Planck*c_light/(357.005*nm), h_Planck*c_light/(356.07*nm), h_Planck*c_light/(354.979*nm), h_Planck*c_light/(354.044*nm), h_Planck*c_light/(352.953*nm), h_Planck*c_light/(352.018*nm), h_Planck*c_light/(350.926*nm), h_Planck*c_light/(349.99*nm), h_Planck*c_light/(349.055*nm), h_Planck*c_light/(347.963*nm), h_Planck*c_light/(347.027*nm), h_Planck*c_light/(345.935*nm), h_Planck*c_light/(344.998*nm), h_Planck*c_light/(344.062*nm), h_Planck*c_light/(342.969*nm), h_Planck*c_light/(342.033*nm), h_Planck*c_light/(340.94*nm), h_Planck*c_light/(340.003*nm), h_Planck*c_light/(339.066*nm), h_Planck*c_light/(337.972*nm), h_Planck*c_light/(337.035*nm), h_Planck*c_light/(335.942*nm), h_Planck*c_light/(335.004*nm), h_Planck*c_light/(334.066*nm), h_Planck*c_light/(332.972*nm), h_Planck*c_light/(332.034*nm), h_Planck*c_light/(330.034*nm), optPhotMaxE_};


    G4double concentration;
    std::vector<G4double> available_concentrations;
    std::map<G4double, std::vector<G4double>> absLength;
    std::map<G4double, std::vector<G4double>> WLS_absLength;

    concentration = 8.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 8.1916*m, 8.404*m, 8.267*m, 8.3063*m, 8.395*m, 8.2733*m, 8.1249*m, 8.2838*m, 8.2509*m, 8.0636*m, 8.0586*m, 7.9863*m, 8.1054*m, 8.2056*m, 8.1715*m, 8.1047*m, 8.1463*m, 8.3699*m, 8.1579*m, 8.0404*m, 8.087*m, 8.3966*m, 7.9784*m, 8.047*m, 8.2303*m, 8.0724*m, 8.1228*m, 8.0335*m, 8.3199*m, 7.8313*m, 7.9673*m, 7.8829*m, 7.9094*m, 7.7826*m, 7.7865*m, 7.9038*m, 7.8015*m, 7.9485*m, 7.851*m, 7.801*m, 7.899*m, 7.9238*m, 7.9347*m, 7.8373*m, 7.7262*m, 7.7338*m, 7.7912*m, 7.7733*m, 7.65*m, 7.7625*m, 7.6698*m, 7.5919*m, 7.7012*m, 7.51*m, 7.5907*m, 7.5763*m, 7.4511*m, 7.4605*m, 7.6127*m, 7.4185*m, 7.4424*m, 7.3628*m, 7.521*m, 7.472*m, 7.4042*m, 7.3167*m, 7.4539*m, 7.449*m, 7.334*m, 7.4076*m, 7.3089*m, 7.4103*m, 7.2427*m, 7.318*m, 7.3341*m, 7.2615*m, 7.2738*m, 7.1922*m, 7.3088*m, 7.1006*m, 7.2539*m, 7.1282*m, 7.1048*m, 7.178*m, 7.0218*m, 7.0991*m, 7.0897*m, 7.0512*m, 7.0086*m, 7.1743*m, 7.0664*m, 6.9257*m, 6.9815*m, 7.0505*m, 6.906*m, 6.84*m, 7.0087*m, 6.7677*m, 6.8387*m, 6.7982*m, 6.8388*m, 6.7578*m, 6.6626*m, 6.7713*m, 6.7381*m, 6.6308*m, 6.5512*m, 6.6281*m, 6.644*m, 6.6424*m, 6.5701*m, 6.5256*m, 6.5549*m, 6.5148*m, 6.4564*m, 6.5716*m, 6.5163*m, 6.4409*m, 6.4623*m, 6.3864*m, 6.4046*m, 6.3118*m, 6.3459*m, 6.3183*m, 6.331*m, 6.2791*m, 6.3088*m, 6.2379*m, 6.1737*m, 6.1324*m, 6.1732*m, 6.1744*m, 6.0149*m, 5.9462*m, 6.0302*m, 5.9644*m, 5.982*m, 5.918*m, 5.8833*m, 5.8408*m, 5.9222*m, 5.7917*m, 5.7442*m, 5.7619*m, 5.6372*m, 5.6313*m, 5.5925*m, 5.5817*m, 5.5416*m, 5.4838*m, 5.3322*m, 5.3111*m, 5.2893*m, 5.2315*m, 5.2032*m, 5.1054*m, 5.0002*m, 4.9245*m, 4.8419*m, 4.7264*m, 4.596*m, 4.4732*m, 4.3267*m, 4.1138*m, 3.8952*m, 3.7231*m, 3.4888*m, 3.2308*m, 3.011*m, 2.7741*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 2.5242*m, 2.291*m, 2.0659*m, 1.8376*m, 1.5903*m, 1.3739*m, 1.1562*m, 0.9358*m, 0.7675*m, 0.6202*m, 0.4925*m, 0.3851*m, 0.2863*m, 0.2178*m, 0.165*m, 0.1265*m, 0.0949*m, 0.0691*m, 0.0528*m, 0.0418*m, 0.0328*m, 0.0261*m, 0.0211*m, 0.0183*m, 0.0167*m, 0.0133*m, 0.0117*m, 0.0105*m, 0.0096891*m, 0.0089716*m, 0.0084991*m, 0.0082002*m, 0.0080207*m, 0.0079387*m, 0.0078843*m, 0.0078883*m, 0.0078813*m, 0.0078778*m, 0.0078473*m, 0.007795*m, 0.007701*m, 0.0075812*m, 0.0074206*m, 0.0072474*m, 0.0070338*m, 0.0067589*m, 0.0064389*m, 0.0061184*m, 0.0058425*m, 0.005552*m, 0.0053203*m, 0.0050843*m, 0.0049818*m, 0.004893*m, 0.0048477*m, 0.0048643*m, 0.0048915*m, 0.0049706*m, 0.0050252*m, 0.0051141*m, 0.0051646*m, 0.0052636*m, 0.0053242*m, 0.005371*m, 0.0054318*m, 0.0054395*m, 0.0054937*m, 0.0054865*m, 0.0055009*m, 0.0054696*m, 0.005451*m, 0.0054451*m, 0.0054307*m, 0.0054562*m, 0.0055016*m, 0.0056283*m, 0.0057409*m, 0.0059369*m, 0.0061275*m, 0.0063551*m, 0.0064948*m, 0.0067181*m, 0.0069183*m, 0.0071367*m, 0.0073426*m, 0.0075311*m, 0.0076923*m, 0.0079742*m, 0.008156*m, 0.0084221*m, 0.0085684*m, 0.008885*m, 0.0091589*m, 0.0094779*m, 0.01*m, 0.0103*m, 0.0108*m, 0.0113*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 16.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 6.9477*m, 7.2589*m, 7.0569*m, 7.1144*m, 7.2455*m, 7.0662*m, 6.8523*m, 7.0815*m, 7.0335*m, 6.7656*m, 6.7586*m, 6.6574*m, 6.8246*m, 6.968*m, 6.9189*m, 6.8236*m, 6.8828*m, 7.2082*m, 6.8994*m, 6.733*m, 6.7986*m, 7.2479*m, 6.6465*m, 6.7423*m, 7.0037*m, 6.7779*m, 6.8493*m, 6.7233*m, 7.1344*m, 6.4447*m, 6.6311*m, 6.515*m, 6.5513*m, 6.3791*m, 6.3844*m, 6.5435*m, 6.4045*m, 6.6051*m, 6.4715*m, 6.4039*m, 6.537*m, 6.571*m, 6.586*m, 6.453*m, 6.3036*m, 6.3138*m, 6.3906*m, 6.3666*m, 6.2029*m, 6.3521*m, 6.229*m, 6.1267*m, 6.2704*m, 6.0208*m, 6.1252*m, 6.1065*m, 5.9454*m, 5.9575*m, 6.1539*m, 5.9041*m, 5.9343*m, 5.8339*m, 6.035*m, 5.9721*m, 5.886*m, 5.7761*m, 5.9491*m, 5.9427*m, 5.7977*m, 5.8903*m, 5.7665*m, 5.8937*m, 5.6844*m, 5.7777*m, 5.7979*m, 5.7076*m, 5.7229*m, 5.6225*m, 5.7663*m, 5.5113*m, 5.6983*m, 5.5447*m, 5.5163*m, 5.6051*m, 5.4169*m, 5.5095*m, 5.4982*m, 5.452*m, 5.4012*m, 5.6006*m, 5.4702*m, 5.3035*m, 5.3691*m, 5.4511*m, 5.2804*m, 5.2036*m, 5.4014*m, 5.1203*m, 5.202*m, 5.1553*m, 5.2022*m, 5.109*m, 5.001*m, 5.1244*m, 5.0865*m, 4.9653*m, 4.8764*m, 4.9622*m, 4.98*m, 4.9783*m, 4.8974*m, 4.8482*m, 4.8805*m, 4.8363*m, 4.7721*m, 4.8991*m, 4.8378*m, 4.7552*m, 4.7786*m, 4.6961*m, 4.7158*m, 4.6159*m, 4.6524*m, 4.6228*m, 4.6364*m, 4.5809*m, 4.6126*m, 4.5373*m, 4.4696*m, 4.4264*m, 4.4691*m, 4.4704*m, 4.305*m, 4.235*m, 4.3208*m, 4.2534*m, 4.2714*m, 4.2064*m, 4.1715*m, 4.1289*m, 4.2107*m, 4.08*m, 4.033*m, 4.0504*m, 3.9283*m, 3.9225*m, 3.885*m, 3.8745*m, 3.836*m, 3.7808*m, 3.6383*m, 3.6186*m, 3.5984*m, 3.5451*m, 3.5192*m, 3.4303*m, 3.336*m, 3.2689*m, 3.1965*m, 3.0966*m, 2.9856*m, 2.8828*m, 2.7623*m, 2.5911*m, 2.4199*m, 2.2885*m, 2.114*m, 1.9275*m, 1.773*m, 1.611*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 1.4449*m, 1.2941*m, 1.1523*m, 1.012*m, 0.864*m, 0.7377*m, 0.6137*m, 0.4909*m, 0.3991*m, 0.32*m, 0.2525*m, 0.1963*m, 0.1452*m, 0.1101*m, 0.0832*m, 0.0637*m, 0.0477*m, 0.0347*m, 0.0264*m, 0.021*m, 0.0164*m, 0.0131*m, 0.0106*m, 0.0091694*m, 0.008354*m, 0.0066686*m, 0.0058575*m, 0.0052449*m, 0.0048445*m, 0.0044858*m, 0.0042495*m, 0.0041001*m, 0.0040103*m, 0.0039693*m, 0.0039422*m, 0.0039442*m, 0.0039406*m, 0.0039389*m, 0.0039237*m, 0.0038975*m, 0.0038505*m, 0.0037906*m, 0.0037103*m, 0.0036237*m, 0.0035169*m, 0.0033795*m, 0.0032195*m, 0.0030592*m, 0.0029212*m, 0.002776*m, 0.0026601*m, 0.0025422*m, 0.0024909*m, 0.0024465*m, 0.0024239*m, 0.0024321*m, 0.0024458*m, 0.0024853*m, 0.0025126*m, 0.002557*m, 0.0025823*m, 0.0026318*m, 0.0026621*m, 0.0026855*m, 0.0027159*m, 0.0027198*m, 0.0027468*m, 0.0027433*m, 0.0027504*m, 0.0027348*m, 0.0027255*m, 0.0027226*m, 0.0027154*m, 0.0027281*m, 0.0027508*m, 0.0028141*m, 0.0028705*m, 0.0029685*m, 0.0030637*m, 0.0031776*m, 0.0032474*m, 0.0033591*m, 0.0034591*m, 0.0035684*m, 0.0036713*m, 0.0037656*m, 0.0038462*m, 0.0039871*m, 0.004078*m, 0.0042111*m, 0.0042842*m, 0.0044425*m, 0.0045795*m, 0.0047389*m, 0.0049784*m, 0.0051419*m, 0.0053926*m, 0.0056263*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 24.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 6.0318*m, 6.3885*m, 6.1558*m, 6.2217*m, 6.3729*m, 6.1664*m, 5.9244*m, 6.1839*m, 6.1292*m, 5.8276*m, 5.8197*m, 5.7077*m, 5.8934*m, 6.0547*m, 5.9993*m, 5.8923*m, 5.9586*m, 6.3296*m, 5.9773*m, 5.7913*m, 5.8644*m, 6.3757*m, 5.6957*m, 5.8016*m, 6.0952*m, 5.8413*m, 5.9211*m, 5.7806*m, 6.2446*m, 5.4753*m, 5.6787*m, 5.5516*m, 5.5912*m, 5.4045*m, 5.4102*m, 5.5827*m, 5.4319*m, 5.6502*m, 5.5044*m, 5.4312*m, 5.5756*m, 5.6127*m, 5.6292*m, 5.4842*m, 5.3234*m, 5.3343*m, 5.4169*m, 5.391*m, 5.2162*m, 5.3754*m, 5.2438*m, 5.1356*m, 5.288*m, 5.0245*m, 5.134*m, 5.1143*m, 4.946*m, 4.9585*m, 5.1643*m, 4.9032*m, 4.9345*m, 4.8307*m, 5.0394*m, 4.9737*m, 4.8844*m, 4.7715*m, 4.9498*m, 4.9432*m, 4.7936*m, 4.8889*m, 4.7616*m, 4.8924*m, 4.6779*m, 4.7731*m, 4.7938*m, 4.7015*m, 4.7171*m, 4.6152*m, 4.7615*m, 4.5033*m, 4.692*m, 4.5368*m, 4.5084*m, 4.5977*m, 4.4092*m, 4.5015*m, 4.4902*m, 4.4441*m, 4.3936*m, 4.5931*m, 4.4622*m, 4.2969*m, 4.3617*m, 4.4432*m, 4.2742*m, 4.199*m, 4.3937*m, 4.1179*m, 4.1975*m, 4.152*m, 4.1976*m, 4.107*m, 4.0027*m, 4.1219*m, 4.0852*m, 3.9684*m, 3.8836*m, 3.9655*m, 3.9826*m, 3.9809*m, 3.9036*m, 3.8568*m, 3.8875*m, 3.8455*m, 3.7848*m, 3.9052*m, 3.847*m, 3.7689*m, 3.7909*m, 3.7133*m, 3.7318*m, 3.6383*m, 3.6724*m, 3.6447*m, 3.6574*m, 3.6058*m, 3.6352*m, 3.5653*m, 3.5028*m, 3.4631*m, 3.5023*m, 3.5034*m, 3.3521*m, 3.2886*m, 3.3664*m, 3.3053*m, 3.3216*m, 3.2628*m, 3.2313*m, 3.193*m, 3.2666*m, 3.1492*m, 3.1073*m, 3.1228*m, 3.0145*m, 3.0094*m, 2.9763*m, 2.9671*m, 2.9332*m, 2.8849*m, 2.7611*m, 2.7441*m, 2.7267*m, 2.6809*m, 2.6587*m, 2.5829*m, 2.5029*m, 2.4464*m, 2.3858*m, 2.3026*m, 2.2109*m, 2.1267*m, 2.0287*m, 1.8911*m, 1.7552*m, 1.6519*m, 1.5164*m, 1.3734*m, 1.2564*m, 1.1351*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 1.0121*m, 0.9017*m, 0.7989*m, 0.6983*m, 0.5931*m, 0.5042*m, 0.4177*m, 0.3327*m, 0.2696*m, 0.2157*m, 0.1697*m, 0.1318*m, 0.0973*m, 0.0737*m, 0.0556*m, 0.0425*m, 0.0318*m, 0.0231*m, 0.0176*m, 0.014*m, 0.011*m, 0.0087088*m, 0.0070573*m, 0.0061148*m, 0.0055709*m, 0.0044457*m, 0.003905*m, 0.0034966*m, 0.0032297*m, 0.0029905*m, 0.002833*m, 0.0027334*m, 0.0026736*m, 0.0026462*m, 0.0026281*m, 0.0026294*m, 0.0026271*m, 0.0026259*m, 0.0026158*m, 0.0025983*m, 0.002567*m, 0.0025271*m, 0.0024735*m, 0.0024158*m, 0.0023446*m, 0.002253*m, 0.0021463*m, 0.0020395*m, 0.0019475*m, 0.0018507*m, 0.0017734*m, 0.0016948*m, 0.0016606*m, 0.001631*m, 0.0016159*m, 0.0016214*m, 0.0016305*m, 0.0016569*m, 0.0016751*m, 0.0017047*m, 0.0017215*m, 0.0017545*m, 0.0017747*m, 0.0017903*m, 0.0018106*m, 0.0018132*m, 0.0018312*m, 0.0018288*m, 0.0018336*m, 0.0018232*m, 0.001817*m, 0.001815*m, 0.0018102*m, 0.0018187*m, 0.0018339*m, 0.0018761*m, 0.0019136*m, 0.001979*m, 0.0020425*m, 0.0021184*m, 0.0021649*m, 0.0022394*m, 0.0023061*m, 0.0023789*m, 0.0024475*m, 0.0025104*m, 0.0025641*m, 0.0026581*m, 0.0027187*m, 0.0028074*m, 0.0028561*m, 0.0029617*m, 0.003053*m, 0.0031593*m, 0.003319*m, 0.0034279*m, 0.0035951*m, 0.0037509*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 32.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 5.3293*m, 5.7045*m, 5.4588*m, 5.528*m, 5.6879*m, 5.4699*m, 5.2178*m, 5.4883*m, 5.4309*m, 5.1179*m, 5.1099*m, 4.9951*m, 5.1858*m, 5.3531*m, 5.2954*m, 5.1846*m, 5.2532*m, 5.642*m, 5.2726*m, 5.0807*m, 5.1558*m, 5.6909*m, 4.9829*m, 5.0913*m, 5.3954*m, 5.1321*m, 5.2144*m, 5.0697*m, 5.5521*m, 4.7594*m, 4.9656*m, 4.8364*m, 4.8765*m, 4.6882*m, 4.6939*m, 4.868*m, 4.7157*m, 4.9365*m, 4.7887*m, 4.715*m, 4.8607*m, 4.8984*m, 4.9151*m, 4.7684*m, 4.6071*m, 4.6179*m, 4.7006*m, 4.6747*m, 4.5003*m, 4.659*m, 4.5278*m, 4.4205*m, 4.5717*m, 4.3111*m, 4.419*m, 4.3994*m, 4.2342*m, 4.2465*m, 4.4489*m, 4.1924*m, 4.223*m, 4.1219*m, 4.3257*m, 4.2613*m, 4.1742*m, 4.0645*m, 4.2379*m, 4.2315*m, 4.0859*m, 4.1785*m, 4.055*m, 4.1819*m, 3.9742*m, 4.0661*m, 4.0861*m, 3.997*m, 4.012*m, 3.914*m, 4.0548*m, 3.8071*m, 3.9878*m, 3.839*m, 3.8119*m, 3.8972*m, 3.7176*m, 3.8053*m, 3.7946*m, 3.7507*m, 3.7028*m, 3.8928*m, 3.7679*m, 3.6115*m, 3.6727*m, 3.7499*m, 3.5901*m, 3.5195*m, 3.703*m, 3.4437*m, 3.5181*m, 3.4755*m, 3.5182*m, 3.4336*m, 3.3366*m, 3.4475*m, 3.4132*m, 3.3049*m, 3.2267*m, 3.3022*m, 3.318*m, 3.3165*m, 3.2451*m, 3.202*m, 3.2303*m, 3.1916*m, 3.136*m, 3.2465*m, 3.193*m, 3.1214*m, 3.1416*m, 3.0706*m, 3.0875*m, 3.0024*m, 3.0334*m, 3.0082*m, 3.0198*m, 2.9729*m, 2.9996*m, 2.9362*m, 2.8798*m, 2.8441*m, 2.8794*m, 2.8804*m, 2.7446*m, 2.6879*m, 2.7574*m, 2.7028*m, 2.7173*m, 2.6649*m, 2.6369*m, 2.603*m, 2.6683*m, 2.5643*m, 2.5272*m, 2.5409*m, 2.4455*m, 2.4411*m, 2.4121*m, 2.404*m, 2.3744*m, 2.3323*m, 2.2247*m, 2.21*m, 2.195*m, 2.1555*m, 2.1363*m, 2.0712*m, 2.0028*m, 1.9546*m, 1.9031*m, 1.8327*m, 1.7554*m, 1.6848*m, 1.603*m, 1.4888*m, 1.3769*m, 1.2925*m, 1.1822*m, 1.0668*m, 0.9729*m, 0.8762*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.7788*m, 0.6919*m, 0.6114*m, 0.533*m, 0.4516*m, 0.383*m, 0.3166*m, 0.2516*m, 0.2036*m, 0.1626*m, 0.1278*m, 0.0991*m, 0.0732*m, 0.0553*m, 0.0418*m, 0.0319*m, 0.0239*m, 0.0174*m, 0.0132*m, 0.0105*m, 0.0082148*m, 0.006533*m, 0.0052939*m, 0.0045868*m, 0.0041787*m, 0.0033343*m, 0.0029288*m, 0.0026224*m, 0.0024223*m, 0.0022429*m, 0.0021248*m, 0.0020501*m, 0.0020052*m, 0.0019847*m, 0.0019711*m, 0.0019721*m, 0.0019703*m, 0.0019694*m, 0.0019618*m, 0.0019487*m, 0.0019252*m, 0.0018953*m, 0.0018551*m, 0.0018119*m, 0.0017584*m, 0.0016897*m, 0.0016097*m, 0.0015296*m, 0.0014606*m, 0.001388*m, 0.0013301*m, 0.0012711*m, 0.0012455*m, 0.0012232*m, 0.0012119*m, 0.0012161*m, 0.0012229*m, 0.0012427*m, 0.0012563*m, 0.0012785*m, 0.0012912*m, 0.0013159*m, 0.001331*m, 0.0013427*m, 0.0013579*m, 0.0013599*m, 0.0013734*m, 0.0013716*m, 0.0013752*m, 0.0013674*m, 0.0013627*m, 0.0013613*m, 0.0013577*m, 0.001364*m, 0.0013754*m, 0.0014071*m, 0.0014352*m, 0.0014842*m, 0.0015319*m, 0.0015888*m, 0.0016237*m, 0.0016795*m, 0.0017296*m, 0.0017842*m, 0.0018356*m, 0.0018828*m, 0.0019231*m, 0.0019936*m, 0.002039*m, 0.0021055*m, 0.0021421*m, 0.0022213*m, 0.0022897*m, 0.0023695*m, 0.0024892*m, 0.0025709*m, 0.0026963*m, 0.0028132*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 40.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 4.7733*m, 5.1527*m, 4.9036*m, 4.9735*m, 5.1358*m, 4.9148*m, 4.6618*m, 4.9334*m, 4.8754*m, 4.5624*m, 4.5544*m, 4.4407*m, 4.6299*m, 4.7973*m, 4.7394*m, 4.6287*m, 4.6972*m, 5.0891*m, 4.7165*m, 4.5254*m, 4.6001*m, 5.1389*m, 4.4286*m, 4.5359*m, 4.8397*m, 4.5764*m, 4.6584*m, 4.5145*m, 4.9979*m, 4.2091*m, 4.4115*m, 4.2845*m, 4.3239*m, 4.1396*m, 4.1451*m, 4.3154*m, 4.1664*m, 4.3829*m, 4.2377*m, 4.1657*m, 4.3084*m, 4.3454*m, 4.3618*m, 4.2179*m, 4.0606*m, 4.0712*m, 4.1517*m, 4.1264*m, 3.9572*m, 4.1111*m, 3.9838*m, 3.8803*m, 4.0263*m, 3.7751*m, 3.8787*m, 3.8599*m, 3.7015*m, 3.7132*m, 3.9076*m, 3.6616*m, 3.6908*m, 3.5945*m, 3.7891*m, 3.7274*m, 3.6442*m, 3.54*m, 3.705*m, 3.6989*m, 3.5603*m, 3.6484*m, 3.531*m, 3.6516*m, 3.4546*m, 3.5415*m, 3.5605*m, 3.4761*m, 3.4903*m, 3.3977*m, 3.5308*m, 3.2973*m, 3.4674*m, 3.3272*m, 3.3018*m, 3.3819*m, 3.2135*m, 3.2956*m, 3.2855*m, 3.2445*m, 3.1997*m, 3.3778*m, 3.2606*m, 3.1147*m, 3.1716*m, 3.2437*m, 3.0948*m, 3.0293*m, 3.1999*m, 2.9593*m, 3.028*m, 2.9886*m, 3.0281*m, 2.9499*m, 2.8606*m, 2.9627*m, 2.9311*m, 2.8315*m, 2.7598*m, 2.829*m, 2.8435*m, 2.8421*m, 2.7767*m, 2.7373*m, 2.7631*m, 2.7278*m, 2.6771*m, 2.778*m, 2.7291*m, 2.6638*m, 2.6822*m, 2.6176*m, 2.633*m, 2.5557*m, 2.5838*m, 2.561*m, 2.5715*m, 2.529*m, 2.5532*m, 2.4959*m, 2.445*m, 2.4128*m, 2.4446*m, 2.4455*m, 2.3235*m, 2.2728*m, 2.335*m, 2.2861*m, 2.2991*m, 2.2522*m, 2.2273*m, 2.197*m, 2.2553*m, 2.1626*m, 2.1297*m, 2.1418*m, 2.0573*m, 2.0533*m, 2.0277*m, 2.0206*m, 1.9945*m, 1.9573*m, 1.8629*m, 1.85*m, 1.8368*m, 1.8023*m, 1.7855*m, 1.7287*m, 1.6692*m, 1.6275*m, 1.5829*m, 1.5221*m, 1.4556*m, 1.3949*m, 1.325*m, 1.2277*m, 1.1328*m, 1.0615*m, 0.9687*m, 0.8721*m, 0.7938*m, 0.7135*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.6329*m, 0.5613*m, 0.4952*m, 0.431*m, 0.3645*m, 0.3088*m, 0.2549*m, 0.2023*m, 0.1636*m, 0.1305*m, 0.1025*m, 0.0795*m, 0.0586*m, 0.0443*m, 0.0334*m, 0.0256*m, 0.0191*m, 0.0139*m, 0.0106*m, 0.0083913*m, 0.006573*m, 0.0052271*m, 0.0042356*m, 0.0036698*m, 0.0033433*m, 0.0026674*m, 0.002343*m, 0.002098*m, 0.0019378*m, 0.0017943*m, 0.0016998*m, 0.00164*m, 0.0016041*m, 0.0015877*m, 0.0015769*m, 0.0015777*m, 0.0015763*m, 0.0015756*m, 0.0015695*m, 0.001559*m, 0.0015402*m, 0.0015162*m, 0.0014841*m, 0.0014495*m, 0.0014068*m, 0.0013518*m, 0.0012878*m, 0.0012237*m, 0.0011685*m, 0.0011104*m, 0.0010641*m, 0.0010169*m, 0.00099637*m, 0.0009786*m, 0.00096954*m, 0.00097286*m, 0.00097831*m, 0.00099412*m, 0.001005*m, 0.0010228*m, 0.0010329*m, 0.0010527*m, 0.0010648*m, 0.0010742*m, 0.0010864*m, 0.0010879*m, 0.0010987*m, 0.0010973*m, 0.0011002*m, 0.0010939*m, 0.0010902*m, 0.001089*m, 0.0010861*m, 0.0010912*m, 0.0011003*m, 0.0011257*m, 0.0011482*m, 0.0011874*m, 0.0012255*m, 0.001271*m, 0.001299*m, 0.0013436*m, 0.0013837*m, 0.0014273*m, 0.0014685*m, 0.0015062*m, 0.0015385*m, 0.0015948*m, 0.0016312*m, 0.0016844*m, 0.0017137*m, 0.001777*m, 0.0018318*m, 0.0018956*m, 0.0019914*m, 0.0020567*m, 0.0021571*m, 0.0022505*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 48.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 4.3224*m, 4.6983*m, 4.4509*m, 4.52*m, 4.6815*m, 4.462*m, 4.2129*m, 4.4803*m, 4.4231*m, 4.1156*m, 4.1078*m, 3.997*m, 4.1816*m, 4.346*m, 4.289*m, 4.1805*m, 4.2476*m, 4.6349*m, 4.2666*m, 4.0795*m, 4.1524*m, 4.6845*m, 3.9853*m, 4.0898*m, 4.3878*m, 4.1293*m, 4.2095*m, 4.0689*m, 4.5443*m, 3.7728*m, 3.9687*m, 3.8456*m, 3.8837*m, 3.7059*m, 3.7112*m, 3.8756*m, 3.7317*m, 3.9409*m, 3.8004*m, 3.731*m, 3.8687*m, 3.9046*m, 3.9205*m, 3.7813*m, 3.6301*m, 3.6402*m, 3.7175*m, 3.6932*m, 3.531*m, 3.6786*m, 3.5565*m, 3.4576*m, 3.5972*m, 3.3577*m, 3.4562*m, 3.4383*m, 3.2879*m, 3.299*m, 3.4837*m, 3.2501*m, 3.2777*m, 3.1868*m, 3.3709*m, 3.3124*m, 3.2337*m, 3.1354*m, 3.2912*m, 3.2854*m, 3.1545*m, 3.2376*m, 3.1269*m, 3.2407*m, 3.0551*m, 3.1368*m, 3.1547*m, 3.0753*m, 3.0886*m, 3.0018*m, 3.1268*m, 2.9079*m, 3.0672*m, 2.9358*m, 2.9121*m, 2.987*m, 2.8298*m, 2.9064*m, 2.8969*m, 2.8587*m, 2.817*m, 2.9832*m, 2.8737*m, 2.738*m, 2.7909*m, 2.8579*m, 2.7196*m, 2.659*m, 2.8171*m, 2.5943*m, 2.6578*m, 2.6214*m, 2.6579*m, 2.5856*m, 2.5035*m, 2.5975*m, 2.5683*m, 2.4767*m, 2.411*m, 2.4744*m, 2.4878*m, 2.4865*m, 2.4264*m, 2.3904*m, 2.414*m, 2.3817*m, 2.3353*m, 2.4277*m, 2.3828*m, 2.3232*m, 2.34*m, 2.2811*m, 2.2951*m, 2.2248*m, 2.2503*m, 2.2295*m, 2.2391*m, 2.2005*m, 2.2225*m, 2.1704*m, 2.1242*m, 2.0951*m, 2.1239*m, 2.1247*m, 2.0144*m, 1.9687*m, 2.0248*m, 1.9807*m, 1.9924*m, 1.9502*m, 1.9278*m, 1.9006*m, 1.953*m, 1.8697*m, 1.8402*m, 1.8511*m, 1.7754*m, 1.7719*m, 1.749*m, 1.7426*m, 1.7193*m, 1.6862*m, 1.6022*m, 1.5908*m, 1.5791*m, 1.5485*m, 1.5337*m, 1.4834*m, 1.4309*m, 1.3941*m, 1.3549*m, 1.3015*m, 1.2432*m, 1.1902*m, 1.1291*m, 1.0445*m, 0.9622*m, 0.9005*m, 0.8206*m, 0.7375*m, 0.6704*m, 0.6018*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.5331*m, 0.4722*m, 0.4161*m, 0.3618*m, 0.3057*m, 0.2587*m, 0.2133*m, 0.1692*m, 0.1367*m, 0.109*m, 0.0856*m, 0.0663*m, 0.0489*m, 0.037*m, 0.0279*m, 0.0213*m, 0.0159*m, 0.0116*m, 0.0088311*m, 0.0069937*m, 0.0054781*m, 0.0043563*m, 0.0035299*m, 0.0030583*m, 0.0027862*m, 0.0022229*m, 0.0019525*m, 0.0017483*m, 0.0016148*m, 0.0014953*m, 0.0014165*m, 0.0013667*m, 0.0013368*m, 0.0013231*m, 0.0013141*m, 0.0013147*m, 0.0013135*m, 0.001313*m, 0.0013079*m, 0.0012992*m, 0.0012835*m, 0.0012635*m, 0.0012368*m, 0.0012079*m, 0.0011723*m, 0.0011265*m, 0.0010732*m, 0.0010197*m, 0.00097374*m, 0.00092533*m, 0.00088671*m, 0.00084739*m, 0.00083031*m, 0.0008155*m, 0.00080795*m, 0.00081071*m, 0.00081526*m, 0.00082844*m, 0.00083753*m, 0.00085235*m, 0.00086077*m, 0.00087727*m, 0.00088736*m, 0.00089516*m, 0.0009053*m, 0.00090659*m, 0.00091562*m, 0.00091442*m, 0.00091681*m, 0.00091159*m, 0.00090849*m, 0.00090752*m, 0.00090512*m, 0.00090937*m, 0.00091693*m, 0.00093804*m, 0.00095682*m, 0.00098948*m, 0.0010212*m, 0.0010592*m, 0.0010825*m, 0.0011197*m, 0.001153*m, 0.0011895*m, 0.0012238*m, 0.0012552*m, 0.0012821*m, 0.001329*m, 0.0013593*m, 0.0014037*m, 0.0014281*m, 0.0014808*m, 0.0015265*m, 0.0015796*m, 0.0016595*m, 0.001714*m, 0.0017975*m, 0.0018754*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 56.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 3.9493*m, 4.3176*m, 4.0747*m, 4.1424*m, 4.301*m, 4.0856*m, 3.8428*m, 4.1035*m, 4.0475*m, 3.7486*m, 3.741*m, 3.634*m, 3.8125*m, 3.9723*m, 3.9168*m, 3.8114*m, 3.8765*m, 4.2552*m, 3.895*m, 3.7136*m, 3.7842*m, 4.304*m, 3.6227*m, 3.7236*m, 4.0131*m, 3.7618*m, 3.8396*m, 3.7034*m, 4.1661*m, 3.4185*m, 3.6067*m, 3.4883*m, 3.5249*m, 3.3544*m, 3.3595*m, 3.5171*m, 3.3791*m, 3.5799*m, 3.445*m, 3.3785*m, 3.5105*m, 3.5449*m, 3.5603*m, 3.4266*m, 3.2821*m, 3.2917*m, 3.3656*m, 3.3423*m, 3.1878*m, 3.3284*m, 3.2119*m, 3.1181*m, 3.2507*m, 3.0233*m, 3.1167*m, 3.0997*m, 2.9574*m, 2.9679*m, 3.1428*m, 2.9218*m, 2.9478*m, 2.8621*m, 3.0359*m, 2.9805*m, 2.9063*m, 2.8138*m, 2.9606*m, 2.9551*m, 2.8318*m, 2.91*m, 2.8058*m, 2.9129*m, 2.7385*m, 2.8152*m, 2.832*m, 2.7574*m, 2.7699*m, 2.6885*m, 2.8057*m, 2.6007*m, 2.7498*m, 2.6268*m, 2.6047*m, 2.6747*m, 2.528*m, 2.5993*m, 2.5905*m, 2.5549*m, 2.516*m, 2.6711*m, 2.5688*m, 2.4426*m, 2.4917*m, 2.5542*m, 2.4255*m, 2.3693*m, 2.5162*m, 2.3095*m, 2.3682*m, 2.3345*m, 2.3683*m, 2.3014*m, 2.2256*m, 2.3124*m, 2.2855*m, 2.201*m, 2.1405*m, 2.1988*m, 2.2111*m, 2.2099*m, 2.1546*m, 2.1215*m, 2.1432*m, 2.1135*m, 2.0709*m, 2.1558*m, 2.1146*m, 2.0598*m, 2.0752*m, 2.0212*m, 2.034*m, 1.9697*m, 1.993*m, 1.9741*m, 1.9828*m, 1.9475*m, 1.9676*m, 1.92*m, 1.8779*m, 1.8513*m, 1.8776*m, 1.8783*m, 1.7779*m, 1.7364*m, 1.7873*m, 1.7473*m, 1.7579*m, 1.7197*m, 1.6993*m, 1.6747*m, 1.7221*m, 1.6467*m, 1.62*m, 1.6298*m, 1.5615*m, 1.5583*m, 1.5376*m, 1.5319*m, 1.5109*m, 1.4811*m, 1.4056*m, 1.3953*m, 1.3849*m, 1.3574*m, 1.3441*m, 1.2991*m, 1.2522*m, 1.2193*m, 1.1843*m, 1.1367*m, 1.0849*m, 1.0378*m, 0.9837*m, 0.9089*m, 0.8363*m, 0.782*m, 0.7117*m, 0.6388*m, 0.5802*m, 0.5203*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.4604*m, 0.4075*m, 0.3588*m, 0.3117*m, 0.2631*m, 0.2225*m, 0.1834*m, 0.1454*m, 0.1174*m, 0.0936*m, 0.0735*m, 0.0569*m, 0.0419*m, 0.0317*m, 0.0239*m, 0.0183*m, 0.0137*m, 0.0099274*m, 0.0075705*m, 0.0059952*m, 0.0046959*m, 0.0037342*m, 0.0030258*m, 0.0026216*m, 0.0023883*m, 0.0019053*m, 0.0016736*m, 0.0014985*m, 0.0013842*m, 0.0012817*m, 0.0012142*m, 0.0011715*m, 0.0011458*m, 0.0011341*m, 0.0011263*m, 0.0011269*m, 0.0011259*m, 0.0011254*m, 0.001121*m, 0.0011136*m, 0.0011001*m, 0.001083*m, 0.0010601*m, 0.0010353*m, 0.0010048*m, 0.00096556*m, 0.00091985*m, 0.00087406*m, 0.00083464*m, 0.00079314*m, 0.00076004*m, 0.00072633*m, 0.00071169*m, 0.000699*m, 0.00069253*m, 0.0006949*m, 0.00069879*m, 0.00071009*m, 0.00071788*m, 0.00073059*m, 0.00073781*m, 0.00075194*m, 0.0007606*m, 0.00076728*m, 0.00077597*m, 0.00077707*m, 0.00078481*m, 0.00078379*m, 0.00078584*m, 0.00078137*m, 0.00077871*m, 0.00077787*m, 0.00077582*m, 0.00077946*m, 0.00078594*m, 0.00080404*m, 0.00082013*m, 0.00084813*m, 0.00087535*m, 0.00090788*m, 0.00092782*m, 0.00095973*m, 0.00098832*m, 0.0010195*m, 0.0010489*m, 0.0010759*m, 0.0010989*m, 0.0011392*m, 0.0011651*m, 0.0012032*m, 0.0012241*m, 0.0012693*m, 0.0013084*m, 0.001354*m, 0.0014224*m, 0.0014691*m, 0.0015408*m, 0.0016075*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 64.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 3.6355*m, 3.9939*m, 3.7571*m, 3.823*m, 3.9777*m, 3.7677*m, 3.5325*m, 3.7852*m, 3.7308*m, 3.4416*m, 3.4343*m, 3.3314*m, 3.5033*m, 3.6578*m, 3.6041*m, 3.5022*m, 3.5651*m, 3.9329*m, 3.583*m, 3.408*m, 3.476*m, 3.9806*m, 3.3206*m, 3.4176*m, 3.6973*m, 3.4544*m, 3.5294*m, 3.3981*m, 3.8461*m, 3.125*m, 3.3052*m, 3.1918*m, 3.2268*m, 3.0639*m, 3.0688*m, 3.2193*m, 3.0874*m, 3.2795*m, 3.1503*m, 3.0868*m, 3.213*m, 3.246*m, 3.2606*m, 3.1328*m, 2.9949*m, 3.0042*m, 3.0745*m, 3.0523*m, 2.9053*m, 3.039*m, 2.9283*m, 2.8392*m, 2.9651*m, 2.7495*m, 2.8379*m, 2.8218*m, 2.6873*m, 2.6972*m, 2.8626*m, 2.6537*m, 2.6782*m, 2.5975*m, 2.7614*m, 2.7091*m, 2.6391*m, 2.552*m, 2.6903*m, 2.6851*m, 2.569*m, 2.6426*m, 2.5445*m, 2.6453*m, 2.4813*m, 2.5533*m, 2.5691*m, 2.499*m, 2.5108*m, 2.4345*m, 2.5444*m, 2.3523*m, 2.4919*m, 2.3767*m, 2.356*m, 2.4215*m, 2.2843*m, 2.351*m, 2.3427*m, 2.3094*m, 2.2732*m, 2.4181*m, 2.3225*m, 2.2048*m, 2.2505*m, 2.3088*m, 2.1888*m, 2.1366*m, 2.2733*m, 2.081*m, 2.1355*m, 2.1042*m, 2.1356*m, 2.0735*m, 2.0033*m, 2.0837*m, 2.0587*m, 1.9805*m, 1.9245*m, 1.9785*m, 1.9899*m, 1.9887*m, 1.9376*m, 1.907*m, 1.9271*m, 1.8996*m, 1.8603*m, 1.9387*m, 1.9006*m, 1.8501*m, 1.8643*m, 1.8145*m, 1.8263*m, 1.7671*m, 1.7886*m, 1.7711*m, 1.7791*m, 1.7467*m, 1.7651*m, 1.7214*m, 1.6827*m, 1.6584*m, 1.6824*m, 1.6832*m, 1.5911*m, 1.5532*m, 1.5997*m, 1.5631*m, 1.5728*m, 1.5378*m, 1.5192*m, 1.4967*m, 1.5401*m, 1.4712*m, 1.4469*m, 1.4558*m, 1.3936*m, 1.3906*m, 1.3718*m, 1.3666*m, 1.3475*m, 1.3205*m, 1.2519*m, 1.2426*m, 1.2331*m, 1.2083*m, 1.1962*m, 1.1555*m, 1.1131*m, 1.0834*m, 1.0519*m, 1.009*m, 0.9624*m, 0.9201*m, 0.8715*m, 0.8044*m, 0.7395*m, 0.691*m, 0.6283*m, 0.5635*m, 0.5114*m, 0.4582*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.4052*m, 0.3584*m, 0.3154*m, 0.2738*m, 0.231*m, 0.1953*m, 0.1608*m, 0.1274*m, 0.1029*m, 0.082*m, 0.0643*m, 0.0498*m, 0.0367*m, 0.0277*m, 0.0209*m, 0.016*m, 0.012*m, 0.0086875*m, 0.0066248*m, 0.0052462*m, 0.0041091*m, 0.0032676*m, 0.0026476*m, 0.0022939*m, 0.0020898*m, 0.0016672*m, 0.0014644*m, 0.0013112*m, 0.0012111*m, 0.0011214*m, 0.0010624*m, 0.001025*m, 0.0010026*m, 0.00099233*m, 0.00098554*m, 0.00098604*m, 0.00098516*m, 0.00098472*m, 0.00098092*m, 0.00097437*m, 0.00096262*m, 0.00094764*m, 0.00092757*m, 0.00090593*m, 0.00087922*m, 0.00084487*m, 0.00080487*m, 0.0007648*m, 0.00073031*m, 0.000694*m, 0.00066504*m, 0.00063554*m, 0.00062273*m, 0.00061162*m, 0.00060596*m, 0.00060804*m, 0.00061144*m, 0.00062133*m, 0.00062814*m, 0.00063926*m, 0.00064558*m, 0.00065795*m, 0.00066552*m, 0.00067137*m, 0.00067897*m, 0.00067994*m, 0.00068671*m, 0.00068581*m, 0.00068761*m, 0.0006837*m, 0.00068137*m, 0.00068064*m, 0.00067884*m, 0.00068202*m, 0.0006877*m, 0.00070353*m, 0.00071762*m, 0.00074211*m, 0.00076593*m, 0.00079439*m, 0.00081185*m, 0.00083977*m, 0.00086478*m, 0.00089209*m, 0.00091782*m, 0.00094139*m, 0.00096154*m, 0.00099678*m, 0.0010195*m, 0.0010528*m, 0.001071*m, 0.0011106*m, 0.0011449*m, 0.0011847*m, 0.0012446*m, 0.0012855*m, 0.0013482*m, 0.0014066*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 72.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 3.3679*m, 3.7154*m, 3.4855*m, 3.5493*m, 3.6996*m, 3.4957*m, 3.2686*m, 3.5126*m, 3.46*m, 3.1811*m, 3.1741*m, 3.0753*m, 3.2404*m, 3.3894*m, 3.3376*m, 3.2394*m, 3.3*m, 3.6561*m, 3.3172*m, 3.1488*m, 3.2142*m, 3.7024*m, 3.0649*m, 3.158*m, 3.4276*m, 3.1934*m, 3.2656*m, 3.1393*m, 3.5717*m, 2.878*m, 3.0502*m, 2.9417*m, 2.9752*m, 2.8197*m, 2.8243*m, 2.968*m, 2.8421*m, 3.0256*m, 2.9021*m, 2.8415*m, 2.9619*m, 2.9935*m, 3.0076*m, 2.8853*m, 2.754*m, 2.7628*m, 2.8298*m, 2.8087*m, 2.6688*m, 2.796*m, 2.6907*m, 2.6061*m, 2.7257*m, 2.5212*m, 2.6049*m, 2.5897*m, 2.4624*m, 2.4717*m, 2.6283*m, 2.4307*m, 2.4538*m, 2.3776*m, 2.5325*m, 2.483*m, 2.4169*m, 2.3348*m, 2.4652*m, 2.4603*m, 2.3508*m, 2.4202*m, 2.3277*m, 2.4228*m, 2.2682*m, 2.336*m, 2.3509*m, 2.2849*m, 2.296*m, 2.2243*m, 2.3277*m, 2.1472*m, 2.2782*m, 2.1701*m, 2.1506*m, 2.2121*m, 2.0835*m, 2.1459*m, 2.1382*m, 2.107*m, 2.0731*m, 2.2089*m, 2.1192*m, 2.0091*m, 2.0519*m, 2.1064*m, 1.9942*m, 1.9455*m, 2.0732*m, 1.8936*m, 1.9445*m, 1.9153*m, 1.9446*m, 1.8867*m, 1.8213*m, 1.8962*m, 1.8729*m, 1.8001*m, 1.7481*m, 1.7983*m, 1.8088*m, 1.8078*m, 1.7603*m, 1.7319*m, 1.7505*m, 1.7251*m, 1.6886*m, 1.7613*m, 1.7259*m, 1.6791*m, 1.6923*m, 1.6462*m, 1.6571*m, 1.6023*m, 1.6221*m, 1.606*m, 1.6134*m, 1.5834*m, 1.6005*m, 1.56*m, 1.5243*m, 1.5018*m, 1.524*m, 1.5247*m, 1.4398*m, 1.4049*m, 1.4478*m, 1.414*m, 1.423*m, 1.3908*m, 1.3737*m, 1.353*m, 1.3929*m, 1.3295*m, 1.3071*m, 1.3154*m, 1.2582*m, 1.2556*m, 1.2383*m, 1.2336*m, 1.2161*m, 1.1913*m, 1.1286*m, 1.1201*m, 1.1114*m, 1.0887*m, 1.0777*m, 1.0405*m, 1.0018*m, 0.9748*m, 0.9461*m, 0.9071*m, 0.8647*m, 0.8263*m, 0.7823*m, 0.7215*m, 0.6628*m, 0.619*m, 0.5625*m, 0.5041*m, 0.4572*m, 0.4094*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.3618*m, 0.3198*m, 0.2813*m, 0.2441*m, 0.2059*m, 0.1739*m, 0.1432*m, 0.1134*m, 0.0915*m, 0.0729*m, 0.0572*m, 0.0443*m, 0.0326*m, 0.0247*m, 0.0186*m, 0.0142*m, 0.0106*m, 0.007723*m, 0.0058892*m, 0.0046636*m, 0.0036527*m, 0.0029046*m, 0.0023535*m, 0.0020391*m, 0.0018576*m, 0.0014819*m, 0.0013017*m, 0.0011655*m, 0.0010766*m, 0.00099684*m, 0.00094434*m, 0.00091114*m, 0.00089118*m, 0.00088207*m, 0.00087604*m, 0.00087648*m, 0.0008757*m, 0.00087531*m, 0.00087192*m, 0.00086611*m, 0.00085566*m, 0.00084235*m, 0.00082451*m, 0.00080527*m, 0.00078153*m, 0.00075099*m, 0.00071544*m, 0.00067982*m, 0.00064916*m, 0.00061689*m, 0.00059114*m, 0.00056492*m, 0.00055354*m, 0.00054367*m, 0.00053863*m, 0.00054048*m, 0.0005435*m, 0.00055229*m, 0.00055835*m, 0.00056823*m, 0.00057385*m, 0.00058484*m, 0.00059158*m, 0.00059677*m, 0.00060353*m, 0.00060439*m, 0.00061041*m, 0.00060961*m, 0.00061121*m, 0.00060773*m, 0.00060566*m, 0.00060501*m, 0.00060341*m, 0.00060624*m, 0.00061129*m, 0.00062536*m, 0.00063788*m, 0.00065966*m, 0.00068083*m, 0.00070613*m, 0.00072164*m, 0.00074646*m, 0.0007687*m, 0.00079297*m, 0.00081584*m, 0.00083679*m, 0.0008547*m, 0.00088603*m, 0.00090622*m, 0.00093579*m, 0.00095204*m, 0.00098722*m, 0.0010177*m, 0.0010531*m, 0.0011063*m, 0.0011426*m, 0.0011984*m, 0.0012503*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 80.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 3.137*m, 3.4732*m, 3.2505*m, 3.3122*m, 3.4578*m, 3.2604*m, 3.0414*m, 3.2767*m, 3.2258*m, 2.9573*m, 2.9506*m, 2.8558*m, 3.0143*m, 3.1577*m, 3.1078*m, 3.0133*m, 3.0716*m, 3.4156*m, 3.0881*m, 2.9263*m, 2.989*m, 3.4606*m, 2.8459*m, 2.9351*m, 3.1946*m, 2.9691*m, 3.0385*m, 2.9172*m, 3.3339*m, 2.6671*m, 2.8318*m, 2.7279*m, 2.7599*m, 2.6115*m, 2.6159*m, 2.7531*m, 2.6329*m, 2.8082*m, 2.6901*m, 2.6323*m, 2.7473*m, 2.7775*m, 2.7909*m, 2.6741*m, 2.549*m, 2.5573*m, 2.6212*m, 2.601*m, 2.468*m, 2.5889*m, 2.4887*m, 2.4084*m, 2.522*m, 2.3279*m, 2.4072*m, 2.3928*m, 2.2722*m, 2.2811*m, 2.4295*m, 2.2422*m, 2.2641*m, 2.1921*m, 2.3386*m, 2.2918*m, 2.2292*m, 2.1517*m, 2.2749*m, 2.2703*m, 2.1667*m, 2.2323*m, 2.145*m, 2.2347*m, 2.0889*m, 2.1528*m, 2.1669*m, 2.1046*m, 2.115*m, 2.0475*m, 2.1449*m, 1.9749*m, 2.0983*m, 1.9965*m, 1.9782*m, 2.036*m, 1.9152*m, 1.9738*m, 1.9665*m, 1.9372*m, 1.9053*m, 2.033*m, 1.9487*m, 1.8454*m, 1.8855*m, 1.9366*m, 1.8314*m, 1.7858*m, 1.9055*m, 1.7373*m, 1.7848*m, 1.7575*m, 1.7849*m, 1.7308*m, 1.6697*m, 1.7396*m, 1.7179*m, 1.6499*m, 1.6014*m, 1.6481*m, 1.658*m, 1.657*m, 1.6127*m, 1.5862*m, 1.6036*m, 1.5799*m, 1.5459*m, 1.6136*m, 1.5807*m, 1.5371*m, 1.5493*m, 1.5064*m, 1.5166*m, 1.4656*m, 1.484*m, 1.469*m, 1.4759*m, 1.448*m, 1.4639*m, 1.4263*m, 1.3932*m, 1.3723*m, 1.3929*m, 1.3935*m, 1.3148*m, 1.2825*m, 1.3222*m, 1.2909*m, 1.2992*m, 1.2694*m, 1.2535*m, 1.2344*m, 1.2713*m, 1.2127*m, 1.1921*m, 1.1997*m, 1.1469*m, 1.1444*m, 1.1285*m, 1.1241*m, 1.108*m, 1.0851*m, 1.0273*m, 1.0195*m, 1.0115*m, 0.9906*m, 0.9805*m, 0.9463*m, 0.9108*m, 0.886*m, 0.8596*m, 0.8239*m, 0.785*m, 0.7499*m, 0.7096*m, 0.6541*m, 0.6005*m, 0.5605*m, 0.5091*m, 0.456*m, 0.4134*m, 0.37*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.3268*m, 0.2888*m, 0.2539*m, 0.2203*m, 0.1857*m, 0.1568*m, 0.1291*m, 0.1022*m, 0.0825*m, 0.0657*m, 0.0515*m, 0.0399*m, 0.0294*m, 0.0222*m, 0.0168*m, 0.0128*m, 0.0095717*m, 0.0069512*m, 0.0053006*m, 0.0041974*m, 0.0032876*m, 0.0026142*m, 0.0021182*m, 0.0018352*m, 0.0016719*m, 0.0013337*m, 0.0011715*m, 0.001049*m, 0.00096891*m, 0.00089716*m, 0.00084991*m, 0.00082002*m, 0.00080207*m, 0.00079387*m, 0.00078843*m, 0.00078883*m, 0.00078813*m, 0.00078778*m, 0.00078473*m, 0.0007795*m, 0.0007701*m, 0.00075812*m, 0.00074206*m, 0.00072474*m, 0.00070338*m, 0.00067589*m, 0.00064389*m, 0.00061184*m, 0.00058425*m, 0.0005552*m, 0.00053203*m, 0.00050843*m, 0.00049818*m, 0.0004893*m, 0.00048477*m, 0.00048643*m, 0.00048915*m, 0.00049706*m, 0.00050252*m, 0.00051141*m, 0.00051646*m, 0.00052636*m, 0.00053242*m, 0.0005371*m, 0.00054318*m, 0.00054395*m, 0.00054937*m, 0.00054865*m, 0.00055009*m, 0.00054696*m, 0.0005451*m, 0.00054451*m, 0.00054307*m, 0.00054562*m, 0.00055016*m, 0.00056283*m, 0.00057409*m, 0.00059369*m, 0.00061275*m, 0.00063551*m, 0.00064948*m, 0.00067181*m, 0.00069183*m, 0.00071367*m, 0.00073426*m, 0.00075311*m, 0.00076923*m, 0.00079742*m, 0.0008156*m, 0.00084221*m, 0.00085684*m, 0.0008885*m, 0.00091589*m, 0.00094779*m, 0.00099569*m, 0.0010284*m, 0.0010785*m, 0.0011253*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 88.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.9357*m, 3.2606*m, 3.0452*m, 3.1048*m, 3.2457*m, 3.0547*m, 2.8437*m, 3.0705*m, 3.0214*m, 2.7629*m, 2.7564*m, 2.6656*m, 2.8176*m, 2.9557*m, 2.9076*m, 2.8167*m, 2.8727*m, 3.2048*m, 2.8887*m, 2.7331*m, 2.7934*m, 3.2484*m, 2.656*m, 2.7416*m, 2.9912*m, 2.7742*m, 2.8409*m, 2.7244*m, 3.1258*m, 2.485*m, 2.6425*m, 2.5431*m, 2.5738*m, 2.432*m, 2.4362*m, 2.5672*m, 2.4523*m, 2.6199*m, 2.507*m, 2.4518*m, 2.5617*m, 2.5905*m, 2.6034*m, 2.4917*m, 2.3723*m, 2.3803*m, 2.4412*m, 2.4219*m, 2.2952*m, 2.4104*m, 2.3149*m, 2.2386*m, 2.3466*m, 2.1622*m, 2.2375*m, 2.2237*m, 2.1093*m, 2.1177*m, 2.2586*m, 2.0809*m, 2.1017*m, 2.0334*m, 2.1723*m, 2.1278*m, 2.0686*m, 1.9952*m, 2.1118*m, 2.1075*m, 2.0094*m, 2.0715*m, 1.9889*m, 2.0738*m, 1.9358*m, 1.9963*m, 2.0096*m, 1.9507*m, 1.9605*m, 1.8967*m, 1.9888*m, 1.8283*m, 1.9447*m, 1.8486*m, 1.8314*m, 1.8859*m, 1.772*m, 1.8272*m, 1.8204*m, 1.7928*m, 1.7627*m, 1.8831*m, 1.8036*m, 1.7063*m, 1.744*m, 1.7922*m, 1.6932*m, 1.6503*m, 1.7629*m, 1.6047*m, 1.6494*m, 1.6237*m, 1.6495*m, 1.5986*m, 1.5413*m, 1.607*m, 1.5865*m, 1.5228*m, 1.4774*m, 1.5212*m, 1.5304*m, 1.5295*m, 1.488*m, 1.4632*m, 1.4794*m, 1.4572*m, 1.4255*m, 1.4888*m, 1.458*m, 1.4172*m, 1.4286*m, 1.3885*m, 1.398*m, 1.3504*m, 1.3676*m, 1.3536*m, 1.36*m, 1.334*m, 1.3488*m, 1.3137*m, 1.2828*m, 1.2634*m, 1.2826*m, 1.2831*m, 1.2098*m, 1.1797*m, 1.2167*m, 1.1875*m, 1.1953*m, 1.1675*m, 1.1528*m, 1.135*m, 1.1693*m, 1.1148*m, 1.0956*m, 1.1027*m, 1.0536*m, 1.0513*m, 1.0366*m, 1.0325*m, 1.0175*m, 0.9963*m, 0.9428*m, 0.9355*m, 0.9281*m, 0.9087*m, 0.8994*m, 0.8678*m, 0.8349*m, 0.812*m, 0.7876*m, 0.7546*m, 0.7188*m, 0.6864*m, 0.6493*m, 0.5982*m, 0.5489*m, 0.5122*m, 0.465*m, 0.4162*m, 0.3772*m, 0.3375*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.298*m, 0.2632*m, 0.2314*m, 0.2006*m, 0.1691*m, 0.1428*m, 0.1175*m, 0.093*m, 0.075*m, 0.0598*m, 0.0469*m, 0.0363*m, 0.0267*m, 0.0202*m, 0.0152*m, 0.0116*m, 0.0087023*m, 0.0063197*m, 0.0048189*m, 0.003816*m, 0.0029888*m, 0.0023766*m, 0.0019257*m, 0.0016684*m, 0.0015199*m, 0.0012125*m, 0.001065*m, 0.00095362*m, 0.00088083*m, 0.0008156*m, 0.00077264*m, 0.00074547*m, 0.00072915*m, 0.0007217*m, 0.00071676*m, 0.00071712*m, 0.00071648*m, 0.00071616*m, 0.00071339*m, 0.00070863*m, 0.00070009*m, 0.0006892*m, 0.0006746*m, 0.00065886*m, 0.00063943*m, 0.00061445*m, 0.00058536*m, 0.00055622*m, 0.00053113*m, 0.00050473*m, 0.00048366*m, 0.00046221*m, 0.00045289*m, 0.00044482*m, 0.0004407*m, 0.00044221*m, 0.00044468*m, 0.00045187*m, 0.00045683*m, 0.00046492*m, 0.00046951*m, 0.00047851*m, 0.00048402*m, 0.00048827*m, 0.0004938*m, 0.0004945*m, 0.00049943*m, 0.00049877*m, 0.00050008*m, 0.00049723*m, 0.00049554*m, 0.00049501*m, 0.0004937*m, 0.00049602*m, 0.00050015*m, 0.00051166*m, 0.0005219*m, 0.00053972*m, 0.00055704*m, 0.00057774*m, 0.00059043*m, 0.00061074*m, 0.00062893*m, 0.00064879*m, 0.00066751*m, 0.00068465*m, 0.0006993*m, 0.00072493*m, 0.00074146*m, 0.00076565*m, 0.00077894*m, 0.00080773*m, 0.00083263*m, 0.00086162*m, 0.00090517*m, 0.00093488*m, 0.00098048*m, 0.001023*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 96.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.7587*m, 3.0726*m, 2.8643*m, 2.9218*m, 3.0581*m, 2.8735*m, 2.6701*m, 2.8887*m, 2.8413*m, 2.5925*m, 2.5863*m, 2.4991*m, 2.6451*m, 2.778*m, 2.7316*m, 2.6442*m, 2.6981*m, 3.0186*m, 2.7134*m, 2.5639*m, 2.6218*m, 3.0608*m, 2.4899*m, 2.572*m, 2.8122*m, 2.6034*m, 2.6674*m, 2.5555*m, 2.9421*m, 2.3262*m, 2.477*m, 2.3818*m, 2.4111*m, 2.2755*m, 2.2795*m, 2.4048*m, 2.295*m, 2.4553*m, 2.3472*m, 2.2945*m, 2.3995*m, 2.4272*m, 2.4395*m, 2.3326*m, 2.2186*m, 2.2262*m, 2.2843*m, 2.2659*m, 2.1451*m, 2.255*m, 2.1639*m, 2.0911*m, 2.1941*m, 2.0184*m, 2.0901*m, 2.077*m, 1.9682*m, 1.9762*m, 2.1102*m, 1.9412*m, 1.9609*m, 1.8962*m, 2.028*m, 1.9858*m, 1.9295*m, 1.8599*m, 1.9706*m, 1.9665*m, 1.8734*m, 1.9323*m, 1.8539*m, 1.9345*m, 1.8037*m, 1.8609*m, 1.8735*m, 1.8178*m, 1.8271*m, 1.7667*m, 1.8539*m, 1.7019*m, 1.8121*m, 1.7211*m, 1.7048*m, 1.7564*m, 1.6487*m, 1.7009*m, 1.6944*m, 1.6683*m, 1.64*m, 1.7537*m, 1.6786*m, 1.5867*m, 1.6223*m, 1.6678*m, 1.5744*m, 1.5339*m, 1.6401*m, 1.491*m, 1.5331*m, 1.5089*m, 1.5331*m, 1.4852*m, 1.4313*m, 1.4931*m, 1.4739*m, 1.4138*m, 1.3712*m, 1.4123*m, 1.421*m, 1.4202*m, 1.3811*m, 1.3578*m, 1.3731*m, 1.3522*m, 1.3224*m, 1.3819*m, 1.353*m, 1.3146*m, 1.3254*m, 1.2877*m, 1.2967*m, 1.252*m, 1.2681*m, 1.255*m, 1.261*m, 1.2366*m, 1.2505*m, 1.2176*m, 1.1887*m, 1.1704*m, 1.1884*m, 1.189*m, 1.1203*m, 1.0921*m, 1.1267*m, 1.0995*m, 1.1067*m, 1.0807*m, 1.067*m, 1.0503*m, 1.0824*m, 1.0315*m, 1.0136*m, 1.0202*m, 0.9744*m, 0.9723*m, 0.9585*m, 0.9547*m, 0.9407*m, 0.9209*m, 0.8711*m, 0.8643*m, 0.8574*m, 0.8394*m, 0.8307*m, 0.8013*m, 0.7707*m, 0.7494*m, 0.7268*m, 0.6961*m, 0.6629*m, 0.6328*m, 0.5984*m, 0.5511*m, 0.5055*m, 0.4715*m, 0.4279*m, 0.3829*m, 0.3469*m, 0.3102*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2739*m, 0.2418*m, 0.2125*m, 0.1842*m, 0.1552*m, 0.131*m, 0.1078*m, 0.0853*m, 0.0688*m, 0.0548*m, 0.043*m, 0.0333*m, 0.0245*m, 0.0185*m, 0.014*m, 0.0107*m, 0.0079777*m, 0.0057934*m, 0.0044175*m, 0.0034981*m, 0.0027398*m, 0.0021786*m, 0.0017653*m, 0.0015294*m, 0.0013933*m, 0.0011114*m, 0.00097626*m, 0.00087415*m, 0.00080742*m, 0.00074763*m, 0.00070826*m, 0.00068335*m, 0.00066839*m, 0.00066155*m, 0.00065703*m, 0.00065736*m, 0.00065677*m, 0.00065648*m, 0.00065394*m, 0.00064958*m, 0.00064175*m, 0.00063176*m, 0.00061838*m, 0.00060395*m, 0.00058615*m, 0.00056325*m, 0.00053658*m, 0.00050987*m, 0.00048687*m, 0.00046266*m, 0.00044336*m, 0.00042369*m, 0.00041515*m, 0.00040775*m, 0.00040398*m, 0.00040536*m, 0.00040763*m, 0.00041422*m, 0.00041876*m, 0.00042617*m, 0.00043039*m, 0.00043863*m, 0.00044368*m, 0.00044758*m, 0.00045265*m, 0.00045329*m, 0.00045781*m, 0.00045721*m, 0.0004584*m, 0.0004558*m, 0.00045425*m, 0.00045376*m, 0.00045256*m, 0.00045468*m, 0.00045847*m, 0.00046902*m, 0.00047841*m, 0.00049474*m, 0.00051062*m, 0.00052959*m, 0.00054123*m, 0.00055984*m, 0.00057652*m, 0.00059473*m, 0.00061188*m, 0.0006276*m, 0.00064103*m, 0.00066452*m, 0.00067967*m, 0.00070184*m, 0.00071403*m, 0.00074042*m, 0.00076324*m, 0.00078982*m, 0.00082974*m, 0.00085698*m, 0.00089877*m, 0.00093772*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 104.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.6018*m, 2.905*m, 2.7036*m, 2.7592*m, 2.8911*m, 2.7125*m, 2.5165*m, 2.7272*m, 2.6815*m, 2.4419*m, 2.4359*m, 2.3522*m, 2.4924*m, 2.6204*m, 2.5757*m, 2.4916*m, 2.5434*m, 2.8527*m, 2.5582*m, 2.4144*m, 2.47*m, 2.8936*m, 2.3434*m, 2.4222*m, 2.6534*m, 2.4523*m, 2.514*m, 2.4063*m, 2.7788*m, 2.1865*m, 2.331*m, 2.2397*m, 2.2678*m, 2.138*m, 2.1418*m, 2.2618*m, 2.1566*m, 2.3102*m, 2.2066*m, 2.1561*m, 2.2567*m, 2.2832*m, 2.295*m, 2.1926*m, 2.0836*m, 2.0908*m, 2.1464*m, 2.1288*m, 2.0134*m, 2.1183*m, 2.0313*m, 1.9619*m, 2.0602*m, 1.8926*m, 1.9609*m, 1.9484*m, 1.8448*m, 1.8524*m, 1.9801*m, 1.8191*m, 1.8379*m, 1.7763*m, 1.9018*m, 1.8616*m, 1.808*m, 1.7418*m, 1.8471*m, 1.8431*m, 1.7546*m, 1.8106*m, 1.7361*m, 1.8127*m, 1.6884*m, 1.7428*m, 1.7548*m, 1.7018*m, 1.7106*m, 1.6533*m, 1.7361*m, 1.5919*m, 1.6964*m, 1.6101*m, 1.5946*m, 1.6436*m, 1.5415*m, 1.5909*m, 1.5848*m, 1.5601*m, 1.5332*m, 1.641*m, 1.5697*m, 1.4828*m, 1.5165*m, 1.5596*m, 1.4711*m, 1.4328*m, 1.5333*m, 1.3923*m, 1.4321*m, 1.4092*m, 1.4321*m, 1.3869*m, 1.3359*m, 1.3943*m, 1.3761*m, 1.3195*m, 1.2792*m, 1.318*m, 1.3263*m, 1.3255*m, 1.2886*m, 1.2666*m, 1.2811*m, 1.2614*m, 1.2333*m, 1.2894*m, 1.2621*m, 1.2259*m, 1.2361*m, 1.2006*m, 1.209*m, 1.1669*m, 1.1821*m, 1.1698*m, 1.1755*m, 1.1525*m, 1.1655*m, 1.1346*m, 1.1074*m, 1.0902*m, 1.1071*m, 1.1077*m, 1.0431*m, 1.0167*m, 1.0492*m, 1.0236*m, 1.0304*m, 1.006*m, 0.9931*m, 0.9775*m, 1.0076*m, 0.9598*m, 0.943*m, 0.9492*m, 0.9063*m, 0.9043*m, 0.8913*m, 0.8878*m, 0.8747*m, 0.8562*m, 0.8095*m, 0.8032*m, 0.7967*m, 0.7799*m, 0.7717*m, 0.7442*m, 0.7157*m, 0.6958*m, 0.6747*m, 0.6461*m, 0.615*m, 0.587*m, 0.555*m, 0.5109*m, 0.4684*m, 0.4369*m, 0.3963*m, 0.3545*m, 0.321*m, 0.2871*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2533*m, 0.2236*m, 0.1965*m, 0.1703*m, 0.1434*m, 0.1211*m, 0.0996*m, 0.0788*m, 0.0635*m, 0.0506*m, 0.0397*m, 0.0307*m, 0.0226*m, 0.0171*m, 0.0129*m, 0.0098466*m, 0.0073645*m, 0.005348*m, 0.0040779*m, 0.0032291*m, 0.0025291*m, 0.0020111*m, 0.0016295*m, 0.0014118*m, 0.0012861*m, 0.0010259*m, 0.00090116*m, 0.00080691*m, 0.00074531*m, 0.00069012*m, 0.00065377*m, 0.00063079*m, 0.00061697*m, 0.00061067*m, 0.00060649*m, 0.00060679*m, 0.00060625*m, 0.00060598*m, 0.00060364*m, 0.00059961*m, 0.00059238*m, 0.00058317*m, 0.00057082*m, 0.00055749*m, 0.00054106*m, 0.00051992*m, 0.0004953*m, 0.00047065*m, 0.00044942*m, 0.00042708*m, 0.00040925*m, 0.0003911*m, 0.00038322*m, 0.00037638*m, 0.0003729*m, 0.00037418*m, 0.00037627*m, 0.00038236*m, 0.00038655*m, 0.00039339*m, 0.00039728*m, 0.00040489*m, 0.00040955*m, 0.00041315*m, 0.00041783*m, 0.00041842*m, 0.00042259*m, 0.00042204*m, 0.00042314*m, 0.00042074*m, 0.0004193*m, 0.00041885*m, 0.00041775*m, 0.00041971*m, 0.0004232*m, 0.00043294*m, 0.00044161*m, 0.00045668*m, 0.00047134*m, 0.00048886*m, 0.0004996*m, 0.00051678*m, 0.00053217*m, 0.00054898*m, 0.00056481*m, 0.00057932*m, 0.00059172*m, 0.0006134*m, 0.00062739*m, 0.00064786*m, 0.00065911*m, 0.00068346*m, 0.00070453*m, 0.00072907*m, 0.00076591*m, 0.00079106*m, 0.00082964*m, 0.00086559*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 112.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.4618*m, 2.7548*m, 2.5601*m, 2.6138*m, 2.7413*m, 2.5687*m, 2.3797*m, 2.5829*m, 2.5387*m, 2.3078*m, 2.3021*m, 2.2216*m, 2.3564*m, 2.4797*m, 2.4367*m, 2.3556*m, 2.4056*m, 2.7042*m, 2.4198*m, 2.2814*m, 2.3349*m, 2.7437*m, 2.2131*m, 2.2889*m, 2.5116*m, 2.3178*m, 2.3772*m, 2.2736*m, 2.6327*m, 2.0626*m, 2.2012*m, 2.1136*m, 2.1406*m, 2.0161*m, 2.0198*m, 2.1348*m, 2.034*m, 2.1812*m, 2.0819*m, 2.0335*m, 2.1299*m, 2.1553*m, 2.1667*m, 2.0685*m, 1.9641*m, 1.971*m, 2.0242*m, 2.0074*m, 1.8969*m, 1.9973*m, 1.914*m, 1.8477*m, 1.9416*m, 1.7816*m, 1.8468*m, 1.8349*m, 1.736*m, 1.7432*m, 1.8651*m, 1.7115*m, 1.7294*m, 1.6707*m, 1.7903*m, 1.7519*m, 1.7009*m, 1.6378*m, 1.7381*m, 1.7344*m, 1.65*m, 1.7034*m, 1.6324*m, 1.7054*m, 1.587*m, 1.6388*m, 1.6502*m, 1.5997*m, 1.6081*m, 1.5536*m, 1.6324*m, 1.4952*m, 1.5946*m, 1.5125*m, 1.4978*m, 1.5443*m, 1.4473*m, 1.4943*m, 1.4885*m, 1.465*m, 1.4395*m, 1.5419*m, 1.4742*m, 1.3917*m, 1.4236*m, 1.4645*m, 1.3806*m, 1.3443*m, 1.4396*m, 1.3059*m, 1.3435*m, 1.3219*m, 1.3436*m, 1.3007*m, 1.2525*m, 1.3077*m, 1.2905*m, 1.2369*m, 1.1988*m, 1.2356*m, 1.2433*m, 1.2426*m, 1.2077*m, 1.1869*m, 1.2006*m, 1.1819*m, 1.1554*m, 1.2084*m, 1.1826*m, 1.1485*m, 1.158*m, 1.1245*m, 1.1324*m, 1.0927*m, 1.1071*m, 1.0954*m, 1.1008*m, 1.079*m, 1.0914*m, 1.0622*m, 1.0365*m, 1.0203*m, 1.0363*m, 1.0368*m, 0.9759*m, 0.951*m, 0.9816*m, 0.9575*m, 0.9639*m, 0.9409*m, 0.9287*m, 0.914*m, 0.9424*m, 0.8974*m, 0.8816*m, 0.8874*m, 0.847*m, 0.8451*m, 0.833*m, 0.8296*m, 0.8173*m, 0.7999*m, 0.7561*m, 0.7501*m, 0.7441*m, 0.7282*m, 0.7206*m, 0.6948*m, 0.668*m, 0.6493*m, 0.6295*m, 0.6027*m, 0.5736*m, 0.5474*m, 0.5174*m, 0.4761*m, 0.4364*m, 0.4069*m, 0.369*m, 0.33*m, 0.2988*m, 0.2671*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2357*m, 0.208*m, 0.1827*m, 0.1583*m, 0.1333*m, 0.1125*m, 0.0925*m, 0.0732*m, 0.059*m, 0.047*m, 0.0369*m, 0.0285*m, 0.021*m, 0.0159*m, 0.012*m, 0.009144*m, 0.0068388*m, 0.0049662*m, 0.0037867*m, 0.0029985*m, 0.0023485*m, 0.0018675*m, 0.0015131*m, 0.0013109*m, 0.0011943*m, 0.00095266*m, 0.00083679*m, 0.00074927*m, 0.00069208*m, 0.00064083*m, 0.00060708*m, 0.00058573*m, 0.0005729*m, 0.00056705*m, 0.00056317*m, 0.00056345*m, 0.00056295*m, 0.0005627*m, 0.00056052*m, 0.00055678*m, 0.00055007*m, 0.00054151*m, 0.00053004*m, 0.00051767*m, 0.00050241*m, 0.00048278*m, 0.00045992*m, 0.00043703*m, 0.00041732*m, 0.00039657*m, 0.00038002*m, 0.00036317*m, 0.00035585*m, 0.0003495*m, 0.00034626*m, 0.00034745*m, 0.0003494*m, 0.00035504*m, 0.00035894*m, 0.00036529*m, 0.0003689*m, 0.00037597*m, 0.0003803*m, 0.00038364*m, 0.00038798*m, 0.00038854*m, 0.00039241*m, 0.00039189*m, 0.00039292*m, 0.00039068*m, 0.00038935*m, 0.00038894*m, 0.00038791*m, 0.00038973*m, 0.00039297*m, 0.00040202*m, 0.00041007*m, 0.00042406*m, 0.00043768*m, 0.00045394*m, 0.00046391*m, 0.00047987*m, 0.00049416*m, 0.00050977*m, 0.00052447*m, 0.00053794*m, 0.00054945*m, 0.00056959*m, 0.00058257*m, 0.00060158*m, 0.00061203*m, 0.00063464*m, 0.00065421*m, 0.00067699*m, 0.00071121*m, 0.00073455*m, 0.00077038*m, 0.00080376*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 120.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.3362*m, 2.6194*m, 2.431*m, 2.4829*m, 2.6063*m, 2.4393*m, 2.2569*m, 2.453*m, 2.4103*m, 2.1877*m, 2.1822*m, 2.1047*m, 2.2345*m, 2.3534*m, 2.3119*m, 2.2337*m, 2.2819*m, 2.5704*m, 2.2956*m, 2.1622*m, 2.2137*m, 2.6087*m, 2.0966*m, 2.1695*m, 2.3842*m, 2.1974*m, 2.2545*m, 2.1548*m, 2.5012*m, 1.952*m, 2.0851*m, 2.001*m, 2.0268*m, 1.9074*m, 1.9109*m, 2.0213*m, 1.9245*m, 2.0659*m, 1.9705*m, 1.9241*m, 2.0166*m, 2.041*m, 2.0519*m, 1.9576*m, 1.8575*m, 1.8641*m, 1.9151*m, 1.899*m, 1.7932*m, 1.8894*m, 1.8096*m, 1.7461*m, 1.836*m, 1.6828*m, 1.7452*m, 1.7338*m, 1.6392*m, 1.6461*m, 1.7627*m, 1.6159*m, 1.6329*m, 1.5769*m, 1.6912*m, 1.6545*m, 1.6057*m, 1.5456*m, 1.6413*m, 1.6377*m, 1.5572*m, 1.6081*m, 1.5404*m, 1.61*m, 1.4971*m, 1.5464*m, 1.5573*m, 1.5092*m, 1.5172*m, 1.4652*m, 1.5403*m, 1.4096*m, 1.5043*m, 1.4261*m, 1.4121*m, 1.4564*m, 1.364*m, 1.4087*m, 1.4032*m, 1.3808*m, 1.3566*m, 1.4541*m, 1.3896*m, 1.3111*m, 1.3415*m, 1.3804*m, 1.3005*m, 1.266*m, 1.3567*m, 1.2295*m, 1.2653*m, 1.2448*m, 1.2654*m, 1.2247*m, 1.1789*m, 1.2313*m, 1.215*m, 1.1641*m, 1.1279*m, 1.1628*m, 1.1702*m, 1.1694*m, 1.1364*m, 1.1166*m, 1.1296*m, 1.1119*m, 1.0867*m, 1.1371*m, 1.1125*m, 1.0802*m, 1.0893*m, 1.0575*m, 1.065*m, 1.0273*m, 1.041*m, 1.0299*m, 1.035*m, 1.0144*m, 1.0261*m, 0.9985*m, 0.9741*m, 0.9588*m, 0.9739*m, 0.9744*m, 0.9168*m, 0.8932*m, 0.9222*m, 0.8994*m, 0.9055*m, 0.8837*m, 0.8722*m, 0.8583*m, 0.8851*m, 0.8426*m, 0.8277*m, 0.8332*m, 0.795*m, 0.7933*m, 0.7818*m, 0.7786*m, 0.767*m, 0.7506*m, 0.7092*m, 0.7036*m, 0.6979*m, 0.683*m, 0.6758*m, 0.6515*m, 0.6263*m, 0.6087*m, 0.59*m, 0.5648*m, 0.5375*m, 0.5128*m, 0.4846*m, 0.4458*m, 0.4085*m, 0.3808*m, 0.3453*m, 0.3087*m, 0.2794*m, 0.2498*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2203*m, 0.1944*m, 0.1707*m, 0.1479*m, 0.1245*m, 0.1051*m, 0.0864*m, 0.0684*m, 0.0551*m, 0.0439*m, 0.0344*m, 0.0266*m, 0.0196*m, 0.0148*m, 0.0112*m, 0.0085349*m, 0.0063832*m, 0.0046352*m, 0.0035343*m, 0.0027987*m, 0.0021919*m, 0.001743*m, 0.0014123*m, 0.0012236*m, 0.0011147*m, 0.00088915*m, 0.00078101*m, 0.00069932*m, 0.00064594*m, 0.00059811*m, 0.0005666*m, 0.00054668*m, 0.00053471*m, 0.00052924*m, 0.00052562*m, 0.00052589*m, 0.00052542*m, 0.00052519*m, 0.00052315*m, 0.00051966*m, 0.0005134*m, 0.00050541*m, 0.00049471*m, 0.00048316*m, 0.00046892*m, 0.0004506*m, 0.00042926*m, 0.00040789*m, 0.0003895*m, 0.00037013*m, 0.00035469*m, 0.00033895*m, 0.00033212*m, 0.0003262*m, 0.00032318*m, 0.00032429*m, 0.0003261*m, 0.00033137*m, 0.00033501*m, 0.00034094*m, 0.00034431*m, 0.00035091*m, 0.00035495*m, 0.00035806*m, 0.00036212*m, 0.00036263*m, 0.00036625*m, 0.00036577*m, 0.00036672*m, 0.00036464*m, 0.0003634*m, 0.00036301*m, 0.00036205*m, 0.00036375*m, 0.00036677*m, 0.00037522*m, 0.00038273*m, 0.00039579*m, 0.0004085*m, 0.00042368*m, 0.00043298*m, 0.00044788*m, 0.00046122*m, 0.00047578*m, 0.00048951*m, 0.00050208*m, 0.00051282*m, 0.00053162*m, 0.00054373*m, 0.00056147*m, 0.00057123*m, 0.00059233*m, 0.0006106*m, 0.00063186*m, 0.00066379*m, 0.00068558*m, 0.00071902*m, 0.00075018*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 128.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.2227*m, 2.4966*m, 2.3143*m, 2.3645*m, 2.4839*m, 2.3223*m, 2.1462*m, 2.3356*m, 2.2943*m, 2.0794*m, 2.0741*m, 1.9995*m, 2.1246*m, 2.2393*m, 2.1992*m, 2.1238*m, 2.1703*m, 2.4491*m, 2.1835*m, 2.0549*m, 2.1046*m, 2.4862*m, 1.9917*m, 2.0619*m, 2.2691*m, 2.0888*m, 2.1439*m, 2.0477*m, 2.3822*m, 1.8526*m, 1.9807*m, 1.8997*m, 1.9246*m, 1.8098*m, 1.8132*m, 1.9192*m, 1.8262*m, 1.9622*m, 1.8704*m, 1.8258*m, 1.9148*m, 1.9383*m, 1.9487*m, 1.8581*m, 1.7619*m, 1.7683*m, 1.8172*m, 1.8018*m, 1.7002*m, 1.7925*m, 1.7159*m, 1.6551*m, 1.7413*m, 1.5945*m, 1.6542*m, 1.6433*m, 1.5527*m, 1.5593*m, 1.671*m, 1.5304*m, 1.5467*m, 1.4931*m, 1.6024*m, 1.5673*m, 1.5206*m, 1.4631*m, 1.5547*m, 1.5513*m, 1.4743*m, 1.523*m, 1.4582*m, 1.5248*m, 1.4168*m, 1.464*m, 1.4744*m, 1.4284*m, 1.4361*m, 1.3864*m, 1.4581*m, 1.3333*m, 1.4237*m, 1.349*m, 1.3357*m, 1.378*m, 1.2898*m, 1.3325*m, 1.3272*m, 1.3058*m, 1.2827*m, 1.3758*m, 1.3142*m, 1.2393*m, 1.2683*m, 1.3054*m, 1.2293*m, 1.1964*m, 1.2828*m, 1.1616*m, 1.1957*m, 1.1761*m, 1.1958*m, 1.157*m, 1.1134*m, 1.1633*m, 1.1478*m, 1.0993*m, 1.065*m, 1.0981*m, 1.1051*m, 1.1044*m, 1.073*m, 1.0542*m, 1.0665*m, 1.0497*m, 1.0258*m, 1.0736*m, 1.0503*m, 1.0196*m, 1.0282*m, 0.998*m, 1.0052*m, 0.9694*m, 0.9823*m, 0.9718*m, 0.9766*m, 0.9571*m, 0.9682*m, 0.9419*m, 0.9189*m, 0.9043*m, 0.9187*m, 0.9191*m, 0.8645*m, 0.8421*m, 0.8696*m, 0.848*m, 0.8537*m, 0.8331*m, 0.8222*m, 0.809*m, 0.8345*m, 0.7941*m, 0.78*m, 0.7852*m, 0.7491*m, 0.7474*m, 0.7366*m, 0.7336*m, 0.7226*m, 0.707*m, 0.6679*m, 0.6626*m, 0.6572*m, 0.6431*m, 0.6363*m, 0.6133*m, 0.5894*m, 0.5728*m, 0.5552*m, 0.5314*m, 0.5056*m, 0.4823*m, 0.4557*m, 0.4191*m, 0.384*m, 0.3579*m, 0.3244*m, 0.2899*m, 0.2624*m, 0.2345*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2068*m, 0.1825*m, 0.1602*m, 0.1388*m, 0.1169*m, 0.0986*m, 0.0811*m, 0.0641*m, 0.0517*m, 0.0412*m, 0.0323*m, 0.025*m, 0.0184*m, 0.0139*m, 0.0105*m, 0.0080019*m, 0.0059845*m, 0.0043457*m, 0.0033135*m, 0.0026238*m, 0.002055*m, 0.0016341*m, 0.001324*m, 0.0011471*m, 0.001045*m, 0.00083358*m, 0.00073219*m, 0.00065561*m, 0.00060557*m, 0.00056072*m, 0.00053119*m, 0.00051251*m, 0.00050129*m, 0.00049617*m, 0.00049277*m, 0.00049302*m, 0.00049258*m, 0.00049236*m, 0.00049046*m, 0.00048719*m, 0.00048131*m, 0.00047382*m, 0.00046379*m, 0.00045296*m, 0.00043961*m, 0.00042243*m, 0.00040243*m, 0.0003824*m, 0.00036515*m, 0.000347*m, 0.00033252*m, 0.00031777*m, 0.00031137*m, 0.00030581*m, 0.00030298*m, 0.00030402*m, 0.00030572*m, 0.00031066*m, 0.00031407*m, 0.00031963*m, 0.00032279*m, 0.00032898*m, 0.00033276*m, 0.00033569*m, 0.00033949*m, 0.00033997*m, 0.00034336*m, 0.00034291*m, 0.0003438*m, 0.00034185*m, 0.00034068*m, 0.00034032*m, 0.00033942*m, 0.00034101*m, 0.00034385*m, 0.00035177*m, 0.00035881*m, 0.00037106*m, 0.00038297*m, 0.0003972*m, 0.00040592*m, 0.00041988*m, 0.00043239*m, 0.00044605*m, 0.00045891*m, 0.0004707*m, 0.00048077*m, 0.00049839*m, 0.00050975*m, 0.00052638*m, 0.00053552*m, 0.00055531*m, 0.00057243*m, 0.00059237*m, 0.00062231*m, 0.00064273*m, 0.00067408*m, 0.00070329*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 136.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.1197*m, 2.3849*m, 2.2083*m, 2.2568*m, 2.3726*m, 2.216*m, 2.0458*m, 2.2289*m, 2.1889*m, 1.9814*m, 1.9763*m, 1.9043*m, 2.025*m, 2.1358*m, 2.097*m, 2.0242*m, 2.0691*m, 2.3388*m, 2.0819*m, 1.9578*m, 2.0057*m, 2.3748*m, 1.8968*m, 1.9645*m, 2.1645*m, 1.9904*m, 2.0436*m, 1.9508*m, 2.274*m, 1.7629*m, 1.8862*m, 1.8082*m, 1.8322*m, 1.7217*m, 1.725*m, 1.827*m, 1.7375*m, 1.8684*m, 1.78*m, 1.7371*m, 1.8227*m, 1.8453*m, 1.8554*m, 1.7681*m, 1.6757*m, 1.6818*m, 1.7289*m, 1.714*m, 1.6164*m, 1.7051*m, 1.6315*m, 1.5731*m, 1.6559*m, 1.5149*m, 1.5722*m, 1.5617*m, 1.4749*m, 1.4812*m, 1.5884*m, 1.4535*m, 1.4691*m, 1.4177*m, 1.5226*m, 1.4889*m, 1.4441*m, 1.389*m, 1.4768*m, 1.4735*m, 1.3997*m, 1.4464*m, 1.3843*m, 1.4481*m, 1.3447*m, 1.3898*m, 1.3998*m, 1.3558*m, 1.3631*m, 1.3156*m, 1.3843*m, 1.2648*m, 1.3513*m, 1.2798*m, 1.2671*m, 1.3075*m, 1.2232*m, 1.264*m, 1.259*m, 1.2386*m, 1.2164*m, 1.3054*m, 1.2465*m, 1.175*m, 1.2027*m, 1.2382*m, 1.1654*m, 1.134*m, 1.2165*m, 1.1008*m, 1.1334*m, 1.1147*m, 1.1334*m, 1.0964*m, 1.0548*m, 1.1025*m, 1.0876*m, 1.0414*m, 1.0087*m, 1.0403*m, 1.0469*m, 1.0463*m, 1.0163*m, 0.9984*m, 1.0101*m, 0.9941*m, 0.9713*m, 1.0169*m, 0.9947*m, 0.9654*m, 0.9736*m, 0.9449*m, 0.9517*m, 0.9176*m, 0.9299*m, 0.9199*m, 0.9245*m, 0.9059*m, 0.9165*m, 0.8915*m, 0.8695*m, 0.8557*m, 0.8693*m, 0.8697*m, 0.8178*m, 0.7965*m, 0.8226*m, 0.8021*m, 0.8076*m, 0.788*m, 0.7776*m, 0.7651*m, 0.7892*m, 0.7509*m, 0.7375*m, 0.7425*m, 0.7081*m, 0.7066*m, 0.6963*m, 0.6934*m, 0.683*m, 0.6682*m, 0.6311*m, 0.626*m, 0.6209*m, 0.6075*m, 0.6011*m, 0.5793*m, 0.5567*m, 0.5409*m, 0.5243*m, 0.5017*m, 0.4773*m, 0.4552*m, 0.43*m, 0.3954*m, 0.3622*m, 0.3375*m, 0.3059*m, 0.2734*m, 0.2474*m, 0.221*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.1949*m, 0.1719*m, 0.1509*m, 0.1308*m, 0.1101*m, 0.0928*m, 0.0763*m, 0.0604*m, 0.0487*m, 0.0387*m, 0.0304*m, 0.0235*m, 0.0173*m, 0.0131*m, 0.0098614*m, 0.0075315*m, 0.0056326*m, 0.0040901*m, 0.0031187*m, 0.0024695*m, 0.0019341*m, 0.001538*m, 0.0012461*m, 0.0010796*m, 0.00098354*m, 0.00078454*m, 0.00068912*m, 0.00061705*m, 0.00056995*m, 0.00052774*m, 0.00049995*m, 0.00048237*m, 0.0004718*m, 0.00046698*m, 0.00046378*m, 0.00046402*m, 0.0004636*m, 0.0004634*m, 0.00046161*m, 0.00045853*m, 0.000453*m, 0.00044595*m, 0.00043651*m, 0.00042632*m, 0.00041375*m, 0.00039759*m, 0.00037876*m, 0.00035991*m, 0.00034367*m, 0.00032659*m, 0.00031296*m, 0.00029908*m, 0.00029305*m, 0.00028782*m, 0.00028516*m, 0.00028613*m, 0.00028774*m, 0.00029239*m, 0.0002956*m, 0.00030083*m, 0.0003038*m, 0.00030962*m, 0.00031319*m, 0.00031594*m, 0.00031952*m, 0.00031997*m, 0.00032316*m, 0.00032274*m, 0.00032358*m, 0.00032174*m, 0.00032064*m, 0.0003203*m, 0.00031945*m, 0.00032095*m, 0.00032362*m, 0.00033107*m, 0.0003377*m, 0.00034923*m, 0.00036044*m, 0.00037383*m, 0.00038205*m, 0.00039518*m, 0.00040696*m, 0.00041981*m, 0.00043192*m, 0.00044301*m, 0.00045249*m, 0.00046907*m, 0.00047977*m, 0.00049542*m, 0.00050402*m, 0.00052265*m, 0.00053876*m, 0.00055752*m, 0.0005857*m, 0.00060493*m, 0.00063443*m, 0.00066192*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 144.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.0258*m, 2.2827*m, 2.1116*m, 2.1586*m, 2.2707*m, 2.1191*m, 1.9544*m, 2.1315*m, 2.0928*m, 1.8922*m, 1.8872*m, 1.8178*m, 1.9343*m, 2.0414*m, 2.0039*m, 1.9336*m, 1.9769*m, 2.2381*m, 1.9893*m, 1.8694*m, 1.9156*m, 2.2729*m, 1.8105*m, 1.8759*m, 2.0692*m, 1.9009*m, 1.9523*m, 1.8627*m, 2.1752*m, 1.6815*m, 1.8003*m, 1.7251*m, 1.7482*m, 1.6418*m, 1.645*m, 1.7433*m, 1.657*m, 1.7832*m, 1.698*m, 1.6567*m, 1.7391*m, 1.7609*m, 1.7706*m, 1.6865*m, 1.5975*m, 1.6034*m, 1.6487*m, 1.6343*m, 1.5404*m, 1.6258*m, 1.555*m, 1.4988*m, 1.5784*m, 1.4429*m, 1.498*m, 1.4879*m, 1.4045*m, 1.4106*m, 1.5135*m, 1.3839*m, 1.3989*m, 1.3496*m, 1.4503*m, 1.4179*m, 1.375*m, 1.3221*m, 1.4063*m, 1.4031*m, 1.3323*m, 1.3771*m, 1.3176*m, 1.3788*m, 1.2796*m, 1.3229*m, 1.3324*m, 1.2902*m, 1.2972*m, 1.2516*m, 1.3175*m, 1.203*m, 1.2859*m, 1.2174*m, 1.2052*m, 1.2439*m, 1.1632*m, 1.2022*m, 1.1974*m, 1.1779*m, 1.1567*m, 1.2419*m, 1.1855*m, 1.117*m, 1.1435*m, 1.1775*m, 1.1078*m, 1.0778*m, 1.1568*m, 1.0461*m, 1.0772*m, 1.0593*m, 1.0773*m, 1.0419*m, 1.0021*m, 1.0476*m, 1.0335*m, 0.9893*m, 0.958*m, 0.9882*m, 0.9946*m, 0.9939*m, 0.9653*m, 0.9482*m, 0.9594*m, 0.9441*m, 0.9224*m, 0.9659*m, 0.9447*m, 0.9167*m, 0.9245*m, 0.8971*m, 0.9036*m, 0.8711*m, 0.8828*m, 0.8733*m, 0.8777*m, 0.8599*m, 0.87*m, 0.8462*m, 0.8252*m, 0.812*m, 0.825*m, 0.8254*m, 0.7759*m, 0.7556*m, 0.7805*m, 0.7609*m, 0.7661*m, 0.7475*m, 0.7376*m, 0.7257*m, 0.7487*m, 0.7122*m, 0.6994*m, 0.7041*m, 0.6715*m, 0.6699*m, 0.6601*m, 0.6574*m, 0.6475*m, 0.6334*m, 0.5981*m, 0.5933*m, 0.5885*m, 0.5757*m, 0.5696*m, 0.5489*m, 0.5274*m, 0.5124*m, 0.4966*m, 0.4751*m, 0.4519*m, 0.431*m, 0.4071*m, 0.3743*m, 0.3428*m, 0.3194*m, 0.2894*m, 0.2586*m, 0.234*m, 0.209*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.1843*m, 0.1625*m, 0.1427*m, 0.1236*m, 0.104*m, 0.0877*m, 0.0721*m, 0.057*m, 0.046*m, 0.0366*m, 0.0287*m, 0.0222*m, 0.0163*m, 0.0124*m, 0.009314*m, 0.0071134*m, 0.0053199*m, 0.003863*m, 0.0029455*m, 0.0023323*m, 0.0018267*m, 0.0014525*m, 0.0011769*m, 0.0010197*m, 0.00092891*m, 0.00074096*m, 0.00065084*m, 0.00058276*m, 0.00053828*m, 0.00049842*m, 0.00047217*m, 0.00045557*m, 0.00044559*m, 0.00044104*m, 0.00043802*m, 0.00043824*m, 0.00043785*m, 0.00043766*m, 0.00043596*m, 0.00043305*m, 0.00042783*m, 0.00042118*m, 0.00041226*m, 0.00040263*m, 0.00039077*m, 0.0003755*m, 0.00035772*m, 0.00033991*m, 0.00032458*m, 0.00030844*m, 0.00029557*m, 0.00028246*m, 0.00027677*m, 0.00027183*m, 0.00026932*m, 0.00027024*m, 0.00027175*m, 0.00027615*m, 0.00027918*m, 0.00028412*m, 0.00028692*m, 0.00029242*m, 0.00029579*m, 0.00029839*m, 0.00030177*m, 0.0003022*m, 0.00030521*m, 0.00030481*m, 0.0003056*m, 0.00030386*m, 0.00030283*m, 0.00030251*m, 0.00030171*m, 0.00030312*m, 0.00030564*m, 0.00031268*m, 0.00031894*m, 0.00032983*m, 0.00034041*m, 0.00035306*m, 0.00036082*m, 0.00037323*m, 0.00038435*m, 0.00039649*m, 0.00040792*m, 0.0004184*m, 0.00042735*m, 0.00044301*m, 0.00045311*m, 0.0004679*m, 0.00047602*m, 0.00049361*m, 0.00050883*m, 0.00052655*m, 0.00055316*m, 0.00057132*m, 0.00059918*m, 0.00062515*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 152.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 1.94*m, 2.1889*m, 2.023*m, 2.0685*m, 2.1773*m, 2.0302*m, 1.8709*m, 2.0423*m, 2.0048*m, 1.8107*m, 1.8059*m, 1.7388*m, 1.8514*m, 1.955*m, 1.9187*m, 1.8507*m, 1.8926*m, 2.1456*m, 1.9046*m, 1.7886*m, 1.8333*m, 2.1794*m, 1.7318*m, 1.7949*m, 1.9819*m, 1.8191*m, 1.8688*m, 1.7822*m, 2.0846*m, 1.6072*m, 1.7219*m, 1.6493*m, 1.6716*m, 1.569*m, 1.572*m, 1.6668*m, 1.5837*m, 1.7053*m, 1.6231*m, 1.5833*m, 1.6628*m, 1.6839*m, 1.6933*m, 1.6121*m, 1.5263*m, 1.5319*m, 1.5756*m, 1.5618*m, 1.4713*m, 1.5535*m, 1.4853*m, 1.4312*m, 1.5079*m, 1.3775*m, 1.4304*m, 1.4208*m, 1.3405*m, 1.3463*m, 1.4454*m, 1.3207*m, 1.3352*m, 1.2878*m, 1.3845*m, 1.3534*m, 1.3121*m, 1.2613*m, 1.3423*m, 1.3392*m, 1.2711*m, 1.3142*m, 1.257*m, 1.3158*m, 1.2205*m, 1.2621*m, 1.2712*m, 1.2307*m, 1.2374*m, 1.1936*m, 1.2569*m, 1.147*m, 1.2266*m, 1.1608*m, 1.1491*m, 1.1862*m, 1.1088*m, 1.1462*m, 1.1416*m, 1.1229*m, 1.1025*m, 1.1843*m, 1.1302*m, 1.0645*m, 1.0899*m, 1.1225*m, 1.0557*m, 1.0269*m, 1.1026*m, 0.9965*m, 1.0264*m, 1.0092*m, 1.0264*m, 0.9925*m, 0.9544*m, 0.998*m, 0.9844*m, 0.9421*m, 0.9122*m, 0.9411*m, 0.9472*m, 0.9466*m, 0.9192*m, 0.9028*m, 0.9136*m, 0.8989*m, 0.8781*m, 0.9197*m, 0.8994*m, 0.8727*m, 0.8802*m, 0.8539*m, 0.8601*m, 0.829*m, 0.8403*m, 0.8311*m, 0.8353*m, 0.8184*m, 0.828*m, 0.8052*m, 0.7852*m, 0.7726*m, 0.785*m, 0.7854*m, 0.7381*m, 0.7187*m, 0.7425*m, 0.7238*m, 0.7288*m, 0.7109*m, 0.7015*m, 0.6901*m, 0.7121*m, 0.6773*m, 0.665*m, 0.6695*m, 0.6384*m, 0.6369*m, 0.6276*m, 0.625*m, 0.6155*m, 0.6021*m, 0.5684*m, 0.5639*m, 0.5592*m, 0.5471*m, 0.5412*m, 0.5215*m, 0.501*m, 0.4868*m, 0.4717*m, 0.4513*m, 0.4292*m, 0.4092*m, 0.3865*m, 0.3553*m, 0.3253*m, 0.3031*m, 0.2746*m, 0.2453*m, 0.2219*m, 0.1982*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.1747*m, 0.1541*m, 0.1353*m, 0.1172*m, 0.0986*m, 0.0832*m, 0.0684*m, 0.0541*m, 0.0436*m, 0.0347*m, 0.0272*m, 0.021*m, 0.0155*m, 0.0117*m, 0.0088243*m, 0.0067393*m, 0.00504*m, 0.0036598*m, 0.0027905*m, 0.0022096*m, 0.0017306*m, 0.0013761*m, 0.001115*m, 0.000966*m, 0.00088002*m, 0.00070196*m, 0.00061658*m, 0.00055209*m, 0.00050995*m, 0.00047219*m, 0.00044732*m, 0.00043159*m, 0.00042214*m, 0.00041782*m, 0.00041496*m, 0.00041517*m, 0.0004148*m, 0.00041462*m, 0.00041302*m, 0.00041026*m, 0.00040531*m, 0.00039901*m, 0.00039056*m, 0.00038144*m, 0.0003702*m, 0.00035573*m, 0.00033889*m, 0.00032202*m, 0.0003075*m, 0.00029221*m, 0.00028001*m, 0.0002676*m, 0.0002622*m, 0.00025753*m, 0.00025514*m, 0.00025602*m, 0.00025745*m, 0.00026161*m, 0.00026448*m, 0.00026916*m, 0.00027182*m, 0.00027703*m, 0.00028022*m, 0.00028268*m, 0.00028588*m, 0.00028629*m, 0.00028914*m, 0.00028876*m, 0.00028952*m, 0.00028787*m, 0.00028689*m, 0.00028658*m, 0.00028583*m, 0.00028717*m, 0.00028956*m, 0.00029622*m, 0.00030215*m, 0.00031247*m, 0.0003225*m, 0.00033448*m, 0.00034183*m, 0.00035359*m, 0.00036412*m, 0.00037562*m, 0.00038645*m, 0.00039638*m, 0.00040486*m, 0.0004197*m, 0.00042926*m, 0.00044327*m, 0.00045097*m, 0.00046763*m, 0.00048205*m, 0.00049883*m, 0.00052405*m, 0.00054125*m, 0.00056765*m, 0.00059224*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 160.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 1.861*m, 2.1025*m, 1.9415*m, 1.9857*m, 2.0913*m, 1.9486*m, 1.7941*m, 1.9602*m, 1.9239*m, 1.7359*m, 1.7313*m, 1.6664*m, 1.7753*m, 1.8757*m, 1.8405*m, 1.7746*m, 1.8152*m, 2.0605*m, 1.8268*m, 1.7146*m, 1.7578*m, 2.0933*m, 1.6596*m, 1.7206*m, 1.9017*m, 1.744*m, 1.7921*m, 1.7083*m, 2.0013*m, 1.5393*m, 1.65*m, 1.5799*m, 1.6015*m, 1.5024*m, 1.5053*m, 1.5968*m, 1.5165*m, 1.634*m, 1.5546*m, 1.5162*m, 1.593*m, 1.6133*m, 1.6224*m, 1.544*m, 1.4611*m, 1.4666*m, 1.5088*m, 1.4954*m, 1.4081*m, 1.4874*m, 1.4216*m, 1.3695*m, 1.4434*m, 1.3177*m, 1.3687*m, 1.3594*m, 1.2821*m, 1.2877*m, 1.3831*m, 1.263*m, 1.2769*m, 1.2313*m, 1.3245*m, 1.2945*m, 1.2548*m, 1.2059*m, 1.2838*m, 1.2808*m, 1.2153*m, 1.2567*m, 1.2017*m, 1.2583*m, 1.1666*m, 1.2066*m, 1.2154*m, 1.1764*m, 1.1829*m, 1.1408*m, 1.2016*m, 1.0959*m, 1.1724*m, 1.1092*m, 1.0979*m, 1.1337*m, 1.0592*m, 1.0952*m, 1.0908*m, 1.0727*m, 1.0532*m, 1.1318*m, 1.0798*m, 1.0167*m, 1.0411*m, 1.0724*m, 1.0082*m, 0.9806*m, 1.0533*m, 0.9515*m, 0.9801*m, 0.9636*m, 0.9801*m, 0.9476*m, 0.9111*m, 0.9529*m, 0.9398*m, 0.8993*m, 0.8706*m, 0.8983*m, 0.9041*m, 0.9036*m, 0.8773*m, 0.8616*m, 0.8719*m, 0.8578*m, 0.8379*m, 0.8778*m, 0.8583*m, 0.8327*m, 0.8399*m, 0.8147*m, 0.8207*m, 0.7909*m, 0.8016*m, 0.7929*m, 0.7969*m, 0.7807*m, 0.7899*m, 0.7681*m, 0.7489*m, 0.7368*m, 0.7487*m, 0.7491*m, 0.7038*m, 0.6853*m, 0.708*m, 0.6901*m, 0.6949*m, 0.6778*m, 0.6688*m, 0.6579*m, 0.6789*m, 0.6456*m, 0.6339*m, 0.6382*m, 0.6084*m, 0.607*m, 0.5981*m, 0.5956*m, 0.5865*m, 0.5737*m, 0.5415*m, 0.5372*m, 0.5328*m, 0.5212*m, 0.5156*m, 0.4967*m, 0.4772*m, 0.4636*m, 0.4492*m, 0.4297*m, 0.4086*m, 0.3896*m, 0.3679*m, 0.3381*m, 0.3096*m, 0.2884*m, 0.2612*m, 0.2333*m, 0.2111*m, 0.1885*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, 0., 0.};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.1661*m, 0.1465*m, 0.1286*m, 0.1114*m, 0.0937*m, 0.079*m, 0.065*m, 0.0514*m, 0.0414*m, 0.033*m, 0.0258*m, 0.02*m, 0.0147*m, 0.0111*m, 0.0083834*m, 0.0064025*m, 0.0047881*m, 0.0034768*m, 0.002651*m, 0.0020991*m, 0.0016441*m, 0.0013073*m, 0.0010592*m, 0.0009177*m, 0.00083603*m, 0.00066686*m, 0.00058575*m, 0.00052449*m, 0.00048445*m, 0.00044858*m, 0.00042495*m, 0.00041001*m, 0.00040103*m, 0.00039693*m, 0.00039422*m, 0.00039442*m, 0.00039406*m, 0.00039389*m, 0.00039237*m, 0.00038975*m, 0.00038505*m, 0.00037906*m, 0.00037103*m, 0.00036237*m, 0.00035169*m, 0.00033795*m, 0.00032195*m, 0.00030592*m, 0.00029212*m, 0.0002776*m, 0.00026601*m, 0.00025422*m, 0.00024909*m, 0.00024465*m, 0.00024239*m, 0.00024321*m, 0.00024458*m, 0.00024853*m, 0.00025126*m, 0.0002557*m, 0.00025823*m, 0.00026318*m, 0.00026621*m, 0.00026855*m, 0.00027159*m, 0.00027198*m, 0.00027468*m, 0.00027433*m, 0.00027504*m, 0.00027348*m, 0.00027255*m, 0.00027226*m, 0.00027154*m, 0.00027281*m, 0.00027508*m, 0.00028141*m, 0.00028705*m, 0.00029685*m, 0.00030637*m, 0.00031776*m, 0.00032474*m, 0.00033591*m, 0.00034591*m, 0.00035684*m, 0.00036713*m, 0.00037656*m, 0.00038462*m, 0.00039871*m, 0.0004078*m, 0.00042111*m, 0.00042842*m, 0.00044425*m, 0.00045795*m, 0.00047389*m, 0.00049784*m, 0.00051419*m, 0.00053926*m, 0.00056263*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    G4double cromophore_concentration_ = FindClosestNumber(cromophore_concentration, &available_concentrations);
    if(verbosity){
      G4cout << "Taking the absorption length spectra for " << cromophore_concentration_ << " mg/kg." << G4endl;
    } 

    mpt->AddProperty("ABSLENGTH", abs_energy, absLength[cromophore_concentration_]);
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength[cromophore_concentration_]);

    //std::vector<G4double> absLength = {
    //  attenuation_length,  attenuation_length,
    //  attenuation_length,  attenuation_length,
    //  attenuation_length,  attenuation_length
    //};

    if(verbosity){
      G4cout << "Attenuation length energy entries = " << abs_energy.size() << G4endl;
      G4cout << "Attenuation length entries = " << absLength[cromophore_concentration_].size() << G4endl;
      G4cout << "WLS absorption length energy entries = " << WLS_abs_energy.size() << G4endl;
      G4cout << "WLS absorption length entries = " << WLS_absLength[cromophore_concentration_].size() << G4endl;
    } 

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

  G4MaterialPropertiesTable* FakeG2P_FB118(G4double cromophore_concentration, G4double rindex, G4bool verbosity)
  {
    // Exactly the same G4MaterialPropertiesTable as that of &G2P_FB118(), but with its absorption and emission
    // spectra shifted 55 nm towards higher wavelengths: the WLS absorption length spectrum, the attenuation 
    // length spectrum and the emission spectrum are shifted 55 nm towards bigger wavelengths. 55 nm is the 
    // spectral width between the absorption peak and the emission peak of G2P_FB118(). Thus, this fake material
    // is designed so as to absorb the emission of G2P_FB118() and shift way further into the visible part of the
    // EM spectrum. This fake material could be useful to implement a proof of concept of a Dual Scintillator
    // XArapuca.

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

    G4cout << "rindex=" << rindex << G4endl;

    std::vector<G4double> rIndex = {
      rindex,
      rindex,  rindex,   // 609 , 589.26 nm
      rindex,  rindex,   // 550 , 530 nm
      rindex,  rindex,   // 500 , 490 nm
      rindex,  rindex,   // 481 , 460 nm
      rindex,  rindex,   // 435 , 425 nm
      rindex
    };
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // abs_energy got 324 entries
    std::vector<G4double> abs_energy = {optPhotMinE_, h_Planck*c_light/(656.073*nm), h_Planck*c_light/(654.073*nm), h_Planck*c_light/(653.034*nm), h_Planck*c_light/(651.994*nm), h_Planck*c_light/(650.954*nm), h_Planck*c_light/(650.063*nm), h_Planck*c_light/(649.022*nm), h_Planck*c_light/(647.982*nm), h_Planck*c_light/(646.941*nm), h_Planck*c_light/(646.048*nm), h_Planck*c_light/(645.007*nm), h_Planck*c_light/(643.965*nm), h_Planck*c_light/(643.072*nm), h_Planck*c_light/(642.029*nm), h_Planck*c_light/(640.987*nm), h_Planck*c_light/(639.944*nm), h_Planck*c_light/(639.05*nm), h_Planck*c_light/(638.006*nm), h_Planck*c_light/(636.963*nm), h_Planck*c_light/(636.068*nm), h_Planck*c_light/(635.024*nm), h_Planck*c_light/(633.979*nm), h_Planck*c_light/(632.935*nm), h_Planck*c_light/(632.039*nm), h_Planck*c_light/(630.994*nm), h_Planck*c_light/(629.948*nm), h_Planck*c_light/(629.052*nm), h_Planck*c_light/(628.006*nm), h_Planck*c_light/(626.96*nm), h_Planck*c_light/(626.063*nm), h_Planck*c_light/(625.016*nm), h_Planck*c_light/(623.969*nm), h_Planck*c_light/(623.071*nm), h_Planck*c_light/(622.024*nm), h_Planck*c_light/(620.976*nm), h_Planck*c_light/(619.928*nm), h_Planck*c_light/(619.029*nm), h_Planck*c_light/(617.98*nm), h_Planck*c_light/(616.932*nm), h_Planck*c_light/(616.032*nm), h_Planck*c_light/(614.983*nm), h_Planck*c_light/(613.934*nm), h_Planck*c_light/(613.034*nm), h_Planck*c_light/(611.984*nm), h_Planck*c_light/(610.933*nm), h_Planck*c_light/(610.033*nm), h_Planck*c_light/(608.982*nm), h_Planck*c_light/(607.931*nm), h_Planck*c_light/(607.03*nm), h_Planck*c_light/(605.978*nm), h_Planck*c_light/(604.926*nm), h_Planck*c_light/(604.024*nm), h_Planck*c_light/(602.972*nm), h_Planck*c_light/(602.07*nm), h_Planck*c_light/(601.017*nm), h_Planck*c_light/(599.964*nm), h_Planck*c_light/(599.061*nm), h_Planck*c_light/(598.008*nm), h_Planck*c_light/(596.954*nm), h_Planck*c_light/(596.05*nm), h_Planck*c_light/(594.996*nm), h_Planck*c_light/(593.942*nm), h_Planck*c_light/(593.038*nm), h_Planck*c_light/(591.982*nm), h_Planck*c_light/(590.927*nm), h_Planck*c_light/(590.022*nm), h_Planck*c_light/(588.967*nm), h_Planck*c_light/(588.062*nm), h_Planck*c_light/(587.006*nm), h_Planck*c_light/(585.949*nm), h_Planck*c_light/(585.043*nm), h_Planck*c_light/(583.986*nm), h_Planck*c_light/(582.929*nm), h_Planck*c_light/(582.023*nm), h_Planck*c_light/(580.965*nm), h_Planck*c_light/(580.059*nm), h_Planck*c_light/(579.0*nm), h_Planck*c_light/(577.942*nm), h_Planck*c_light/(577.035*nm), h_Planck*c_light/(575.976*nm), h_Planck*c_light/(575.068*nm), h_Planck*c_light/(574.009*nm), h_Planck*c_light/(572.95*nm), h_Planck*c_light/(572.041*nm), h_Planck*c_light/(570.981*nm), h_Planck*c_light/(570.073*nm), h_Planck*c_light/(569.012*nm), h_Planck*c_light/(567.952*nm), h_Planck*c_light/(567.042*nm), h_Planck*c_light/(565.981*nm), h_Planck*c_light/(565.071*nm), h_Planck*c_light/(564.01*nm), h_Planck*c_light/(562.948*nm), h_Planck*c_light/(562.038*nm), h_Planck*c_light/(560.976*nm), h_Planck*c_light/(560.065*nm), h_Planck*c_light/(559.002*nm), h_Planck*c_light/(557.939*nm), h_Planck*c_light/(557.028*nm), h_Planck*c_light/(555.965*nm), h_Planck*c_light/(555.053*nm), h_Planck*c_light/(553.989*nm), h_Planck*c_light/(552.925*nm), h_Planck*c_light/(552.013*nm), h_Planck*c_light/(550.949*nm), h_Planck*c_light/(550.036*nm), h_Planck*c_light/(548.971*nm), h_Planck*c_light/(548.058*nm), h_Planck*c_light/(546.993*nm), h_Planck*c_light/(545.927*nm), h_Planck*c_light/(545.014*nm), h_Planck*c_light/(543.948*nm), h_Planck*c_light/(543.034*nm), h_Planck*c_light/(541.968*nm), h_Planck*c_light/(541.053*nm), h_Planck*c_light/(539.986*nm), h_Planck*c_light/(539.072*nm), h_Planck*c_light/(538.004*nm), h_Planck*c_light/(536.937*nm), h_Planck*c_light/(536.022*nm), h_Planck*c_light/(534.954*nm), h_Planck*c_light/(534.038*nm), h_Planck*c_light/(532.97*nm), h_Planck*c_light/(532.054*nm), h_Planck*c_light/(530.985*nm), h_Planck*c_light/(530.069*nm), h_Planck*c_light/(529.0*nm), h_Planck*c_light/(527.93*nm), h_Planck*c_light/(527.014*nm), h_Planck*c_light/(525.944*nm), h_Planck*c_light/(525.027*nm), h_Planck*c_light/(523.956*nm), h_Planck*c_light/(523.039*nm), h_Planck*c_light/(521.968*nm), h_Planck*c_light/(521.05*nm), h_Planck*c_light/(519.979*nm), h_Planck*c_light/(519.061*nm), h_Planck*c_light/(517.99*nm), h_Planck*c_light/(517.071*nm), h_Planck*c_light/(515.999*nm), h_Planck*c_light/(514.927*nm), h_Planck*c_light/(514.008*nm), h_Planck*c_light/(512.936*nm), h_Planck*c_light/(512.016*nm), h_Planck*c_light/(510.943*nm), h_Planck*c_light/(510.023*nm), h_Planck*c_light/(508.95*nm), h_Planck*c_light/(508.03*nm), h_Planck*c_light/(506.956*nm), h_Planck*c_light/(506.036*nm), h_Planck*c_light/(504.962*nm), h_Planck*c_light/(504.041*nm), h_Planck*c_light/(502.966*nm), h_Planck*c_light/(502.045*nm), h_Planck*c_light/(500.97*nm), h_Planck*c_light/(500.049*nm), h_Planck*c_light/(498.974*nm), h_Planck*c_light/(498.052*nm), h_Planck*c_light/(496.976*nm), h_Planck*c_light/(496.054*nm), h_Planck*c_light/(494.978*nm), h_Planck*c_light/(494.056*nm), h_Planck*c_light/(492.979*nm), h_Planck*c_light/(492.056*nm), h_Planck*c_light/(490.98*nm), h_Planck*c_light/(490.056*nm), h_Planck*c_light/(488.979*nm), h_Planck*c_light/(488.056*nm), h_Planck*c_light/(486.978*nm), h_Planck*c_light/(486.054*nm), h_Planck*c_light/(484.976*nm), h_Planck*c_light/(484.044*nm), h_Planck*c_light/(482.967*nm), h_Planck*c_light/(482.043*nm), h_Planck*c_light/(480.965*nm), h_Planck*c_light/(480.041*nm), h_Planck*c_light/(478.963*nm), h_Planck*c_light/(478.039*nm), h_Planck*c_light/(476.96*nm), h_Planck*c_light/(476.036*nm), h_Planck*c_light/(474.957*nm), h_Planck*c_light/(474.032*nm), h_Planck*c_light/(472.952*nm), h_Planck*c_light/(472.027*nm), h_Planck*c_light/(470.947*nm), h_Planck*c_light/(470.022*nm), h_Planck*c_light/(468.942*nm), h_Planck*c_light/(468.016*nm), h_Planck*c_light/(466.936*nm), h_Planck*c_light/(466.009*nm), h_Planck*c_light/(464.928*nm), h_Planck*c_light/(464.002*nm), h_Planck*c_light/(463.075*nm), h_Planck*c_light/(461.994*nm), h_Planck*c_light/(461.067*nm), h_Planck*c_light/(459.985*nm), h_Planck*c_light/(459.058*nm), h_Planck*c_light/(457.976*nm), h_Planck*c_light/(457.049*nm), h_Planck*c_light/(455.966*nm), h_Planck*c_light/(455.038*nm), h_Planck*c_light/(453.956*nm), h_Planck*c_light/(453.027*nm), h_Planck*c_light/(451.944*nm), h_Planck*c_light/(451.016*nm), h_Planck*c_light/(449.932*nm), h_Planck*c_light/(449.004*nm), h_Planck*c_light/(448.075*nm), h_Planck*c_light/(446.991*nm), h_Planck*c_light/(446.062*nm), h_Planck*c_light/(444.977*nm), h_Planck*c_light/(444.048*nm), h_Planck*c_light/(442.963*nm), h_Planck*c_light/(442.034*nm), h_Planck*c_light/(440.948*nm), h_Planck*c_light/(440.018*nm), h_Planck*c_light/(438.933*nm), h_Planck*c_light/(438.003*nm), h_Planck*c_light/(437.072*nm), h_Planck*c_light/(435.986*nm), h_Planck*c_light/(435.056*nm), h_Planck*c_light/(433.97*nm), h_Planck*c_light/(433.038*nm), h_Planck*c_light/(431.952*nm), h_Planck*c_light/(431.021*nm), h_Planck*c_light/(429.934*nm), h_Planck*c_light/(429.002*nm), h_Planck*c_light/(428.071*nm), h_Planck*c_light/(426.983*nm), h_Planck*c_light/(426.051*nm), h_Planck*c_light/(424.964*nm), h_Planck*c_light/(424.031*nm), h_Planck*c_light/(422.944*nm), h_Planck*c_light/(422.011*nm), h_Planck*c_light/(420.923*nm), h_Planck*c_light/(419.99*nm), h_Planck*c_light/(419.057*nm), h_Planck*c_light/(417.968*nm), h_Planck*c_light/(417.035*nm), h_Planck*c_light/(415.946*nm), h_Planck*c_light/(415.012*nm), h_Planck*c_light/(413.923*nm), h_Planck*c_light/(412.989*nm), h_Planck*c_light/(412.055*nm), h_Planck*c_light/(410.965*nm), h_Planck*c_light/(410.031*nm), h_Planck*c_light/(408.941*nm), h_Planck*c_light/(408.006*nm), h_Planck*c_light/(407.072*nm), h_Planck*c_light/(405.981*nm), h_Planck*c_light/(405.046*nm), h_Planck*c_light/(403.956*nm), h_Planck*c_light/(403.02*nm), h_Planck*c_light/(401.929*nm), h_Planck*c_light/(400.994*nm), h_Planck*c_light/(400.058*nm), h_Planck*c_light/(398.967*nm), h_Planck*c_light/(398.031*nm), h_Planck*c_light/(396.939*nm), h_Planck*c_light/(396.003*nm), h_Planck*c_light/(395.067*nm), h_Planck*c_light/(393.974*nm), h_Planck*c_light/(393.038*nm), h_Planck*c_light/(391.945*nm), h_Planck*c_light/(391.009*nm), h_Planck*c_light/(390.072*nm), h_Planck*c_light/(388.979*nm), h_Planck*c_light/(388.042*nm), h_Planck*c_light/(386.948*nm), h_Planck*c_light/(386.011*nm), h_Planck*c_light/(385.074*nm), h_Planck*c_light/(383.98*nm), h_Planck*c_light/(383.042*nm), h_Planck*c_light/(381.948*nm), h_Planck*c_light/(381.01*nm), h_Planck*c_light/(380.072*nm), h_Planck*c_light/(378.978*nm), h_Planck*c_light/(378.04*nm), h_Planck*c_light/(376.945*nm), h_Planck*c_light/(376.006*nm), h_Planck*c_light/(375.068*nm), h_Planck*c_light/(373.972*nm), h_Planck*c_light/(373.034*nm), h_Planck*c_light/(371.938*nm), h_Planck*c_light/(370.999*nm), h_Planck*c_light/(370.06*nm), h_Planck*c_light/(368.964*nm), h_Planck*c_light/(368.025*nm), h_Planck*c_light/(366.928*nm), h_Planck*c_light/(365.989*nm), h_Planck*c_light/(365.049*nm), h_Planck*c_light/(363.952*nm), h_Planck*c_light/(363.012*nm), h_Planck*c_light/(362.072*nm), h_Planck*c_light/(360.975*nm), h_Planck*c_light/(360.035*nm), h_Planck*c_light/(358.938*nm), h_Planck*c_light/(357.997*nm), h_Planck*c_light/(357.056*nm), h_Planck*c_light/(355.959*nm), h_Planck*c_light/(355.018*nm), h_Planck*c_light/(354.077*nm), h_Planck*c_light/(352.979*nm), h_Planck*c_light/(352.038*nm), h_Planck*c_light/(350.94*nm), h_Planck*c_light/(349.998*nm), h_Planck*c_light/(349.056*nm), h_Planck*c_light/(347.958*nm), h_Planck*c_light/(347.016*nm), h_Planck*c_light/(346.074*nm), h_Planck*c_light/(344.975*nm), h_Planck*c_light/(344.033*nm), h_Planck*c_light/(342.934*nm), h_Planck*c_light/(341.992*nm), h_Planck*c_light/(341.049*nm), h_Planck*c_light/(339.95*nm), h_Planck*c_light/(339.007*nm), h_Planck*c_light/(338.064*nm), h_Planck*c_light/(336.964*nm), h_Planck*c_light/(336.021*nm), h_Planck*c_light/(335.078*nm), h_Planck*c_light/(333.078*nm), optPhotMaxE_};

    
    // WLS_abs_energy got 102 entries
    std::vector<G4double> WLS_abs_energy = {optPhotMinE_, h_Planck*c_light/(486.052*nm), h_Planck*c_light/(484.052*nm), h_Planck*c_light/(482.974*nm), h_Planck*c_light/(482.05*nm), h_Planck*c_light/(480.971*nm), h_Planck*c_light/(480.046*nm), h_Planck*c_light/(478.967*nm), h_Planck*c_light/(478.042*nm), h_Planck*c_light/(476.962*nm), h_Planck*c_light/(476.037*nm), h_Planck*c_light/(474.957*nm), h_Planck*c_light/(474.032*nm), h_Planck*c_light/(472.951*nm), h_Planck*c_light/(472.025*nm), h_Planck*c_light/(470.945*nm), h_Planck*c_light/(470.018*nm), h_Planck*c_light/(468.937*nm), h_Planck*c_light/(468.011*nm), h_Planck*c_light/(466.93*nm), h_Planck*c_light/(466.002*nm), h_Planck*c_light/(465.075*nm), h_Planck*c_light/(463.994*nm), h_Planck*c_light/(463.066*nm), h_Planck*c_light/(461.984*nm), h_Planck*c_light/(461.056*nm), h_Planck*c_light/(459.974*nm), h_Planck*c_light/(459.046*nm), h_Planck*c_light/(457.963*nm), h_Planck*c_light/(457.034*nm), h_Planck*c_light/(455.951*nm), h_Planck*c_light/(455.023*nm), h_Planck*c_light/(453.939*nm), h_Planck*c_light/(453.01*nm), h_Planck*c_light/(451.926*nm), h_Planck*c_light/(450.997*nm), h_Planck*c_light/(450.068*nm), h_Planck*c_light/(448.983*nm), h_Planck*c_light/(448.054*nm), h_Planck*c_light/(446.969*nm), h_Planck*c_light/(446.039*nm), h_Planck*c_light/(444.954*nm), h_Planck*c_light/(444.023*nm), h_Planck*c_light/(442.938*nm), h_Planck*c_light/(442.007*nm), h_Planck*c_light/(441.077*nm), h_Planck*c_light/(439.991*nm), h_Planck*c_light/(439.06*nm), h_Planck*c_light/(437.973*nm), h_Planck*c_light/(437.042*nm), h_Planck*c_light/(435.956*nm), h_Planck*c_light/(435.024*nm), h_Planck*c_light/(433.937*nm), h_Planck*c_light/(433.005*nm), h_Planck*c_light/(432.073*nm), h_Planck*c_light/(430.986*nm), h_Planck*c_light/(430.053*nm), h_Planck*c_light/(428.966*nm), h_Planck*c_light/(428.033*nm), h_Planck*c_light/(426.945*nm), h_Planck*c_light/(426.012*nm), h_Planck*c_light/(424.924*nm), h_Planck*c_light/(423.991*nm), h_Planck*c_light/(423.057*nm), h_Planck*c_light/(421.968*nm), h_Planck*c_light/(421.035*nm), h_Planck*c_light/(419.946*nm), h_Planck*c_light/(419.012*nm), h_Planck*c_light/(417.922*nm), h_Planck*c_light/(416.988*nm), h_Planck*c_light/(416.054*nm), h_Planck*c_light/(414.964*nm), h_Planck*c_light/(414.03*nm), h_Planck*c_light/(412.939*nm), h_Planck*c_light/(412.005*nm), h_Planck*c_light/(411.07*nm), h_Planck*c_light/(409.979*nm), h_Planck*c_light/(409.044*nm), h_Planck*c_light/(407.953*nm), h_Planck*c_light/(407.018*nm), h_Planck*c_light/(405.926*nm), h_Planck*c_light/(404.99*nm), h_Planck*c_light/(404.055*nm), h_Planck*c_light/(402.963*nm), h_Planck*c_light/(402.027*nm), h_Planck*c_light/(400.935*nm), h_Planck*c_light/(399.998*nm), h_Planck*c_light/(399.062*nm), h_Planck*c_light/(397.969*nm), h_Planck*c_light/(397.033*nm), h_Planck*c_light/(395.94*nm), h_Planck*c_light/(395.003*nm), h_Planck*c_light/(394.066*nm), h_Planck*c_light/(392.972*nm), h_Planck*c_light/(392.035*nm), h_Planck*c_light/(390.942*nm), h_Planck*c_light/(390.004*nm), h_Planck*c_light/(389.066*nm), h_Planck*c_light/(387.972*nm), h_Planck*c_light/(387.034*nm), h_Planck*c_light/(385.034*nm), optPhotMaxE_};


    G4double concentration;
    std::vector<G4double> available_concentrations;
    std::map<G4double, std::vector<G4double>> absLength;
    std::map<G4double, std::vector<G4double>> WLS_absLength;

    concentration = 8.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 8.1916*m, 8.404*m, 8.267*m, 8.3063*m, 8.395*m, 8.2733*m, 8.1249*m, 8.2838*m, 8.2509*m, 8.0636*m, 8.0586*m, 7.9863*m, 8.1054*m, 8.2056*m, 8.1715*m, 8.1047*m, 8.1463*m, 8.3699*m, 8.1579*m, 8.0404*m, 8.087*m, 8.3966*m, 7.9784*m, 8.047*m, 8.2303*m, 8.0724*m, 8.1228*m, 8.0335*m, 8.3199*m, 7.8313*m, 7.9673*m, 7.8829*m, 7.9094*m, 7.7826*m, 7.7865*m, 7.9038*m, 7.8015*m, 7.9485*m, 7.851*m, 7.801*m, 7.899*m, 7.9238*m, 7.9347*m, 7.8373*m, 7.7262*m, 7.7338*m, 7.7912*m, 7.7733*m, 7.65*m, 7.7625*m, 7.6698*m, 7.5919*m, 7.7012*m, 7.51*m, 7.5907*m, 7.5763*m, 7.4511*m, 7.4605*m, 7.6127*m, 7.4185*m, 7.4424*m, 7.3628*m, 7.521*m, 7.472*m, 7.4042*m, 7.3167*m, 7.4539*m, 7.449*m, 7.334*m, 7.4076*m, 7.3089*m, 7.4103*m, 7.2427*m, 7.318*m, 7.3341*m, 7.2615*m, 7.2738*m, 7.1922*m, 7.3088*m, 7.1006*m, 7.2539*m, 7.1282*m, 7.1048*m, 7.178*m, 7.0218*m, 7.0991*m, 7.0897*m, 7.0512*m, 7.0086*m, 7.1743*m, 7.0664*m, 6.9257*m, 6.9815*m, 7.0505*m, 6.906*m, 6.84*m, 7.0087*m, 6.7677*m, 6.8387*m, 6.7982*m, 6.8388*m, 6.7578*m, 6.6626*m, 6.7713*m, 6.7381*m, 6.6308*m, 6.5512*m, 6.6281*m, 6.644*m, 6.6424*m, 6.5701*m, 6.5256*m, 6.5549*m, 6.5148*m, 6.4564*m, 6.5716*m, 6.5163*m, 6.4409*m, 6.4623*m, 6.3864*m, 6.4046*m, 6.3118*m, 6.3459*m, 6.3183*m, 6.331*m, 6.2791*m, 6.3088*m, 6.2379*m, 6.1737*m, 6.1324*m, 6.1732*m, 6.1744*m, 6.0149*m, 5.9462*m, 6.0302*m, 5.9644*m, 5.982*m, 5.918*m, 5.8833*m, 5.8408*m, 5.9222*m, 5.7917*m, 5.7442*m, 5.7619*m, 5.6372*m, 5.6313*m, 5.5925*m, 5.5817*m, 5.5416*m, 5.4838*m, 5.3322*m, 5.3111*m, 5.2893*m, 5.2315*m, 5.2032*m, 5.1054*m, 5.0002*m, 4.9245*m, 4.8419*m, 4.7264*m, 4.596*m, 4.4732*m, 4.3267*m, 4.1138*m, 3.8952*m, 3.7231*m, 3.4888*m, 3.2308*m, 3.011*m, 2.7741*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 2.5242*m, 2.291*m, 2.0659*m, 1.8376*m, 1.5903*m, 1.3739*m, 1.1562*m, 0.9358*m, 0.7675*m, 0.6202*m, 0.4925*m, 0.3851*m, 0.2863*m, 0.2178*m, 0.165*m, 0.1265*m, 0.0949*m, 0.0691*m, 0.0528*m, 0.0418*m, 0.0328*m, 0.0261*m, 0.0211*m, 0.0183*m, 0.0167*m, 0.0133*m, 0.0117*m, 0.0105*m, 0.0096891*m, 0.0089716*m, 0.0084991*m, 0.0082002*m, 0.0080207*m, 0.0079387*m, 0.0078843*m, 0.0078883*m, 0.0078813*m, 0.0078778*m, 0.0078473*m, 0.007795*m, 0.007701*m, 0.0075812*m, 0.0074206*m, 0.0072474*m, 0.0070338*m, 0.0067589*m, 0.0064389*m, 0.0061184*m, 0.0058425*m, 0.005552*m, 0.0053203*m, 0.0050843*m, 0.0049818*m, 0.004893*m, 0.0048477*m, 0.0048643*m, 0.0048915*m, 0.0049706*m, 0.0050252*m, 0.0051141*m, 0.0051646*m, 0.0052636*m, 0.0053242*m, 0.005371*m, 0.0054318*m, 0.0054395*m, 0.0054937*m, 0.0054865*m, 0.0055009*m, 0.0054696*m, 0.005451*m, 0.0054451*m, 0.0054307*m, 0.0054562*m, 0.0055016*m, 0.0056283*m, 0.0057409*m, 0.0059369*m, 0.0061275*m, 0.0063551*m, 0.0064948*m, 0.0067181*m, 0.0069183*m, 0.0071367*m, 0.0073426*m, 0.0075311*m, 0.0076923*m, 0.0079742*m, 0.008156*m, 0.0084221*m, 0.0085684*m, 0.008885*m, 0.0091589*m, 0.0094779*m, 0.01*m, 0.0103*m, 0.0108*m, 0.0113*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 16.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 6.9477*m, 7.2589*m, 7.0569*m, 7.1144*m, 7.2455*m, 7.0662*m, 6.8523*m, 7.0815*m, 7.0335*m, 6.7656*m, 6.7586*m, 6.6574*m, 6.8246*m, 6.968*m, 6.9189*m, 6.8236*m, 6.8828*m, 7.2082*m, 6.8994*m, 6.733*m, 6.7986*m, 7.2479*m, 6.6465*m, 6.7423*m, 7.0037*m, 6.7779*m, 6.8493*m, 6.7233*m, 7.1344*m, 6.4447*m, 6.6311*m, 6.515*m, 6.5513*m, 6.3791*m, 6.3844*m, 6.5435*m, 6.4045*m, 6.6051*m, 6.4715*m, 6.4039*m, 6.537*m, 6.571*m, 6.586*m, 6.453*m, 6.3036*m, 6.3138*m, 6.3906*m, 6.3666*m, 6.2029*m, 6.3521*m, 6.229*m, 6.1267*m, 6.2704*m, 6.0208*m, 6.1252*m, 6.1065*m, 5.9454*m, 5.9575*m, 6.1539*m, 5.9041*m, 5.9343*m, 5.8339*m, 6.035*m, 5.9721*m, 5.886*m, 5.7761*m, 5.9491*m, 5.9427*m, 5.7977*m, 5.8903*m, 5.7665*m, 5.8937*m, 5.6844*m, 5.7777*m, 5.7979*m, 5.7076*m, 5.7229*m, 5.6225*m, 5.7663*m, 5.5113*m, 5.6983*m, 5.5447*m, 5.5163*m, 5.6051*m, 5.4169*m, 5.5095*m, 5.4982*m, 5.452*m, 5.4012*m, 5.6006*m, 5.4702*m, 5.3035*m, 5.3691*m, 5.4511*m, 5.2804*m, 5.2036*m, 5.4014*m, 5.1203*m, 5.202*m, 5.1553*m, 5.2022*m, 5.109*m, 5.001*m, 5.1244*m, 5.0865*m, 4.9653*m, 4.8764*m, 4.9622*m, 4.98*m, 4.9783*m, 4.8974*m, 4.8482*m, 4.8805*m, 4.8363*m, 4.7721*m, 4.8991*m, 4.8378*m, 4.7552*m, 4.7786*m, 4.6961*m, 4.7158*m, 4.6159*m, 4.6524*m, 4.6228*m, 4.6364*m, 4.5809*m, 4.6126*m, 4.5373*m, 4.4696*m, 4.4264*m, 4.4691*m, 4.4704*m, 4.305*m, 4.235*m, 4.3208*m, 4.2534*m, 4.2714*m, 4.2064*m, 4.1715*m, 4.1289*m, 4.2107*m, 4.08*m, 4.033*m, 4.0504*m, 3.9283*m, 3.9225*m, 3.885*m, 3.8745*m, 3.836*m, 3.7808*m, 3.6383*m, 3.6186*m, 3.5984*m, 3.5451*m, 3.5192*m, 3.4303*m, 3.336*m, 3.2689*m, 3.1965*m, 3.0966*m, 2.9856*m, 2.8828*m, 2.7623*m, 2.5911*m, 2.4199*m, 2.2885*m, 2.114*m, 1.9275*m, 1.773*m, 1.611*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 1.4449*m, 1.2941*m, 1.1523*m, 1.012*m, 0.864*m, 0.7377*m, 0.6137*m, 0.4909*m, 0.3991*m, 0.32*m, 0.2525*m, 0.1963*m, 0.1452*m, 0.1101*m, 0.0832*m, 0.0637*m, 0.0477*m, 0.0347*m, 0.0264*m, 0.021*m, 0.0164*m, 0.0131*m, 0.0106*m, 0.0091694*m, 0.008354*m, 0.0066686*m, 0.0058575*m, 0.0052449*m, 0.0048445*m, 0.0044858*m, 0.0042495*m, 0.0041001*m, 0.0040103*m, 0.0039693*m, 0.0039422*m, 0.0039442*m, 0.0039406*m, 0.0039389*m, 0.0039237*m, 0.0038975*m, 0.0038505*m, 0.0037906*m, 0.0037103*m, 0.0036237*m, 0.0035169*m, 0.0033795*m, 0.0032195*m, 0.0030592*m, 0.0029212*m, 0.002776*m, 0.0026601*m, 0.0025422*m, 0.0024909*m, 0.0024465*m, 0.0024239*m, 0.0024321*m, 0.0024458*m, 0.0024853*m, 0.0025126*m, 0.002557*m, 0.0025823*m, 0.0026318*m, 0.0026621*m, 0.0026855*m, 0.0027159*m, 0.0027198*m, 0.0027468*m, 0.0027433*m, 0.0027504*m, 0.0027348*m, 0.0027255*m, 0.0027226*m, 0.0027154*m, 0.0027281*m, 0.0027508*m, 0.0028141*m, 0.0028705*m, 0.0029685*m, 0.0030637*m, 0.0031776*m, 0.0032474*m, 0.0033591*m, 0.0034591*m, 0.0035684*m, 0.0036713*m, 0.0037656*m, 0.0038462*m, 0.0039871*m, 0.004078*m, 0.0042111*m, 0.0042842*m, 0.0044425*m, 0.0045795*m, 0.0047389*m, 0.0049784*m, 0.0051419*m, 0.0053926*m, 0.0056263*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 24.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 6.0318*m, 6.3885*m, 6.1558*m, 6.2217*m, 6.3729*m, 6.1664*m, 5.9244*m, 6.1839*m, 6.1292*m, 5.8276*m, 5.8197*m, 5.7077*m, 5.8934*m, 6.0547*m, 5.9993*m, 5.8923*m, 5.9586*m, 6.3296*m, 5.9773*m, 5.7913*m, 5.8644*m, 6.3757*m, 5.6957*m, 5.8016*m, 6.0952*m, 5.8413*m, 5.9211*m, 5.7806*m, 6.2446*m, 5.4753*m, 5.6787*m, 5.5516*m, 5.5912*m, 5.4045*m, 5.4102*m, 5.5827*m, 5.4319*m, 5.6502*m, 5.5044*m, 5.4312*m, 5.5756*m, 5.6127*m, 5.6292*m, 5.4842*m, 5.3234*m, 5.3343*m, 5.4169*m, 5.391*m, 5.2162*m, 5.3754*m, 5.2438*m, 5.1356*m, 5.288*m, 5.0245*m, 5.134*m, 5.1143*m, 4.946*m, 4.9585*m, 5.1643*m, 4.9032*m, 4.9345*m, 4.8307*m, 5.0394*m, 4.9737*m, 4.8844*m, 4.7715*m, 4.9498*m, 4.9432*m, 4.7936*m, 4.8889*m, 4.7616*m, 4.8924*m, 4.6779*m, 4.7731*m, 4.7938*m, 4.7015*m, 4.7171*m, 4.6152*m, 4.7615*m, 4.5033*m, 4.692*m, 4.5368*m, 4.5084*m, 4.5977*m, 4.4092*m, 4.5015*m, 4.4902*m, 4.4441*m, 4.3936*m, 4.5931*m, 4.4622*m, 4.2969*m, 4.3617*m, 4.4432*m, 4.2742*m, 4.199*m, 4.3937*m, 4.1179*m, 4.1975*m, 4.152*m, 4.1976*m, 4.107*m, 4.0027*m, 4.1219*m, 4.0852*m, 3.9684*m, 3.8836*m, 3.9655*m, 3.9826*m, 3.9809*m, 3.9036*m, 3.8568*m, 3.8875*m, 3.8455*m, 3.7848*m, 3.9052*m, 3.847*m, 3.7689*m, 3.7909*m, 3.7133*m, 3.7318*m, 3.6383*m, 3.6724*m, 3.6447*m, 3.6574*m, 3.6058*m, 3.6352*m, 3.5653*m, 3.5028*m, 3.4631*m, 3.5023*m, 3.5034*m, 3.3521*m, 3.2886*m, 3.3664*m, 3.3053*m, 3.3216*m, 3.2628*m, 3.2313*m, 3.193*m, 3.2666*m, 3.1492*m, 3.1073*m, 3.1228*m, 3.0145*m, 3.0094*m, 2.9763*m, 2.9671*m, 2.9332*m, 2.8849*m, 2.7611*m, 2.7441*m, 2.7267*m, 2.6809*m, 2.6587*m, 2.5829*m, 2.5029*m, 2.4464*m, 2.3858*m, 2.3026*m, 2.2109*m, 2.1267*m, 2.0287*m, 1.8911*m, 1.7552*m, 1.6519*m, 1.5164*m, 1.3734*m, 1.2564*m, 1.1351*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 1.0121*m, 0.9017*m, 0.7989*m, 0.6983*m, 0.5931*m, 0.5042*m, 0.4177*m, 0.3327*m, 0.2696*m, 0.2157*m, 0.1697*m, 0.1318*m, 0.0973*m, 0.0737*m, 0.0556*m, 0.0425*m, 0.0318*m, 0.0231*m, 0.0176*m, 0.014*m, 0.011*m, 0.0087088*m, 0.0070573*m, 0.0061148*m, 0.0055709*m, 0.0044457*m, 0.003905*m, 0.0034966*m, 0.0032297*m, 0.0029905*m, 0.002833*m, 0.0027334*m, 0.0026736*m, 0.0026462*m, 0.0026281*m, 0.0026294*m, 0.0026271*m, 0.0026259*m, 0.0026158*m, 0.0025983*m, 0.002567*m, 0.0025271*m, 0.0024735*m, 0.0024158*m, 0.0023446*m, 0.002253*m, 0.0021463*m, 0.0020395*m, 0.0019475*m, 0.0018507*m, 0.0017734*m, 0.0016948*m, 0.0016606*m, 0.001631*m, 0.0016159*m, 0.0016214*m, 0.0016305*m, 0.0016569*m, 0.0016751*m, 0.0017047*m, 0.0017215*m, 0.0017545*m, 0.0017747*m, 0.0017903*m, 0.0018106*m, 0.0018132*m, 0.0018312*m, 0.0018288*m, 0.0018336*m, 0.0018232*m, 0.001817*m, 0.001815*m, 0.0018102*m, 0.0018187*m, 0.0018339*m, 0.0018761*m, 0.0019136*m, 0.001979*m, 0.0020425*m, 0.0021184*m, 0.0021649*m, 0.0022394*m, 0.0023061*m, 0.0023789*m, 0.0024475*m, 0.0025104*m, 0.0025641*m, 0.0026581*m, 0.0027187*m, 0.0028074*m, 0.0028561*m, 0.0029617*m, 0.003053*m, 0.0031593*m, 0.003319*m, 0.0034279*m, 0.0035951*m, 0.0037509*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 32.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 5.3293*m, 5.7045*m, 5.4588*m, 5.528*m, 5.6879*m, 5.4699*m, 5.2178*m, 5.4883*m, 5.4309*m, 5.1179*m, 5.1099*m, 4.9951*m, 5.1858*m, 5.3531*m, 5.2954*m, 5.1846*m, 5.2532*m, 5.642*m, 5.2726*m, 5.0807*m, 5.1558*m, 5.6909*m, 4.9829*m, 5.0913*m, 5.3954*m, 5.1321*m, 5.2144*m, 5.0697*m, 5.5521*m, 4.7594*m, 4.9656*m, 4.8364*m, 4.8765*m, 4.6882*m, 4.6939*m, 4.868*m, 4.7157*m, 4.9365*m, 4.7887*m, 4.715*m, 4.8607*m, 4.8984*m, 4.9151*m, 4.7684*m, 4.6071*m, 4.6179*m, 4.7006*m, 4.6747*m, 4.5003*m, 4.659*m, 4.5278*m, 4.4205*m, 4.5717*m, 4.3111*m, 4.419*m, 4.3994*m, 4.2342*m, 4.2465*m, 4.4489*m, 4.1924*m, 4.223*m, 4.1219*m, 4.3257*m, 4.2613*m, 4.1742*m, 4.0645*m, 4.2379*m, 4.2315*m, 4.0859*m, 4.1785*m, 4.055*m, 4.1819*m, 3.9742*m, 4.0661*m, 4.0861*m, 3.997*m, 4.012*m, 3.914*m, 4.0548*m, 3.8071*m, 3.9878*m, 3.839*m, 3.8119*m, 3.8972*m, 3.7176*m, 3.8053*m, 3.7946*m, 3.7507*m, 3.7028*m, 3.8928*m, 3.7679*m, 3.6115*m, 3.6727*m, 3.7499*m, 3.5901*m, 3.5195*m, 3.703*m, 3.4437*m, 3.5181*m, 3.4755*m, 3.5182*m, 3.4336*m, 3.3366*m, 3.4475*m, 3.4132*m, 3.3049*m, 3.2267*m, 3.3022*m, 3.318*m, 3.3165*m, 3.2451*m, 3.202*m, 3.2303*m, 3.1916*m, 3.136*m, 3.2465*m, 3.193*m, 3.1214*m, 3.1416*m, 3.0706*m, 3.0875*m, 3.0024*m, 3.0334*m, 3.0082*m, 3.0198*m, 2.9729*m, 2.9996*m, 2.9362*m, 2.8798*m, 2.8441*m, 2.8794*m, 2.8804*m, 2.7446*m, 2.6879*m, 2.7574*m, 2.7028*m, 2.7173*m, 2.6649*m, 2.6369*m, 2.603*m, 2.6683*m, 2.5643*m, 2.5272*m, 2.5409*m, 2.4455*m, 2.4411*m, 2.4121*m, 2.404*m, 2.3744*m, 2.3323*m, 2.2247*m, 2.21*m, 2.195*m, 2.1555*m, 2.1363*m, 2.0712*m, 2.0028*m, 1.9546*m, 1.9031*m, 1.8327*m, 1.7554*m, 1.6848*m, 1.603*m, 1.4888*m, 1.3769*m, 1.2925*m, 1.1822*m, 1.0668*m, 0.9729*m, 0.8762*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.7788*m, 0.6919*m, 0.6114*m, 0.533*m, 0.4516*m, 0.383*m, 0.3166*m, 0.2516*m, 0.2036*m, 0.1626*m, 0.1278*m, 0.0991*m, 0.0732*m, 0.0553*m, 0.0418*m, 0.0319*m, 0.0239*m, 0.0174*m, 0.0132*m, 0.0105*m, 0.0082148*m, 0.006533*m, 0.0052939*m, 0.0045868*m, 0.0041787*m, 0.0033343*m, 0.0029288*m, 0.0026224*m, 0.0024223*m, 0.0022429*m, 0.0021248*m, 0.0020501*m, 0.0020052*m, 0.0019847*m, 0.0019711*m, 0.0019721*m, 0.0019703*m, 0.0019694*m, 0.0019618*m, 0.0019487*m, 0.0019252*m, 0.0018953*m, 0.0018551*m, 0.0018119*m, 0.0017584*m, 0.0016897*m, 0.0016097*m, 0.0015296*m, 0.0014606*m, 0.001388*m, 0.0013301*m, 0.0012711*m, 0.0012455*m, 0.0012232*m, 0.0012119*m, 0.0012161*m, 0.0012229*m, 0.0012427*m, 0.0012563*m, 0.0012785*m, 0.0012912*m, 0.0013159*m, 0.001331*m, 0.0013427*m, 0.0013579*m, 0.0013599*m, 0.0013734*m, 0.0013716*m, 0.0013752*m, 0.0013674*m, 0.0013627*m, 0.0013613*m, 0.0013577*m, 0.001364*m, 0.0013754*m, 0.0014071*m, 0.0014352*m, 0.0014842*m, 0.0015319*m, 0.0015888*m, 0.0016237*m, 0.0016795*m, 0.0017296*m, 0.0017842*m, 0.0018356*m, 0.0018828*m, 0.0019231*m, 0.0019936*m, 0.002039*m, 0.0021055*m, 0.0021421*m, 0.0022213*m, 0.0022897*m, 0.0023695*m, 0.0024892*m, 0.0025709*m, 0.0026963*m, 0.0028132*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 40.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 4.7733*m, 5.1527*m, 4.9036*m, 4.9735*m, 5.1358*m, 4.9148*m, 4.6618*m, 4.9334*m, 4.8754*m, 4.5624*m, 4.5544*m, 4.4407*m, 4.6299*m, 4.7973*m, 4.7394*m, 4.6287*m, 4.6972*m, 5.0891*m, 4.7165*m, 4.5254*m, 4.6001*m, 5.1389*m, 4.4286*m, 4.5359*m, 4.8397*m, 4.5764*m, 4.6584*m, 4.5145*m, 4.9979*m, 4.2091*m, 4.4115*m, 4.2845*m, 4.3239*m, 4.1396*m, 4.1451*m, 4.3154*m, 4.1664*m, 4.3829*m, 4.2377*m, 4.1657*m, 4.3084*m, 4.3454*m, 4.3618*m, 4.2179*m, 4.0606*m, 4.0712*m, 4.1517*m, 4.1264*m, 3.9572*m, 4.1111*m, 3.9838*m, 3.8803*m, 4.0263*m, 3.7751*m, 3.8787*m, 3.8599*m, 3.7015*m, 3.7132*m, 3.9076*m, 3.6616*m, 3.6908*m, 3.5945*m, 3.7891*m, 3.7274*m, 3.6442*m, 3.54*m, 3.705*m, 3.6989*m, 3.5603*m, 3.6484*m, 3.531*m, 3.6516*m, 3.4546*m, 3.5415*m, 3.5605*m, 3.4761*m, 3.4903*m, 3.3977*m, 3.5308*m, 3.2973*m, 3.4674*m, 3.3272*m, 3.3018*m, 3.3819*m, 3.2135*m, 3.2956*m, 3.2855*m, 3.2445*m, 3.1997*m, 3.3778*m, 3.2606*m, 3.1147*m, 3.1716*m, 3.2437*m, 3.0948*m, 3.0293*m, 3.1999*m, 2.9593*m, 3.028*m, 2.9886*m, 3.0281*m, 2.9499*m, 2.8606*m, 2.9627*m, 2.9311*m, 2.8315*m, 2.7598*m, 2.829*m, 2.8435*m, 2.8421*m, 2.7767*m, 2.7373*m, 2.7631*m, 2.7278*m, 2.6771*m, 2.778*m, 2.7291*m, 2.6638*m, 2.6822*m, 2.6176*m, 2.633*m, 2.5557*m, 2.5838*m, 2.561*m, 2.5715*m, 2.529*m, 2.5532*m, 2.4959*m, 2.445*m, 2.4128*m, 2.4446*m, 2.4455*m, 2.3235*m, 2.2728*m, 2.335*m, 2.2861*m, 2.2991*m, 2.2522*m, 2.2273*m, 2.197*m, 2.2553*m, 2.1626*m, 2.1297*m, 2.1418*m, 2.0573*m, 2.0533*m, 2.0277*m, 2.0206*m, 1.9945*m, 1.9573*m, 1.8629*m, 1.85*m, 1.8368*m, 1.8023*m, 1.7855*m, 1.7287*m, 1.6692*m, 1.6275*m, 1.5829*m, 1.5221*m, 1.4556*m, 1.3949*m, 1.325*m, 1.2277*m, 1.1328*m, 1.0615*m, 0.9687*m, 0.8721*m, 0.7938*m, 0.7135*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.6329*m, 0.5613*m, 0.4952*m, 0.431*m, 0.3645*m, 0.3088*m, 0.2549*m, 0.2023*m, 0.1636*m, 0.1305*m, 0.1025*m, 0.0795*m, 0.0586*m, 0.0443*m, 0.0334*m, 0.0256*m, 0.0191*m, 0.0139*m, 0.0106*m, 0.0083913*m, 0.006573*m, 0.0052271*m, 0.0042356*m, 0.0036698*m, 0.0033433*m, 0.0026674*m, 0.002343*m, 0.002098*m, 0.0019378*m, 0.0017943*m, 0.0016998*m, 0.00164*m, 0.0016041*m, 0.0015877*m, 0.0015769*m, 0.0015777*m, 0.0015763*m, 0.0015756*m, 0.0015695*m, 0.001559*m, 0.0015402*m, 0.0015162*m, 0.0014841*m, 0.0014495*m, 0.0014068*m, 0.0013518*m, 0.0012878*m, 0.0012237*m, 0.0011685*m, 0.0011104*m, 0.0010641*m, 0.0010169*m, 0.00099637*m, 0.0009786*m, 0.00096954*m, 0.00097286*m, 0.00097831*m, 0.00099412*m, 0.001005*m, 0.0010228*m, 0.0010329*m, 0.0010527*m, 0.0010648*m, 0.0010742*m, 0.0010864*m, 0.0010879*m, 0.0010987*m, 0.0010973*m, 0.0011002*m, 0.0010939*m, 0.0010902*m, 0.001089*m, 0.0010861*m, 0.0010912*m, 0.0011003*m, 0.0011257*m, 0.0011482*m, 0.0011874*m, 0.0012255*m, 0.001271*m, 0.001299*m, 0.0013436*m, 0.0013837*m, 0.0014273*m, 0.0014685*m, 0.0015062*m, 0.0015385*m, 0.0015948*m, 0.0016312*m, 0.0016844*m, 0.0017137*m, 0.001777*m, 0.0018318*m, 0.0018956*m, 0.0019914*m, 0.0020567*m, 0.0021571*m, 0.0022505*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 48.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 4.3224*m, 4.6983*m, 4.4509*m, 4.52*m, 4.6815*m, 4.462*m, 4.2129*m, 4.4803*m, 4.4231*m, 4.1156*m, 4.1078*m, 3.997*m, 4.1816*m, 4.346*m, 4.289*m, 4.1805*m, 4.2476*m, 4.6349*m, 4.2666*m, 4.0795*m, 4.1524*m, 4.6845*m, 3.9853*m, 4.0898*m, 4.3878*m, 4.1293*m, 4.2095*m, 4.0689*m, 4.5443*m, 3.7728*m, 3.9687*m, 3.8456*m, 3.8837*m, 3.7059*m, 3.7112*m, 3.8756*m, 3.7317*m, 3.9409*m, 3.8004*m, 3.731*m, 3.8687*m, 3.9046*m, 3.9205*m, 3.7813*m, 3.6301*m, 3.6402*m, 3.7175*m, 3.6932*m, 3.531*m, 3.6786*m, 3.5565*m, 3.4576*m, 3.5972*m, 3.3577*m, 3.4562*m, 3.4383*m, 3.2879*m, 3.299*m, 3.4837*m, 3.2501*m, 3.2777*m, 3.1868*m, 3.3709*m, 3.3124*m, 3.2337*m, 3.1354*m, 3.2912*m, 3.2854*m, 3.1545*m, 3.2376*m, 3.1269*m, 3.2407*m, 3.0551*m, 3.1368*m, 3.1547*m, 3.0753*m, 3.0886*m, 3.0018*m, 3.1268*m, 2.9079*m, 3.0672*m, 2.9358*m, 2.9121*m, 2.987*m, 2.8298*m, 2.9064*m, 2.8969*m, 2.8587*m, 2.817*m, 2.9832*m, 2.8737*m, 2.738*m, 2.7909*m, 2.8579*m, 2.7196*m, 2.659*m, 2.8171*m, 2.5943*m, 2.6578*m, 2.6214*m, 2.6579*m, 2.5856*m, 2.5035*m, 2.5975*m, 2.5683*m, 2.4767*m, 2.411*m, 2.4744*m, 2.4878*m, 2.4865*m, 2.4264*m, 2.3904*m, 2.414*m, 2.3817*m, 2.3353*m, 2.4277*m, 2.3828*m, 2.3232*m, 2.34*m, 2.2811*m, 2.2951*m, 2.2248*m, 2.2503*m, 2.2295*m, 2.2391*m, 2.2005*m, 2.2225*m, 2.1704*m, 2.1242*m, 2.0951*m, 2.1239*m, 2.1247*m, 2.0144*m, 1.9687*m, 2.0248*m, 1.9807*m, 1.9924*m, 1.9502*m, 1.9278*m, 1.9006*m, 1.953*m, 1.8697*m, 1.8402*m, 1.8511*m, 1.7754*m, 1.7719*m, 1.749*m, 1.7426*m, 1.7193*m, 1.6862*m, 1.6022*m, 1.5908*m, 1.5791*m, 1.5485*m, 1.5337*m, 1.4834*m, 1.4309*m, 1.3941*m, 1.3549*m, 1.3015*m, 1.2432*m, 1.1902*m, 1.1291*m, 1.0445*m, 0.9622*m, 0.9005*m, 0.8206*m, 0.7375*m, 0.6704*m, 0.6018*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.5331*m, 0.4722*m, 0.4161*m, 0.3618*m, 0.3057*m, 0.2587*m, 0.2133*m, 0.1692*m, 0.1367*m, 0.109*m, 0.0856*m, 0.0663*m, 0.0489*m, 0.037*m, 0.0279*m, 0.0213*m, 0.0159*m, 0.0116*m, 0.0088311*m, 0.0069937*m, 0.0054781*m, 0.0043563*m, 0.0035299*m, 0.0030583*m, 0.0027862*m, 0.0022229*m, 0.0019525*m, 0.0017483*m, 0.0016148*m, 0.0014953*m, 0.0014165*m, 0.0013667*m, 0.0013368*m, 0.0013231*m, 0.0013141*m, 0.0013147*m, 0.0013135*m, 0.001313*m, 0.0013079*m, 0.0012992*m, 0.0012835*m, 0.0012635*m, 0.0012368*m, 0.0012079*m, 0.0011723*m, 0.0011265*m, 0.0010732*m, 0.0010197*m, 0.00097374*m, 0.00092533*m, 0.00088671*m, 0.00084739*m, 0.00083031*m, 0.0008155*m, 0.00080795*m, 0.00081071*m, 0.00081526*m, 0.00082844*m, 0.00083753*m, 0.00085235*m, 0.00086077*m, 0.00087727*m, 0.00088736*m, 0.00089516*m, 0.0009053*m, 0.00090659*m, 0.00091562*m, 0.00091442*m, 0.00091681*m, 0.00091159*m, 0.00090849*m, 0.00090752*m, 0.00090512*m, 0.00090937*m, 0.00091693*m, 0.00093804*m, 0.00095682*m, 0.00098948*m, 0.0010212*m, 0.0010592*m, 0.0010825*m, 0.0011197*m, 0.001153*m, 0.0011895*m, 0.0012238*m, 0.0012552*m, 0.0012821*m, 0.001329*m, 0.0013593*m, 0.0014037*m, 0.0014281*m, 0.0014808*m, 0.0015265*m, 0.0015796*m, 0.0016595*m, 0.001714*m, 0.0017975*m, 0.0018754*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 56.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 3.9493*m, 4.3176*m, 4.0747*m, 4.1424*m, 4.301*m, 4.0856*m, 3.8428*m, 4.1035*m, 4.0475*m, 3.7486*m, 3.741*m, 3.634*m, 3.8125*m, 3.9723*m, 3.9168*m, 3.8114*m, 3.8765*m, 4.2552*m, 3.895*m, 3.7136*m, 3.7842*m, 4.304*m, 3.6227*m, 3.7236*m, 4.0131*m, 3.7618*m, 3.8396*m, 3.7034*m, 4.1661*m, 3.4185*m, 3.6067*m, 3.4883*m, 3.5249*m, 3.3544*m, 3.3595*m, 3.5171*m, 3.3791*m, 3.5799*m, 3.445*m, 3.3785*m, 3.5105*m, 3.5449*m, 3.5603*m, 3.4266*m, 3.2821*m, 3.2917*m, 3.3656*m, 3.3423*m, 3.1878*m, 3.3284*m, 3.2119*m, 3.1181*m, 3.2507*m, 3.0233*m, 3.1167*m, 3.0997*m, 2.9574*m, 2.9679*m, 3.1428*m, 2.9218*m, 2.9478*m, 2.8621*m, 3.0359*m, 2.9805*m, 2.9063*m, 2.8138*m, 2.9606*m, 2.9551*m, 2.8318*m, 2.91*m, 2.8058*m, 2.9129*m, 2.7385*m, 2.8152*m, 2.832*m, 2.7574*m, 2.7699*m, 2.6885*m, 2.8057*m, 2.6007*m, 2.7498*m, 2.6268*m, 2.6047*m, 2.6747*m, 2.528*m, 2.5993*m, 2.5905*m, 2.5549*m, 2.516*m, 2.6711*m, 2.5688*m, 2.4426*m, 2.4917*m, 2.5542*m, 2.4255*m, 2.3693*m, 2.5162*m, 2.3095*m, 2.3682*m, 2.3345*m, 2.3683*m, 2.3014*m, 2.2256*m, 2.3124*m, 2.2855*m, 2.201*m, 2.1405*m, 2.1988*m, 2.2111*m, 2.2099*m, 2.1546*m, 2.1215*m, 2.1432*m, 2.1135*m, 2.0709*m, 2.1558*m, 2.1146*m, 2.0598*m, 2.0752*m, 2.0212*m, 2.034*m, 1.9697*m, 1.993*m, 1.9741*m, 1.9828*m, 1.9475*m, 1.9676*m, 1.92*m, 1.8779*m, 1.8513*m, 1.8776*m, 1.8783*m, 1.7779*m, 1.7364*m, 1.7873*m, 1.7473*m, 1.7579*m, 1.7197*m, 1.6993*m, 1.6747*m, 1.7221*m, 1.6467*m, 1.62*m, 1.6298*m, 1.5615*m, 1.5583*m, 1.5376*m, 1.5319*m, 1.5109*m, 1.4811*m, 1.4056*m, 1.3953*m, 1.3849*m, 1.3574*m, 1.3441*m, 1.2991*m, 1.2522*m, 1.2193*m, 1.1843*m, 1.1367*m, 1.0849*m, 1.0378*m, 0.9837*m, 0.9089*m, 0.8363*m, 0.782*m, 0.7117*m, 0.6388*m, 0.5802*m, 0.5203*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.4604*m, 0.4075*m, 0.3588*m, 0.3117*m, 0.2631*m, 0.2225*m, 0.1834*m, 0.1454*m, 0.1174*m, 0.0936*m, 0.0735*m, 0.0569*m, 0.0419*m, 0.0317*m, 0.0239*m, 0.0183*m, 0.0137*m, 0.0099274*m, 0.0075705*m, 0.0059952*m, 0.0046959*m, 0.0037342*m, 0.0030258*m, 0.0026216*m, 0.0023883*m, 0.0019053*m, 0.0016736*m, 0.0014985*m, 0.0013842*m, 0.0012817*m, 0.0012142*m, 0.0011715*m, 0.0011458*m, 0.0011341*m, 0.0011263*m, 0.0011269*m, 0.0011259*m, 0.0011254*m, 0.001121*m, 0.0011136*m, 0.0011001*m, 0.001083*m, 0.0010601*m, 0.0010353*m, 0.0010048*m, 0.00096556*m, 0.00091985*m, 0.00087406*m, 0.00083464*m, 0.00079314*m, 0.00076004*m, 0.00072633*m, 0.00071169*m, 0.000699*m, 0.00069253*m, 0.0006949*m, 0.00069879*m, 0.00071009*m, 0.00071788*m, 0.00073059*m, 0.00073781*m, 0.00075194*m, 0.0007606*m, 0.00076728*m, 0.00077597*m, 0.00077707*m, 0.00078481*m, 0.00078379*m, 0.00078584*m, 0.00078137*m, 0.00077871*m, 0.00077787*m, 0.00077582*m, 0.00077946*m, 0.00078594*m, 0.00080404*m, 0.00082013*m, 0.00084813*m, 0.00087535*m, 0.00090788*m, 0.00092782*m, 0.00095973*m, 0.00098832*m, 0.0010195*m, 0.0010489*m, 0.0010759*m, 0.0010989*m, 0.0011392*m, 0.0011651*m, 0.0012032*m, 0.0012241*m, 0.0012693*m, 0.0013084*m, 0.001354*m, 0.0014224*m, 0.0014691*m, 0.0015408*m, 0.0016075*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 64.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 3.6355*m, 3.9939*m, 3.7571*m, 3.823*m, 3.9777*m, 3.7677*m, 3.5325*m, 3.7852*m, 3.7308*m, 3.4416*m, 3.4343*m, 3.3314*m, 3.5033*m, 3.6578*m, 3.6041*m, 3.5022*m, 3.5651*m, 3.9329*m, 3.583*m, 3.408*m, 3.476*m, 3.9806*m, 3.3206*m, 3.4176*m, 3.6973*m, 3.4544*m, 3.5294*m, 3.3981*m, 3.8461*m, 3.125*m, 3.3052*m, 3.1918*m, 3.2268*m, 3.0639*m, 3.0688*m, 3.2193*m, 3.0874*m, 3.2795*m, 3.1503*m, 3.0868*m, 3.213*m, 3.246*m, 3.2606*m, 3.1328*m, 2.9949*m, 3.0042*m, 3.0745*m, 3.0523*m, 2.9053*m, 3.039*m, 2.9283*m, 2.8392*m, 2.9651*m, 2.7495*m, 2.8379*m, 2.8218*m, 2.6873*m, 2.6972*m, 2.8626*m, 2.6537*m, 2.6782*m, 2.5975*m, 2.7614*m, 2.7091*m, 2.6391*m, 2.552*m, 2.6903*m, 2.6851*m, 2.569*m, 2.6426*m, 2.5445*m, 2.6453*m, 2.4813*m, 2.5533*m, 2.5691*m, 2.499*m, 2.5108*m, 2.4345*m, 2.5444*m, 2.3523*m, 2.4919*m, 2.3767*m, 2.356*m, 2.4215*m, 2.2843*m, 2.351*m, 2.3427*m, 2.3094*m, 2.2732*m, 2.4181*m, 2.3225*m, 2.2048*m, 2.2505*m, 2.3088*m, 2.1888*m, 2.1366*m, 2.2733*m, 2.081*m, 2.1355*m, 2.1042*m, 2.1356*m, 2.0735*m, 2.0033*m, 2.0837*m, 2.0587*m, 1.9805*m, 1.9245*m, 1.9785*m, 1.9899*m, 1.9887*m, 1.9376*m, 1.907*m, 1.9271*m, 1.8996*m, 1.8603*m, 1.9387*m, 1.9006*m, 1.8501*m, 1.8643*m, 1.8145*m, 1.8263*m, 1.7671*m, 1.7886*m, 1.7711*m, 1.7791*m, 1.7467*m, 1.7651*m, 1.7214*m, 1.6827*m, 1.6584*m, 1.6824*m, 1.6832*m, 1.5911*m, 1.5532*m, 1.5997*m, 1.5631*m, 1.5728*m, 1.5378*m, 1.5192*m, 1.4967*m, 1.5401*m, 1.4712*m, 1.4469*m, 1.4558*m, 1.3936*m, 1.3906*m, 1.3718*m, 1.3666*m, 1.3475*m, 1.3205*m, 1.2519*m, 1.2426*m, 1.2331*m, 1.2083*m, 1.1962*m, 1.1555*m, 1.1131*m, 1.0834*m, 1.0519*m, 1.009*m, 0.9624*m, 0.9201*m, 0.8715*m, 0.8044*m, 0.7395*m, 0.691*m, 0.6283*m, 0.5635*m, 0.5114*m, 0.4582*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.4052*m, 0.3584*m, 0.3154*m, 0.2738*m, 0.231*m, 0.1953*m, 0.1608*m, 0.1274*m, 0.1029*m, 0.082*m, 0.0643*m, 0.0498*m, 0.0367*m, 0.0277*m, 0.0209*m, 0.016*m, 0.012*m, 0.0086875*m, 0.0066248*m, 0.0052462*m, 0.0041091*m, 0.0032676*m, 0.0026476*m, 0.0022939*m, 0.0020898*m, 0.0016672*m, 0.0014644*m, 0.0013112*m, 0.0012111*m, 0.0011214*m, 0.0010624*m, 0.001025*m, 0.0010026*m, 0.00099233*m, 0.00098554*m, 0.00098604*m, 0.00098516*m, 0.00098472*m, 0.00098092*m, 0.00097437*m, 0.00096262*m, 0.00094764*m, 0.00092757*m, 0.00090593*m, 0.00087922*m, 0.00084487*m, 0.00080487*m, 0.0007648*m, 0.00073031*m, 0.000694*m, 0.00066504*m, 0.00063554*m, 0.00062273*m, 0.00061162*m, 0.00060596*m, 0.00060804*m, 0.00061144*m, 0.00062133*m, 0.00062814*m, 0.00063926*m, 0.00064558*m, 0.00065795*m, 0.00066552*m, 0.00067137*m, 0.00067897*m, 0.00067994*m, 0.00068671*m, 0.00068581*m, 0.00068761*m, 0.0006837*m, 0.00068137*m, 0.00068064*m, 0.00067884*m, 0.00068202*m, 0.0006877*m, 0.00070353*m, 0.00071762*m, 0.00074211*m, 0.00076593*m, 0.00079439*m, 0.00081185*m, 0.00083977*m, 0.00086478*m, 0.00089209*m, 0.00091782*m, 0.00094139*m, 0.00096154*m, 0.00099678*m, 0.0010195*m, 0.0010528*m, 0.001071*m, 0.0011106*m, 0.0011449*m, 0.0011847*m, 0.0012446*m, 0.0012855*m, 0.0013482*m, 0.0014066*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 72.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 3.3679*m, 3.7154*m, 3.4855*m, 3.5493*m, 3.6996*m, 3.4957*m, 3.2686*m, 3.5126*m, 3.46*m, 3.1811*m, 3.1741*m, 3.0753*m, 3.2404*m, 3.3894*m, 3.3376*m, 3.2394*m, 3.3*m, 3.6561*m, 3.3172*m, 3.1488*m, 3.2142*m, 3.7024*m, 3.0649*m, 3.158*m, 3.4276*m, 3.1934*m, 3.2656*m, 3.1393*m, 3.5717*m, 2.878*m, 3.0502*m, 2.9417*m, 2.9752*m, 2.8197*m, 2.8243*m, 2.968*m, 2.8421*m, 3.0256*m, 2.9021*m, 2.8415*m, 2.9619*m, 2.9935*m, 3.0076*m, 2.8853*m, 2.754*m, 2.7628*m, 2.8298*m, 2.8087*m, 2.6688*m, 2.796*m, 2.6907*m, 2.6061*m, 2.7257*m, 2.5212*m, 2.6049*m, 2.5897*m, 2.4624*m, 2.4717*m, 2.6283*m, 2.4307*m, 2.4538*m, 2.3776*m, 2.5325*m, 2.483*m, 2.4169*m, 2.3348*m, 2.4652*m, 2.4603*m, 2.3508*m, 2.4202*m, 2.3277*m, 2.4228*m, 2.2682*m, 2.336*m, 2.3509*m, 2.2849*m, 2.296*m, 2.2243*m, 2.3277*m, 2.1472*m, 2.2782*m, 2.1701*m, 2.1506*m, 2.2121*m, 2.0835*m, 2.1459*m, 2.1382*m, 2.107*m, 2.0731*m, 2.2089*m, 2.1192*m, 2.0091*m, 2.0519*m, 2.1064*m, 1.9942*m, 1.9455*m, 2.0732*m, 1.8936*m, 1.9445*m, 1.9153*m, 1.9446*m, 1.8867*m, 1.8213*m, 1.8962*m, 1.8729*m, 1.8001*m, 1.7481*m, 1.7983*m, 1.8088*m, 1.8078*m, 1.7603*m, 1.7319*m, 1.7505*m, 1.7251*m, 1.6886*m, 1.7613*m, 1.7259*m, 1.6791*m, 1.6923*m, 1.6462*m, 1.6571*m, 1.6023*m, 1.6221*m, 1.606*m, 1.6134*m, 1.5834*m, 1.6005*m, 1.56*m, 1.5243*m, 1.5018*m, 1.524*m, 1.5247*m, 1.4398*m, 1.4049*m, 1.4478*m, 1.414*m, 1.423*m, 1.3908*m, 1.3737*m, 1.353*m, 1.3929*m, 1.3295*m, 1.3071*m, 1.3154*m, 1.2582*m, 1.2556*m, 1.2383*m, 1.2336*m, 1.2161*m, 1.1913*m, 1.1286*m, 1.1201*m, 1.1114*m, 1.0887*m, 1.0777*m, 1.0405*m, 1.0018*m, 0.9748*m, 0.9461*m, 0.9071*m, 0.8647*m, 0.8263*m, 0.7823*m, 0.7215*m, 0.6628*m, 0.619*m, 0.5625*m, 0.5041*m, 0.4572*m, 0.4094*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.3618*m, 0.3198*m, 0.2813*m, 0.2441*m, 0.2059*m, 0.1739*m, 0.1432*m, 0.1134*m, 0.0915*m, 0.0729*m, 0.0572*m, 0.0443*m, 0.0326*m, 0.0247*m, 0.0186*m, 0.0142*m, 0.0106*m, 0.007723*m, 0.0058892*m, 0.0046636*m, 0.0036527*m, 0.0029046*m, 0.0023535*m, 0.0020391*m, 0.0018576*m, 0.0014819*m, 0.0013017*m, 0.0011655*m, 0.0010766*m, 0.00099684*m, 0.00094434*m, 0.00091114*m, 0.00089118*m, 0.00088207*m, 0.00087604*m, 0.00087648*m, 0.0008757*m, 0.00087531*m, 0.00087192*m, 0.00086611*m, 0.00085566*m, 0.00084235*m, 0.00082451*m, 0.00080527*m, 0.00078153*m, 0.00075099*m, 0.00071544*m, 0.00067982*m, 0.00064916*m, 0.00061689*m, 0.00059114*m, 0.00056492*m, 0.00055354*m, 0.00054367*m, 0.00053863*m, 0.00054048*m, 0.0005435*m, 0.00055229*m, 0.00055835*m, 0.00056823*m, 0.00057385*m, 0.00058484*m, 0.00059158*m, 0.00059677*m, 0.00060353*m, 0.00060439*m, 0.00061041*m, 0.00060961*m, 0.00061121*m, 0.00060773*m, 0.00060566*m, 0.00060501*m, 0.00060341*m, 0.00060624*m, 0.00061129*m, 0.00062536*m, 0.00063788*m, 0.00065966*m, 0.00068083*m, 0.00070613*m, 0.00072164*m, 0.00074646*m, 0.0007687*m, 0.00079297*m, 0.00081584*m, 0.00083679*m, 0.0008547*m, 0.00088603*m, 0.00090622*m, 0.00093579*m, 0.00095204*m, 0.00098722*m, 0.0010177*m, 0.0010531*m, 0.0011063*m, 0.0011426*m, 0.0011984*m, 0.0012503*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 80.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 3.137*m, 3.4732*m, 3.2505*m, 3.3122*m, 3.4578*m, 3.2604*m, 3.0414*m, 3.2767*m, 3.2258*m, 2.9573*m, 2.9506*m, 2.8558*m, 3.0143*m, 3.1577*m, 3.1078*m, 3.0133*m, 3.0716*m, 3.4156*m, 3.0881*m, 2.9263*m, 2.989*m, 3.4606*m, 2.8459*m, 2.9351*m, 3.1946*m, 2.9691*m, 3.0385*m, 2.9172*m, 3.3339*m, 2.6671*m, 2.8318*m, 2.7279*m, 2.7599*m, 2.6115*m, 2.6159*m, 2.7531*m, 2.6329*m, 2.8082*m, 2.6901*m, 2.6323*m, 2.7473*m, 2.7775*m, 2.7909*m, 2.6741*m, 2.549*m, 2.5573*m, 2.6212*m, 2.601*m, 2.468*m, 2.5889*m, 2.4887*m, 2.4084*m, 2.522*m, 2.3279*m, 2.4072*m, 2.3928*m, 2.2722*m, 2.2811*m, 2.4295*m, 2.2422*m, 2.2641*m, 2.1921*m, 2.3386*m, 2.2918*m, 2.2292*m, 2.1517*m, 2.2749*m, 2.2703*m, 2.1667*m, 2.2323*m, 2.145*m, 2.2347*m, 2.0889*m, 2.1528*m, 2.1669*m, 2.1046*m, 2.115*m, 2.0475*m, 2.1449*m, 1.9749*m, 2.0983*m, 1.9965*m, 1.9782*m, 2.036*m, 1.9152*m, 1.9738*m, 1.9665*m, 1.9372*m, 1.9053*m, 2.033*m, 1.9487*m, 1.8454*m, 1.8855*m, 1.9366*m, 1.8314*m, 1.7858*m, 1.9055*m, 1.7373*m, 1.7848*m, 1.7575*m, 1.7849*m, 1.7308*m, 1.6697*m, 1.7396*m, 1.7179*m, 1.6499*m, 1.6014*m, 1.6481*m, 1.658*m, 1.657*m, 1.6127*m, 1.5862*m, 1.6036*m, 1.5799*m, 1.5459*m, 1.6136*m, 1.5807*m, 1.5371*m, 1.5493*m, 1.5064*m, 1.5166*m, 1.4656*m, 1.484*m, 1.469*m, 1.4759*m, 1.448*m, 1.4639*m, 1.4263*m, 1.3932*m, 1.3723*m, 1.3929*m, 1.3935*m, 1.3148*m, 1.2825*m, 1.3222*m, 1.2909*m, 1.2992*m, 1.2694*m, 1.2535*m, 1.2344*m, 1.2713*m, 1.2127*m, 1.1921*m, 1.1997*m, 1.1469*m, 1.1444*m, 1.1285*m, 1.1241*m, 1.108*m, 1.0851*m, 1.0273*m, 1.0195*m, 1.0115*m, 0.9906*m, 0.9805*m, 0.9463*m, 0.9108*m, 0.886*m, 0.8596*m, 0.8239*m, 0.785*m, 0.7499*m, 0.7096*m, 0.6541*m, 0.6005*m, 0.5605*m, 0.5091*m, 0.456*m, 0.4134*m, 0.37*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.3268*m, 0.2888*m, 0.2539*m, 0.2203*m, 0.1857*m, 0.1568*m, 0.1291*m, 0.1022*m, 0.0825*m, 0.0657*m, 0.0515*m, 0.0399*m, 0.0294*m, 0.0222*m, 0.0168*m, 0.0128*m, 0.0095717*m, 0.0069512*m, 0.0053006*m, 0.0041974*m, 0.0032876*m, 0.0026142*m, 0.0021182*m, 0.0018352*m, 0.0016719*m, 0.0013337*m, 0.0011715*m, 0.001049*m, 0.00096891*m, 0.00089716*m, 0.00084991*m, 0.00082002*m, 0.00080207*m, 0.00079387*m, 0.00078843*m, 0.00078883*m, 0.00078813*m, 0.00078778*m, 0.00078473*m, 0.0007795*m, 0.0007701*m, 0.00075812*m, 0.00074206*m, 0.00072474*m, 0.00070338*m, 0.00067589*m, 0.00064389*m, 0.00061184*m, 0.00058425*m, 0.0005552*m, 0.00053203*m, 0.00050843*m, 0.00049818*m, 0.0004893*m, 0.00048477*m, 0.00048643*m, 0.00048915*m, 0.00049706*m, 0.00050252*m, 0.00051141*m, 0.00051646*m, 0.00052636*m, 0.00053242*m, 0.0005371*m, 0.00054318*m, 0.00054395*m, 0.00054937*m, 0.00054865*m, 0.00055009*m, 0.00054696*m, 0.0005451*m, 0.00054451*m, 0.00054307*m, 0.00054562*m, 0.00055016*m, 0.00056283*m, 0.00057409*m, 0.00059369*m, 0.00061275*m, 0.00063551*m, 0.00064948*m, 0.00067181*m, 0.00069183*m, 0.00071367*m, 0.00073426*m, 0.00075311*m, 0.00076923*m, 0.00079742*m, 0.0008156*m, 0.00084221*m, 0.00085684*m, 0.0008885*m, 0.00091589*m, 0.00094779*m, 0.00099569*m, 0.0010284*m, 0.0010785*m, 0.0011253*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 88.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.9357*m, 3.2606*m, 3.0452*m, 3.1048*m, 3.2457*m, 3.0547*m, 2.8437*m, 3.0705*m, 3.0214*m, 2.7629*m, 2.7564*m, 2.6656*m, 2.8176*m, 2.9557*m, 2.9076*m, 2.8167*m, 2.8727*m, 3.2048*m, 2.8887*m, 2.7331*m, 2.7934*m, 3.2484*m, 2.656*m, 2.7416*m, 2.9912*m, 2.7742*m, 2.8409*m, 2.7244*m, 3.1258*m, 2.485*m, 2.6425*m, 2.5431*m, 2.5738*m, 2.432*m, 2.4362*m, 2.5672*m, 2.4523*m, 2.6199*m, 2.507*m, 2.4518*m, 2.5617*m, 2.5905*m, 2.6034*m, 2.4917*m, 2.3723*m, 2.3803*m, 2.4412*m, 2.4219*m, 2.2952*m, 2.4104*m, 2.3149*m, 2.2386*m, 2.3466*m, 2.1622*m, 2.2375*m, 2.2237*m, 2.1093*m, 2.1177*m, 2.2586*m, 2.0809*m, 2.1017*m, 2.0334*m, 2.1723*m, 2.1278*m, 2.0686*m, 1.9952*m, 2.1118*m, 2.1075*m, 2.0094*m, 2.0715*m, 1.9889*m, 2.0738*m, 1.9358*m, 1.9963*m, 2.0096*m, 1.9507*m, 1.9605*m, 1.8967*m, 1.9888*m, 1.8283*m, 1.9447*m, 1.8486*m, 1.8314*m, 1.8859*m, 1.772*m, 1.8272*m, 1.8204*m, 1.7928*m, 1.7627*m, 1.8831*m, 1.8036*m, 1.7063*m, 1.744*m, 1.7922*m, 1.6932*m, 1.6503*m, 1.7629*m, 1.6047*m, 1.6494*m, 1.6237*m, 1.6495*m, 1.5986*m, 1.5413*m, 1.607*m, 1.5865*m, 1.5228*m, 1.4774*m, 1.5212*m, 1.5304*m, 1.5295*m, 1.488*m, 1.4632*m, 1.4794*m, 1.4572*m, 1.4255*m, 1.4888*m, 1.458*m, 1.4172*m, 1.4286*m, 1.3885*m, 1.398*m, 1.3504*m, 1.3676*m, 1.3536*m, 1.36*m, 1.334*m, 1.3488*m, 1.3137*m, 1.2828*m, 1.2634*m, 1.2826*m, 1.2831*m, 1.2098*m, 1.1797*m, 1.2167*m, 1.1875*m, 1.1953*m, 1.1675*m, 1.1528*m, 1.135*m, 1.1693*m, 1.1148*m, 1.0956*m, 1.1027*m, 1.0536*m, 1.0513*m, 1.0366*m, 1.0325*m, 1.0175*m, 0.9963*m, 0.9428*m, 0.9355*m, 0.9281*m, 0.9087*m, 0.8994*m, 0.8678*m, 0.8349*m, 0.812*m, 0.7876*m, 0.7546*m, 0.7188*m, 0.6864*m, 0.6493*m, 0.5982*m, 0.5489*m, 0.5122*m, 0.465*m, 0.4162*m, 0.3772*m, 0.3375*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.298*m, 0.2632*m, 0.2314*m, 0.2006*m, 0.1691*m, 0.1428*m, 0.1175*m, 0.093*m, 0.075*m, 0.0598*m, 0.0469*m, 0.0363*m, 0.0267*m, 0.0202*m, 0.0152*m, 0.0116*m, 0.0087023*m, 0.0063197*m, 0.0048189*m, 0.003816*m, 0.0029888*m, 0.0023766*m, 0.0019257*m, 0.0016684*m, 0.0015199*m, 0.0012125*m, 0.001065*m, 0.00095362*m, 0.00088083*m, 0.0008156*m, 0.00077264*m, 0.00074547*m, 0.00072915*m, 0.0007217*m, 0.00071676*m, 0.00071712*m, 0.00071648*m, 0.00071616*m, 0.00071339*m, 0.00070863*m, 0.00070009*m, 0.0006892*m, 0.0006746*m, 0.00065886*m, 0.00063943*m, 0.00061445*m, 0.00058536*m, 0.00055622*m, 0.00053113*m, 0.00050473*m, 0.00048366*m, 0.00046221*m, 0.00045289*m, 0.00044482*m, 0.0004407*m, 0.00044221*m, 0.00044468*m, 0.00045187*m, 0.00045683*m, 0.00046492*m, 0.00046951*m, 0.00047851*m, 0.00048402*m, 0.00048827*m, 0.0004938*m, 0.0004945*m, 0.00049943*m, 0.00049877*m, 0.00050008*m, 0.00049723*m, 0.00049554*m, 0.00049501*m, 0.0004937*m, 0.00049602*m, 0.00050015*m, 0.00051166*m, 0.0005219*m, 0.00053972*m, 0.00055704*m, 0.00057774*m, 0.00059043*m, 0.00061074*m, 0.00062893*m, 0.00064879*m, 0.00066751*m, 0.00068465*m, 0.0006993*m, 0.00072493*m, 0.00074146*m, 0.00076565*m, 0.00077894*m, 0.00080773*m, 0.00083263*m, 0.00086162*m, 0.00090517*m, 0.00093488*m, 0.00098048*m, 0.001023*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 96.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.7587*m, 3.0726*m, 2.8643*m, 2.9218*m, 3.0581*m, 2.8735*m, 2.6701*m, 2.8887*m, 2.8413*m, 2.5925*m, 2.5863*m, 2.4991*m, 2.6451*m, 2.778*m, 2.7316*m, 2.6442*m, 2.6981*m, 3.0186*m, 2.7134*m, 2.5639*m, 2.6218*m, 3.0608*m, 2.4899*m, 2.572*m, 2.8122*m, 2.6034*m, 2.6674*m, 2.5555*m, 2.9421*m, 2.3262*m, 2.477*m, 2.3818*m, 2.4111*m, 2.2755*m, 2.2795*m, 2.4048*m, 2.295*m, 2.4553*m, 2.3472*m, 2.2945*m, 2.3995*m, 2.4272*m, 2.4395*m, 2.3326*m, 2.2186*m, 2.2262*m, 2.2843*m, 2.2659*m, 2.1451*m, 2.255*m, 2.1639*m, 2.0911*m, 2.1941*m, 2.0184*m, 2.0901*m, 2.077*m, 1.9682*m, 1.9762*m, 2.1102*m, 1.9412*m, 1.9609*m, 1.8962*m, 2.028*m, 1.9858*m, 1.9295*m, 1.8599*m, 1.9706*m, 1.9665*m, 1.8734*m, 1.9323*m, 1.8539*m, 1.9345*m, 1.8037*m, 1.8609*m, 1.8735*m, 1.8178*m, 1.8271*m, 1.7667*m, 1.8539*m, 1.7019*m, 1.8121*m, 1.7211*m, 1.7048*m, 1.7564*m, 1.6487*m, 1.7009*m, 1.6944*m, 1.6683*m, 1.64*m, 1.7537*m, 1.6786*m, 1.5867*m, 1.6223*m, 1.6678*m, 1.5744*m, 1.5339*m, 1.6401*m, 1.491*m, 1.5331*m, 1.5089*m, 1.5331*m, 1.4852*m, 1.4313*m, 1.4931*m, 1.4739*m, 1.4138*m, 1.3712*m, 1.4123*m, 1.421*m, 1.4202*m, 1.3811*m, 1.3578*m, 1.3731*m, 1.3522*m, 1.3224*m, 1.3819*m, 1.353*m, 1.3146*m, 1.3254*m, 1.2877*m, 1.2967*m, 1.252*m, 1.2681*m, 1.255*m, 1.261*m, 1.2366*m, 1.2505*m, 1.2176*m, 1.1887*m, 1.1704*m, 1.1884*m, 1.189*m, 1.1203*m, 1.0921*m, 1.1267*m, 1.0995*m, 1.1067*m, 1.0807*m, 1.067*m, 1.0503*m, 1.0824*m, 1.0315*m, 1.0136*m, 1.0202*m, 0.9744*m, 0.9723*m, 0.9585*m, 0.9547*m, 0.9407*m, 0.9209*m, 0.8711*m, 0.8643*m, 0.8574*m, 0.8394*m, 0.8307*m, 0.8013*m, 0.7707*m, 0.7494*m, 0.7268*m, 0.6961*m, 0.6629*m, 0.6328*m, 0.5984*m, 0.5511*m, 0.5055*m, 0.4715*m, 0.4279*m, 0.3829*m, 0.3469*m, 0.3102*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2739*m, 0.2418*m, 0.2125*m, 0.1842*m, 0.1552*m, 0.131*m, 0.1078*m, 0.0853*m, 0.0688*m, 0.0548*m, 0.043*m, 0.0333*m, 0.0245*m, 0.0185*m, 0.014*m, 0.0107*m, 0.0079777*m, 0.0057934*m, 0.0044175*m, 0.0034981*m, 0.0027398*m, 0.0021786*m, 0.0017653*m, 0.0015294*m, 0.0013933*m, 0.0011114*m, 0.00097626*m, 0.00087415*m, 0.00080742*m, 0.00074763*m, 0.00070826*m, 0.00068335*m, 0.00066839*m, 0.00066155*m, 0.00065703*m, 0.00065736*m, 0.00065677*m, 0.00065648*m, 0.00065394*m, 0.00064958*m, 0.00064175*m, 0.00063176*m, 0.00061838*m, 0.00060395*m, 0.00058615*m, 0.00056325*m, 0.00053658*m, 0.00050987*m, 0.00048687*m, 0.00046266*m, 0.00044336*m, 0.00042369*m, 0.00041515*m, 0.00040775*m, 0.00040398*m, 0.00040536*m, 0.00040763*m, 0.00041422*m, 0.00041876*m, 0.00042617*m, 0.00043039*m, 0.00043863*m, 0.00044368*m, 0.00044758*m, 0.00045265*m, 0.00045329*m, 0.00045781*m, 0.00045721*m, 0.0004584*m, 0.0004558*m, 0.00045425*m, 0.00045376*m, 0.00045256*m, 0.00045468*m, 0.00045847*m, 0.00046902*m, 0.00047841*m, 0.00049474*m, 0.00051062*m, 0.00052959*m, 0.00054123*m, 0.00055984*m, 0.00057652*m, 0.00059473*m, 0.00061188*m, 0.0006276*m, 0.00064103*m, 0.00066452*m, 0.00067967*m, 0.00070184*m, 0.00071403*m, 0.00074042*m, 0.00076324*m, 0.00078982*m, 0.00082974*m, 0.00085698*m, 0.00089877*m, 0.00093772*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 104.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.6018*m, 2.905*m, 2.7036*m, 2.7592*m, 2.8911*m, 2.7125*m, 2.5165*m, 2.7272*m, 2.6815*m, 2.4419*m, 2.4359*m, 2.3522*m, 2.4924*m, 2.6204*m, 2.5757*m, 2.4916*m, 2.5434*m, 2.8527*m, 2.5582*m, 2.4144*m, 2.47*m, 2.8936*m, 2.3434*m, 2.4222*m, 2.6534*m, 2.4523*m, 2.514*m, 2.4063*m, 2.7788*m, 2.1865*m, 2.331*m, 2.2397*m, 2.2678*m, 2.138*m, 2.1418*m, 2.2618*m, 2.1566*m, 2.3102*m, 2.2066*m, 2.1561*m, 2.2567*m, 2.2832*m, 2.295*m, 2.1926*m, 2.0836*m, 2.0908*m, 2.1464*m, 2.1288*m, 2.0134*m, 2.1183*m, 2.0313*m, 1.9619*m, 2.0602*m, 1.8926*m, 1.9609*m, 1.9484*m, 1.8448*m, 1.8524*m, 1.9801*m, 1.8191*m, 1.8379*m, 1.7763*m, 1.9018*m, 1.8616*m, 1.808*m, 1.7418*m, 1.8471*m, 1.8431*m, 1.7546*m, 1.8106*m, 1.7361*m, 1.8127*m, 1.6884*m, 1.7428*m, 1.7548*m, 1.7018*m, 1.7106*m, 1.6533*m, 1.7361*m, 1.5919*m, 1.6964*m, 1.6101*m, 1.5946*m, 1.6436*m, 1.5415*m, 1.5909*m, 1.5848*m, 1.5601*m, 1.5332*m, 1.641*m, 1.5697*m, 1.4828*m, 1.5165*m, 1.5596*m, 1.4711*m, 1.4328*m, 1.5333*m, 1.3923*m, 1.4321*m, 1.4092*m, 1.4321*m, 1.3869*m, 1.3359*m, 1.3943*m, 1.3761*m, 1.3195*m, 1.2792*m, 1.318*m, 1.3263*m, 1.3255*m, 1.2886*m, 1.2666*m, 1.2811*m, 1.2614*m, 1.2333*m, 1.2894*m, 1.2621*m, 1.2259*m, 1.2361*m, 1.2006*m, 1.209*m, 1.1669*m, 1.1821*m, 1.1698*m, 1.1755*m, 1.1525*m, 1.1655*m, 1.1346*m, 1.1074*m, 1.0902*m, 1.1071*m, 1.1077*m, 1.0431*m, 1.0167*m, 1.0492*m, 1.0236*m, 1.0304*m, 1.006*m, 0.9931*m, 0.9775*m, 1.0076*m, 0.9598*m, 0.943*m, 0.9492*m, 0.9063*m, 0.9043*m, 0.8913*m, 0.8878*m, 0.8747*m, 0.8562*m, 0.8095*m, 0.8032*m, 0.7967*m, 0.7799*m, 0.7717*m, 0.7442*m, 0.7157*m, 0.6958*m, 0.6747*m, 0.6461*m, 0.615*m, 0.587*m, 0.555*m, 0.5109*m, 0.4684*m, 0.4369*m, 0.3963*m, 0.3545*m, 0.321*m, 0.2871*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2533*m, 0.2236*m, 0.1965*m, 0.1703*m, 0.1434*m, 0.1211*m, 0.0996*m, 0.0788*m, 0.0635*m, 0.0506*m, 0.0397*m, 0.0307*m, 0.0226*m, 0.0171*m, 0.0129*m, 0.0098466*m, 0.0073645*m, 0.005348*m, 0.0040779*m, 0.0032291*m, 0.0025291*m, 0.0020111*m, 0.0016295*m, 0.0014118*m, 0.0012861*m, 0.0010259*m, 0.00090116*m, 0.00080691*m, 0.00074531*m, 0.00069012*m, 0.00065377*m, 0.00063079*m, 0.00061697*m, 0.00061067*m, 0.00060649*m, 0.00060679*m, 0.00060625*m, 0.00060598*m, 0.00060364*m, 0.00059961*m, 0.00059238*m, 0.00058317*m, 0.00057082*m, 0.00055749*m, 0.00054106*m, 0.00051992*m, 0.0004953*m, 0.00047065*m, 0.00044942*m, 0.00042708*m, 0.00040925*m, 0.0003911*m, 0.00038322*m, 0.00037638*m, 0.0003729*m, 0.00037418*m, 0.00037627*m, 0.00038236*m, 0.00038655*m, 0.00039339*m, 0.00039728*m, 0.00040489*m, 0.00040955*m, 0.00041315*m, 0.00041783*m, 0.00041842*m, 0.00042259*m, 0.00042204*m, 0.00042314*m, 0.00042074*m, 0.0004193*m, 0.00041885*m, 0.00041775*m, 0.00041971*m, 0.0004232*m, 0.00043294*m, 0.00044161*m, 0.00045668*m, 0.00047134*m, 0.00048886*m, 0.0004996*m, 0.00051678*m, 0.00053217*m, 0.00054898*m, 0.00056481*m, 0.00057932*m, 0.00059172*m, 0.0006134*m, 0.00062739*m, 0.00064786*m, 0.00065911*m, 0.00068346*m, 0.00070453*m, 0.00072907*m, 0.00076591*m, 0.00079106*m, 0.00082964*m, 0.00086559*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 112.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.4618*m, 2.7548*m, 2.5601*m, 2.6138*m, 2.7413*m, 2.5687*m, 2.3797*m, 2.5829*m, 2.5387*m, 2.3078*m, 2.3021*m, 2.2216*m, 2.3564*m, 2.4797*m, 2.4367*m, 2.3556*m, 2.4056*m, 2.7042*m, 2.4198*m, 2.2814*m, 2.3349*m, 2.7437*m, 2.2131*m, 2.2889*m, 2.5116*m, 2.3178*m, 2.3772*m, 2.2736*m, 2.6327*m, 2.0626*m, 2.2012*m, 2.1136*m, 2.1406*m, 2.0161*m, 2.0198*m, 2.1348*m, 2.034*m, 2.1812*m, 2.0819*m, 2.0335*m, 2.1299*m, 2.1553*m, 2.1667*m, 2.0685*m, 1.9641*m, 1.971*m, 2.0242*m, 2.0074*m, 1.8969*m, 1.9973*m, 1.914*m, 1.8477*m, 1.9416*m, 1.7816*m, 1.8468*m, 1.8349*m, 1.736*m, 1.7432*m, 1.8651*m, 1.7115*m, 1.7294*m, 1.6707*m, 1.7903*m, 1.7519*m, 1.7009*m, 1.6378*m, 1.7381*m, 1.7344*m, 1.65*m, 1.7034*m, 1.6324*m, 1.7054*m, 1.587*m, 1.6388*m, 1.6502*m, 1.5997*m, 1.6081*m, 1.5536*m, 1.6324*m, 1.4952*m, 1.5946*m, 1.5125*m, 1.4978*m, 1.5443*m, 1.4473*m, 1.4943*m, 1.4885*m, 1.465*m, 1.4395*m, 1.5419*m, 1.4742*m, 1.3917*m, 1.4236*m, 1.4645*m, 1.3806*m, 1.3443*m, 1.4396*m, 1.3059*m, 1.3435*m, 1.3219*m, 1.3436*m, 1.3007*m, 1.2525*m, 1.3077*m, 1.2905*m, 1.2369*m, 1.1988*m, 1.2356*m, 1.2433*m, 1.2426*m, 1.2077*m, 1.1869*m, 1.2006*m, 1.1819*m, 1.1554*m, 1.2084*m, 1.1826*m, 1.1485*m, 1.158*m, 1.1245*m, 1.1324*m, 1.0927*m, 1.1071*m, 1.0954*m, 1.1008*m, 1.079*m, 1.0914*m, 1.0622*m, 1.0365*m, 1.0203*m, 1.0363*m, 1.0368*m, 0.9759*m, 0.951*m, 0.9816*m, 0.9575*m, 0.9639*m, 0.9409*m, 0.9287*m, 0.914*m, 0.9424*m, 0.8974*m, 0.8816*m, 0.8874*m, 0.847*m, 0.8451*m, 0.833*m, 0.8296*m, 0.8173*m, 0.7999*m, 0.7561*m, 0.7501*m, 0.7441*m, 0.7282*m, 0.7206*m, 0.6948*m, 0.668*m, 0.6493*m, 0.6295*m, 0.6027*m, 0.5736*m, 0.5474*m, 0.5174*m, 0.4761*m, 0.4364*m, 0.4069*m, 0.369*m, 0.33*m, 0.2988*m, 0.2671*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2357*m, 0.208*m, 0.1827*m, 0.1583*m, 0.1333*m, 0.1125*m, 0.0925*m, 0.0732*m, 0.059*m, 0.047*m, 0.0369*m, 0.0285*m, 0.021*m, 0.0159*m, 0.012*m, 0.009144*m, 0.0068388*m, 0.0049662*m, 0.0037867*m, 0.0029985*m, 0.0023485*m, 0.0018675*m, 0.0015131*m, 0.0013109*m, 0.0011943*m, 0.00095266*m, 0.00083679*m, 0.00074927*m, 0.00069208*m, 0.00064083*m, 0.00060708*m, 0.00058573*m, 0.0005729*m, 0.00056705*m, 0.00056317*m, 0.00056345*m, 0.00056295*m, 0.0005627*m, 0.00056052*m, 0.00055678*m, 0.00055007*m, 0.00054151*m, 0.00053004*m, 0.00051767*m, 0.00050241*m, 0.00048278*m, 0.00045992*m, 0.00043703*m, 0.00041732*m, 0.00039657*m, 0.00038002*m, 0.00036317*m, 0.00035585*m, 0.0003495*m, 0.00034626*m, 0.00034745*m, 0.0003494*m, 0.00035504*m, 0.00035894*m, 0.00036529*m, 0.0003689*m, 0.00037597*m, 0.0003803*m, 0.00038364*m, 0.00038798*m, 0.00038854*m, 0.00039241*m, 0.00039189*m, 0.00039292*m, 0.00039068*m, 0.00038935*m, 0.00038894*m, 0.00038791*m, 0.00038973*m, 0.00039297*m, 0.00040202*m, 0.00041007*m, 0.00042406*m, 0.00043768*m, 0.00045394*m, 0.00046391*m, 0.00047987*m, 0.00049416*m, 0.00050977*m, 0.00052447*m, 0.00053794*m, 0.00054945*m, 0.00056959*m, 0.00058257*m, 0.00060158*m, 0.00061203*m, 0.00063464*m, 0.00065421*m, 0.00067699*m, 0.00071121*m, 0.00073455*m, 0.00077038*m, 0.00080376*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 120.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.3362*m, 2.6194*m, 2.431*m, 2.4829*m, 2.6063*m, 2.4393*m, 2.2569*m, 2.453*m, 2.4103*m, 2.1877*m, 2.1822*m, 2.1047*m, 2.2345*m, 2.3534*m, 2.3119*m, 2.2337*m, 2.2819*m, 2.5704*m, 2.2956*m, 2.1622*m, 2.2137*m, 2.6087*m, 2.0966*m, 2.1695*m, 2.3842*m, 2.1974*m, 2.2545*m, 2.1548*m, 2.5012*m, 1.952*m, 2.0851*m, 2.001*m, 2.0268*m, 1.9074*m, 1.9109*m, 2.0213*m, 1.9245*m, 2.0659*m, 1.9705*m, 1.9241*m, 2.0166*m, 2.041*m, 2.0519*m, 1.9576*m, 1.8575*m, 1.8641*m, 1.9151*m, 1.899*m, 1.7932*m, 1.8894*m, 1.8096*m, 1.7461*m, 1.836*m, 1.6828*m, 1.7452*m, 1.7338*m, 1.6392*m, 1.6461*m, 1.7627*m, 1.6159*m, 1.6329*m, 1.5769*m, 1.6912*m, 1.6545*m, 1.6057*m, 1.5456*m, 1.6413*m, 1.6377*m, 1.5572*m, 1.6081*m, 1.5404*m, 1.61*m, 1.4971*m, 1.5464*m, 1.5573*m, 1.5092*m, 1.5172*m, 1.4652*m, 1.5403*m, 1.4096*m, 1.5043*m, 1.4261*m, 1.4121*m, 1.4564*m, 1.364*m, 1.4087*m, 1.4032*m, 1.3808*m, 1.3566*m, 1.4541*m, 1.3896*m, 1.3111*m, 1.3415*m, 1.3804*m, 1.3005*m, 1.266*m, 1.3567*m, 1.2295*m, 1.2653*m, 1.2448*m, 1.2654*m, 1.2247*m, 1.1789*m, 1.2313*m, 1.215*m, 1.1641*m, 1.1279*m, 1.1628*m, 1.1702*m, 1.1694*m, 1.1364*m, 1.1166*m, 1.1296*m, 1.1119*m, 1.0867*m, 1.1371*m, 1.1125*m, 1.0802*m, 1.0893*m, 1.0575*m, 1.065*m, 1.0273*m, 1.041*m, 1.0299*m, 1.035*m, 1.0144*m, 1.0261*m, 0.9985*m, 0.9741*m, 0.9588*m, 0.9739*m, 0.9744*m, 0.9168*m, 0.8932*m, 0.9222*m, 0.8994*m, 0.9055*m, 0.8837*m, 0.8722*m, 0.8583*m, 0.8851*m, 0.8426*m, 0.8277*m, 0.8332*m, 0.795*m, 0.7933*m, 0.7818*m, 0.7786*m, 0.767*m, 0.7506*m, 0.7092*m, 0.7036*m, 0.6979*m, 0.683*m, 0.6758*m, 0.6515*m, 0.6263*m, 0.6087*m, 0.59*m, 0.5648*m, 0.5375*m, 0.5128*m, 0.4846*m, 0.4458*m, 0.4085*m, 0.3808*m, 0.3453*m, 0.3087*m, 0.2794*m, 0.2498*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2203*m, 0.1944*m, 0.1707*m, 0.1479*m, 0.1245*m, 0.1051*m, 0.0864*m, 0.0684*m, 0.0551*m, 0.0439*m, 0.0344*m, 0.0266*m, 0.0196*m, 0.0148*m, 0.0112*m, 0.0085349*m, 0.0063832*m, 0.0046352*m, 0.0035343*m, 0.0027987*m, 0.0021919*m, 0.001743*m, 0.0014123*m, 0.0012236*m, 0.0011147*m, 0.00088915*m, 0.00078101*m, 0.00069932*m, 0.00064594*m, 0.00059811*m, 0.0005666*m, 0.00054668*m, 0.00053471*m, 0.00052924*m, 0.00052562*m, 0.00052589*m, 0.00052542*m, 0.00052519*m, 0.00052315*m, 0.00051966*m, 0.0005134*m, 0.00050541*m, 0.00049471*m, 0.00048316*m, 0.00046892*m, 0.0004506*m, 0.00042926*m, 0.00040789*m, 0.0003895*m, 0.00037013*m, 0.00035469*m, 0.00033895*m, 0.00033212*m, 0.0003262*m, 0.00032318*m, 0.00032429*m, 0.0003261*m, 0.00033137*m, 0.00033501*m, 0.00034094*m, 0.00034431*m, 0.00035091*m, 0.00035495*m, 0.00035806*m, 0.00036212*m, 0.00036263*m, 0.00036625*m, 0.00036577*m, 0.00036672*m, 0.00036464*m, 0.0003634*m, 0.00036301*m, 0.00036205*m, 0.00036375*m, 0.00036677*m, 0.00037522*m, 0.00038273*m, 0.00039579*m, 0.0004085*m, 0.00042368*m, 0.00043298*m, 0.00044788*m, 0.00046122*m, 0.00047578*m, 0.00048951*m, 0.00050208*m, 0.00051282*m, 0.00053162*m, 0.00054373*m, 0.00056147*m, 0.00057123*m, 0.00059233*m, 0.0006106*m, 0.00063186*m, 0.00066379*m, 0.00068558*m, 0.00071902*m, 0.00075018*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 128.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.2227*m, 2.4966*m, 2.3143*m, 2.3645*m, 2.4839*m, 2.3223*m, 2.1462*m, 2.3356*m, 2.2943*m, 2.0794*m, 2.0741*m, 1.9995*m, 2.1246*m, 2.2393*m, 2.1992*m, 2.1238*m, 2.1703*m, 2.4491*m, 2.1835*m, 2.0549*m, 2.1046*m, 2.4862*m, 1.9917*m, 2.0619*m, 2.2691*m, 2.0888*m, 2.1439*m, 2.0477*m, 2.3822*m, 1.8526*m, 1.9807*m, 1.8997*m, 1.9246*m, 1.8098*m, 1.8132*m, 1.9192*m, 1.8262*m, 1.9622*m, 1.8704*m, 1.8258*m, 1.9148*m, 1.9383*m, 1.9487*m, 1.8581*m, 1.7619*m, 1.7683*m, 1.8172*m, 1.8018*m, 1.7002*m, 1.7925*m, 1.7159*m, 1.6551*m, 1.7413*m, 1.5945*m, 1.6542*m, 1.6433*m, 1.5527*m, 1.5593*m, 1.671*m, 1.5304*m, 1.5467*m, 1.4931*m, 1.6024*m, 1.5673*m, 1.5206*m, 1.4631*m, 1.5547*m, 1.5513*m, 1.4743*m, 1.523*m, 1.4582*m, 1.5248*m, 1.4168*m, 1.464*m, 1.4744*m, 1.4284*m, 1.4361*m, 1.3864*m, 1.4581*m, 1.3333*m, 1.4237*m, 1.349*m, 1.3357*m, 1.378*m, 1.2898*m, 1.3325*m, 1.3272*m, 1.3058*m, 1.2827*m, 1.3758*m, 1.3142*m, 1.2393*m, 1.2683*m, 1.3054*m, 1.2293*m, 1.1964*m, 1.2828*m, 1.1616*m, 1.1957*m, 1.1761*m, 1.1958*m, 1.157*m, 1.1134*m, 1.1633*m, 1.1478*m, 1.0993*m, 1.065*m, 1.0981*m, 1.1051*m, 1.1044*m, 1.073*m, 1.0542*m, 1.0665*m, 1.0497*m, 1.0258*m, 1.0736*m, 1.0503*m, 1.0196*m, 1.0282*m, 0.998*m, 1.0052*m, 0.9694*m, 0.9823*m, 0.9718*m, 0.9766*m, 0.9571*m, 0.9682*m, 0.9419*m, 0.9189*m, 0.9043*m, 0.9187*m, 0.9191*m, 0.8645*m, 0.8421*m, 0.8696*m, 0.848*m, 0.8537*m, 0.8331*m, 0.8222*m, 0.809*m, 0.8345*m, 0.7941*m, 0.78*m, 0.7852*m, 0.7491*m, 0.7474*m, 0.7366*m, 0.7336*m, 0.7226*m, 0.707*m, 0.6679*m, 0.6626*m, 0.6572*m, 0.6431*m, 0.6363*m, 0.6133*m, 0.5894*m, 0.5728*m, 0.5552*m, 0.5314*m, 0.5056*m, 0.4823*m, 0.4557*m, 0.4191*m, 0.384*m, 0.3579*m, 0.3244*m, 0.2899*m, 0.2624*m, 0.2345*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.2068*m, 0.1825*m, 0.1602*m, 0.1388*m, 0.1169*m, 0.0986*m, 0.0811*m, 0.0641*m, 0.0517*m, 0.0412*m, 0.0323*m, 0.025*m, 0.0184*m, 0.0139*m, 0.0105*m, 0.0080019*m, 0.0059845*m, 0.0043457*m, 0.0033135*m, 0.0026238*m, 0.002055*m, 0.0016341*m, 0.001324*m, 0.0011471*m, 0.001045*m, 0.00083358*m, 0.00073219*m, 0.00065561*m, 0.00060557*m, 0.00056072*m, 0.00053119*m, 0.00051251*m, 0.00050129*m, 0.00049617*m, 0.00049277*m, 0.00049302*m, 0.00049258*m, 0.00049236*m, 0.00049046*m, 0.00048719*m, 0.00048131*m, 0.00047382*m, 0.00046379*m, 0.00045296*m, 0.00043961*m, 0.00042243*m, 0.00040243*m, 0.0003824*m, 0.00036515*m, 0.000347*m, 0.00033252*m, 0.00031777*m, 0.00031137*m, 0.00030581*m, 0.00030298*m, 0.00030402*m, 0.00030572*m, 0.00031066*m, 0.00031407*m, 0.00031963*m, 0.00032279*m, 0.00032898*m, 0.00033276*m, 0.00033569*m, 0.00033949*m, 0.00033997*m, 0.00034336*m, 0.00034291*m, 0.0003438*m, 0.00034185*m, 0.00034068*m, 0.00034032*m, 0.00033942*m, 0.00034101*m, 0.00034385*m, 0.00035177*m, 0.00035881*m, 0.00037106*m, 0.00038297*m, 0.0003972*m, 0.00040592*m, 0.00041988*m, 0.00043239*m, 0.00044605*m, 0.00045891*m, 0.0004707*m, 0.00048077*m, 0.00049839*m, 0.00050975*m, 0.00052638*m, 0.00053552*m, 0.00055531*m, 0.00057243*m, 0.00059237*m, 0.00062231*m, 0.00064273*m, 0.00067408*m, 0.00070329*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 136.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.1197*m, 2.3849*m, 2.2083*m, 2.2568*m, 2.3726*m, 2.216*m, 2.0458*m, 2.2289*m, 2.1889*m, 1.9814*m, 1.9763*m, 1.9043*m, 2.025*m, 2.1358*m, 2.097*m, 2.0242*m, 2.0691*m, 2.3388*m, 2.0819*m, 1.9578*m, 2.0057*m, 2.3748*m, 1.8968*m, 1.9645*m, 2.1645*m, 1.9904*m, 2.0436*m, 1.9508*m, 2.274*m, 1.7629*m, 1.8862*m, 1.8082*m, 1.8322*m, 1.7217*m, 1.725*m, 1.827*m, 1.7375*m, 1.8684*m, 1.78*m, 1.7371*m, 1.8227*m, 1.8453*m, 1.8554*m, 1.7681*m, 1.6757*m, 1.6818*m, 1.7289*m, 1.714*m, 1.6164*m, 1.7051*m, 1.6315*m, 1.5731*m, 1.6559*m, 1.5149*m, 1.5722*m, 1.5617*m, 1.4749*m, 1.4812*m, 1.5884*m, 1.4535*m, 1.4691*m, 1.4177*m, 1.5226*m, 1.4889*m, 1.4441*m, 1.389*m, 1.4768*m, 1.4735*m, 1.3997*m, 1.4464*m, 1.3843*m, 1.4481*m, 1.3447*m, 1.3898*m, 1.3998*m, 1.3558*m, 1.3631*m, 1.3156*m, 1.3843*m, 1.2648*m, 1.3513*m, 1.2798*m, 1.2671*m, 1.3075*m, 1.2232*m, 1.264*m, 1.259*m, 1.2386*m, 1.2164*m, 1.3054*m, 1.2465*m, 1.175*m, 1.2027*m, 1.2382*m, 1.1654*m, 1.134*m, 1.2165*m, 1.1008*m, 1.1334*m, 1.1147*m, 1.1334*m, 1.0964*m, 1.0548*m, 1.1025*m, 1.0876*m, 1.0414*m, 1.0087*m, 1.0403*m, 1.0469*m, 1.0463*m, 1.0163*m, 0.9984*m, 1.0101*m, 0.9941*m, 0.9713*m, 1.0169*m, 0.9947*m, 0.9654*m, 0.9736*m, 0.9449*m, 0.9517*m, 0.9176*m, 0.9299*m, 0.9199*m, 0.9245*m, 0.9059*m, 0.9165*m, 0.8915*m, 0.8695*m, 0.8557*m, 0.8693*m, 0.8697*m, 0.8178*m, 0.7965*m, 0.8226*m, 0.8021*m, 0.8076*m, 0.788*m, 0.7776*m, 0.7651*m, 0.7892*m, 0.7509*m, 0.7375*m, 0.7425*m, 0.7081*m, 0.7066*m, 0.6963*m, 0.6934*m, 0.683*m, 0.6682*m, 0.6311*m, 0.626*m, 0.6209*m, 0.6075*m, 0.6011*m, 0.5793*m, 0.5567*m, 0.5409*m, 0.5243*m, 0.5017*m, 0.4773*m, 0.4552*m, 0.43*m, 0.3954*m, 0.3622*m, 0.3375*m, 0.3059*m, 0.2734*m, 0.2474*m, 0.221*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.1949*m, 0.1719*m, 0.1509*m, 0.1308*m, 0.1101*m, 0.0928*m, 0.0763*m, 0.0604*m, 0.0487*m, 0.0387*m, 0.0304*m, 0.0235*m, 0.0173*m, 0.0131*m, 0.0098614*m, 0.0075315*m, 0.0056326*m, 0.0040901*m, 0.0031187*m, 0.0024695*m, 0.0019341*m, 0.001538*m, 0.0012461*m, 0.0010796*m, 0.00098354*m, 0.00078454*m, 0.00068912*m, 0.00061705*m, 0.00056995*m, 0.00052774*m, 0.00049995*m, 0.00048237*m, 0.0004718*m, 0.00046698*m, 0.00046378*m, 0.00046402*m, 0.0004636*m, 0.0004634*m, 0.00046161*m, 0.00045853*m, 0.000453*m, 0.00044595*m, 0.00043651*m, 0.00042632*m, 0.00041375*m, 0.00039759*m, 0.00037876*m, 0.00035991*m, 0.00034367*m, 0.00032659*m, 0.00031296*m, 0.00029908*m, 0.00029305*m, 0.00028782*m, 0.00028516*m, 0.00028613*m, 0.00028774*m, 0.00029239*m, 0.0002956*m, 0.00030083*m, 0.0003038*m, 0.00030962*m, 0.00031319*m, 0.00031594*m, 0.00031952*m, 0.00031997*m, 0.00032316*m, 0.00032274*m, 0.00032358*m, 0.00032174*m, 0.00032064*m, 0.0003203*m, 0.00031945*m, 0.00032095*m, 0.00032362*m, 0.00033107*m, 0.0003377*m, 0.00034923*m, 0.00036044*m, 0.00037383*m, 0.00038205*m, 0.00039518*m, 0.00040696*m, 0.00041981*m, 0.00043192*m, 0.00044301*m, 0.00045249*m, 0.00046907*m, 0.00047977*m, 0.00049542*m, 0.00050402*m, 0.00052265*m, 0.00053876*m, 0.00055752*m, 0.0005857*m, 0.00060493*m, 0.00063443*m, 0.00066192*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 144.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 2.0258*m, 2.2827*m, 2.1116*m, 2.1586*m, 2.2707*m, 2.1191*m, 1.9544*m, 2.1315*m, 2.0928*m, 1.8922*m, 1.8872*m, 1.8178*m, 1.9343*m, 2.0414*m, 2.0039*m, 1.9336*m, 1.9769*m, 2.2381*m, 1.9893*m, 1.8694*m, 1.9156*m, 2.2729*m, 1.8105*m, 1.8759*m, 2.0692*m, 1.9009*m, 1.9523*m, 1.8627*m, 2.1752*m, 1.6815*m, 1.8003*m, 1.7251*m, 1.7482*m, 1.6418*m, 1.645*m, 1.7433*m, 1.657*m, 1.7832*m, 1.698*m, 1.6567*m, 1.7391*m, 1.7609*m, 1.7706*m, 1.6865*m, 1.5975*m, 1.6034*m, 1.6487*m, 1.6343*m, 1.5404*m, 1.6258*m, 1.555*m, 1.4988*m, 1.5784*m, 1.4429*m, 1.498*m, 1.4879*m, 1.4045*m, 1.4106*m, 1.5135*m, 1.3839*m, 1.3989*m, 1.3496*m, 1.4503*m, 1.4179*m, 1.375*m, 1.3221*m, 1.4063*m, 1.4031*m, 1.3323*m, 1.3771*m, 1.3176*m, 1.3788*m, 1.2796*m, 1.3229*m, 1.3324*m, 1.2902*m, 1.2972*m, 1.2516*m, 1.3175*m, 1.203*m, 1.2859*m, 1.2174*m, 1.2052*m, 1.2439*m, 1.1632*m, 1.2022*m, 1.1974*m, 1.1779*m, 1.1567*m, 1.2419*m, 1.1855*m, 1.117*m, 1.1435*m, 1.1775*m, 1.1078*m, 1.0778*m, 1.1568*m, 1.0461*m, 1.0772*m, 1.0593*m, 1.0773*m, 1.0419*m, 1.0021*m, 1.0476*m, 1.0335*m, 0.9893*m, 0.958*m, 0.9882*m, 0.9946*m, 0.9939*m, 0.9653*m, 0.9482*m, 0.9594*m, 0.9441*m, 0.9224*m, 0.9659*m, 0.9447*m, 0.9167*m, 0.9245*m, 0.8971*m, 0.9036*m, 0.8711*m, 0.8828*m, 0.8733*m, 0.8777*m, 0.8599*m, 0.87*m, 0.8462*m, 0.8252*m, 0.812*m, 0.825*m, 0.8254*m, 0.7759*m, 0.7556*m, 0.7805*m, 0.7609*m, 0.7661*m, 0.7475*m, 0.7376*m, 0.7257*m, 0.7487*m, 0.7122*m, 0.6994*m, 0.7041*m, 0.6715*m, 0.6699*m, 0.6601*m, 0.6574*m, 0.6475*m, 0.6334*m, 0.5981*m, 0.5933*m, 0.5885*m, 0.5757*m, 0.5696*m, 0.5489*m, 0.5274*m, 0.5124*m, 0.4966*m, 0.4751*m, 0.4519*m, 0.431*m, 0.4071*m, 0.3743*m, 0.3428*m, 0.3194*m, 0.2894*m, 0.2586*m, 0.234*m, 0.209*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.1843*m, 0.1625*m, 0.1427*m, 0.1236*m, 0.104*m, 0.0877*m, 0.0721*m, 0.057*m, 0.046*m, 0.0366*m, 0.0287*m, 0.0222*m, 0.0163*m, 0.0124*m, 0.009314*m, 0.0071134*m, 0.0053199*m, 0.003863*m, 0.0029455*m, 0.0023323*m, 0.0018267*m, 0.0014525*m, 0.0011769*m, 0.0010197*m, 0.00092891*m, 0.00074096*m, 0.00065084*m, 0.00058276*m, 0.00053828*m, 0.00049842*m, 0.00047217*m, 0.00045557*m, 0.00044559*m, 0.00044104*m, 0.00043802*m, 0.00043824*m, 0.00043785*m, 0.00043766*m, 0.00043596*m, 0.00043305*m, 0.00042783*m, 0.00042118*m, 0.00041226*m, 0.00040263*m, 0.00039077*m, 0.0003755*m, 0.00035772*m, 0.00033991*m, 0.00032458*m, 0.00030844*m, 0.00029557*m, 0.00028246*m, 0.00027677*m, 0.00027183*m, 0.00026932*m, 0.00027024*m, 0.00027175*m, 0.00027615*m, 0.00027918*m, 0.00028412*m, 0.00028692*m, 0.00029242*m, 0.00029579*m, 0.00029839*m, 0.00030177*m, 0.0003022*m, 0.00030521*m, 0.00030481*m, 0.0003056*m, 0.00030386*m, 0.00030283*m, 0.00030251*m, 0.00030171*m, 0.00030312*m, 0.00030564*m, 0.00031268*m, 0.00031894*m, 0.00032983*m, 0.00034041*m, 0.00035306*m, 0.00036082*m, 0.00037323*m, 0.00038435*m, 0.00039649*m, 0.00040792*m, 0.0004184*m, 0.00042735*m, 0.00044301*m, 0.00045311*m, 0.0004679*m, 0.00047602*m, 0.00049361*m, 0.00050883*m, 0.00052655*m, 0.00055316*m, 0.00057132*m, 0.00059918*m, 0.00062515*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 152.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 1.94*m, 2.1889*m, 2.023*m, 2.0685*m, 2.1773*m, 2.0302*m, 1.8709*m, 2.0423*m, 2.0048*m, 1.8107*m, 1.8059*m, 1.7388*m, 1.8514*m, 1.955*m, 1.9187*m, 1.8507*m, 1.8926*m, 2.1456*m, 1.9046*m, 1.7886*m, 1.8333*m, 2.1794*m, 1.7318*m, 1.7949*m, 1.9819*m, 1.8191*m, 1.8688*m, 1.7822*m, 2.0846*m, 1.6072*m, 1.7219*m, 1.6493*m, 1.6716*m, 1.569*m, 1.572*m, 1.6668*m, 1.5837*m, 1.7053*m, 1.6231*m, 1.5833*m, 1.6628*m, 1.6839*m, 1.6933*m, 1.6121*m, 1.5263*m, 1.5319*m, 1.5756*m, 1.5618*m, 1.4713*m, 1.5535*m, 1.4853*m, 1.4312*m, 1.5079*m, 1.3775*m, 1.4304*m, 1.4208*m, 1.3405*m, 1.3463*m, 1.4454*m, 1.3207*m, 1.3352*m, 1.2878*m, 1.3845*m, 1.3534*m, 1.3121*m, 1.2613*m, 1.3423*m, 1.3392*m, 1.2711*m, 1.3142*m, 1.257*m, 1.3158*m, 1.2205*m, 1.2621*m, 1.2712*m, 1.2307*m, 1.2374*m, 1.1936*m, 1.2569*m, 1.147*m, 1.2266*m, 1.1608*m, 1.1491*m, 1.1862*m, 1.1088*m, 1.1462*m, 1.1416*m, 1.1229*m, 1.1025*m, 1.1843*m, 1.1302*m, 1.0645*m, 1.0899*m, 1.1225*m, 1.0557*m, 1.0269*m, 1.1026*m, 0.9965*m, 1.0264*m, 1.0092*m, 1.0264*m, 0.9925*m, 0.9544*m, 0.998*m, 0.9844*m, 0.9421*m, 0.9122*m, 0.9411*m, 0.9472*m, 0.9466*m, 0.9192*m, 0.9028*m, 0.9136*m, 0.8989*m, 0.8781*m, 0.9197*m, 0.8994*m, 0.8727*m, 0.8802*m, 0.8539*m, 0.8601*m, 0.829*m, 0.8403*m, 0.8311*m, 0.8353*m, 0.8184*m, 0.828*m, 0.8052*m, 0.7852*m, 0.7726*m, 0.785*m, 0.7854*m, 0.7381*m, 0.7187*m, 0.7425*m, 0.7238*m, 0.7288*m, 0.7109*m, 0.7015*m, 0.6901*m, 0.7121*m, 0.6773*m, 0.665*m, 0.6695*m, 0.6384*m, 0.6369*m, 0.6276*m, 0.625*m, 0.6155*m, 0.6021*m, 0.5684*m, 0.5639*m, 0.5592*m, 0.5471*m, 0.5412*m, 0.5215*m, 0.501*m, 0.4868*m, 0.4717*m, 0.4513*m, 0.4292*m, 0.4092*m, 0.3865*m, 0.3553*m, 0.3253*m, 0.3031*m, 0.2746*m, 0.2453*m, 0.2219*m, 0.1982*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.1747*m, 0.1541*m, 0.1353*m, 0.1172*m, 0.0986*m, 0.0832*m, 0.0684*m, 0.0541*m, 0.0436*m, 0.0347*m, 0.0272*m, 0.021*m, 0.0155*m, 0.0117*m, 0.0088243*m, 0.0067393*m, 0.00504*m, 0.0036598*m, 0.0027905*m, 0.0022096*m, 0.0017306*m, 0.0013761*m, 0.001115*m, 0.000966*m, 0.00088002*m, 0.00070196*m, 0.00061658*m, 0.00055209*m, 0.00050995*m, 0.00047219*m, 0.00044732*m, 0.00043159*m, 0.00042214*m, 0.00041782*m, 0.00041496*m, 0.00041517*m, 0.0004148*m, 0.00041462*m, 0.00041302*m, 0.00041026*m, 0.00040531*m, 0.00039901*m, 0.00039056*m, 0.00038144*m, 0.0003702*m, 0.00035573*m, 0.00033889*m, 0.00032202*m, 0.0003075*m, 0.00029221*m, 0.00028001*m, 0.0002676*m, 0.0002622*m, 0.00025753*m, 0.00025514*m, 0.00025602*m, 0.00025745*m, 0.00026161*m, 0.00026448*m, 0.00026916*m, 0.00027182*m, 0.00027703*m, 0.00028022*m, 0.00028268*m, 0.00028588*m, 0.00028629*m, 0.00028914*m, 0.00028876*m, 0.00028952*m, 0.00028787*m, 0.00028689*m, 0.00028658*m, 0.00028583*m, 0.00028717*m, 0.00028956*m, 0.00029622*m, 0.00030215*m, 0.00031247*m, 0.0003225*m, 0.00033448*m, 0.00034183*m, 0.00035359*m, 0.00036412*m, 0.00037562*m, 0.00038645*m, 0.00039638*m, 0.00040486*m, 0.0004197*m, 0.00042926*m, 0.00044327*m, 0.00045097*m, 0.00046763*m, 0.00048205*m, 0.00049883*m, 0.00052405*m, 0.00054125*m, 0.00056765*m, 0.00059224*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    concentration = 160.0;
    absLength[concentration] = {noAbsLength_, noAbsLength_, 1.861*m, 2.1025*m, 1.9415*m, 1.9857*m, 2.0913*m, 1.9486*m, 1.7941*m, 1.9602*m, 1.9239*m, 1.7359*m, 1.7313*m, 1.6664*m, 1.7753*m, 1.8757*m, 1.8405*m, 1.7746*m, 1.8152*m, 2.0605*m, 1.8268*m, 1.7146*m, 1.7578*m, 2.0933*m, 1.6596*m, 1.7206*m, 1.9017*m, 1.744*m, 1.7921*m, 1.7083*m, 2.0013*m, 1.5393*m, 1.65*m, 1.5799*m, 1.6015*m, 1.5024*m, 1.5053*m, 1.5968*m, 1.5165*m, 1.634*m, 1.5546*m, 1.5162*m, 1.593*m, 1.6133*m, 1.6224*m, 1.544*m, 1.4611*m, 1.4666*m, 1.5088*m, 1.4954*m, 1.4081*m, 1.4874*m, 1.4216*m, 1.3695*m, 1.4434*m, 1.3177*m, 1.3687*m, 1.3594*m, 1.2821*m, 1.2877*m, 1.3831*m, 1.263*m, 1.2769*m, 1.2313*m, 1.3245*m, 1.2945*m, 1.2548*m, 1.2059*m, 1.2838*m, 1.2808*m, 1.2153*m, 1.2567*m, 1.2017*m, 1.2583*m, 1.1666*m, 1.2066*m, 1.2154*m, 1.1764*m, 1.1829*m, 1.1408*m, 1.2016*m, 1.0959*m, 1.1724*m, 1.1092*m, 1.0979*m, 1.1337*m, 1.0592*m, 1.0952*m, 1.0908*m, 1.0727*m, 1.0532*m, 1.1318*m, 1.0798*m, 1.0167*m, 1.0411*m, 1.0724*m, 1.0082*m, 0.9806*m, 1.0533*m, 0.9515*m, 0.9801*m, 0.9636*m, 0.9801*m, 0.9476*m, 0.9111*m, 0.9529*m, 0.9398*m, 0.8993*m, 0.8706*m, 0.8983*m, 0.9041*m, 0.9036*m, 0.8773*m, 0.8616*m, 0.8719*m, 0.8578*m, 0.8379*m, 0.8778*m, 0.8583*m, 0.8327*m, 0.8399*m, 0.8147*m, 0.8207*m, 0.7909*m, 0.8016*m, 0.7929*m, 0.7969*m, 0.7807*m, 0.7899*m, 0.7681*m, 0.7489*m, 0.7368*m, 0.7487*m, 0.7491*m, 0.7038*m, 0.6853*m, 0.708*m, 0.6901*m, 0.6949*m, 0.6778*m, 0.6688*m, 0.6579*m, 0.6789*m, 0.6456*m, 0.6339*m, 0.6382*m, 0.6084*m, 0.607*m, 0.5981*m, 0.5956*m, 0.5865*m, 0.5737*m, 0.5415*m, 0.5372*m, 0.5328*m, 0.5212*m, 0.5156*m, 0.4967*m, 0.4772*m, 0.4636*m, 0.4492*m, 0.4297*m, 0.4086*m, 0.3896*m, 0.3679*m, 0.3381*m, 0.3096*m, 0.2884*m, 0.2612*m, 0.2333*m, 0.2111*m, 0.1885*m, 2.6105*m, 2.6169*m, 2.4706*m, 2.4975*m, 2.2821*m, 2.427*m, 2.3882*m, 2.198*m, 2.1299*m, 2.0846*m, 2.1567*m, 2.0835*m, 2.0423*m, 1.9375*m, 1.96*m, 1.963*m, 1.8889*m, 1.8522*m, 1.7337*m, 1.8264*m, 1.7205*m, 1.7155*m, 1.627*m, 1.5588*m, 1.5526*m, 1.5337*m, 1.4664*m, 1.4388*m, 1.3298*m, 1.3794*m, 1.2917*m, 1.2517*m, 1.2382*m, 1.1879*m, 1.1149*m, 1.0593*m, 1.0037*m, 0.9401*m, 0.8733*m, 0.8266*m, 0.7783*m, 0.7448*m, 0.6882*m, 0.6273*m, 0.5982*m, 0.5628*m, 0.5343*m, 0.4993*m, 0.471*m, 0.4464*m, 0.4263*m, 0.407*m, 0.3878*m, 0.3689*m, 0.3487*m, 0.3354*m, 0.3213*m, 0.303*m, 0.2915*m, 0.2756*m, 0.2617*m, 0.2464*m, 0.2359*m, 0.2242*m, 0.2155*m, 0.2043*m, 0.1952*m, 0.1885*m, 0.1803*m, 0.1732*m, 0.1644*m, 0.1596*m, 0.1549*m, 0.1488*m, 0.1429*m, 0.1376*m, 0.1332*m, 0.1288*m, 0.1239*m, 0.1204*m, 0.1159*m, 0.1128*m, 0.1089*m, 0.1058*m, 0.1031*m, 0.1001*m, 0.0971*m, 0.094*m, 0.0909*m, 0.0883*m, 0.0849*m, 0.0825*m, 0.08*m, 0.0779*m, 0.0759*m, 0.0739*m, 0.0722*m, 0.0703*m, 0.0687*m, 0.0674*m, 0.0657*m, 0.0642*m, 0.063*m, 0.0617*m, 0.0604*m, 0.0593*m, 0.0582*m, 0.057*m, 0.0561*m, 0.0553*m, 0.0544*m, 0.0537*m, 0.0529*m, 0.0521*m, 0.0513*m, 0.0504*m, 0.0494*m, 0.0485*m, 0.0477*m, 0.0471*m, 0.0462*m, 0.0454*m, 0.0447*m, 0.0438*m, 0.043*m, 0.042*m, 0.041*m, 0.0402*m, 0.039*m, 0.0379*m, 0.0366*m, 0.0349*m, 0.0334*m, 0.0311*m, 0.0289*m, 0.0267*m, 0.0236*m, 0.0212*m, 0.0189*m, 0.0163*m, 0.0142*m, 0.0121*m, 0.0105*m, 0.0089402*m, 0.0073743*m, 0.0063534*m, 0.0053878*m, 0.0046087*m, 0.0039625*m, 0.0033818*m, noAbsLength_, noAbsLength_};
    WLS_absLength[concentration] = {noAbsLength_, noAbsLength_, 0.1661*m, 0.1465*m, 0.1286*m, 0.1114*m, 0.0937*m, 0.079*m, 0.065*m, 0.0514*m, 0.0414*m, 0.033*m, 0.0258*m, 0.02*m, 0.0147*m, 0.0111*m, 0.0083834*m, 0.0064025*m, 0.0047881*m, 0.0034768*m, 0.002651*m, 0.0020991*m, 0.0016441*m, 0.0013073*m, 0.0010592*m, 0.0009177*m, 0.00083603*m, 0.00066686*m, 0.00058575*m, 0.00052449*m, 0.00048445*m, 0.00044858*m, 0.00042495*m, 0.00041001*m, 0.00040103*m, 0.00039693*m, 0.00039422*m, 0.00039442*m, 0.00039406*m, 0.00039389*m, 0.00039237*m, 0.00038975*m, 0.00038505*m, 0.00037906*m, 0.00037103*m, 0.00036237*m, 0.00035169*m, 0.00033795*m, 0.00032195*m, 0.00030592*m, 0.00029212*m, 0.0002776*m, 0.00026601*m, 0.00025422*m, 0.00024909*m, 0.00024465*m, 0.00024239*m, 0.00024321*m, 0.00024458*m, 0.00024853*m, 0.00025126*m, 0.0002557*m, 0.00025823*m, 0.00026318*m, 0.00026621*m, 0.00026855*m, 0.00027159*m, 0.00027198*m, 0.00027468*m, 0.00027433*m, 0.00027504*m, 0.00027348*m, 0.00027255*m, 0.00027226*m, 0.00027154*m, 0.00027281*m, 0.00027508*m, 0.00028141*m, 0.00028705*m, 0.00029685*m, 0.00030637*m, 0.00031776*m, 0.00032474*m, 0.00033591*m, 0.00034591*m, 0.00035684*m, 0.00036713*m, 0.00037656*m, 0.00038462*m, 0.00039871*m, 0.0004078*m, 0.00042111*m, 0.00042842*m, 0.00044425*m, 0.00045795*m, 0.00047389*m, 0.00049784*m, 0.00051419*m, 0.00053926*m, 0.00056263*m, noAbsLength_, noAbsLength_};
    available_concentrations.push_back(concentration);

    G4double cromophore_concentration_ = FindClosestNumber(cromophore_concentration, &available_concentrations);
    if(verbosity){
      G4cout << "Taking the absorption length spectra for " << cromophore_concentration_ << " mg/kg." << G4endl;
    } 

    mpt->AddProperty("ABSLENGTH", abs_energy, absLength[cromophore_concentration_]);
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength[cromophore_concentration_]);

    //std::vector<G4double> absLength = {
    //  attenuation_length,  attenuation_length,
    //  attenuation_length,  attenuation_length,
    //  attenuation_length,  attenuation_length
    //};

    if(verbosity){
      G4cout << "Attenuation length energy entries = " << abs_energy.size() << G4endl;
      G4cout << "Attenuation length entries = " << absLength[cromophore_concentration_].size() << G4endl;
      G4cout << "WLS absorption length energy entries = " << WLS_abs_energy.size() << G4endl;
      G4cout << "WLS absorption length entries = " << WLS_absLength[cromophore_concentration_].size() << G4endl;
    } 

    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_, h_Planck*c_light/(657.0*nm), h_Planck*c_light/(656.0*nm), 
                    h_Planck*c_light/(654.44*nm), h_Planck*c_light/(634.1*nm), 
                    h_Planck*c_light/(602.19*nm), h_Planck*c_light/(577.45*nm), 
                    h_Planck*c_light/(557.87*nm), h_Planck*c_light/(545.87*nm), 
                    h_Planck*c_light/(536.27*nm), h_Planck*c_light/(527.86*nm), 
                    h_Planck*c_light/(522.23*nm), h_Planck*c_light/(519.41*nm), 
                    h_Planck*c_light/(517.0*nm), h_Planck*c_light/(512.2*nm), 
                    h_Planck*c_light/(506.61*nm), h_Planck*c_light/(501.42*nm), 
                    h_Planck*c_light/(498.61*nm), h_Planck*c_light/(495.39*nm), 
                    h_Planck*c_light/(492.96*nm), h_Planck*c_light/(491.34*nm), 
                    h_Planck*c_light/(490.13*nm), h_Planck*c_light/(488.52*nm), 
                    h_Planck*c_light/(485.73*nm), h_Planck*c_light/(484.14*nm), 
                    h_Planck*c_light/(482.57*nm), h_Planck*c_light/(480.21*nm), 
                    h_Planck*c_light/(477.08*nm), h_Planck*c_light/(473.54*nm), 
                    h_Planck*c_light/(469.61*nm), h_Planck*c_light/(466.48*nm), 
                    h_Planck*c_light/(464.54*nm), h_Planck*c_light/(461.81*nm), 
                    h_Planck*c_light/(459.45*nm), h_Planck*c_light/(456.67*nm), 
                    h_Planck*c_light/(449.5*nm), h_Planck*c_light/(435.54*nm), 
                    h_Planck*c_light/(434.0*nm), h_Planck*c_light/(433.0*nm),
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

  std::vector<G4double> abs_energy =  {optPhotMinE_, 2.7588*eV, 2.7733*eV, 4.7943*eV, 4.8604*eV, optPhotMaxE_};
  std::vector<G4double> absLength =   {1.*m, 1.*m, 1.*m, 1.*m, 1.*m, 1.*m};

  // The following absorption length is the result of a defective computation which I have to fix. Until fixing that computation
  // I am going to assume a long absorption length for this material  
  //  std::vector<G4double> abs_energy =  {optPhotMinE_, 2.7588*eV, 2.7733*eV, 2.7879*eV, 2.8027*eV, 2.8227*eV, 2.8407*eV, 2.8532*eV, 
  //                                      2.8687*eV, 2.8818*eV, 2.8949*eV, 2.9149*eV, 2.9257*eV, 2.9366*eV, 2.9488*eV, 2.9641*eV, 
  //                                      2.9822*eV, 2.9949*eV, 3.0068*eV, 3.0175*eV, 3.0282*eV, 3.0423*eV, 3.0614*eV, 3.0763*eV, 
  //                                      3.0928*eV, 3.111*eV, 3.1282*eV, 3.1398*eV, 3.1514*eV, 3.1631*eV, 3.1748*eV, 3.1974*eV, 
  //                                      3.2199*eV, 3.2449*eV, 3.2671*eV, 3.2787*eV, 3.3024*eV, 3.3232*eV, 3.3451*eV, 3.3709*eV, 
  //                                      3.3895*eV, 3.4163*eV, 3.442*eV, 3.4641*eV, 3.4818*eV, 3.5056*eV, 3.5165*eV, 3.5393*eV, 
  //                                      3.5568*eV, 3.5747*eV, 3.5944*eV, 3.6118*eV, 3.6314*eV, 3.6495*eV, 3.6612*eV, 3.6782*eV, 
  //                                      3.6997*eV, 3.7171*eV, 3.7345*eV, 3.7485*eV, 3.7637*eV, 3.7804*eV, 3.8014*eV, 3.8142*eV, 
  //                                      3.8328*eV, 3.8472*eV, 3.8574*eV, 3.866*eV, 3.8852*eV, 3.8888*eV, 3.9015*eV, 3.9104*eV, 
  //                                      3.9169*eV, 3.9267*eV, 3.9389*eV, 3.9497*eV, 3.9634*eV, 3.9742*eV, 3.9851*eV, 3.9945*eV, 
  //                                      4.0007*eV, 4.0085*eV, 4.018*eV, 4.029*eV, 4.0402*eV, 4.0563*eV, 4.0679*eV, 4.0756*eV, 
  //                                      4.0886*eV, 4.0993*eV, 4.1133*eV, 4.1299*eV, 4.1433*eV, 4.1534*eV, 4.1635*eV, 4.1805*eV, 
  //                                      4.2011*eV, 4.222*eV, 4.2421*eV, 4.2655*eV, 4.3001*eV, 4.338*eV, 4.38*eV, 4.4166*eV, 
  //                                      4.4665*eV, 4.511*eV, 4.5629*eV, 4.6161*eV, 4.6774*eV, 4.7404*eV, 4.7943*eV, 4.8604*eV, optPhotMaxE_};
  //  std::vector<G4double> absLength =   {23.3403*mm, 23.3403*mm, 22.8518*mm, 22.8518*mm, 22.8518*mm, 23.3403*mm, 23.0327*mm, 22.8518*mm, 22.8518*mm, 
  //                                      22.8518*mm, 23.3403*mm, 22.8518*mm, 22.3851*mm, 23.3403*mm, 23.6499*mm, 23.6499*mm, 23.3403*mm, 23.3403*mm, 
  //                                      23.6499*mm, 23.6499*mm, 23.6499*mm, 22.8518*mm, 23.3403*mm, 23.3403*mm, 23.3403*mm, 23.3403*mm, 23.0327*mm, 
  //                                      23.0327*mm, 23.0327*mm, 23.3374*mm, 23.6499*mm, 23.3403*mm, 22.6792*mm, 22.8518*mm, 21.889*mm, 22.2598*mm, 
  //                                      21.9337*mm, 21.9337*mm, 21.5324*mm, 21.0836*mm, 20.9742*mm, 20.6806*mm, 20.4822*mm, 19.8994*mm, 19.4537*mm, 
  //                                      18.6188*mm, 18.6188*mm, 17.8475*mm, 17.136*mm, 16.4761*mm, 16.0121*mm, 15.1974*mm, 14.1325*mm, 13.3397*mm, 
  //                                      12.8225*mm, 12.0863*mm, 11.1549*mm, 10.4281*mm, 9.8144*mm, 9.205*mm, 8.3988*mm, 7.7092*mm, 7.1125*mm, 
  //                                      6.7153*mm, 6.1858*mm, 5.8705*mm, 5.4229*mm, 5.1275*mm, 4.6954*mm, 4.571*mm, 4.257*mm, 4.0262*mm, 
  //                                      3.8241*mm, 3.6253*mm, 3.354*mm, 3.1466*mm, 2.864*mm, 2.6756*mm, 2.5024*mm, 2.3898*mm, 2.2897*mm, 
  //                                      2.188*mm, 2.0444*mm, 1.9268*mm, 1.8156*mm, 1.6477*mm, 1.5535*mm, 1.4743*mm, 1.3886*mm, 1.3037*mm, 
  //                                      1.1916*mm, 1.1031*mm, 1.0222*mm, 0.9416*mm, 0.8991*mm, 0.8175*mm, 0.734*mm, 0.6535*mm, 0.5682*mm, 
  //                                      0.5059*mm, 0.4079*mm, 0.3696*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 
  //                                      0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm, 0.3401*mm};
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
      optPhotMinE_,
      hc_ / (750. * nm), hc_ / (740. * nm), hc_ / (380. * nm), hc_ / (370. * nm),
      optPhotMaxE_
    };
    std::vector<G4double> absLength = {
      noAbsLength_,
      noAbsLength_, 3.5 * m, 3.5 * m, noAbsLength_,
      noAbsLength_
    };
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,
      hc_ / (490. * nm), hc_ / (485. * nm), hc_ / (475. * nm), hc_ / (454. * nm),
      hc_ / (443. * nm), hc_ / (430. * nm), hc_ / (410. * nm), hc_ / (405. * nm),
      hc_ / (359. * nm), hc_ / (350. * nm), hc_ / (345. * nm),
      optPhotMaxE_
    };
    std::vector<G4double> WLS_absLength = {
      noAbsLength_,
      noAbsLength_, 44.2 * mm, 5.39 * mm, 0.395 * mm, // 490, 485, 475, 454 nm
      0.462 * mm, 0.354 * mm, 0.571 * mm, 0.612 * mm, // 443, 430, 410, 405 nm
      4.51 * mm,  4.81  * mm, noAbsLength_,           // 359, 350, 345  nm
      noAbsLength_
    };
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);
    //for (int i=0; i<WLS_abs_entries; i++)
    //  G4cout << "* Y11 WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (hc_ / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / mm << " mm" << G4endl;

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_,
      hc_ / (580. * nm), hc_ / (550. * nm), hc_ / (530. * nm),
      hc_ / (525. * nm), hc_ / (520. * nm), hc_ / (515. * nm),
      hc_ / (510. * nm), hc_ / (505. * nm), hc_ / (500. * nm),
      hc_ / (495. * nm), hc_ / (490. * nm), hc_ / (485. * nm),
      hc_ / (480. * nm), hc_ / (475. * nm), hc_ / (470. * nm),
      hc_ / (465. * nm), hc_ / (460. * nm), hc_ / (455. * nm),
      hc_ / (450. * nm), hc_ / (445. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.000,
      0.000, 0.200, 0.300, // 580, 550, 530 nm
      0.400, 0.600, 0.750, // 525, 520, 515 nm
      0.750, 0.720, 0.700, // 510, 505, 500 nm
      0.680, 0.650, 0.700, // 495, 490, 485 nm
      0.900, 1.000, 0.950, // 480, 475, 470 nm
      0.500, 0.300, 0.100, // 465, 460, 455 nm
      0.050, 0.000,        // 450, 445 nm
      0.000
    };
    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 8.5 * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.87);

    return mpt;
  }



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
      optPhotMinE_,
      hc_ / (750. * nm), hc_ / (740. * nm), hc_ / (380. * nm), hc_ / (370. * nm),
      optPhotMaxE_
    };
    std::vector<G4double> absLength = {
      noAbsLength_,
      noAbsLength_, 3.5 * m, 3.5 * m, noAbsLength_,
      noAbsLength_
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
      hc_ / (418. * nm), hc_ / (412. * nm), hc_ / (405. * nm), hc_ / (400. * nm),
      hc_ / (394. * nm), hc_ / (387. * nm), hc_ / (384. * nm), hc_ / (382. * nm),
      hc_ / (378. * nm), hc_ / (370. * nm), hc_ / (361. * nm), hc_ / (353. * nm),
      hc_ / (345. * nm), hc_ / (341. * nm), hc_ / (336. * nm), hc_ / (331. * nm),
      hc_ / (316. * nm), hc_ / (301. * nm), hc_ / (280. * nm),
      optPhotMaxE_
    };

    float minAbsLength = 0.395 * mm;

    std::vector<float> B2_absorption {
      -0.01, -0.06, -0.26, -0.44, // 418, 412, 405, 400
      -0.59, -0.59, -0.64, -0.77, // 394, 387, 384, 382
      -0.92, -1.00, -0.93, -0.85, // 378, 370, 361, 353
      -0.87, -0.87, -0.77, -0.56, // 345, 341, 336, 331
      -0.35, -0.22, -0.12         // 316, 301, 280
    };

    std::vector<G4double> WLS_absLength {noAbsLength_};

    for (auto &abs_value : B2_absorption)
      WLS_absLength.push_back(- minAbsLength / abs_value);

    WLS_absLength.push_back(noAbsLength_);

    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    std::vector<G4double> WLS_emi_energy = {
      optPhotMinE_,
      hc_ / (542 * nm), hc_ / (525 * nm), hc_ / (508 * nm), hc_ / (497 * nm),
      hc_ / (488 * nm), hc_ / (479 * nm), hc_ / (473 * nm), hc_ / (467 * nm),
      hc_ / (463 * nm), hc_ / (458 * nm), hc_ / (454 * nm), hc_ / (449 * nm),
      hc_ / (445 * nm), hc_ / (442 * nm), hc_ / (440 * nm), hc_ / (438 * nm),
      hc_ / (433 * nm), hc_ / (429 * nm), hc_ / (424 * nm), hc_ / (420 * nm),
      hc_ / (418 * nm), hc_ / (416 * nm), hc_ / (411 * nm), hc_ / (404 * nm),
      hc_ / (402 * nm), hc_ / (399 * nm), hc_ / (398 * nm), hc_ / (396 * nm),
      hc_ / (395 * nm), hc_ / (394 * nm), hc_ / (392 * nm), hc_ / (391 * nm),
      hc_ / (386 * nm), hc_ / (380 * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_emiSpectrum = {
      0.000,
      0.053, 0.070, 0.109, 0.143, // 542, 525, 508, 497
      0.199, 0.270, 0.337, 0.423, // 488, 479, 473, 467
      0.497, 0.582, 0.615, 0.645, // 463, 458, 454, 449
      0.679, 0.750, 0.801, 0.857, // 445, 442, 440, 438
      0.957, 0.999, 0.949, 0.906, // 433, 429, 424, 420
      0.855, 0.809, 0.750, 0.750, // 418, 416, 411, 404
      0.719, 0.671, 0.590, 0.500, // 402, 399, 398, 396
      0.421, 0.327, 0.217, 0.138, // 395, 394, 392, 391
      0.065, 0.023,               // 386, 380
      0.000
    };

    mpt->AddProperty("WLSCOMPONENT",  WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 8.5 * ns);

    // WLS Quantum Efficiency
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.87);

    return mpt;
  }



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

  // Copper Optical Properties Table
  G4MaterialPropertiesTable * Copper()
  {
       G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

      // Reflectivity
      std::vector<G4double> refl_energies = {
              optPhotMinE_, hc_ / (700. * nm), hc_ / (650. * nm),
                            hc_ / (600. * nm), hc_ / (500. * nm),
                            hc_ / (400. * nm), hc_ / (300. * nm),
                            hc_ / (170. * nm), optPhotMaxE_ };

      // We assume a reflectivity of the copper of 20% at VUV
      // Measurements may be required to update these values
      // Visible spectrum taken from: E. Hecht, Optics, 5th edn., Pg 142, Fig 4.69
      std::vector<G4double> reflectivities = { 0.90, 0.90, 0.70, 0.60, 0.40, 0.35, 0.30, 0.20, 0.20};

      // We assume that the reflectivity is mostly specular.
      // Measurements may be required to update these values
      // Add Properties
      mpt->AddProperty("SPECULARLOBECONSTANT", {optPhotMinE_, optPhotMaxE_}, {0., 0.});
      mpt->AddProperty("SPECULARSPIKECONSTANT",{optPhotMinE_, optPhotMaxE_}, {0.75, 0.75});
      mpt->AddProperty("BACKSCATTERCONSTANT",  {optPhotMinE_, optPhotMaxE_}, {0., 0.});
      mpt->AddProperty("REFLECTIVITY", refl_energies, reflectivities);
      return mpt;
  }

  // Stainles Steel Optical Properties Table
  G4MaterialPropertiesTable * Steel()
  {
      G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

      // Reflectivity
      std::vector<G4double> refl_energies = {
              optPhotMinE_, hc_ / (500. * nm), hc_ / (350. * nm),
                            hc_ / (300. * nm), hc_ / (170. * nm),  optPhotMaxE_ };

      // We assume a reflectivity of the stainless steel of 20% at VUV
      // Measurements may be required to update these values
      // Visible spectrum taken from: https://doi.org/10.1063/1.331503
      std::vector<G4double> reflectivities = { 0.60, 0.60, 0.50, 0.40, 0.20, 0.20};

      // We assume that the reflectivity is mostly specular.
      // Measurements may be required to update these values
      // Add properties
      mpt->AddProperty("SPECULARLOBECONSTANT", {optPhotMinE_, optPhotMaxE_}, {0., 0.});
      mpt->AddProperty("SPECULARSPIKECONSTANT",{optPhotMinE_, optPhotMaxE_}, {0.75, 0.75});
      mpt->AddProperty("BACKSCATTERCONSTANT",  {optPhotMinE_, optPhotMaxE_}, {0., 0.});
      mpt->AddProperty("REFLECTIVITY", refl_energies, reflectivities);
      return mpt;
  }

  /// Generic material, to be modifed by the user ///
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
