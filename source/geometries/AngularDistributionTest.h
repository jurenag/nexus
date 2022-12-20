#ifndef ANGULAR_DISTRIBUTION_TEST
#define ANGULAR_DISTRIBUTION_TEST

#include "GeometryBase.h"

#include <G4ThreeVector.hh>

class G4MaterialPropertiesTable;
class G4GenericMessenger;
class G4VPhysicalVolume;

namespace nexus {

  /// The aim of this geometry is to evaluate the angular distribution of the light,
  /// that is conveniently generated, after interacting with a certain test sample.
  /// The geometry consists of a sample and a light-detection plane (LDP), which is
  /// basically a 100% efficient thin disk. The sample, which is arranged parallel 
  /// to the LDP, is placed on top of the LDP. The Construct() method calls
  /// ConstructSample() and ConstructLDP(), which take care of the construction of 
  /// the sample and the LDP respectively. For each different sample you want to 
  /// test, you might write one method. Say, for a sample called X, you might write
  /// the method ConstructX(). Depending on which sample you want to test, you might
  /// call one method or another from ConstructSample(). 

class AngularDistributionTest: public GeometryBase
  {
  public:
    ///Default constructor
    AngularDistributionTest();
    ///Destructor
    ~AngularDistributionTest();

    void Construct();
    void ConstructSample(G4VPhysicalVolume*);
    void ConstructWorldTightDichroicFilter(G4VPhysicalVolume*);
    void ConstructFramedDichroicFilter(G4VPhysicalVolume*);
    void ConstructFramedPTPAlone(G4VPhysicalVolume*);
    void ConstructLDP(G4VPhysicalVolume*);
    G4bool geometry_is_ill_formed();
    G4ThreeVector GenerateVertex(const G4String&) const;

  private:
    G4int config_code_;                         /// < 1 - World-tight dichroic filter.
                                                ///   2 - Framed dichroic filter.
                                                ///   3 - Same as config_code_==2, but 
                                                ///       removing the DF physical volume
                                                ///       and leaving the PTP floating.
    G4int filter_code_;                         /// < 1 - Use perfectly transparent filter.
                                                ///   2 - Use OPTO filter (substrate MPT is set to
                                                ///       that of OPTO filter, i.e. SCHOTT B270)
    G4double world_radius_, world_height_;      /// < The world is a cylinder with these dimensions.
    G4double sample_LDP_distance_;              /// < The distance between the sample and the LDP.
    G4double lgv_zpos_wrt_world_center_;        /// < Light-generation-vertex z coordinate wrt the center of the world.
    G4MaterialPropertiesTable* sm_mpt_;         /// < The optical properties of the media which fills the cylinder.

    // config_code_==1 AND 2 AND 3 parameters
    G4double DF_thickn_;
    G4double ptp_thickn_; 

    // config_code_==2 AND 3 parameters
    G4double DF_x_size_, DF_y_size_;  /// < DF dimensions
    G4double rib_height_;             /// < Rib height (wrt to the rib that surround the filter)
    G4double shallowness_;            /// < Shallowness parameter of DF within rib.
                                      /// < It ranges from 0. to 1., 0. meaning that the DF is the 
                                      /// < farthest away from the light source it can be, while 
                                      /// < still being fully contained within the rib.
    G4bool randomly_sample_lgv_;      /// < Whether the generation vertex is randomly sampled over
                                      ///   the coating filter surface.

    G4GenericMessenger* msg_;
  };

}

#endif

//    The world is a cylinder whose height is equal to world_height_
//    and its diameter is equal to 2*world_radius_. The LDP is a cylinder
//    that fits tightly at the bottom of the world. Thus, LDP radius is 
//    that of the world radius, and the LDP height is set to be 1/1000
//    of the world_height
//    _______________________________
//    |                             |    ^
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    | world_height_
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |                             |    |
//    |_____________________________|    v
//    |             LDP             |    ^
//    |_____________________________|    v  // < The LDP height is 1/1000 of the world height
//
//    <----------------------------->
//             2*world_radius_

