#ifndef X_ARAPUCA_RODS
#define X_ARAPUCA_RODS

#include "GeometryBase.h"

class G4VPhysicalVolume;
class G4MaterialPropertiesTable;
class G4GenericMessenger;
class G4UserLimits;

class SiPMMPPC;

namespace nexus {

  class SiPMMPPC; // Declaring it outside of the namespace environment gives a compile-time error

  /// X-Arapuca but replacing the WLS plate with WLS rods

  class XArapucaRods: public GeometryBase
  {

  public:
    ///Constructor
    XArapucaRods();
    ///Destructor
    ~XArapucaRods();

    void Construct();                   ///< Constructs the geometry

  private:

    void ConstructRods(G4VPhysicalVolume*) const;                   ///< Called by Construct(). 
    void ConstructBoards(G4VPhysicalVolume*) const;                 ///< Called by Construct(). Adds the board-mounted SiPMs.
    void ConstructReflectiveCase(G4VPhysicalVolume*) const;         ///< Called by Construct(). Adds the reflective case that encloses the XArapucaRods geometry.
    void ConstructDichroicAssemblies(G4VPhysicalVolume*) const;     ///< Called by Construct(). Adds the dichroic filters plane, i.e. dichroic filters embedded in a frame
                                                                    ///< If double_sided_, the DFAs are placed on both sides of the XArapucaRods. If not, one reflective cover is added 
                                                                    ///< on the dead side.
    void ConstructReflectiveSeparators(G4VPhysicalVolume*) const;   ///< Called by Construct().
    
    G4ThreeVector GenerateVertex(const G4String&) const;
    G4bool geometry_is_ill_formed();                                    ///<Check whether the WLS plate fits within the XArapuca cavity and the given number of photosensors fit along the XArapuca long side
    SiPMMPPC * ChooseSiPM(G4int) const;                                 ///< Returns a pointer to a new object of a class which inherits from the SiPMMPPC abstract class.
                                                                        ///< Such class depends on the integer value given to ChooseSiPM according to the following enumeration
                                                                        ///< 1: HamamatsuS133606050VE
                                                                        ///< 2: HamamatsuS133605075HQR
                                                                        ///< 3: FbkNuvHdCryoTT
                                                                        ///< 4: BroadcomAFBRS4N44P044M
                                                                        ///< any other integer: PerfectSiPMMPPC

  private:

    ///---- General attributes ----///
    G4String surrounding_media_;                                    ///< Which media to place the XArapucaRods in
                                                                    ///< 'lar'  - The XArapucaRods is placed in Liquid ARgon
                                                                    ///< 'gar'  - The XArapucaRods is placed in Gaseous ARgon
                                                                    ///< 'air'  - The XArapucaRods is placed in air - ¡Note that the implemented air has no bulk-absorption length yet!
                                                                    ///< Default behaviour is that of surrounding_media_=='lar'.
    G4double DFA_thickn_;                                           ///< Frame thickness 
    G4double DF_thickn_;                                            ///< Overall thickness of the dichroic filters (MLS+substrate). Must be smaller than the frame thickness.
    G4double DF_substrate_thickn_;                                  ///< Thickness of the dichroic filters substrate. Must be smaller than DF_thickn_.
    G4MaterialPropertiesTable* DF_substrate_mpt_;                   ///< Material Properties Table of the dichroic filters substrate.
    G4double DF_pos_wrt_DFA_pos_;                                   ///< Position (height) of the dichroic filters with respect to the DFA position (height). This parameter can take values from 0 to 1. See * below
    G4bool DF_are_coated_;                                          ///< Whether the filters are coated with PTP or not
    G4double coating_thickn_;                                       ///< Thickness of the coating layer that is deposited over the dichroic filters
    G4double coating_rindex_;                                       ///< Refractive index of the coating layer
    G4double outter_frame_width_along_wlsplength_;                  ///< Outter frame dimension along the length of the X-ARAPUCA
    G4double outter_frame_width_along_wlspwidth_;                   ///< Outter frame dimension along the width of the X-ARAPUCA
    G4double inner_frames_width_along_wlsplength_;                  ///< Inner frames dimension along the length of the X-ARAPUCA
    G4double inner_frames_width_along_wlspwidth_;                   ///< Inner frames dimension along the width of the X-ARAPUCA
    G4int df_no_along_wlsplength_;                                  ///< Number of dichroic filters along the length of the X-ARAPUCA
    G4int df_no_along_wlspwidth_;                                   ///< Number of dichroic filters along the width of the X-ARAPUCA
    G4int DFA_config_code_;                                         ///< Value in (1,2) which labels which DFA configuration will be simulated
                                                                    ///< 1 -    regular DFA whose ribs have a rectangular cross sections with dimensions (DFA_thickn_ x inner_frames_width_along_x_)
                                                                    ///< 2 -    DFA whose ribs have a triangular cross sections. The dimension of the height (resp. base) of such triangle is DFA_thickn_ (resp. inner_frames_width_along_x_)
    G4bool DFA_frame_is_reflective_;                                ///< Whether the FR4 DFA frame is vikuiti-coated or not
    G4bool DFA_frame_is_specular_;                                  ///< Whether the vikuiti coating of the DFA frame is specular-spikely reflective or diffusively reflective. Only makes a difference if DFA_frame_is_reflective_==True.
    G4bool remove_DFs_;                                             ///< Whether to remove the dichroic filters or not
    G4bool remove_DFA_frame_;                                       ///< Whether to remove the dichroic filters assembly frame
    G4double secondary_wls_attlength_;                              ///< Attenuation length of the secondary WLShifter, in case EJ286 is used. Attenuation length for the WLS rods.
    G4double cromophore_concentration_;                             ///< Cromophore concentration (in miligrams of cromophore per kilogram of PMMA) of the secondary WLShifter, in case G2P_FB118 is used. For config_code_==1 (resp. 2), this is the cromophore concentration for the WLS plate (WLS fibers).
    G4double case_thickn_;                                          ///< Thickness of the reflective case that encloses the XArapucaRods geometry
    G4int SiPM_code_;                                               ///< Integer signalling which SiPM to construct
                                                                    ///< 1                  -> Hamamatsu S13360-6050VE
                                                                    ///< 2                  -> Hamamatsu S13360-5075HD-HQR
                                                                    ///< 3                  -> FBK-NUV-HD-CRYO-TT
                                                                    ///< 4                  -> Broadcom AFBR-S4N44P044M (2x2 SiPM array)
                                                                    ///< Any other integer  -> PerfectSiPMMPPC (100% efficiency)
    G4int num_phsensors_;                                           ///< The number of photosensors per side which is perpendicular to the rods
    G4double rod_sipm_gap_;                                         ///< Gap between photosensors and the secondary-WLS material. A negative gap can help modelate the immersion of the SiPMs into the flat dimples. Be careful not to collide the SiPMs into the secondary-WLS material.
    G4bool ref_phsensors_supports_;                                 ///< Whether photosensors supports are reflective (the FR4 box that supports the SiPM, regardless PS_config_code_)
    G4bool double_sided_;                                           ///< Whether there are dichroic filters on both sides of the secondary-WLS material
    G4String generation_region_;                                    ///< Where to place the generation vertex (GV).
                                                                    ///< 'random'   - The GV is randomly sampled over the DFA (including the frame)
                                                                    ///< 'dichroic' - The GV is randomly sampled over the DFs (not including the frame)
                                                                    ///< 'custom'   - The GV is randomly sampled within a circle centered in (gen_x_, -, gen_z_)
                                                                    ///<              whose diameter is equal to gen_diameter_.
                                                                    ///< Default behaviour is that of generation_region_=='random'.
    G4double gen_x_, gen_z_;                                        ///< Average GV coordinates. It is only used if generation_region_=='custom' is True.
    G4double gen_diameter_;                                         ///< Diameter of the circle where the GV could be randomly sampled. It is only used if generation_region_=='custom' is True.
    G4String path_to_inwards_dichroic_data_;                        ///< Absolute path to dichroic data file that is to be sampled for the light trying to enter the XA cavity
    G4String path_to_outwards_dichroic_data_;                       ///< Absolute path to dichroic data file that is to be sampled for the light trying to escape the XA cavity
    G4double world_extra_thickn_;                                   ///< Extra thickness for the surrounding box world to wrap the XArapuca
    G4double tunneling_probability_;                                ///< Probability that a photon tunnels through the surface of the WLSPlate (inwards or outwards), meaning that it is 
                                                                    ///< straightfoward-ly transmitted (no deflection due to Frensel refraction). The goal of this parameter is to let 
                                                                    ///< us model surface imperfections of the WLS plate.
    G4bool add_blocks_between_sipms_;                               ///< Whether to add reflective blocks filling the gaps in between the board-mounted SiPMs. This option is only available 
                                                                    ///< if with_boards_ is true. These blocks are nearly coplanar with the SiPM surface: Their surface is 0.1 millimeters
                                                                    ///< behind the SiPMs surface - i.e. the SiPMs protrude 0.1 mm from these blocks. There are two reasons for this offset:
                                                                    ///<    1)  These blocks are part of a SiPMBoard object, whose logical volume is wrapped with a G4LogicalSkinSurface.
                                                                    ///<        On the other hand, if tunneling_probability_ is different from 0.0, then the secondary-WLS material logical 
                                                                    ///<        volume will also be wrapped with a G4LogicalSkinSurface. If that's the case, and gap_ is 0.0, then two 
                                                                    ///<        G4LogicalSkinSurface's match in space, giving an undefined behaviour which is prone to bugs. P.e. last time 
                                                                    ///<        this surfaces-match was observed, photons arriving to that surface from within a WLSPlate were absorbed.
                                                                    ///<    2)  Setting the blocks to be perfectly coplanar with the SiPM surface is not realistic either way, so we are
                                                                    ///<        safe introducing this 0.1 mm gap in between the secondary-WLS material and the SiPMs surface.                                                                    
    G4double rod_length_;                                           ///< Length of each rod (X)
    G4double rod_height_;                                           ///< Height of each rod (Y)
    G4double rod_width_;                                            ///< Width of each rod (Z)
    G4int rods_no_;                                                 ///< Number of rods that are installed within the XArapucaRods
    G4bool add_reflective_separators_;                              ///< This parameter only makes a difference if rods_no_ is >=2. Whether to add reflective separators in between each two rods.
                                                                    ///< The width (along Y) of each separator matches rods_lateral_gap_/2.
    G4double rods_lateral_gap_;                                     ///< This parameter is used for
                                                                    ///<  1)  the gap (along Y) in between one extremal rod and the surrouding reflective case, and
                                                                    ///<  2)  if rods_no_>=2 is true, then it is also the gap (along Y) among any two adjacent rods. Note that
                                                                    ///<      if add_reflective_separators_ is True, then the reflective separators are inserted in these gaps.
    G4bool along_long_side_;                                        ///< Whether the rods are parallel to the long side of the XArapucaRods or not.
                                                                    ///< The two SiPMs board are always installed perpendicular to the rods

    G4bool cut_rods_;                                               ///< Whether to split up each rod into two pieces
    G4double cut_angle_;                                            ///< This parameter only makes a difference if cut_rods_ is True. Angle of the cut, with respect to the Z axis (the axis along which the rods are rod_width_ long)
    G4double cut_thickness_;                                        ///< This parameter only makes a difference if cut_rods_ is True. Thickness of the cut/crack that is carved from the rod.
    G4bool place_separator_at_the_cut_;                             ///< This parameter only makes a difference if cut_rods_ is true. Whether to insert a reflective separator in the rod cut. If that's the case,
                                                                    ///< then the material properties table of such foils is opticalprops::specularspikeVIKUITI() and its thickness match cut_thickness_/2.

    ///---- General internal attributes ----///                     ///< These attributes are internal. They must not be set by the user.
    G4double internal_length_, internal_width_, internal_thickn_;   ///< Internal dimensions of the reflector cavity
    G4double overall_length_, overall_width_, overall_thickn_;      ///< Overall XArapucaRods dimensions
    G4double DFA_length_, DFA_width_;                               ///< Dichroic filters assembly transverse dimensions.
    G4double DF_length_, DF_width_;                                 ///< Dichroic filter transverse dimensions. All of the filters of the DFA have the same dimensions.
    G4double residual_displacement_;                                ///< This parameter only makes a difference if cut_rods_ is true. It is meant to be a small fraction of cut_thickness_.
                                                                    ///< For more information, check the documentation wihtin XArapucaRods::ConstructRods.

    G4GenericMessenger* msg_;                                       ///< Messenger for the definition of control commands
    G4UserLimits* ul_;                                              ///< Useful to set a maximum track length in the plate
  };
}

#endif

/*
________ ___________________ ________
        |                   |  
        |       DF          | <--- DF_pos_wrt_DFA_pos_ = 1.  (Shallow)
        |___________________|
        |                   |
        |                   |
        |                   |       
        |                   |
  DFA   |                   |  DFA
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |___________________|
        |                   |
        |        DF         | <--- DF_pos_wrt_DFA_pos_ = 0. (Deepest position within DFA)
________|___________________|________


_____________________________________
                WLSPLATE
_____________________________________

________ ___________________ ________
        |                   |  
        |       DF          | <--- DF_pos_wrt_DFA_pos_ = 0. (Deepest position within DFA)
        |___________________|
        |                   |
        |                   |
        |                   |       
        |                   |
  DFA   |                   |  DFA
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |                   |
        |___________________|
        |                   |
        |        DF         | <--- DF_pos_wrt_DFA_pos_ = 1.  (Shallow)
________|___________________|________



*/
