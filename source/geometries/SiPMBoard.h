#ifndef SIPM_BOARD
#define SIPM_BOARD

#include "GeometryBase.h"

namespace nexus {

  class SiPMBoard: public GeometryBase
  {
  public:
    // Tries to imitate the concatenation of four SiPM-mounting-boards for FD1 supercells 
    // See DUNE FD SP TDR vol. IV, Figure 5.9 of section 5.4.:Light Collectors
    SiPMBoard();

    // Destructor
    ~SiPMBoard();

    inline void SetBaseID(G4int input)                      { base_id_                  = input; }
    inline void SetNumPhsensors(G4int input)                { num_phsensors_            = input; }
    inline void SetSiPMCode(G4int input)                    { SiPM_code_                = input; }
    inline void SetBoardLength(G4double input)              { board_length_             = input; }
    inline void SetBoardHeight(G4double input)              { board_height_             = input; }
    inline void SetBoardThickness(G4double input)           { board_thickn_             = input; }
    inline void SetReflectiveSupports(G4bool input)         { ref_phsensors_supports_   = input; }
    inline void SetReflectivityScaleFactor(G4double input)  { reflectivity_scale_factor_ = input; }
    inline void SetAddBlocks(G4bool input)                  { add_blocks_between_sipms_ = input; }
    inline void SetSiPMProtrusion(G4double input)           { sipm_protrusion_          = input; }

    G4int GetNumPhsensors()         const;  ///< Returns the number of photosensors attached to the board 
    G4int GetSiPMCode()             const;  ///< Returns the SiPM code
    G4double GetBoardLength()       const;  ///< Returns the length of the board
    G4double GetBoardHeight()       const;  ///< Returns the height of the board
    G4double GetOverallHeight()     const;  ///< Returns the span of the geometry along the height axis 
                                            ///< (i.e. the sipm height if the sipm heighth is greater 
                                            ///< than the board height, and the board height otherwise)
    G4double GetBoardThickness()    const;  ///< Returns the thickness of the board alone
    G4double GetOverallThickness()  const;  ///< Returns the span of the geometry along the thickness axis 
                                            ///< (i.e. the sum of the thickness of the board plus the SiPM thickness)
    G4double GetHasBlocks()         const;  ///< Returns the add_blocks_between_sipms_ attribute
    G4double GetSiPMProtrusion()    const;  ///< Returns the sipm_protrusion_ attribute

    G4bool GeometryIsIllFormed() const;     ///< Whether the provided parameters describe a feasible geometry

    void Construct();
    void ConstructBoard(G4LogicalVolume*);  ///< Constructs the board where the SiPMs are placed. By default, this 
                                            ///< board is wrapped by a G4LogicalSkinSurface with 
                                            ///< opticalprops::Vikuiti() G4MaterialPropertiesTable.
    void ConstructSiPMs(G4LogicalVolume*);

  private:

    G4int base_id_;                     ///< ID of the first SiPM that is mounted. The i-th SiPM after the first 
                                        ///< one has an ID equal to base_id_+i
    G4int num_phsensors_;               ///< Number of SiPMs that are mounted on the board
    G4int SiPM_code_;                   ///< Integer signalling which SiPM to mount on the board
                                        ///< 1                  -> Hamamatsu S13360-6050VE
                                        ///< 2                  -> Hamamatsu S13360-5075HD-HQR
                                        ///< 3                  -> FBK-NUV-HD-CRYO-TT
                                        ///< 4                  -> Broadcom AFBR-S4N44P044M
                                        ///< Any other integer  -> PerfectSiPMMPPC (100% efficiency)
    G4double board_thickn_;             ///< Thickness of the mounting board
    G4double board_length_;             ///< Board length
    G4double board_height_;             ///< Board height (width)
    G4bool   ref_phsensors_supports_;   ///< Whether the FR4 supports of the SiPMs are vikuiti-coated    
    G4double reflectivity_scale_factor_;///< Scale factor for the vikuiti reflectivity curve which is used for both,
                                        ///< the SiPMs supports and the board. It must belong to the [0., 1.] range.
    G4bool   add_blocks_between_sipms_; ///< Whether to add blocks filling the space in between each two
                                        ///< two sipms. Blocks are also added at both ends of the board.
                                        ///< Analogously to no-blocks case, when the blocks are placed
                                        ///< they are wrapped by a reflective vikuiti-like G4LogicalSkinSurface.
    G4double sipm_protrusion_;          ///< This parameter only makes a difference if add_blocks_between_sipms_ 
                                        ///< is True. In that case, this parameter is a length which must have 
                                        ///< a value in between 0.0 and the sipm thickness. It gives the distance 
                                        ///< which the SiPM protrudes from the surrounding blocks. The case where
                                        ///< sipm_protrusion_ is 0.0 means that the blocks are coplanar with the
                                        ///< SiPMs surface. If it matches the sipm thickness, no blocks are placed.
  };

} // namespace nexus

#endif