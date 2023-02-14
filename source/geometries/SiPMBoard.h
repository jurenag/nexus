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

    inline void SetBaseID(G4int input)              { base_id_                  = input; }
    inline void SetNumPhsensors(G4int input)        { num_phsensors_            = input; }
    inline void SetSiPMCode(G4int input)            { SiPM_code_                = input; }
    inline void SetBoardLength(G4double input)      { board_length_             = input; }
    inline void SetBoardHeight(G4double input)      { board_height_             = input; }
    inline void SetBoardThickness(G4double input)   { board_thickn_             = input; }
    inline void SetReflectiveSupports(G4bool input) { ref_phsensors_supports_   = input; }

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

    G4bool GeometryIsIllFormed() const;     ///< Whether the provided parameters describe a feasible geometry

    void Construct();
    void ConstructBoard(G4LogicalVolume*);
    void ConstructSiPMs(G4LogicalVolume*);

  private:

    G4int base_id_;                     ///< ID of the first SiPM that is mounted. The i-th SiPM after the first 
                                        ///< one has an ID equal to base_id_+i
    G4int num_phsensors_;               ///< Number of SiPMs that are mounted on the board
    G4int SiPM_code_;                   ///< Integer signalling which SiPM to mount on the board
                                        ///< 1                  -> Hamamatsu S13360-6050VE
                                        ///< 2                  -> Hamamatsu S13360-5075HD-HQR
                                        ///< Any other integer  -> FBK-NUV-HD-CRYO-TT
    G4double board_thickn_;             ///< Thickness of the mounting board
    G4double board_length_;             ///< Board length
    G4double board_height_;             ///< Board height (width)
    G4bool   ref_phsensors_supports_;   ///< Whether the FR4 supports of the SiPMs are vikuiti-coated    
  };

} // namespace nexus

#endif