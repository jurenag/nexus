#include "WLSPlate.h"

#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "FactoryBase.h"  
#include "LArSphere.h"

#include <G4VSolid.hh>
#include <G4MultiUnion.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4SubtractionSolid.hh>
#include <G4Transform3D.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Orb.hh>
#include <G4Para.hh>
#include <G4RotationMatrix.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4GenericMessenger.hh>

using namespace nexus;

REGISTER_CLASS(WLSPlate, GeometryBase)

namespace nexus{

  WLSPlate::WLSPlate( G4bool with_LAr, 
                      G4bool dimples_at_x_plus, 
                      G4bool dimples_at_x_minus,
                      G4bool dimples_at_z_plus, 
                      G4bool dimples_at_z_minus, 
                      G4String dimple_type,
                      G4int how_many_dimples, 
                      G4double flat_dimple_width, 
                      G4double flat_dimple_depth, 
                      G4double curvy_dimple_radius,
                      G4bool cut_plate,
                      G4double cut_angle,
                      G4double cut_thickness,
                      G4double tunneling_probability):
  GeometryBase(),
  dx_(487.*mm),
  dy_(3.5*mm),
  dz_(93.*mm), //<Default WLS plate dimensions taken from FD1 TDR vol. IX
  with_LAr_env_(with_LAr),
  dimples_at_x_plus_(dimples_at_x_plus),
  dimples_at_x_minus_(dimples_at_x_minus),
  dimples_at_z_plus_(dimples_at_z_plus),
  dimples_at_z_minus_(dimples_at_z_minus),
  dimple_type_(dimple_type),
  how_many_dimples_(how_many_dimples),
  flat_dimple_width_(flat_dimple_width),
  flat_dimple_depth_(flat_dimple_depth),
  curvy_dimple_radius_(curvy_dimple_radius),
  cut_plate_(cut_plate),
  cut_angle_(cut_angle),
  cut_thickness_(cut_thickness),
  tunneling_probability_(tunneling_probability),
  generation_y_pos_(0.0), // Generating the photons inside the WLSPlate
  generation_mode_("random"),
  wrap_with_collector_(true),
  mpt_(opticalprops::G2P_FB118( 16., 
                                1.502, 
                                true))
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/WLSPlate/",
				"Control commands of geometry WLSPlate.");
    
    G4GenericMessenger::Command& dx_cmd =
      msg_->DeclareProperty("depth", dx_,
			    "Depth of the WLS plate.");
    dx_cmd.SetUnitCategory("Length");
    dx_cmd.SetParameterName("depth", false);
    dx_cmd.SetRange("depth>0.");

    G4GenericMessenger::Command& dy_cmd =
      msg_->DeclareProperty("height", dy_,
			    "Height of the WLS plate.");
    dy_cmd.SetUnitCategory("Length");
    dy_cmd.SetParameterName("height", false);
    dy_cmd.SetRange("height>0.");

    G4GenericMessenger::Command& dz_cmd =
      msg_->DeclareProperty("width", dz_,
			    "Width of the WLS plate.");
    dz_cmd.SetUnitCategory("Length");
    dz_cmd.SetParameterName("width", false);
    dz_cmd.SetRange("width>0.");

    msg_->DeclareProperty("LAr", with_LAr_env_,
    "Build a LAr sphere surrounding the WLSPlate.");

  }

  WLSPlate::WLSPlate( G4double dx, 
                      G4double dy, 
                      G4double dz, 
                      G4bool with_LAr, 
                      G4bool dimples_at_x_plus, 
                      G4bool dimples_at_x_minus, 
                      G4bool dimples_at_z_plus, 
                      G4bool dimples_at_z_minus, 
                      G4String dimple_type, 
                      G4int how_many_dimples,
                      G4double flat_dimple_width, 
                      G4double flat_dimple_depth, 
                      G4double curvy_dimple_radius,
                      G4bool cut_plate,
                      G4double cut_angle,
                      G4double cut_thickness,
                      G4double tunneling_probability):
  GeometryBase(),
  dx_(dx),
  dy_(dy),
  dz_(dz),
  with_LAr_env_(with_LAr),
  dimples_at_x_plus_(dimples_at_x_plus),
  dimples_at_x_minus_(dimples_at_x_minus),
  dimples_at_z_plus_(dimples_at_z_plus),
  dimples_at_z_minus_(dimples_at_z_minus),
  dimple_type_(dimple_type),
  how_many_dimples_(how_many_dimples),
  flat_dimple_width_(flat_dimple_width),
  flat_dimple_depth_(flat_dimple_depth),
  curvy_dimple_radius_(curvy_dimple_radius),
  cut_plate_(cut_plate),
  cut_angle_(cut_angle),
  cut_thickness_(cut_thickness),
  tunneling_probability_(tunneling_probability),
  generation_y_pos_(1.*cm),
  generation_mode_("random"),
  wrap_with_collector_(false),
  mpt_(opticalprops::G2P_FB118( 16., 
                                1.502, 
                                true))
  {
  }

  WLSPlate::WLSPlate( G4double dx, 
                      G4double dy, 
                      G4double dz, 
                      G4MaterialPropertiesTable* mpt, 
                      G4bool with_LAr, 
                      G4bool dimples_at_x_plus,
                      G4bool dimples_at_x_minus,
                      G4bool dimples_at_z_plus,
                      G4bool dimples_at_z_minus,
                      G4String dimple_type, 
                      G4int how_many_dimples,
                      G4double flat_dimple_width, 
                      G4double flat_dimple_depth, 
                      G4double curvy_dimple_radius,
                      G4bool cut_plate,
                      G4double cut_angle,
                      G4double cut_thickness,
                      G4double tunneling_probability):
  GeometryBase(),
  dx_(dx),
  dy_(dy),
  dz_(dz),
  mpt_(mpt),
  with_LAr_env_(with_LAr),
  dimples_at_x_plus_(dimples_at_x_plus),
  dimples_at_x_minus_(dimples_at_x_minus),
  dimples_at_z_plus_(dimples_at_z_plus),
  dimples_at_z_minus_(dimples_at_z_minus),
  dimple_type_(dimple_type),
  how_many_dimples_(how_many_dimples),
  flat_dimple_width_(flat_dimple_width),
  flat_dimple_depth_(flat_dimple_depth),
  curvy_dimple_radius_(curvy_dimple_radius),
  cut_plate_(cut_plate),
  cut_angle_(cut_angle),
  cut_thickness_(cut_thickness),
  tunneling_probability_(tunneling_probability),
  generation_y_pos_(1.*cm),
  generation_mode_("random"),
  wrap_with_collector_(false)
  {
  }

  WLSPlate::~WLSPlate()
  {
    if(msg_) delete msg_;
  }

  void WLSPlate::Construct()
  {

    G4LogicalVolume* sphere_logic = nullptr;
    if(with_LAr_env_){
        G4String name = "LAR_SPHERE";

        // Define solid volume as a sphere
        G4Orb* sphere_solid = new G4Orb(name, dx_+dy_+dz_+(5.*m));

        // Define the material
        G4Material* lAr = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");

        // Set the optical material properties
        lAr->SetMaterialPropertiesTable(opticalprops::LAr());

        // Define the logical volume of the sphere using the material
        // and the solid volume defined above
        sphere_logic = new G4LogicalVolume(sphere_solid, lAr, name);
        sphere_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
        this->SetLogicalVolume(sphere_logic);
    }
    ConstructWLSPlate(sphere_logic);
    if(wrap_with_collector_) ConstructCollector(sphere_logic);
    return;

  }

  void WLSPlate::ConstructWLSPlate(G4LogicalVolume* world_logic_vol)
  {

    const G4String plate_name = "WLS_PLATE";

    G4VSolid* geometry_solid = nullptr;
    G4Box* whole_plate_solid = new G4Box(plate_name, dx_/2., dy_/2., dz_/2.);

    G4bool has_at_least_one_dimple;
    has_at_least_one_dimple = dimples_at_x_plus_ || dimples_at_x_minus_;
    has_at_least_one_dimple = has_at_least_one_dimple || dimples_at_z_plus_ || dimples_at_z_minus_;
    has_at_least_one_dimple = has_at_least_one_dimple && (how_many_dimples_>=1);

    if(has_at_least_one_dimple)
    {
        G4double tolerance = 0.5*mm; // To avoid boolean subtraction of matching surfaces
        G4VSolid* carving_solid = nullptr;
        if(dimple_type_=="flat"){
            // Yes, using 2*flat_dimple_depth_ as the whole dimension of the carving along the z-axis is ok 
            // (the carvings are subtracted from the very edge of the plate)
            carving_solid = dynamic_cast<G4VSolid*>(new G4Box("AUX", flat_dimple_width_/2., dy_/2. +tolerance, flat_dimple_depth_));
        }
        else if(dimple_type_=="spherical"){
            carving_solid = dynamic_cast<G4VSolid*>(new G4Orb("AUX", curvy_dimple_radius_));
        }
        else{ //Default is cylindrical dimples
            carving_solid = dynamic_cast<G4VSolid*>(
                                        new G4Tubs("AUX", 0., curvy_dimple_radius_, dy_/2. +tolerance, 0., twopi));
        }

        G4RotationMatrix* rot = new G4RotationMatrix();
        if(dimple_type_=="cylindrical"){
            rot->rotateX(+90.*deg);
        }

        G4Transform3D* transform_ptr = nullptr;
        G4MultiUnion* carvings_multiunion_solid = new G4MultiUnion("CARVINGS");

        if(dimples_at_z_plus_){
          for(G4int i=0; i<how_many_dimples_; i++){
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector((-1.*dx_/2.)+(1.*(0.5+i)*dx_/(1.*how_many_dimples_)), 0., +dz_/2.));
            carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
          }
        }

        if(dimples_at_z_minus_){
          for(G4int i=0; i<how_many_dimples_; i++){
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector((-1.*dx_/2.)+(1.*(0.5+i)*dx_/(1.*how_many_dimples_)), 0., -1.*dz_/2.));
            carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
          }
        }

        if(dimple_type_=="flat"){
            rot->rotateY(+90.*deg);
        }

        if(dimples_at_x_plus_){
          for(G4int i=0; i<how_many_dimples_; i++){
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector(+dx_/2., 0., (-1.*dz_/2.)+(1.*(0.5+i)*dz_/(1.*how_many_dimples_))));
            carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
          }
        }

        if(dimples_at_x_minus_){
          for(G4int i=0; i<how_many_dimples_; i++){
            transform_ptr = new G4Transform3D(*rot, G4ThreeVector(-1.*dx_/2., 0., (-1.*dz_/2.)+(1.*(0.5+i)*dz_/(1.*how_many_dimples_))));
            carvings_multiunion_solid->AddNode(*carving_solid, *transform_ptr);
          }
        }

        carvings_multiunion_solid->Voxelize();

        G4SubtractionSolid* dimpled_plate_solid = new G4SubtractionSolid(plate_name, 
                                                                        whole_plate_solid, carvings_multiunion_solid);
        geometry_solid = dynamic_cast<G4VSolid*>(dimpled_plate_solid);
    }
    else{
        geometry_solid = dynamic_cast<G4VSolid*>(whole_plate_solid);
    }

    if(cut_plate_){

      G4Para* subtrahend_solid = new G4Para("SUBTRAHEND",
                                            (dx_/2.)+(cut_thickness_/2.),
                                            2.*(dy_/2.), 
                                            2.*(dz_/2.),
                                            0.0,
                                            cut_angle_,                   // Parallelepiped angle with respect to 
                                                                          // the WLS-plate side which is dz_ long
                                            0.0);
      
      G4SubtractionSolid* half_plate_solid = new G4SubtractionSolid(plate_name, 
                                                                    geometry_solid, 
                                                                    subtrahend_solid, 
                                                                    nullptr, 
                                                                    G4ThreeVector(dx_/2., 0., 0.)); // Placing the geometric center of the parallelepiped
                                                                                                    // (subtrahend) right onto the edge of the WLS plate, 
                                                                                                    // so that we are left with one half of the plate minus 
                                                                                                    // half of the cut. The other half of the cut is carved 
                                                                                                    // from the other half of the plate.
      G4MultiUnion* multiunion_geometry_solid = new G4MultiUnion(plate_name);

      G4RotationMatrix* rot_2 = new G4RotationMatrix();
      rot_2->rotateY(180.*deg);

      multiunion_geometry_solid->AddNode( half_plate_solid, 
                                          G4Transform3D(  G4RotationMatrix{}, 
                                                          G4ThreeVector(0., 0., 0.)));  // Note that the geometric center of half_plate_solid 
                                                                                        // is that of the original (uncut) WLS plate, which is 
                                                                                        // geometrically outside the half_plate_solid volume.
      multiunion_geometry_solid->AddNode( half_plate_solid, 
                                          G4Transform3D(  *rot_2,                               // G4Transform3D is a typedef of HepGeom::Transform3D, whose .cc
                                                          G4ThreeVector((1e-6*dx_), 0., 0.)));  // source code is here: 
                                                                                                // apc.u-paris.fr/~franco/g4doxy4.10/html/_transform3_d_8cc_source.html#l00142
                                                                                                // I had to add a residual displacement (1e-6 of dx_), because otherwise, 
                                                                                                // apparently the affine transformation has no inverse (i.e. its inverse 
                                                                                                // transformation has zero determinant) and the error at line 157 of 
                                                                                                // Transform3D.cc pops up. I.e. an affine transformation which consists 
                                                                                                // just of a 180º rotation about the Y axis apparently has no inverse.
      multiunion_geometry_solid->Voxelize();                                                                                  
      geometry_solid = dynamic_cast<G4VSolid*>(multiunion_geometry_solid);
    }                          

    G4Material* pvt = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    if(mpt_ && pvt){
      pvt->SetMaterialPropertiesTable(mpt_);
    }
    else{
      G4Exception("[WLSPlate]", "Construct()", JustWarning,
      "The optical properties of the WLS Plate were not set");
    }

    G4LogicalVolume* geometry_logic =
      new G4LogicalVolume(geometry_solid, pvt, plate_name);

    if(tunneling_probability_!=0.0)
    {
      G4String surface_name = "IMPERFECT_SURFACE";
      G4OpticalSurface* imperfect_surface =
                new G4OpticalSurface( surface_name, 
                                      unified, 
                                      polished, 
                                      dielectric_dielectric);
      imperfect_surface->SetMaterialPropertiesTable(opticalprops::ImperfectDielectricDielectricSurface(tunneling_probability_));
      new G4LogicalSkinSurface(surface_name, geometry_logic, imperfect_surface);
    }

    if(world_logic_vol){
        new G4PVPlacement(nullptr, G4ThreeVector{}, geometry_logic, plate_name, 
                        world_logic_vol, false, 0, true);
    }
    else{
        this->SetLogicalVolume(geometry_logic);
    }
    return;
  }

  void WLSPlate::ConstructCollector(G4LogicalVolume* world_logic_vol){

    G4double collector_thickn = 1.*mm;
    G4double plate_collector_gap = 0.1*mm;

    G4String collector_name = "SURROUNDING_COLLECTOR";

    G4Box* aux = new G4Box( "AUX",    (dx_/2.)+plate_collector_gap+collector_thickn, 
                                      (dy_/2.)+plate_collector_gap+collector_thickn, 
                                      (dz_/2.)+plate_collector_gap+collector_thickn);

    G4Box* subtrahend = new G4Box("AUX",  (dx_/2.)+plate_collector_gap,   
                                          (dy_/2.)+plate_collector_gap,
                                          (dz_/2.)+plate_collector_gap);

    G4SubtractionSolid* collector_solid = new G4SubtractionSolid( collector_name, 
                                                                  aux, subtrahend);

    G4Material* concrete = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");

    G4LogicalVolume* collector_logic = new G4LogicalVolume(collector_solid, concrete, collector_name);

    new G4PVPlacement(nullptr, G4ThreeVector{}, collector_logic, collector_name, 
                      world_logic_vol, false, 0, true);
    return;
  }

  G4ThreeVector WLSPlate::GenerateVertex(const G4String&) const
  {
    if(generation_mode_=="random"){
      G4double tolerance = 1.*mm;
      G4double x_pos =  UniformRandomInRange( (dx_/2.)-tolerance, 
                                              (-dx_/2.)+tolerance);
      G4double z_pos =  UniformRandomInRange( (dz_/2.)-tolerance, 
                                              (-dz_/2)+tolerance); 
      return G4ThreeVector(x_pos, generation_y_pos_, z_pos);

    }
    else{
      return G4ThreeVector(0., generation_y_pos_, 0.);
    }
  }
}