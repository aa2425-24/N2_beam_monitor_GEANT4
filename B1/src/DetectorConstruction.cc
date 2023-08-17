//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"

#include "G4Tubs.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

#include "G4VSensitiveDetector.hh"
#include "SensitiveDetector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PhysicalConstants.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Material.hh"


namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 1*cm;
  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_N");


  //***Elements
  G4double z_fh;
  G4double z_fc;
  G4double a_fh;
  G4double a_fc;
  auto fH = new G4Element("H", "H", z_fh = 1., a_fh = 1.01 * g / mole);
  auto fC = new G4Element("C", "C", z_fc = 6., a_fc = 12.01 * g / mole);

  //my Nitrogen
  G4double a;  // atomic mass
  G4double z_nbr;  // atomic number
  G4double density;
  auto my_N = new G4Material("Nitrogen", z_nbr=7., a=28.01*g/mole,density=4.61*kg/m3 ,
                  kStateGas, 293*kelvin, 400000.*pascal);

  // Glass
  auto detector_glass = new G4Material("Glass", density = 1.032 * g / cm3, 2);
  detector_glass->AddElement(fC, 91.533 * perCent);
  detector_glass->AddElement(fH, 8.467 * perCent);

  //Aluminum 
  G4Material* accAl = nist->FindOrBuildMaterial("G4_Al");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = false;
  std::vector<G4double> N_Energy = {3.10*eV,3.13*eV,3.17*eV,3.19*eV,3.21*eV,3.25*eV,3.28*eV,3.30*eV,3.33*eV,3.36*eV,3.45*eV,3.49*eV,3.53*eV,3.65*eV,3.70*eV,3.83*eV,
  3.94*eV,4.02*eV,4.08*eV,4.14*eV,4.29*eV,4.35*eV,4.49*eV,4.58*eV,4.71*eV,4.80*eV,4.94*eV,5.02*eV,
  5.16*eV,5.25*eV,5.38*eV,5.47*eV,5.66*eV,5.77*eV,5.89*eV};

  std::vector<G4double> N_SCINT = {0.,0.027398206,0.,0.02733212,0.005209789,0.001489691,0.150673609,0.02722473,0.06773277,0.00507211,0.,0.504182056,
  0.065717143,0.001139986,0.688229158,0.,0.275416676,0.046940423,0.000850859,0.240291901,0.,0.107538646,0.000622311,0.21611264,0.,
  0.40391013,0.,0.646980402,0.,1.,0.,0.609969419,0.,0.092228693,0.001927512};

  std::vector<G4double> N_RIND  = {1.0003, 1.0003, 1.0003 , 1.0003, 1.0003,1.0003, 1.0003, 1.0003 , 1.0003, 1.0003,1.0003, 1.0003, 1.0003 , 1.0003, 1.0003,
  1.0003, 1.0003, 1.0003 , 1.0003, 1.0003,1.0003, 1.0003, 1.0003 , 1.0003, 1.0003,1.0003, 1.0003, 1.0003 , 1.0003, 1.0003,
  1.0003, 1.0003, 1.0003 , 1.0003, 1.0003};

  std::vector<G4double> N_ABSL  = {10. * m,10. * m,10. * m ,10. * m,10. * m,10. * m,10. * m,10. * m ,10. * m,10. * m,10. * m,10. * m,10. * m ,10. * m,10. * m
  ,10. * m,10. * m,10. * m ,10. * m,10. * m,10. * m,10. * m,10. * m ,10. * m,10. * m,10. * m,10. * m,10. * m ,10. * m,10. * m
  ,10. * m,10. * m,10. * m ,10. * m,10. * m};

//Properties to make the scintillation work
   G4MaterialPropertiesTable* N_mt = new G4MaterialPropertiesTable(); 
  N_mt->AddProperty("SCINTILLATIONCOMPONENT1", N_Energy, N_SCINT);
  N_mt->AddProperty("RINDEX", N_Energy, N_RIND);
  N_mt->AddProperty("ABSLENGTH", N_Energy, N_ABSL);
  N_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 0 * ns);
  N_mt->AddConstProperty("SCINTILLATIONYIELD", 10000. / MeV);
  N_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  N_mt->AddConstProperty("RESOLUTIONSCALE", 1);
  //env_mat->SetMaterialPropertiesTable(N_mt);
  my_N->SetMaterialPropertiesTable(N_mt); 


  //Optical properties of the detector glass 
  std::vector<G4double> glass_AbsLength = { 420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm, 420. * cm,
  420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm, 420. * cm,420. * cm, 420. * cm};
  G4MaterialPropertiesTable* glass_mt   = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH", N_Energy, glass_AbsLength);
  glass_mt->AddProperty("RINDEX", "Fused Silica");
  detector_glass->SetMaterialPropertiesTable(glass_mt); 
  
 
  // World
  //
  G4double world_sizeXY = 0.5*env_sizeXY;
  G4double world_sizeZ  = 0.5*env_sizeXY;
  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_void");

  G4Material* world_mat = new G4Material("Vacuum", 1.,1.01 * g / mole,
                           density = universe_mean_density, kStateGas,
                           0.1 * kelvin, 1.e-19 * pascal);

  auto solidWorld = new G4Box("World",                           // its name
  world_sizeXY,world_sizeXY,world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    22,                                         // copy number
    checkOverlaps);                            // overlaps checking


//The Nitrogen volume cylinder slab 
G4double inner1 = 0.*cm;
G4double outer1 =2.5*cm; 
G4double hz1 = 1*cm;
G4double sAngle1 = 0.* deg; 
G4double sweep1 = 360.*deg;
//solid
auto solidCyl = new G4Tubs("solidCyl", inner1, outer1,hz1,sAngle1,sweep1);
//logic shit 
auto logicCyl = new G4LogicalVolume(solidCyl, my_N, "Logical Beam Detector");

new G4PVPlacement(0,
G4ThreeVector(),
logicCyl,
"logicCyl",
logicWorld, 
false,
21,
checkOverlaps);


///// Ring of sensitive detectors
G4double ri = 0.*cm;
G4double ry = 4*mm;
G4double z = 0.1*cm;
G4double m = 0.*deg;
G4double sw = 360.*deg;

auto solFiberopticglass = new G4Tubs("Ofiber",ri,ry,z,m,sw);
auto logFiberopticglass = new G4LogicalVolume(solFiberopticglass, detector_glass, "Ofiber");

G4int nb_det = 12;
G4int nb_ring = 1; 

 G4double dPhi = twopi/nb_det, half_dPhi = 0.5*dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);


  // Ring to which the detectors are attached to
auto solidRing = new G4Tubs("Ring", 0, 2.5*cm, 1*cm, 0., twopi);

auto logicRing = new G4LogicalVolume(solidRing,  // its solid
   my_N,                                   // its material
    "Ring");         //byt tbx till envmat!!!!

    new G4PVPlacement(0,
G4ThreeVector(),logicRing,
"logicRing",
logicWorld, 
false,
20,
checkOverlaps);

 
//Place detectors in the ring 
                                 // its name
 for (G4int icrys = 0; icrys < nb_det ; icrys++) {
    G4double phi = icrys*dPhi;
    G4double r = 2.5*cm;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg);
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
    G4ThreeVector position = r*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);

    new G4PVPlacement(transform,             //rotation,position
                      logFiberopticglass,            //its logical volume
                      "Ofiber",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      checkOverlaps);       // checking overlaps
  } 

  //The accelerator

  /* G4double zAcc = 5. *cm; 
  G4double inner3 =2.5*cm; 
  G4double outer3 =2.6*cm;
  auto solidAcc = new G4Tubs("solidAcc",inner3,outer3,zAcc,sAngle1,sweep1);
  auto logicAcc = new G4LogicalVolume(solidAcc,accAl,"solidAcc");

  new G4PVPlacement(0,G4ThreeVector(0,0,-6*cm),
  logicAcc,"logicAcc",logicWorld,false,30); 
 */

//Aluminum hull 
/*   G4double hz2 = 1*mm;

//solid
auto solidAlHull = new G4Tubs("solidCyl", inner1, outer1,hz2,sAngle1,sweep1);
//logic shit 
auto logicAlHull = new G4LogicalVolume(solidAlHull, accAl,"Al hull");

new G4PVPlacement(0,
G4ThreeVector(0,0,-2*cm),
logicCyl,
"logicCyl",
logicWorld, 
false,
33,
checkOverlaps);

new G4PVPlacement(0,
G4ThreeVector(0,0,2*cm),
logicCyl,
"logicCyl",
logicWorld, 
false,
33,
checkOverlaps);  

 */
 

  // Set logicOfiber as scoring volume
//fScoringVolume = logFiberopticglass;

  
G4OpticalSurface* OpSurface= new G4OpticalSurface("detectorsurface");
OpSurface -> SetType(dielectric_metal);
OpSurface -> SetFinish(ground);
OpSurface -> SetModel(unified);
G4LogicalSkinSurface* reflsurf = new G4LogicalSkinSurface("reflsurf",logicCyl,OpSurface); 
G4LogicalSkinSurface* reflsurf1 = new G4LogicalSkinSurface("reflsurf",logicRing,OpSurface);   

  //
  //always return the physical World

  //
  return physWorld;




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

void DetectorConstruction::ConstructSDandField(){
  // its name
  G4String trackerChamberSDname = "/SD";
  SensitiveDetector* envSD = new SensitiveDetector(trackerChamberSDname, "SDHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(envSD);
  SetSensitiveDetector("Ofiber",envSD,true); 
}
//Comment graveyard 
// N_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45. * ns);
  // Set the Birks Constant for the N scintillator
 // env_mat->GetIonisation()->SetBirksConstant(0.126 * mm / MeV); //????
  // N_mt->AddProperty("SCINTILLATIONCOMPONENT2", N_Energy, N_SCINT);
  // N_mt->AddConstProperty("SCINTILLATIONYIELD2", 1.0); 
//G4Material* detector_glass = nist->FindOrBuildMaterial("GLASS");
} 
  /*   
   new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    logicEnv,                 // its logical volume
    "Envelope",               // its name
    logicWorld,               // its mother  volume
    false,                    // no boolean operation
    1,                        // copy number
    checkOverlaps);           // overlaps checking */ 

 //

/*   G4MaterialPropertiesTable* glass_mt = new G4MaterialPropertiesTable();

std::vector<G4double> glass_energy {1.0*eV, 4.0*eV, 6.0*eV};
std::vector<G4double> glass_ref{1.0,1.0,1.0};
std::vector<G4double> glass_abs {100*cm,100*cm,100*cm};


glass_mt->AddProperty("RINDEX",glass_energy,glass_ref);
glass_mt->AddProperty("ABSLENGTH",glass_energy,glass_ref);

detector_glass->SetMaterialPropertiesTable(glass_mt);  */

 //
  // Envelope
  //
/*   auto solidEnv = new G4Box("Envelope",                    // its name
    0.5*env_sizeXY, 0.5 * env_sizeXY,env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
    env_mat,                                     // its material
    "Envelope");  

 */
 
/* 
///Sensitive detector volume// the fiber optic "cable"
G4double inner = 0.*cm;
G4double outer = 0.00001*mm; 
G4double hz = 0.2*cm;
G4double sAngle = 0.* deg; 
G4double sweep = 360.*deg;
//solid
G4Tubs* solidOfiber = new G4Tubs("Ofiber",inner,outer,hz,sAngle,sweep); 
//logic
G4LogicalVolume* logicOfiber = new G4LogicalVolume(solidOfiber, env_mat, "Ofiber");
//place
new G4PVPlacement(0,G4ThreeVector(0,0,1*cm),logicOfiber,"Ofiber", logicWorld, false,2,checkOverlaps); 

 */

  //std::vector<G4double> N_Energy = {3.1 *eV,3.75*eV, 4.1*eV, 5.4*eV, 6.2*eV};
  //std::vector<G4double> N_SCINT = {0.02, 0.7, 0.21, 1.0, 0.001};


  /* 
  // Conical section shape
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.5*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 0.1*cm;
  G4double shape1_hz = 0.5*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
   G4ThreeVector pos1 = G4ThreeVector(0, 0*cm, 1.5*cm);


  auto solidShape1 = new G4Cons("Shape1", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb,
    shape1_hz, shape1_phimin, shape1_phimax);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
    my_N,                                        // its material
    "Shape1");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos1,                     // at position
    logicShape1,              // its logical volume
    "Shape1",                 // its name
    logicCyl,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

 */

/* 
//Rectangular detector instead 
G4double xdet = 1. * cm;
G4double ydet = 3. *mm;
G4double zdet = 0.1*cm; 
G4int Nb_of_rec_det = 24;
G4double dphi = twopi/Nb_of_rec_det; 

auto solidRecDet = new G4Box("RecDet",xdet,ydet,zdet);
auto logicRecDet = new G4LogicalVolume(solidRecDet, detector_glass, "Ofiber");




//Place detectors in the ring 
                                 // its name
for (G4int icrys = 0; icrys < Nb_of_rec_det ; icrys++) {
    G4double phi = icrys*dphi;
    G4double r = 2.5*cm;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg);
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
    G4ThreeVector position = r*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);

    new G4PVPlacement(transform,             //rotation,position
                      logicRecDet,            //its logical volume
                      "Ofiber",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      checkOverlaps);       // checking overlaps
  } */