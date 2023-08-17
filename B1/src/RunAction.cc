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
/// \file B1/src/RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingAction.hh"
#include "G4AnalysisManager.hh"

#include<iostream>
#include<fstream>


namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // add new units for dose
  //
  /* const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); */

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);
  /* accumulableManager->RegisterAccumulable(position_x);
  accumulableManager->RegisterAccumulable(position_y);
  accumulableManager->RegisterAccumulable(position_z);
  accumulableManager->RegisterAccumulable(postEkin);
  accumulableManager->RegisterAccumulable(preEkin); */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{

/*   long seeds[2];
  time_t systime = time(NULL);
  seeds[0] = systime;
  seeds[1] = (systime*G4UniformRand());
  G4Random::setTheSeed(seeds);

  auto seedvec = G4Random::getTheSeed();
G4cout << "seed first position " << seedvec[0] << " and the second " << seedvec[1] << G4endl; */
G4cout << "seed from run " << G4Random::getTheSeed() << G4endl;
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();


  
/* if(fPrimary){
  sG4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
  G4bool polarized = fPrimary->GetPolarized();
  G4double polarization = fPrimary->GetPolarization();
  fRun->SetPrimary(particle,energy,polarized,polarization); 
}
auto analysisManager = G4AnalysisManager::Instance();
if(analysisManager->IsActive()){
analysisManager->OpenFile();
}
 */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
//write to file

 G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();


  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
  // Print
  //
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";

/* 
     results
     << G4endl
     << "--------------------End of Global Run-----------------------"; */
  }
  else {
    /* results
     << G4endl
     << "--------------------End of Local Run------------------------"; */
  }

 /*  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : "
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl; */


     //write to file
   /*  results << G4endl
     << " Position of the particle " << G4endl <<
     " x: "
     << position_x 
     << G4endl 
     << " y: "<< position_y
     << G4endl
     <<" z: "<< position_z << G4endl;

    results << "Kinetic Energy in preStep of the particle: " << preEkin << G4endl; 
    results << "Kinetic Energy in postStep of the particle: " << postEkin << G4endl; 
    results << "Total Energy Deposited: " << fEdep.GetValue() << " MeV" << G4endl;


 */

/* 
     G4cout 
     << G4endl
     << " Position of the particle " << " x: "
     << position_x 
     << G4endl 
     << " y: "<< position_y
<< G4endl
 <<" z: "<< position_z << G4endl
 << "------------------------------------------------------------"
<< G4endl;  

G4cout << "Kinetic Energy in preStep of the particle: " << preEkin << G4endl; 
G4cout << "Kinetic Energy in postStep of the particle: " << postEkin << G4endl; 

G4cout << "Total Energy Deposited: " << fEdep.GetValue() << " MeV" << G4endl;  */
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

void RunAction::AddPosX(G4double pos_x){
  position_x = pos_x; 
}
void RunAction::AddPosY(G4double pos_y){
  position_y = pos_y; 
}
void RunAction::AddPosZ(G4double pos_z){
  position_z = pos_z; 
}
void RunAction::AddPostEkin(G4double ekin){
  postEkin = ekin; 
}
void RunAction::AddPreEkin(G4double pre){
  preEkin = pre; 
}

/*void AddParentID(G4int pID){
 prntID = pID; 
}
void AddTrackID(G4int tID){
tID = track;  
}
void AddParticleName(G4String name){
name = pname; 
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
