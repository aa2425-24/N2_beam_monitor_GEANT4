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
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4StepPoint.hh"
#include "G4UnitsTable.hh"

#include<iostream>
#include<fstream>



namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double ener = 0.;
G4double tEner =0.;
G4double x_pos;
G4double y_pos;
G4double z_pos; 
G4double s_length;
G4int countPh = 0; 
G4int countelec = 0;
G4double countedep = 0;
G4String pname;
G4int tID;
G4int pID; 

void SteppingAction::UserSteppingAction(const G4Step* step)
{

std::ofstream resultsVerify; 
resultsVerify.open("resultsVerify.csv",std::ios::app); 

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
  countedep += edepStep;

  //Collect position of this step
  G4Track* ftrack = step->GetTrack(); 
  G4double x_pos = ftrack->GetPosition().x();
  G4double y_pos = ftrack->GetPosition().y();
  G4double z_pos = ftrack->GetPosition().z();

//These fEventAction-> xyz functions is not used anywhere and can be removed 
   fEventAction->AddPosX(x_pos);
   fEventAction->AddPosY(y_pos);
   fEventAction->AddPosZ(z_pos);
   //Kinetic energy post step 
   G4double Ekin = step->GetPostStepPoint()->GetKineticEnergy();
   fEventAction->AddPostEkin(Ekin);

   // -- pre step 
   G4double preEkin = step->GetPreStepPoint()->GetKineticEnergy();
   fEventAction->AddPreEkin(preEkin); 

   //name the process
  G4String process = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
 // G4cout << process << G4endl; 

  //Trying to find all secondaries 
  //G4int parentID = step->GetTrackID();
  //secondary = step->GetfSecondaries();
  //G4cout << step->GetTrack()->GetDefinition() << G4endl; 
  G4Track* ftrack1 = step->GetTrack();
  G4int pID = ftrack1->GetParentID();
  G4int trackID = ftrack1->GetTrackID();
  G4String pName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

  G4double locTime = ftrack1->GetLocalTime(); 
  G4double globTime = ftrack1->GetGlobalTime();
  G4double propTime = ftrack1->GetProperTime();
  
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  

   //get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();  

       auto pressure
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume()->GetMaterial();

     //G4cout << pressure << G4endl; 

  G4double seed = G4Random::getTheSeed(); 
// check if we are in scoring volume
  /* if (volume == fScoringVolume){
    if(pName == "opticalphoton"){
      countPh++;
     G4cout << "HERE IS THE photon COUNT" << G4endl<< "------------------------------------" << G4endl << countPh <<G4endl; 
    }
    else if(pName == "e-"){
      countelec++; 
     G4cout << "HERE IS THE ELECTRON COUNT" << G4endl<< "------------------------------------" << G4endl << countelec << G4endl; 
    }
  } */
 
resultsVerify << pID << "," << trackID << "," << pName << "," << process << "," << x_pos/CLHEP::mm << "," << y_pos/CLHEP::mm << "," << z_pos/CLHEP::mm
  << "," <<  edepStep/CLHEP::eV << "," <<  preEkin/CLHEP::eV << "," << Ekin/CLHEP::eV << "," << globTime/CLHEP::ns << "," << locTime/CLHEP::ns << "," << propTime/CLHEP::ns << "," << seed << G4endl;

resultsVerify.close();

}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


