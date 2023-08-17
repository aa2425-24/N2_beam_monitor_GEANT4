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
/// \file B1/src/EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
//#include "MyAnalysis.hh"
#include "G4AnalysisManager.hh"
#include "G4PrimaryVertex.hh"
#include "G4EventManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
  pos_x=0.;
  pos_y=0.;
  pos_z=0.;
  postEkin = 0.;
  preEkin = 0.; 
  name = " "; 
  trackID = 0; 
  pID = 0; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{

  G4PrimaryVertex* vertex = new G4PrimaryVertex();
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);

  //Accumulate Positions 
  fRunAction->AddPosX(pos_x);
  fRunAction->AddPosY(pos_y);
  fRunAction->AddPosZ(pos_z);

  //Accumulate Kinetic energy 
  fRunAction->AddPostEkin(postEkin); 


  //Acc Ekin pre step 
  fRunAction->AddPreEkin(preEkin); 


  
  //Send information about the particle to runevent
  
  /*
fRunAction->AddTrackID(trackID);
  fRunAction->AddParentID(pID); 
  fRunAction->AddParticleName(name); 
*/
G4int eventID = event->GetEventID();
  if ( eventID < 100 || eventID % 100 == 0) {
    G4cout << ">>> Event: " << eventID  << G4endl;
    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);

    if(hc->GetSize() == 0){
      G4cout << "No hits" << G4endl; 
    }

    G4cout << "    "
           << hc->GetSize() << " hits stored in this event" << G4endl;
  }
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


