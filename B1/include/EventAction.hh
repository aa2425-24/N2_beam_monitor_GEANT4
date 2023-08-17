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
/// \file B1/include/EventAction.hh
/// \brief Definition of the B1::EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4PrimaryVertex.hh"

#include<iostream>
#include<fstream>



/// Event action class
///

namespace B1
{

class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;
    
    void AddEdep(G4double edep) { fEdep += edep; }
    
    void AddPosX(G4double posx){pos_x = posx;}
    void AddPosY(G4double posy){pos_y = posy;}
    void AddPosZ(G4double posz){pos_z = posz;}
    
    void AddPostEkin(G4double ekin){postEkin = ekin;}
    void AddPreEkin(G4double pre){preEkin = pre;}
    
    void AddParentID(G4int parent){pID = parent;}
    void AddTrackID(G4int track){trackID = track;}
    void AddParticleName(G4String pname){name = pname;}

  private:
    RunAction* fRunAction = nullptr;
    G4double   fEdep = 0.;
    G4double pos_x = 0.; 
    G4double pos_y = 0.; 
    G4double pos_z = 0.; 
    G4double postEkin = 0.; 
    G4double preEkin = 0.; 

    G4String name = " "; 
    G4int trackID = 0; 
    G4int pID = 0; 
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


