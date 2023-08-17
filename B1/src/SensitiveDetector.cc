
#include <iostream> 
#include "G4VSensitiveDetector.hh"
#include "SensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

#include<iostream>
#include<fstream>



SensitiveDetector::SensitiveDetector(const G4String& sdname,
                                    const G4String& hitsCollectionName) : G4VSensitiveDetector(sdname){
G4cout<<"Creating SD with name: "<<sdname<<G4endl;
collectionName.insert(hitsCollectionName);

}

void SensitiveDetector::Initialize(G4HCofThisEvent* hce){

    //create
    fHitsCollection = new SDHitsCollection(SensitiveDetectorName, collectionName[0]);
    //Add Collection
    G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID,fHitsCollection); 

    

    /* const auto detConstruction = static_cast<const B1::DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume(); */
}

G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory* hist){

    G4double edep = step->GetTotalEnergyDeposit(); 
//G4double edep = step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy();
   G4int volume
    = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();  
    seednbr = CLHEP::HepRandom::getTheSeed();

 
if(volume == 0 || volume == 1 ||volume == 2 ||volume == 3 ||volume == 4 ||volume == 5 ||
volume == 6 ||volume == 7 ||volume == 8 ||volume == 9 ||volume == 10 ||volume == 11){
    //if(edep == 0.) return false;
    
    auto newHit = new SDHit();

    newHit->SetTrackID(step->GetTrack()->GetTrackID());
    newHit->SetPos(step->GetPostStepPoint()->GetPosition());
    newHit->SetEdep(edep);
    newHit->SetParentID(step->GetTrack()->GetParentID());
    newHit->SetParticleName(step->GetTrack()->GetParticleDefinition()->GetParticleName());
    newHit->SetGlobalTime(step->GetTrack()->GetGlobalTime());
    newHit->SetLocalTime(step->GetTrack()->GetLocalTime());
    newHit->SetCopyNo(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
    newHit->SetSeedNbr(seednbr);

    fHitsCollection->insert(newHit);
    //step->GetTrack()->SetTrackStatus(fStopAndKill);
}
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*){
    G4cout << "seed " << G4Random::getTheSeed() << G4endl;

    G4int nofHits = fHitsCollection->entries(); 
    G4cout << G4endl
    << "Number of hits in the target volume is " << nofHits 
    << " Hits: " << G4endl;
    if(verboseLevel > 1){
    for(G4int i = 0; i<nofHits; i++) (*fHitsCollection)[i]->Print();
    }
    for(G4int i = 0; i<nofHits; i++) (*fHitsCollection)[i]->toFile();
}
