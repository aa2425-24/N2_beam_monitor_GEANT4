
#ifndef SENSTIVEDETECTOR_H
#define SENSTIVEDETECTOR_H

#include "G4VSensitiveDetector.hh"
#include "SDHit.hh"
#include <vector>

#include<iostream>
#include<fstream>



class G4Step; 
class G4G4HCofThisEvent; 

class SensitiveDetector : public G4VSensitiveDetector {
 public : 
  SensitiveDetector(const G4String& SDname, const G4String& hitsCollectionName);
 ~SensitiveDetector() override = default;

public:
void Initialize(G4HCofThisEvent* hitCollection) override; 
G4bool ProcessHits(G4Step* step, G4TouchableHistory* hist) override;
void EndOfEvent(G4HCofThisEvent* hitCollection) override;

private: 
SDHitsCollection* fHitsCollection = nullptr; 
G4long seednbr = 0; 
};

#endif