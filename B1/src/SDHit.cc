#include "SDHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iomanip>

#include<iostream>
#include<fstream>



G4ThreadLocal G4Allocator<SDHit>* SDHitAllocator = nullptr; 

G4bool SDHit::operator==(const SDHit& right) const
{
  return ( this == &right ) ? true : false;
}

void SDHit::Print(){
    G4cout
    << " Particle Name " << fName
    << " Parent ID: " << fParentID 
    << " Track ID: " << fTrackID
    << " Edep: " << std::setw(7) << G4BestUnit(fEdep,"Energy")
    << " Position " << std::setw(7) << G4BestUnit(fPos,"Length")
    << " Gobal Time " << std::setw(7) << G4BestUnit(fGtime, "Time")
    << " Local Time " << std::setw(7) << G4BestUnit(fLtime, "Time")
    << G4endl;  
}

void SDHit::toFile(){
  std::ofstream results; 
  results.open("results.csv", std::ios::app); 
    results 
    << fParentID << "," 
    << fTrackID << "," 
    << fName << "," 
    << fPos/CLHEP::mm << ","
    << fEdep/CLHEP::eV << "," 
    << fGtime/CLHEP::ns << "," 
    << fLtime/CLHEP::ns << ","
    << fCopyNo << ","
    << fSeed
    << G4endl;
  results.close();
}

