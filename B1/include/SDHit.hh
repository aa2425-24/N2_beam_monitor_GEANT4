#ifndef SDHIT_h
#define SDHIT_h

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

#include<iostream>
#include<fstream>



class SDHit : public G4VHit{
    public:

    SDHit() = default; 
    SDHit(const SDHit&) = default; 
    ~SDHit() override = default; 

    SDHit& operator=(const SDHit&) = default;
    G4bool operator==(const SDHit&) const;

    inline void* operator new(size_t); 
    inline void operator delete(void*);

    void Print() override; 

    //Set methods 
    void SetEdep(G4double ed){fEdep = ed;};
    void SetPos(G4ThreeVector xyz){fPos = xyz;};
    void SetTrackID(G4int track){fTrackID = track;};
    void SetParentID(G4int pid){fParentID = pid;};
    void SetParticleName(G4String name){fName = name;};
    void SetGlobalTime(G4double gt){fGtime = gt;};
    void SetLocalTime(G4double lt){fLtime = lt;};
    void SetCopyNo(G4int nr){fCopyNo = nr;};
    void SetSeedNbr(G4double seed){fSeed = seed;};
    void toFile();
    //get methods 
    G4int GetTrackID() const {return fTrackID;};
    G4double GetEdep() const {return fEdep; }; 
    G4ThreeVector GetPos() const {return fPos; };
    G4int GetParentID() const{return fParentID; };
    G4String GetParticleName() const {return fName; }; 
    G4double GetGlobalTime() const {return fGtime; }; 
    G4double GetLocalTime() const {return fLtime; };
    G4int GetCopyNo() const {return fCopyNo;}; 
    G4double GetSeedNbr() const {return fSeed;};
    private: 
    G4int fTrackID = -1;
    G4int fParentID = -1;
    G4double fEdep = 0.;
    G4ThreeVector fPos; 
    G4String fName = " ";
    G4double fGtime = 0.;
    G4double fLtime = 0.; 
    G4int fCopyNo = 0.;
    G4int count = 0;
    G4double fSeed = 0;

};

using SDHitsCollection = G4THitsCollection<SDHit>;

extern G4ThreadLocal G4Allocator<SDHit>* SDHitAllocator; 

inline void* SDHit::operator new(size_t){
    if(!SDHitAllocator)
        SDHitAllocator = new G4Allocator<SDHit>;
    return (void*) SDHitAllocator->MallocSingle();
}

inline void SDHit::operator delete(void* hit){
    SDHitAllocator->FreeSingle((SDHit*) hit);
}

#endif