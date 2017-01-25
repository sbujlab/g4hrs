#ifndef __REMOLLGENERICDETECTORHIT_HH
#define __REMOLLGENERICDETECTORHIT_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class g4hrsGenericDetectorHit : public G4VHit {
    public:
	g4hrsGenericDetectorHit(G4int, G4int);
	~g4hrsGenericDetectorHit();

	g4hrsGenericDetectorHit(const g4hrsGenericDetectorHit &right);
	const g4hrsGenericDetectorHit& operator=(const g4hrsGenericDetectorHit &right);
	G4int operator==(const g4hrsGenericDetectorHit &right) const;

	inline void *operator new(size_t);
	inline void operator delete(void *aHit);
	void *operator new(size_t,void*p){return p;}

    private:

    public:
	G4int fDetID;
	G4int fCopyID;

	// Position and momentum in lab coordinates
	G4ThreeVector f3X;
	G4ThreeVector f3P;
	// Total momentum, energy, mass
	G4double fP, fE, fM;
	// Origin
	G4ThreeVector f3V;
	// Geant4 track ID, particle type, and mother ID
	G4int    fTrID, fPID, fmTrID;
	// Process generator type
	G4int    fGen;

};


typedef G4THitsCollection<g4hrsGenericDetectorHit> g4hrsGenericDetectorHitsCollection;

extern G4Allocator<g4hrsGenericDetectorHit> g4hrsGenericDetectorHitAllocator;

inline void* g4hrsGenericDetectorHit::operator new(size_t){
    void *aHit;
    aHit = (void *) g4hrsGenericDetectorHitAllocator.MallocSingle();
    return aHit;
}

inline void g4hrsGenericDetectorHit::operator delete(void *aHit){
    g4hrsGenericDetectorHitAllocator.FreeSingle( (g4hrsGenericDetectorHit*) aHit);
}

#endif//__REMOLLGENERICDETECTORHIT_HH
