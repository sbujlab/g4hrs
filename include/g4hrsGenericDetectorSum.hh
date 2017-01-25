#ifndef __REMOLLGENERICDETECTORSUM_HH
#define __REMOLLGENERICDETECTORSUM_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class g4hrsGenericDetectorSum : public G4VHit {
    public:
	g4hrsGenericDetectorSum(G4int, G4int);
	~g4hrsGenericDetectorSum();

	g4hrsGenericDetectorSum(const g4hrsGenericDetectorSum &right);
	const g4hrsGenericDetectorSum& operator=(const g4hrsGenericDetectorSum &right);
	G4int operator==(const g4hrsGenericDetectorSum &right) const;

	inline void *operator new(size_t);
	inline void operator delete(void *aHit);
	void *operator new(size_t,void*p){return p;}

    public:
	G4int    fDetID;
	G4int    fCopyID;
	G4double fEdep;
};


typedef G4THitsCollection<g4hrsGenericDetectorSum> g4hrsGenericDetectorSumCollection;

extern G4Allocator<g4hrsGenericDetectorSum> g4hrsGenericDetectorSumAllocator;

inline void* g4hrsGenericDetectorSum::operator new(size_t){
    void *aHit;
    aHit = (void *) g4hrsGenericDetectorSumAllocator.MallocSingle();
    return aHit;
}

inline void g4hrsGenericDetectorSum::operator delete(void *aHit){
    g4hrsGenericDetectorSumAllocator.FreeSingle( (g4hrsGenericDetectorSum*) aHit);
}

#endif//__REMOLLGENERICDETECTORSUM_HH
