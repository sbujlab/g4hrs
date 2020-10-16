#ifndef __REMOLLEVENT_HH
#define __REMOLLEVENT_HH

/*!
   Event information class.  This needs to
   contain all the information we need to
   generate particles and provide sufficient
   output.
*/

#include <vector>
#include "G4ThreeVector.hh"

class G4ParticleDefinition;

class g4hrsEvent {
    public:
	 g4hrsEvent();
	~g4hrsEvent();

	void ProduceNewParticle( G4ThreeVector, G4ThreeVector, G4String );
	void SetEffCrossSection( G4double xs ){ fEffXs = xs; }
	void SetRate( G4double rate ){ fRate = rate; }
	G4double GetRate(){ return fRate; }
	void SetAsymmetry( G4double A ){ fAsym = A; }
	void SetSensitivity(G4double S) { fSens = S; }
	void SetWeight(G4double V) {fWeight = V;}

	void SetQ2( G4double q2 ){ fQ2 = q2; }
	void SetW2( G4double w2 ){ fW2 = w2; }
	void SetXbj( G4double x ){ fXbj = x; }
	void SetThCoM( G4double th ){ fThCoM = th; }
        void SetNuclA( G4int val){ fNuclA  = val; }


	void Reset();
	void UndoLastParticle();

	G4bool EventIsSane();
	void   Print();

    private:

    public:
	// Interaction information
	G4ThreeVector fBeamMomentum;
	G4ThreeVector fVertexPos;

	// Particles to be produced
	std::vector<G4ThreeVector>    fPartPos;
	std::vector<G4ThreeVector>    fPartMom;  // Generated direction (no ms)
	std::vector<G4ThreeVector>    fPartRealMom; // Direction to go to Geant4
	std::vector<G4ParticleDefinition *> fPartType;

	G4double fBeamE;
	G4double fRate;
	G4double fEffXs;
	G4double fAsym, fmAsym;
	G4double fSens;
	G4double fWeight;

	G4double fQ2;
	G4double fW2;
	G4double fXbj;
	G4double fThCoM;
 	G4double fVEP;   // scattered momentum right after the vertex
        G4int fNuclA; //Nuclear A 
        G4double fEvThSamp; //vertex theta
};

#endif//__REMOLLEVENT_HH
