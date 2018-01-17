#ifndef g4hrsTune_h
#define g4hrsTune_h 

#include "globals.hh"	//for units and g4io 
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4String.hh"

class g4hrsTune {

public:
	g4hrsTune();
	virtual ~g4hrsTune();
	static g4hrsTune* GetTune();
	void SetTune(G4String);

private:
	static g4hrsTune* fTune;

public:
	double snakeMagnet;
	double kappaQuad1;
	double kappaQuad2;
	double kappaQuad3;
	double bendDipole;
	double septumAngle;
	double septumMomentum;
	double septumCurrent;
	double HRSMomentum;
	double HRSAngle;
	int quadsOn;
	double sosQuad;	
};

#endif
