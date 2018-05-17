#ifndef __REMOLLVEVENTGEN_HH
#define __REMOLLVEVENTGEN_HH

#include "g4hrstypes.hh"
#include "g4hrsglobs.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "g4hrsVertex.hh"

/*!
   Generic base class for event generators
   This provides an interface for everyone to
   derive from.

   Ultimately this returns a g4hrsEvent which is
   what the PrimaryGeneratorAction is going to use and
   contains information that will go in the output.

   It needs to be aware of g4hrsBeamTarget and g4hrsRunData,
   take a generically generated event assuming ideal beam
   and transform it into what is going to be simulated.
*/

class g4hrsEvent;
class g4hrsBeamTarget;
class g4hrsRunData;

class g4hrsVEventGen {
    public:
	g4hrsVEventGen();
	virtual ~g4hrsVEventGen();

	g4hrsEvent *GenerateEvent();

	G4String GetName(){ return fName; }

	void SetSampType( SampType_t st ) { fSampType = st; }
	void SetDoMultScatt( G4bool ms ){ fApplyMultScatt = ms; }

	virtual void SetThMin( double th ){ fTh_min = th; }
	virtual void SetThMax( double th ){ fTh_max = th; }
	virtual void SetPhMin( double ph ){ fPh_min = ph; }
	virtual void SetPhMax( double ph ){ fPh_max = ph; }
	virtual void SetEmin( double ){ G4cerr << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__ << " :  Generator does not respond to this command" << G4endl; }
	virtual void SetEmax( double ){ G4cerr << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__ << " :  Generator does not respond to this command" << G4endl; }
	virtual void SetThCoM_min( double ){ G4cerr << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__ << " :  Generator does not respond to this command" << G4endl; }
	virtual void SetThCoM_max( double ){ G4cerr << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__ << " :  Generator does not respond to this command" << G4endl; }

	G4double fThCoM_min, fThCoM_max;
	G4double fTh_min, fTh_max;
	G4double fPh_min, fPh_max;
	G4double fE_min, fE_max;
	G4double fSeptumAngle;

	G4ThreeVector fSetVPos;
	G4bool fIsVPosSet;
	G4ThreeVector fSetVMom;
	G4bool fIsVMomSet;
	G4double fSetVTheta;
	G4bool fIsVThetaSet;

    private:
	const G4String fName;

	g4hrsBeamTarget *fBeamTarg;
	g4hrsRunData    *fRunData;
	void PolishEvent(g4hrsEvent *);
	
	// Pure virtual function that needs to be filled out
	virtual void SamplePhysics(g4hrsVertex *, g4hrsEvent *) = 0;

    protected:

	SampType_t fSampType;
	G4bool     fApplyMultScatt;

};


#endif//__REMOLLVEVENTGEN_HH
