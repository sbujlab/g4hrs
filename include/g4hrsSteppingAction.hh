
#ifndef __REMOLLSTEPPINGACTION_HH
#define __REMOLLSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class g4hrsTransportFunction; 

class g4hrsSteppingAction : public G4UserSteppingAction
{
  public:
    g4hrsSteppingAction();
    virtual ~g4hrsSteppingAction(){};

    virtual void UserSteppingAction(const G4Step*);

    void SetEnableKryptonite(G4bool k){ fEnableKryptonite = k; }

  private:
    G4bool drawFlag;
	G4double rad;
    G4bool fEnableKryptonite;
	g4hrsTransportFunction* fTransportFunction;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };

	int nelements;

	G4double fSeptumAngle;		// septum angle obtained from messenger (constant)
	G4double fHRSAngle;		// HRS angle obtained from messenger (constant)
	G4double septum_angle;		// local septum angle (changes sign depending on L/R HRS)
	G4double hrs_angle;		// local HRS angle (changes sign depending on L/R HRS)
	G4double fHRSMomentum; 		// HRS central momentum
	bool goodParticle;
	double sign; 	// y/phi sign flip for using transport function on left arm

	G4int fLHRS;
	G4int fRHRS;
		
	float r0[5];
	G4double x_tf[12];
	G4double y_tf[12];
	G4double th_tf[12];
	G4double ph_tf[12];

	int numVB;
	int numVar;

	G4double VBdata[14][12];

	G4String VBnames[14];


};

#endif//__REMOLLSTEPPINGACTION_HH
