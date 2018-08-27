#ifndef g4hrsMessenger_HH
#define g4hrsMessenger_HH

#include "globals.hh"
#include "g4hrstypes.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4VModularPhysicsList.hh"

/*!
 *   Global messenger class
 */

class g4hrsIO;
class g4hrsDetectorConstruction;
class g4hrsEventAction;
class g4hrsPrimaryGeneratorAction;
class g4hrsBeamTarget;
class g4hrsSteppingAction;
class g4hrsEMFieldSetup;
class g4hrsEMField;
class g4hrsTune;

class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;

class g4hrsMessenger : public G4UImessenger {
    public:
       	g4hrsMessenger();
       	~g4hrsMessenger();

	void SetIO( g4hrsIO *io ){ fIO = io; }
	void SetPriGen( g4hrsPrimaryGeneratorAction *pg ){ fprigen = pg; }
	void SetDetCon( g4hrsDetectorConstruction *dc ){ fdetcon= dc; }
	void SetEvAct( g4hrsEventAction *ev ){ fevact = ev; }
	void SetStepAct( g4hrsSteppingAction *st ){ fStepAct = st; }
	void SetPhysList( G4VModularPhysicsList *l ){ fPhysicsList = l; }

	void SetNewValue(G4UIcommand* cmd, G4String newValue);

    private:
	g4hrsIO *fIO;
	g4hrsDetectorConstruction *fdetcon;
	g4hrsEventAction *fevact;
	g4hrsPrimaryGeneratorAction *fprigen;
	g4hrsBeamTarget *fBeamTarg;
	g4hrsSteppingAction *fStepAct;
	g4hrsTune* ftune;
	G4VModularPhysicsList *fPhysicsList;

	G4UIcmdWithAnInteger *seedCmd;
	G4UIcmdWithABool     *kryptCmd;
	G4UIcmdWithABool     *opticalCmd;
        G4UIcmdWithABool     *dumpGeometryCmd;

	G4UIcmdWithAString   *detfilesCmd;

	G4UIcmdWithAString   *newfieldCmd;
	G4UIcmdWithAString   *fieldScaleCmd;
	G4UIcmdWithAString   *fieldCurrCmd;
	G4UIcmdWithAString   *genSelectCmd;

	G4UIcmdWithAnInteger		*snakeModCmd;
	G4UIcmdWithADoubleAndUnit	*sepAngCmd;
	G4UIcmdWithADoubleAndUnit	*hrsAngCmd;
	G4UIcmdWithADoubleAndUnit	*hrsMomCmd;
	G4UIcmdWithAnInteger		*hrsSetupCmd;
	G4UIcmdWithADoubleAndUnit	*sepMomCmd;
	G4UIcmdWithADoubleAndUnit	*sepCurCmd;
	G4UIcmdWithADoubleAndUnit	*q1kappaCmd;
	G4UIcmdWithADoubleAndUnit	*q2kappaCmd;
	G4UIcmdWithADoubleAndUnit	*dBendCmd;
	G4UIcmdWithADoubleAndUnit	*q3kappaCmd;
	G4UIcmdWithAString		*fTuneCmd;

	G4UIcmdWithADoubleAndUnit *minEKillCmd;

	G4UIcmdWithADoubleAndUnit *tgtLenCmd;
	G4UIcmdWithADoubleAndUnit *tgtPosCmd;
	G4UIcmdWithAString 	  *tgtMatCmd;

	G4UIcmdWithADoubleAndUnit *beamCurrCmd;
	G4UIcmdWithADoubleAndUnit *beamECmd;

	G4UIcmdWithABool       *rasTypeCmd;

	G4UIcmdWithADoubleAndUnit *rasXCmd;
	G4UIcmdWithADoubleAndUnit *rasYCmd;

	G4UIcmdWithADoubleAndUnit *beamX0Cmd;
	G4UIcmdWithADoubleAndUnit *beamY0Cmd;

	G4UIcmdWithADoubleAndUnit *beamth0Cmd;
	G4UIcmdWithADoubleAndUnit *beamph0Cmd;

	G4UIcmdWithADoubleAndUnit *beamCorrThCmd;
	G4UIcmdWithADoubleAndUnit *beamCorrPhCmd;

	G4UIcmdWithADoubleAndUnit *beamdthCmd;
	G4UIcmdWithADoubleAndUnit *beamdphCmd;

	G4UIcmdWithAString   *fileCmd;
	G4UIcmdWithAString   *pionCmd;
	G4UIcmdWithAString   *LUNDFileCmd;//Dominic Lunde linking GenLUND

	////////////////////////////////////////////////
	// To general event generators
	G4UIcmdWithADoubleAndUnit *thminCmd;
	G4UIcmdWithADoubleAndUnit *thmaxCmd;
	G4UIcmdWithADoubleAndUnit *phminCmd;
	G4UIcmdWithADoubleAndUnit *phmaxCmd;
	G4UIcmdWithADoubleAndUnit *thCoMminCmd;
	G4UIcmdWithADoubleAndUnit *thCoMmaxCmd;
	G4UIcmdWithADoubleAndUnit *EminCmd;
	G4UIcmdWithADoubleAndUnit *EmaxCmd;

	G4UIcmdWith3VectorAndUnit* vPosHCSCmd;
	G4UIcmdWith3VectorAndUnit* vPosTCSCmd;
	G4UIcmdWith3Vector* vMomHCSCmd;	
	G4UIcmdWith3Vector* vMomTCSCmd;	
};

#endif//g4hrsMessenger_HH























