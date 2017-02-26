#include "g4hrsMessenger.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

#include "g4hrsDetectorConstruction.hh"
#include "g4hrsIO.hh"
#include "g4hrsEventAction.hh"
#include "g4hrsVEventGen.hh"
#include "g4hrsPrimaryGeneratorAction.hh"
#include "g4hrsBeamTarget.hh"
#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"
#include "g4hrsSteppingAction.hh"
#include "g4hrsGenFlat.hh"

#include "G4UImanager.hh"
#include "G4RunManager.hh"

#include "G4GDMLParser.hh"
#include "G4VPhysicalVolume.hh"

#include <iostream>

g4hrsMessenger::g4hrsMessenger(){
    /*  Initialize all the things it talks to to NULL */

    fIO           = NULL;
    fdetcon       = NULL;
    fevact        = NULL;
    fprigen       = NULL;
    fBeamTarg     = NULL;
    fStepAct      = NULL;
    fPhysicsList  = NULL;

    // Grab singleton beam/target
    fBeamTarg = g4hrsBeamTarget::GetBeamTarget();

    detfilesCmd = new G4UIcmdWithAString("/g4hrs/setgeofile",this);
    detfilesCmd->SetGuidance("Set geometry GDML files");
    detfilesCmd->SetParameterName("geofilename", false);
    detfilesCmd->AvailableForStates(G4State_PreInit); // Only have this in pre-init or GDML freaks out

    seedCmd = new G4UIcmdWithAnInteger("/g4hrs/seed",this);
    seedCmd->SetGuidance("Set random engine seed");
    seedCmd->SetParameterName("seed", false);

    kryptCmd = new G4UIcmdWithABool("/g4hrs/kryptonite",this);
    kryptCmd->SetGuidance("Treat W, Pb, Cu as kryptonite");
    kryptCmd->SetParameterName("krypt", false);

    opticalCmd = new G4UIcmdWithABool("/g4hrs/optical",this);
    opticalCmd->SetGuidance("Enable optical physics");
    opticalCmd->SetParameterName("optical", false);
    opticalCmd->AvailableForStates(G4State_Idle); // Only have this AFTER we've initalized geometry

    dumpGeometryCmd = new G4UIcmdWithABool("/g4hrs/dumpgeometry",this);
    dumpGeometryCmd->SetGuidance("Dump the geometry tree");
    dumpGeometryCmd->SetParameterName("overlap_check", true);
    dumpGeometryCmd->SetDefaultValue(false);
    dumpGeometryCmd->AvailableForStates(G4State_Idle); // Only have this AFTER we've initalized geometry

    newfieldCmd = new G4UIcmdWithAString("/g4hrs/addfield",this);
    newfieldCmd->SetGuidance("Add magnetic field");
    newfieldCmd->SetParameterName("filename", false);

    fieldScaleCmd = new G4UIcmdWithAString("/g4hrs/scalefield",this);
    fieldScaleCmd->SetGuidance("Scale magnetic field");
    fieldScaleCmd->SetParameterName("filename", false);

    fieldCurrCmd = new G4UIcmdWithAString("/g4hrs/magcurrent",this);
    fieldCurrCmd->SetGuidance("Scale magnetic field by current");
    fieldCurrCmd->SetParameterName("filename", false);

    tgtLenCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/targlen",this);
    tgtLenCmd->SetGuidance("Target length");
    tgtLenCmd->SetParameterName("targlen", false);
    tgtLenCmd->AvailableForStates(G4State_Idle); // Only have this AFTER we've initialized geometry

    tgtPosCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/targpos",this);
    tgtPosCmd->SetGuidance("Target position");
    tgtPosCmd->SetParameterName("targpos", false);
    tgtPosCmd->AvailableForStates(G4State_Idle); // Only have this AFTER we've initialized geometry

    tgtMatCmd = new G4UIcmdWithAString("/g4hrs/targmat", this);
    tgtMatCmd->SetGuidance("Target material");
    tgtMatCmd->SetParameterName("targmat", false);
    tgtMatCmd->AvailableForStates(G4State_Idle); // Only have this AFTER we've initialized geometry

    beamCurrCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beamcurr",this);
    beamCurrCmd->SetGuidance("Beam current");
    beamCurrCmd->SetParameterName("beamcurr", false);

    beamECmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beamene",this);
    beamECmd->SetGuidance("Beam energy");
    beamECmd->SetParameterName("beamene", false);

    genSelectCmd = new G4UIcmdWithAString("/g4hrs/gen",this);
    genSelectCmd->SetGuidance("Select physics generator");
    genSelectCmd->SetParameterName("generator", false);

    fileCmd = new G4UIcmdWithAString("/g4hrs/filename",this);
    fileCmd->SetGuidance("Output filename");
    fileCmd->SetParameterName("filename", false);


    thminCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/thmin",this);
    thminCmd->SetGuidance("Minimum generation angle");
    thminCmd->SetParameterName("thmin", false);

    thmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/thmax",this);
    thmaxCmd->SetGuidance("Minimum generation angle");
    thmaxCmd->SetParameterName("thmax", false);

    thCoMminCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/thcommin",this);
    thCoMminCmd->SetGuidance("Minimum CoM generation angle");
    thCoMminCmd->SetParameterName("thcommin", false);

    thCoMmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/thcommax",this);
    thCoMmaxCmd->SetGuidance("Minimum CoM generation angle");
    thCoMmaxCmd->SetParameterName("thcommax", false);

    EminCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/emin",this);
    EminCmd->SetGuidance("Minimum generation energy");
    EminCmd->SetParameterName("emin", false);

    EmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/emax",this);
    EmaxCmd->SetGuidance("Maximum generation energy");
    EmaxCmd->SetParameterName("emax", false);


    //////////////////////////////////////////////////
    // beam info

    rasTypeCmd = new G4UIcmdWithABool("/g4hrs/oldras", this);
    rasTypeCmd->SetGuidance("Old (no ang corln) or new (ang corl) raster");
    rasTypeCmd->SetParameterName("oldras", true);

    rasXCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/rasx",this);
    rasXCmd->SetGuidance("Square raster width in x (horizontal)");
    rasXCmd->SetParameterName("rasx", false);

    rasYCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/rasy",this);
    rasYCmd->SetGuidance("Square raster width y (vertical)");
    rasYCmd->SetParameterName("rasy", false);

    beamX0Cmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beam_x0",this);
    beamX0Cmd->SetGuidance("beam initial position in x (horizontal)");
    beamX0Cmd->SetParameterName("beamX0", false);

    beamY0Cmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beam_y0",this);
    beamY0Cmd->SetGuidance("beam initial position in y (vertical)");
    beamY0Cmd->SetParameterName("beamY0", false);

    beamth0Cmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beam_th0",this);
    beamth0Cmd->SetGuidance("beam initial direction in x (horizontal)");
    beamth0Cmd->SetParameterName("beamth0", false);

    beamph0Cmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beam_ph0",this);
    beamph0Cmd->SetGuidance("beam initial direction in y (vertical)");
    beamph0Cmd->SetParameterName("beamph0", false);

    beamCorrThCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beam_corrth",this);
    beamCorrThCmd->SetGuidance("beam correlated angle (horizontal)");
    beamCorrThCmd->SetParameterName("beam_corrth", false);

    beamCorrPhCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beam_corrph",this);
    beamCorrPhCmd->SetGuidance("beam correlated angle (vertical)");
    beamCorrPhCmd->SetParameterName("beam_corrph", false);

    beamdthCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beam_dth",this);
    beamdthCmd->SetGuidance("beam gaussian spread in direction x (horizontal)");
    beamdthCmd->SetParameterName("beamdth", false);

    beamdphCmd = new G4UIcmdWithADoubleAndUnit("/g4hrs/beam_dph",this);
    beamdphCmd->SetGuidance("beam gaussian spread in direction y (vertical)");
    beamdphCmd->SetParameterName("beamdph", false);



    /*
       fExpType = kNeutronExp;

       runCmd = new G4UIcmdWithAnInteger("/g4sbs/run",this);
       runCmd->SetGuidance("Run simulation with x events");
       runCmd->SetParameterName("nevt", false);

       gemconfigCmd = new G4UIcmdWithAnInteger("/g4sbs/gemconfig",this);
       gemconfigCmd->SetGuidance("Change between GEM configurations");
       gemconfigCmd->SetParameterName("gemconfig", false);


       sigfileCmd = new G4UIcmdWithAString("/g4sbs/sigmafile",this);
       sigfileCmd->SetGuidance("GEM Sigma filename");
       sigfileCmd->SetParameterName("sigmafile", false);

       tgtCmd = new G4UIcmdWithAString("/g4sbs/target",this);
       tgtCmd->SetGuidance("Target type from LH2, LD2, H2, 3He");
       tgtCmd->SetParameterName("targtype", false);

       kineCmd = new G4UIcmdWithAString("/g4sbs/kine",this);
       kineCmd->SetGuidance("Kinematic type");
       kineCmd->SetParameterName("kinetype", false);

       expCmd = new G4UIcmdWithAString("/g4sbs/exp",this);
       expCmd->SetGuidance("Experiment type");
       expCmd->SetParameterName("exptype", false);

       geantinoCmd = new G4UIcmdWithABool("/g4sbs/shootgeantino", this);
       geantinoCmd->SetGuidance("Shoot a geantino instead of e-");
       geantinoCmd->SetParameterName("shootgeantino", false);


       tgtDenCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targden",this);
       tgtDenCmd->SetGuidance("Target density");
       tgtDenCmd->SetParameterName("targden", false);

       tgtPresCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targpres",this);
       tgtPresCmd->SetGuidance("Gaseous Target pressure");
       tgtPresCmd->SetParameterName("targpres", false);

       beamcurCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamcur",this);
       beamcurCmd->SetGuidance("Beam current");
       beamcurCmd->SetParameterName("beamcur", false);

       runtimeCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/runtime",this);
       runtimeCmd->SetGuidance("Run time");
       runtimeCmd->SetParameterName("runtime", false);

       rasterxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/rasterx",this);
       rasterxCmd->SetGuidance("Raster x size");
       rasterxCmd->SetParameterName("size", false);

       rasteryCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/rastery",this);
       rasteryCmd->SetGuidance("Raster y size");
       rasteryCmd->SetParameterName("size", false);

       beamECmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamE",this);
       beamECmd->SetGuidance("Beam Energy");
       beamECmd->SetParameterName("energy", false);

       bbangCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/bbang",this);
       bbangCmd->SetGuidance("BigBite angle");
       bbangCmd->SetParameterName("angle", false);

       bbdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/bbdist",this);
       bbdistCmd->SetGuidance("BigBite distance");
       bbdistCmd->SetParameterName("dist", false);

       hcalangCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcalang",this);
       hcalangCmd->SetGuidance("HCAL angle");
    hcalangCmd->SetParameterName("angle", false);

    hcaldistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcaldist",this);
    hcaldistCmd->SetGuidance("HCAL distance");
    hcaldistCmd->SetParameterName("dist", false);

    hmagdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/48D48dist",this);
    hmagdistCmd->SetGuidance("48D48 distance");
    hmagdistCmd->SetParameterName("dist", false);

    gemresCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/gemres",this);
    gemresCmd->SetGuidance("GEM resolution");
    gemresCmd->SetParameterName("dist", false);

    // Detector position commands

    cerDisCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/cerdist",this);
    cerDisCmd->SetGuidance("Cerenkov distance from front GEM");
    cerDisCmd->SetParameterName("dist", false);

    cerDepCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/cerdepth",this);
    cerDepCmd->SetGuidance("Cerenkov gas depth");
    cerDepCmd->SetParameterName("dist", false);

    gemSepCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/gemsep",this);
    gemSepCmd->SetGuidance("GEM separation from front to back set");
    gemSepCmd->SetParameterName("dist", false);

    bbCalDistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/bbcaldist",this);
    bbCalDistCmd->SetGuidance("BigBite caloriter distance from front GEM");
    bbCalDistCmd->SetParameterName("dist", false);

    thminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/thmin",this);
    thminCmd->SetGuidance("Minimum electron generation polar angle");
    thminCmd->SetParameterName("angle", false);

    thmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/thmax",this);
    thmaxCmd->SetGuidance("Maximum electron generation polar angle");
    thmaxCmd->SetParameterName("angle", false);

    phminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/phmin",this);
    phminCmd->SetGuidance("Minimum electron generation azimuthal angle");
    phminCmd->SetParameterName("angle", false);

    phmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/phmax",this);
    phmaxCmd->SetGuidance("Maximum electron generation azimuthal angle");
    phmaxCmd->SetParameterName("angle", false);
    */

}

g4hrsMessenger::~g4hrsMessenger(){
}


void g4hrsMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue){
    if( cmd == seedCmd ){
	G4int seed = seedCmd->GetNewIntValue(newValue);
	CLHEP::HepRandom::setTheSeed(seed);
	g4hrsRun::GetRun()->GetData()->SetSeed(seed);
    }

    if( cmd == kryptCmd ){
	G4bool krypt = kryptCmd->GetNewBoolValue(newValue);
	fStepAct->SetEnableKryptonite(krypt);
    }


    if( cmd == newfieldCmd ){
//	fField->AddNewField( newValue );
    }

    if( cmd == fieldScaleCmd ){
	std::istringstream iss(newValue);

	G4String scalefile, scalestr;
	G4double scaleval;

	iss >> scalefile;
	iss >> scalestr;

	scaleval = atof(scalestr.data());
//	fField->SetFieldScale( scalefile, scaleval );
    }

    if( cmd == fieldCurrCmd ){
	std::istringstream iss(newValue);

	G4String scalefile, scalestr, scaleunit;
	G4double scaleval;

	iss >> scalefile;
	iss >> scalestr;
	iss >> scaleunit;

	if( scaleunit != "A" ){
	    // FIXME: less snark and more functionality?
	    G4cerr << __FILE__ << " line " << __LINE__ <<  ":\n\tGraaaah - just put the current for " <<  scalefile <<  " in amps..." << G4endl;
	    exit(1);
	}

	scaleval = atof(scalestr.data());
//	fField->SetMagnetCurrent( scalefile, scaleval );
    }

    if( cmd == tgtLenCmd ){
	G4double len = tgtLenCmd->GetNewDoubleValue(newValue);
	fBeamTarg->SetTargetLen(len);
    }

    if( cmd == tgtPosCmd ){
	G4double pos = tgtPosCmd->GetNewDoubleValue(newValue);
	fBeamTarg->SetTargetPos(pos);
    }

    if ( cmd == tgtMatCmd ){
	fBeamTarg->fTargetMaterial = newValue;
    }

    if( cmd == genSelectCmd ){
	fprigen->SetGenerator( newValue );
    }

    if( cmd == beamCurrCmd ){
	G4double cur = beamCurrCmd->GetNewDoubleValue(newValue);
	fBeamTarg->SetBeamCurrent(cur);
    }

    if( cmd == beamECmd ){
	G4double ene = beamECmd->GetNewDoubleValue(newValue);
	fBeamTarg->fBeamE = ene;
    }

    if( cmd == fileCmd ){
	fIO->SetFilename(newValue);
    }
    

    if( cmd == EminCmd ){
	G4double en = EminCmd->GetNewDoubleValue(newValue);
	g4hrsVEventGen *agen = fprigen->GetGenerator();
	if( agen ){
	    agen->fE_min = en;
	}
    }

    if( cmd == EmaxCmd ){
	G4double en = EmaxCmd->GetNewDoubleValue(newValue);
	g4hrsGenFlat *agen = dynamic_cast<g4hrsGenFlat *>(fprigen->GetGenerator());
	if( agen ){
	    agen->fE_max = en;
	}
    }

    if( cmd == thminCmd ){
	G4double th = thminCmd->GetNewDoubleValue(newValue);
	g4hrsVEventGen *agen = fprigen->GetGenerator();
	if( agen ){
	    agen->fTh_min = th;
	}
    }

    if( cmd == thmaxCmd ){
	G4double th = thminCmd->GetNewDoubleValue(newValue);
	g4hrsVEventGen *agen = fprigen->GetGenerator();
	if( agen ){
	    agen->fTh_max = th;
	}
    }

    if( cmd == thCoMminCmd ){
	G4double th = thCoMminCmd->GetNewDoubleValue(newValue);
	g4hrsVEventGen *agen = fprigen->GetGenerator();
	if( agen ){
	    agen->fThCoM_min = th;
	}
    }

    if( cmd == thCoMmaxCmd ){
	G4double th = thCoMminCmd->GetNewDoubleValue(newValue);
	g4hrsVEventGen *agen = fprigen->GetGenerator();
	if( agen ){
	    agen->fThCoM_max = th;
	}
    }

    if( cmd == rasTypeCmd ){
	G4bool rasType = rasTypeCmd->GetNewBoolValue(newValue);
	fBeamTarg->fOldRaster = rasType;
    }

    if( cmd == rasXCmd ){
	G4double x = rasXCmd->GetNewDoubleValue(newValue);
	fBeamTarg->fRasterX = x;
    }

    if( cmd == rasYCmd ){
	G4double y = rasYCmd->GetNewDoubleValue(newValue);
	fBeamTarg->fRasterY = y;
    }

    if( cmd == beamX0Cmd ){
	G4double x = beamX0Cmd->GetNewDoubleValue(newValue);
	fBeamTarg->fX0 = x;
    }

    if( cmd == beamY0Cmd ){
	G4double y = beamY0Cmd->GetNewDoubleValue(newValue);
	fBeamTarg->fY0 = y;
    }

    if( cmd == beamth0Cmd ){
	G4double x = beamth0Cmd->GetNewDoubleValue(newValue);
	fBeamTarg->fTh0 = x;
    }

    if( cmd == beamph0Cmd ){
	G4double y = beamph0Cmd->GetNewDoubleValue(newValue);
	fBeamTarg->fPh0 = y;
    }

    if( cmd == beamCorrThCmd ){
	G4double x = beamCorrThCmd->GetNewDoubleValue(newValue);
	fBeamTarg->fCorrTh = tan(x);
    }

    if( cmd == beamCorrPhCmd ){
	G4double y = beamCorrPhCmd->GetNewDoubleValue(newValue);
	fBeamTarg->fCorrPh = tan(y);
    }

    if( cmd == beamdthCmd ){
	G4double x = beamdthCmd->GetNewDoubleValue(newValue);
	fBeamTarg->fdTh = x;
    }

    if( cmd == beamdphCmd ){
	G4double y = beamdphCmd->GetNewDoubleValue(newValue);
	fBeamTarg->fdPh = y;
    }
}
