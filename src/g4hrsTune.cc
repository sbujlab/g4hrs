#include "g4hrsTune.hh"

g4hrsTune* g4hrsTune::fTune=0;
g4hrsTune* g4hrsTune::GetTune() {
	if (!fTune) {
		fTune = new g4hrsTune();
	}
	return fTune;
}

g4hrsTune::g4hrsTune() {

	fTune = this;
	
	// This value should not be changed
	momentum_scale = 1.*GeV;


	// Default values
	kappaQuad1 = (-0.8476 * tesla) / (-3.7632047 * 1.063);
	kappaQuad2 = (0.93528 * tesla) / (-3.7632047 * 1.063);
	bendDipole = -0.4205 * tesla / 1.063;
	kappaQuad3 = (1.15762 * tesla) / (-3.7632047 * 1.063);
	septumAngle = 5.*deg;
	septumMomentum = 1.063*GeV;
	septumCurrent = 488.5*ampere;
	HRSAngle = 12.5*deg;
	HRSMomentum = 1.063*GeV;
	quadsOn = 1;
	sosQuad = 0;
}

g4hrsTune::~g4hrsTune() {
	//Nothing to see here
}

void g4hrsTune::SetTune(G4String mTune) {
	
	if(mTune == "PREX") {
		//Do nothing, this is the default tune
	} else if(mTune == "B") {

		kappaQuad1 = 0.140065 * tesla;
		kappaQuad2 = -0.212361 * tesla;
		bendDipole = -0.4042 * tesla;
		kappaQuad3 = -0.330965 * tesla;
		septumCurrent = 488.5*1.05*ampere;

	} else {
		G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		":  Invalid tune specified" << G4endl; 
		exit(1);
	}

}

double g4hrsTune::GetMomentumScale() {

	return momentum_scale;

}
