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

	
	momentum_scale = 1.*GeV;
	// Default values
	kappaQuad1 = -0.8476 * tesla;
	kappaQuad2 = 0.8680 * tesla;
	bendDipole = -0.4192 * tesla;
	kappaQuad3 = 1.1748 * tesla;
	septumAngle = 5.*deg;
	septumMomentum = 1.063*GeV;
	septumCurrent = 528*ampere;
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

		// Put Tune B here

	} else {
		G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		":  Invalid tune specified" << G4endl; 
		exit(1);
	}

}

double g4hrsTune::GetMomentumScale() {

	return momentum_scale;

}
