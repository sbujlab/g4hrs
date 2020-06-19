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
	sosQuad = 1;
}

g4hrsTune::~g4hrsTune() {
	//Nothing to see here
}

void g4hrsTune::SetTune(G4String mTune) {
	
	if(mTune == "B") {
		// This is the tune with PREX-I target position
		kappaQuad1 = 0.140065 * tesla;
		kappaQuad2 = -0.212361 * tesla;
		bendDipole = -0.4042 * tesla;
		kappaQuad3 = -0.330965 * tesla;
		septumCurrent = 488.5*1.05*ampere;

	} else if(mTune == "PREXII") {
		// This tune is for PREX-II with a -5 cm target shift
		kappaQuad1 = 0.167238 * tesla;
		kappaQuad2 = -0.212596 * tesla;
		bendDipole = -0.4042 * tesla;
		kappaQuad3 = -0.333733 * tesla;
		septumCurrent = 497.6*ampere;

	} else if(mTune == "TuneP_CREX") {
		// This tune is for CREX with a -10 cm target shift
		// Septum current is 5% higher than for the PREX-II tune
		kappaQuad1 = 0.22453 * tesla;
		kappaQuad2 = -0.23284 * tesla;
		bendDipole = -0.394539 * tesla;
		kappaQuad3 = -0.264908 * tesla;
		septumCurrent = 505.683*1.05*ampere;

	} else if(mTune == "TuneP_PREX"){
                //RR Tune P-- needs to be worked on
                kappaQuad1 = 0.22453* tesla;
                kappaQuad2 = -0.23284 * tesla;
                bendDipole = -0.394539 * tesla;
                kappaQuad3 = -0.264908 * tesla;
                septumCurrent = 505.683*ampere;

        } else if(mTune == "TuneE"){  
          //From Weibin with CREX map 
          //Central angle is 4.85121, Septum scaled up 5%
               kappaQuad1 = 0.22553 * tesla;
               kappaQuad2 = -0.233232 * tesla;
               bendDipole = -0.394655 * tesla;
               kappaQuad3 = -0.264745 * tesla;
               septumCurrent = 492.489*1.05*ampere;

        } else if(mTune == "QuadsOff"){ 
              //No quad tune to test central ray
               kappaQuad1 = 0.22553 * tesla;
               kappaQuad2 = 0.0 * tesla;
               bendDipole = -0.394543 * tesla;
               kappaQuad3 = 0.0 * tesla;
               septumCurrent = 505.683 * ampere;
        } else {
		G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		":  Invalid tune specified" << G4endl; 
		exit(1);
	}

}

double g4hrsTune::GetMomentumScale() {

	return momentum_scale;

}
