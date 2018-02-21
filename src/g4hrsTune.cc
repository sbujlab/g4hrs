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
	
	// Default values
	snakeMagnet = -3.76320471828;
	kappaQuad1 = -0.8476 * tesla / snakeMagnet;
	kappaQuad2 = 0.8680 * tesla / snakeMagnet;
	bendDipole = -0.4192 * tesla;
	kappaQuad3 = 1.1748 * tesla / snakeMagnet;
	septumAngle = 5.*deg;
	septumMomentum = 1.063*GeV;
	septumCurrent = 528*ampere;
	HRSAngle = 12.5*deg;
	HRSMomentum = 1.063*GeV;

	//Copy+pasted from EMFieldSetup, in case these values are needed down the road	
/*
	//snakeMagnet = 1. / -4.77577734;//unitless, STD - correct, I believe
	//snakeMagnet = 1. / -4.77577734 * ( 4.00 / 1.063 ) ;//unitless, PREX tune B - correct, I believe
	G4double snakeMagnet = -4.77577734 / 1.063;//Nickie's calculation//this is just 4.77577734 / 1.063 //THIS CURRENTLY BEING USED, TTK 01/16/2018
	//snakeMagnet *= 0.83756;
	snakeMagnet *= 0.83762;  //THIS CURRENTLY BEING USED, TTK 01/16/2018
	//snakeMagnet *= 0.7;//This is tuned with hrstrans #1
	//snakeMagnet *= 1.33;//This is tuned with hrstrans when hrstrans is tuned to JLR
	//snakeMagnet *= 1.24;//This is my new tune, with snake d.dat problems fixed #2
	//snakeMagnet *= 0.965; #3
*/

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

void g4hrsTune::SetQ1(double k1) {
	kappaQuad1 = k1/snakeMagnet;
}
 
void g4hrsTune::SetQ2(double k2) {
	kappaQuad2 = k2/snakeMagnet;
}
 
void g4hrsTune::SetQ3(double k3) {
	kappaQuad3 = k3/snakeMagnet;
} 
