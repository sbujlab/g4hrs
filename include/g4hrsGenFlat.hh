#ifndef __REMOLLGENFLAT_HH 
#define __REMOLLGENFLAT_HH 
/*!
 * Flat event generator
 *
 * Seamus Riordan
 * February 5, 2014
 *
*/

#include "g4hrsVEventGen.hh"

class g4hrsGenFlat : public g4hrsVEventGen {
    public:
	 g4hrsGenFlat();
	~g4hrsGenFlat();

	double fE_max;

    private:
	void SamplePhysics(g4hrsVertex *, g4hrsEvent *);


};

#endif//__REMOLLGENMOLLER_HH 
