#ifndef __REMOLLGENBEAM_HH 
#define __REMOLLGENBEAM_HH 
/*!
 * Boring beam event generator
 *
 * Seamus Riordan
 * July 9, 2013
 *
*/

#include "g4hrsVEventGen.hh"

class g4hrsBeamTarget;

class g4hrsGenBeam : public g4hrsVEventGen {
    public:
	g4hrsGenBeam();
	~g4hrsGenBeam();

    private:
	void SamplePhysics(g4hrsVertex *, g4hrsEvent *);

	g4hrsBeamTarget *fBeamTarg;

	double fZpos;
};

#endif//__REMOLLGENBEAM_HH 
