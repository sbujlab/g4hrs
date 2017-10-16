#ifndef __G4HRSTRANSPORTFUNCTION_HH
#define __G4HRSTRANSPORTFUNCTION_HH

enum {sen, sm, sex, col, q1ex, q2ex, den, dex, q3en, q3m, q3ex, fp};

class g4hrsTransportFunction
{
public:
	g4hrsTransportFunction();
	~g4hrsTransportFunction();

	bool CallTransportFunction(float*, double*, double*, double*, double*);	

};

#endif
