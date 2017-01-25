#ifndef __REMOLLRUNDATA_HH
#define __REMOLLRUNDATA_HH

#include "TObject.h"

#include <vector>
#include <string>

#include "g4hrstypes.hh"
#include "g4hrsTextFile.hh"

/*!
 * All the information on the run
 * This will get put into the output
 * stream
*/

class TGeoManager;

class g4hrsRunData : public TObject {
  using TObject::Print;
    public:
	 g4hrsRunData();
	~g4hrsRunData();

	unsigned long long int GetNthrown(){ return fNthrown; }
	void SetNthrown(unsigned long long int n){ fNthrown = n; }

	void Init();

	void SetGenName(const char *n){ strcpy(fGenName, n); }
	const char *GetGenName(){ return fGenName; }

	void SetBeamE(double E){ fBeamE = E; }
	void SetSeed(unsigned int seed){ fSeed = seed; }

	void AddMagData(filedata_t d){fMagData.push_back(d);}
	void SetMacroFile(const char *fn){ fMacro = g4hrsTextFile(fn); }
	void AddGDMLFile(const char *fn);
	void ClearGDMLFiles(){ fGDMLFiles.clear(); }

	void RecreateGDML(const char *adir = NULL, bool clobber = false);

	g4hrsTextFile GetGDMLFile(int i){ return fGDMLFiles[i]; }

	void Print();

	TTimeStamp fRunTime;

	long int  fNthrown;
	unsigned int  fSeed;
	double fBeamE;
	char fGenName[__RUNSTR_LEN];

	char fHostName[__RUNSTR_LEN];
	char fRunPath[__RUNSTR_LEN];

	g4hrsTextFile              fMacro;
	std::vector<g4hrsTextFile> fGDMLFiles;

	std::vector<filedata_t> fMagData;

	ClassDef(g4hrsRunData, 1);
};

#endif//__REMOLLRUNDATA_HH
