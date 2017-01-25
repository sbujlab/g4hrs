#ifndef __REMOLLTEXTFILE_HH
#define __REMOLLTEXTFILE_HH

#define __STRLEN 1024

#include "TObject.h"

class g4hrsTextFile : public TObject {
  using TObject::Print;
     public:       
	 g4hrsTextFile();
	 g4hrsTextFile(const g4hrsTextFile &);
	 const g4hrsTextFile& operator=(const g4hrsTextFile &);
	 g4hrsTextFile(const char *);
	~g4hrsTextFile();

	 void copyFileIn(const char *);

	void Print();

	const char *GetFilename(){ return fFilename; }
	unsigned long long int GetBufferSize(){ return fBufferSize; }
	
	void Recreate(const char *fn = NULL, bool clobber = false);
	void RecreateInDir(const char *path, bool clobber = false);

    private:
	int fFilenameSize;
	char *fFilename;

	unsigned long long int fBufferSize;
	char *fBuffer;

	const char *GetBaseFile(const char *fp = NULL);

	ClassDef(g4hrsTextFile, 1);
};

#endif//__REMOLLTEXTFILE_HH
