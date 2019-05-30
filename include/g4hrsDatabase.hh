#ifndef G4HRSDATABASE_HH
#define G4HRSDATABASE_HH

#include <string>
#include <iostream>
#include <vector>
#include "G4String.hh"
using std::string;
using std::vector;
using std::ifstream;

class g4hrsDatabase{

public:

	g4hrsDatabase(G4String);
	~g4hrsDatabase();
	double Interpolate(double,double,int,int);
        bool Table; //RR adding a flag here
	
private:

	int experiment;
	void LoadTable(string, int);
	vector< vector<double> > XS;
	vector< vector<double> > asym;
	vector< vector<double> > XS_str;
	vector< vector<double> > asym_str;
	vector<double> energy;
	vector<double> angle;
	int n_E;
	int n_Th;
	double E_min, E_step;
		
};

#endif
