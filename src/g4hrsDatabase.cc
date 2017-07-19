#include "g4hrsDatabase.hh"
#include "G4SystemOfUnits.hh"

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
using std::cout;
using std::endl;

g4hrsDatabase::g4hrsDatabase(G4String targmat) {


	if(targmat == "Pb208") {
		E_min = 0.55*GeV;
		E_step = 0.05*GeV;
		n_E = 14;
		n_Th = 66;
		LoadTable("horpb.dat",0);
		LoadTable("horpb1.dat",1);
	}
	if(targmat == "Ca48") {
		E_min = 0.5*GeV;
		E_step = 0.05*GeV;
		n_E = 62;
		n_Th = 141;
		LoadTable("ca48_fsu.dat",0);
		LoadTable("ca48_fsu_stretched.dat",1);
	}

	if(targmat == "Ca40") {
		E_min = 0.5*GeV;
		E_step = 0.05*GeV;
		n_E = 62;
		n_Th = 141;
		LoadTable("ca40_fsu.dat",0);
		LoadTable("ca40_fsu_stretched.dat",1);
	}
// end constructer 
} 

g4hrsDatabase::~g4hrsDatabase()	{
}

void g4hrsDatabase::LoadTable(string filename, int stretch) {
	// stretch = {0, 1} = {R_n not stretched, R_n stretched 1%)
	
	ifstream datafile(filename, std::ifstream::in);

        if( !datafile ){
            fprintf(stderr, "ERROR %s line %d:  %s opening %s failed\n", __FILE__, __LINE__, __FUNCTION__, filename.c_str());
            exit(1);
        }
	
	double thisEnergy, thisAngle, thisXs, thisAsym, ignore;
	string dummy;

	for(int i=0; i<n_E; i++) {

		datafile >> dummy;
		thisEnergy = E_min + double(i)*E_step;
		energy.push_back(thisEnergy);

		vector<double> row_xs;
		vector<double> row_asym;

		if(filename.find("pb") != string::npos) {	
			for(int j=0; j<n_Th; j++) {
				datafile >> thisAngle >> thisXs >> ignore >> thisAsym >> ignore >> ignore;
				angle.push_back(thisAngle*deg);			
				row_xs.push_back(thisXs*millibarn);
				row_asym.push_back(thisAsym);
			} // end j (angle) for loop
		} else if(filename.find("ca") != string::npos) {
			for(int j=0; j<n_Th; j++) {
				datafile >> thisAngle >> thisXs >> thisAsym;
				angle.push_back(thisAngle*deg);			
				row_xs.push_back(thisXs*millibarn);
				row_asym.push_back(thisAsym);
			} // end j (angle) for loop
		}
	
		if(stretch == 0) {
			xs.push_back(row_xs);
			asym.push_back(row_asym);
		} 
		
		if(stretch == 1) {
			xs_str.push_back(row_xs);
			asym_str.push_back(row_asym);
		}

	} // end i (energy) for loop 

	cout << "Table " << filename << " loaded" << endl;

// end LoadTable
}

double g4hrsDatabase::Interpolate(double thisE, double thisTh, int stretch, int value) { 
	// stretch = {0, 1} = {R_n not stretched, R_n stretched 1%)
	// value = {0, 1} = {cross section, asymmetry}	
	double th0, th1, e0, e1;
	int i0, i1, j0, j1;

        // Check ranges
        if( thisTh < angle[0] ) return 0;
        if( thisTh > angle[angle.size()-1] ) return 0;
        if( thisE < energy[0] ) return 0;
        if( thisE > energy[energy.size()-1] ) return 0;


	for(int i=1; i<n_E; i++) {
		if(energy[i-1] < thisE && energy[i] > thisE) {
			e0 = energy[i-1];
			e1 = energy[i];
			i0 = i-1;
			i1 = i;
			break; 
		}  
	} 	
	for(int j=1; j<n_Th; j++) {
		if(angle[j-1] < thisTh && angle[j] > thisTh) {
			th0 = angle[j-1];
			th1 = angle[j];
			j0 = j-1;
			j1 = j;
			break;
		}
	}	
	double xs_e0_th0 = xs[i0][j0];
	double xs_e0_th1 = xs[i0][j1];
	double xs_e1_th0 = xs[i1][j0];
	double xs_e1_th1 = xs[i1][j1];
	double asym_e0_th0 = asym[i0][j0];
	double asym_e0_th1 = asym[i0][j1];
	double asym_e1_th0 = asym[i1][j0];
	double asym_e1_th1 = asym[i1][j1];
	double answer;
	if(value == 0) {
		double m_e0 = (xs_e0_th1 - xs_e0_th0)/(th1 - th0);
		double m_e1 = (xs_e1_th1 - xs_e1_th0)/(th1 - th0);
		double xs0 = xs_e0_th1 + m_e0*(thisTh - th0); 
		double xs1 = xs_e1_th1 + m_e1*(thisTh - th0);
		double m_th = (xs1 - xs0)/(e1 - e0);
		answer = xs0 + m_th*(thisE - e0); 
	}	
	if(value == 1) {
		double m_e0 = (asym_e0_th1 - asym_e0_th0)/(th1 - th0);
		double m_e1 = (asym_e1_th1 - asym_e1_th0)/(th1 - th0);
		double asym0 = asym_e0_th1 + m_e0*(thisTh - th0); 
		double asym1 = asym_e1_th1 + m_e1*(thisTh - th0);
		double m_th = (asym1 - asym0)/(e1 - e0);
		answer = asym0 + m_th*(thisE - e0); 
	}	
	
	return answer;

}












