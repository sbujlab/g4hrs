#include "g4hrsTransportFunction.hh"
#include "prex_forward.hh"
#include <cmath>
#include "G4ios.hh"

g4hrsTransportFunction::g4hrsTransportFunction() {
}

g4hrsTransportFunction::~g4hrsTransportFunction() {
}

bool g4hrsTransportFunction::CallTransportFunction(float* R0, double* x_tf, double* t_tf, double* y_tf, double* p_tf) {

	int nelement = 5;
	bool goodParticle = true;	
	
//	G4cout << "GEANT4 R0" << G4endl;
//	for(int i = 0; i<5; i++){G4cout << R0[i] << G4endl;}
	
	// [sen] = sen	
	x_tf[sen] = x_sp_sen_(R0, &nelement)/1000.;
	t_tf[sen] = atan(t_sp_sen_(R0, &nelement));
	y_tf[sen] = y_sp_sen_(R0, &nelement)/1000.;
	p_tf[sen] = atan(p_sp_sen_(R0, &nelement));
	
//	G4cout << "GEANT4 sen " << x_tf[sen] << "\n";
	
	// [sm] = sm
	x_tf[sm] = x_sp_sm_(R0, &nelement)/1000.;
	t_tf[sm] = atan(t_sp_sm_(R0, &nelement));
	y_tf[sm] = y_sp_sm_(R0, &nelement)/1000.;
	p_tf[sm] = atan(p_sp_sm_(R0, &nelement));
	
	// [sex] = sen
	x_tf[sex] = x_sp_sex_(R0, &nelement)/1000.;
	t_tf[sex] = atan(t_sp_sex_(R0, &nelement));
	y_tf[sex] = y_sp_sex_(R0, &nelement)/1000.;
	p_tf[sex] = atan(p_sp_sex_(R0, &nelement));

	// [col] = col
	x_tf[col] = x_sp_col_(R0, &nelement)/1000.;
	t_tf[col] = atan(t_sp_col_(R0, &nelement));
	y_tf[col] = y_sp_col_(R0, &nelement)/1000.;
	p_tf[col] = atan(p_sp_col_(R0, &nelement));

	// [q1ex] = q1ex
	x_tf[q1ex] = x_sp_q1ex_(R0, &nelement)/1000.;
	t_tf[q1ex] = atan(t_sp_q1ex_(R0, &nelement));
	y_tf[q1ex] = y_sp_q1ex_(R0, &nelement)/1000.;
	p_tf[q1ex] = atan(p_sp_q1ex_(R0, &nelement));

	// [q2ex] = q2ex
	x_tf[q2ex] = x_sp_q2ex_(R0, &nelement)/1000.;
	t_tf[q2ex] = atan(t_sp_q2ex_(R0, &nelement));
	y_tf[q2ex] = y_sp_q2ex_(R0, &nelement)/1000.;
	p_tf[q2ex] = atan(p_sp_q2ex_(R0, &nelement));
	
	// [den] = den
	x_tf[den] = x_sp_den_(R0, &nelement)/1000.;
	t_tf[den] = atan(t_sp_den_(R0, &nelement));
	y_tf[den] = y_sp_den_(R0, &nelement)/1000.;
	p_tf[den] = atan(p_sp_den_(R0, &nelement));

	// [dex] = dex
	x_tf[dex] = x_sp_dex_(R0, &nelement)/1000.;
	t_tf[dex] = atan(t_sp_dex_(R0, &nelement));
	y_tf[dex] = y_sp_dex_(R0, &nelement)/1000.;
	p_tf[dex] = atan(p_sp_dex_(R0, &nelement));

	// [q3en] = q3en
	x_tf[q3en] = x_sp_q3en_(R0, &nelement)/1000.;
	t_tf[q3en] = atan(t_sp_q3en_(R0, &nelement));
	y_tf[q3en] = y_sp_q3en_(R0, &nelement)/1000.;
	p_tf[q3en] = atan(p_sp_q3en_(R0, &nelement));

	// [q3m] = q3m
	x_tf[q3m] = x_sp_q3m_(R0, &nelement)/1000.;
	t_tf[q3m] = atan(t_sp_q3m_(R0, &nelement));
	y_tf[q3m] = y_sp_q3m_(R0, &nelement)/1000.;
	p_tf[q3m] = atan(p_sp_q3m_(R0, &nelement));

	// [q3ex] = q3ex
	x_tf[q3ex] = x_sp_q3ex_(R0, &nelement)/1000.;
	t_tf[q3ex] = atan(t_sp_q3ex_(R0, &nelement));
	y_tf[q3ex] = y_sp_q3ex_(R0, &nelement)/1000.;
	p_tf[q3ex] = atan(p_sp_q3ex_(R0, &nelement));

	// [fp] = fp
	x_tf[fp] = x_sp_fp_(R0, &nelement)/1000.;
	t_tf[fp] = atan(t_sp_fp_(R0, &nelement));
	y_tf[fp] = y_sp_fp_(R0, &nelement)/1000.;
	p_tf[fp] = atan(p_sp_fp_(R0, &nelement));


	// Check apertures
/*	
	if( pow(x_tf[q1ex],2.) + pow(y_tf[q1ex],2.) > pow(0.1492,2.)) {
		goodParticle = false;
	}
	if( pow(x_tf[q2ex],2.) + pow(y_tf[q2ex],2.) > pow(0.3,2.)) {
		goodParticle = false;
	}
	if( pow(x_tf[q3en],2.) + pow(y_tf[q3en],2.) > pow(0.3,2.)) {
		goodParticle = false;
	}
*/
	return goodParticle;
	
}

























