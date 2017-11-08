
#ifndef __REMOLLSTEPPINGACTION_HH
#define __REMOLLSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class g4hrsTransportFunction; 

class g4hrsSteppingAction : public G4UserSteppingAction
{
  public:
    g4hrsSteppingAction();
    virtual ~g4hrsSteppingAction(){};

    virtual void UserSteppingAction(const G4Step*);

    void SetEnableKryptonite(G4bool k){ fEnableKryptonite = k; }

  private:
    G4bool drawFlag;
	G4double rad;
    G4bool fEnableKryptonite;
	g4hrsTransportFunction* fTransportFunction;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };

	int nelements;

	G4double fSeptumAngle;		// septum angle obtained from messenger (constant)
	G4double fHRSAngle;		// HRS angle obtained from messenger (constant)
	G4double septum_angle;		// local septum angle (changes sign depending on L/R HRS)
	G4double hrs_angle;		// local HRS angle (changes sign depending on L/R HRS)
	G4double fHRSMomentum; 		// HRS central momentum
	bool goodParticle;
	double sign; 	// y/phi sign flip for using transport function on left arm

	G4int fLHRS;
	G4int fRHRS;
		
	float r0[5];
	G4double x_tf[12];
	G4double y_tf[12];
	G4double th_tf[12];
	G4double ph_tf[12];

	G4double fX0;
	G4double fY0;
	G4double fZ0;
	G4double fP0;
	G4double fTh0;
	G4double fPh0;

	G4double fX0_tr;
	G4double fY0_tr;
	G4double fZ0_tr;
	G4double fP0_tr;
	G4double fTh0_tr;
	G4double fPh0_tr;
	
	G4double fX_sen;
	G4double fY_sen;
	G4double fZ_sen;
	G4double fTh_sen;
	G4double fPh_sen;

	G4double fX_sen_tr;
	G4double fY_sen_tr;
	G4double fZ_sen_tr;
	G4double fTh_sen_tr;
	G4double fPh_sen_tr;

	G4double fX_sen_tf;
	G4double fY_sen_tf;
	G4double fZ_sen_tf;
	G4double fTh_sen_tf;
	G4double fPh_sen_tf;
	
	G4double fX_sm;
	G4double fY_sm;
	G4double fZ_sm;
	G4double fTh_sm;
	G4double fPh_sm;

	G4double fX_sm_tr;
	G4double fY_sm_tr;
	G4double fZ_sm_tr;
	G4double fTh_sm_tr;
	G4double fPh_sm_tr;

	G4double fX_sm_tf;
	G4double fY_sm_tf;
	G4double fZ_sm_tf;
	G4double fTh_sm_tf;
	G4double fPh_sm_tf;
	
	G4double fX_sex;
	G4double fY_sex;
	G4double fZ_sex;
	G4double fTh_sex;
	G4double fPh_sex;

	G4double fX_sex_tr;
	G4double fY_sex_tr;
	G4double fZ_sex_tr;
	G4double fTh_sex_tr;
	G4double fPh_sex_tr;

	G4double fX_sex_tf;
	G4double fY_sex_tf;
	G4double fZ_sex_tf;
	G4double fTh_sex_tf;
	G4double fPh_sex_tf;
	
	G4double fX_col;
	G4double fY_col;
	G4double fZ_col;
	G4double fTh_col;
	G4double fPh_col;

	G4double fX_col_tr;
	G4double fY_col_tr;
	G4double fZ_col_tr;
	G4double fTh_col_tr;
	G4double fPh_col_tr;

	G4double fX_col_tf;
	G4double fY_col_tf;
	G4double fZ_col_tf;
	G4double fTh_col_tf;
	G4double fPh_col_tf;
	
	G4double fX_q1en;
	G4double fY_q1en;
	G4double fZ_q1en;
	G4double fTh_q1en;
	G4double fPh_q1en;

	G4double fX_q1en_tr;
	G4double fY_q1en_tr;
	G4double fZ_q1en_tr;
	G4double fTh_q1en_tr;
	G4double fPh_q1en_tr;

	G4double fX_q1ex;
	G4double fY_q1ex;
	G4double fZ_q1ex;
	G4double fTh_q1ex;
	G4double fPh_q1ex;

	G4double fX_q1ex_tr;
	G4double fY_q1ex_tr;
	G4double fZ_q1ex_tr;
	G4double fTh_q1ex_tr;
	G4double fPh_q1ex_tr;

	G4double fX_q1ex_tf;
	G4double fY_q1ex_tf;
	G4double fZ_q1ex_tf;
	G4double fTh_q1ex_tf;
	G4double fPh_q1ex_tf;
	
	G4double fX_q2en;
	G4double fY_q2en;
	G4double fZ_q2en;
	G4double fTh_q2en;
	G4double fPh_q2en;

	G4double fX_q2en_tr;
	G4double fY_q2en_tr;
	G4double fZ_q2en_tr;
	G4double fTh_q2en_tr;
	G4double fPh_q2en_tr;

	G4double fX_q2ex;
	G4double fY_q2ex;
	G4double fZ_q2ex;
	G4double fTh_q2ex;
	G4double fPh_q2ex;

	G4double fX_q2ex_tr;
	G4double fY_q2ex_tr;
	G4double fZ_q2ex_tr;
	G4double fTh_q2ex_tr;
	G4double fPh_q2ex_tr;

	G4double fX_q2ex_tf;
	G4double fY_q2ex_tf;
	G4double fZ_q2ex_tf;
	G4double fTh_q2ex_tf;
	G4double fPh_q2ex_tf;
	
	G4double fX_den;
	G4double fY_den;
	G4double fZ_den;
	G4double fTh_den;
	G4double fPh_den;

	G4double fX_den_tr;
	G4double fY_den_tr;
	G4double fZ_den_tr;
	G4double fTh_den_tr;
	G4double fPh_den_tr;

	G4double fX_den_tf;
	G4double fY_den_tf;
	G4double fZ_den_tf;
	G4double fTh_den_tf;
	G4double fPh_den_tf;
	
	G4double fX_dex;
	G4double fY_dex;
	G4double fZ_dex;
	G4double fTh_dex;
	G4double fPh_dex;

	G4double fX_dex_tr;
	G4double fY_dex_tr;
	G4double fZ_dex_tr;
	G4double fTh_dex_tr;
	G4double fPh_dex_tr;

	G4double fX_dex_tf;
	G4double fY_dex_tf;
	G4double fZ_dex_tf;
	G4double fTh_dex_tf;
	G4double fPh_dex_tf;
	
	G4double fX_q3en;
	G4double fY_q3en;
	G4double fZ_q3en;
	G4double fTh_q3en;
	G4double fPh_q3en;

	G4double fX_q3en_tr;
	G4double fY_q3en_tr;
	G4double fZ_q3en_tr;
	G4double fTh_q3en_tr;
	G4double fPh_q3en_tr;

	G4double fX_q3en_tf;
	G4double fY_q3en_tf;
	G4double fZ_q3en_tf;
	G4double fTh_q3en_tf;
	G4double fPh_q3en_tf;
	
	G4double fX_q3ex;
	G4double fY_q3ex;
	G4double fZ_q3ex;
	G4double fTh_q3ex;
	G4double fPh_q3ex;

	G4double fX_q3ex_tr;
	G4double fY_q3ex_tr;
	G4double fZ_q3ex_tr;
	G4double fTh_q3ex_tr;
	G4double fPh_q3ex_tr;

	G4double fX_q3ex_tf;
	G4double fY_q3ex_tf;
	G4double fZ_q3ex_tf;
	G4double fTh_q3ex_tf;
	G4double fPh_q3ex_tf;
	
	G4double fX_vdc;
	G4double fY_vdc;
	G4double fZ_vdc;
	G4double fTh_vdc;
	G4double fPh_vdc;

	G4double fX_vdc_tr;
	G4double fY_vdc_tr;
	G4double fZ_vdc_tr;
	G4double fTh_vdc_tr;
	G4double fPh_vdc_tr;

	G4double fX_fp;
	G4double fY_fp;
	G4double fZ_fp;
	G4double fTh_fp;
	G4double fPh_fp;

	G4double fX_fp_tr;
	G4double fY_fp_tr;
	G4double fZ_fp_tr;
	G4double fTh_fp_tr;
	G4double fPh_fp_tr;

	G4double fX_fp_tf;
	G4double fY_fp_tf;
	G4double fZ_fp_tf;
	G4double fTh_fp_tf;
	G4double fPh_fp_tf;

};

#endif//__REMOLLSTEPPINGACTION_HH
