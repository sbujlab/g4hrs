
#ifndef __REMOLLSTEPPINGACTION_HH
#define __REMOLLSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class hamcPREXTrans;

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
	hamcPREXTrans* fTransportFunction;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };

	int* np;

	G4double fSeptumAngle;
	G4double fHRSAngle;
	G4double fHRSMom;
		
	float r0[5];
	G4double x_tf[10];
	G4double y_tf[10];
	G4double z_tf[10];
	G4double th_tf[10];
	G4double ph_tf[10];

	G4double fX0;
	G4double fY0;
	G4double fZ0;
	G4double fP0;
	G4double fTh0;
	G4double fPh0;
	
	G4double fX0_tr;
	G4double fY0_tr;
	G4double fZ0_tr;
	G4double gP0_tr;
	G4double fTh0_tr;
	G4double fPh0_tr;

	G4double fX_sen;
	G4double fY_sen;
	G4double fZ_sen;
	G4double fP_sen;
	G4double fTheta_sen;
	G4double fPhi_sen;

	G4double fX_sen_tr;
	G4double fY_sen_tr;
	G4double fZ_sen_tr;
	G4double fP_sen_tr;
	G4double fTheta_sen_tr;
	G4double fPhi_sen_tr;

	G4double fX_sen_tf;
	G4double fY_sen_tf;
	G4double fZ_sen_tf;
	G4double fP_sen_tf;
	G4double fTheta_sen_tf;
	G4double fPhi_sen_tf;

	G4double fX_sm;
	G4double fY_sm;
	G4double fZ_sm;
	G4double fP_sm;
	G4double fTheta_sm;
	G4double fPhi_sm;

	G4double fX_sm_tr;
	G4double fY_sm_tr;
	G4double fZ_sm_tr;
	G4double fP_sm_tr;
	G4double fTheta_sm_tr;
	G4double fPhi_sm_tr;

	G4double fX_sm_tf;
	G4double fY_sm_tf;
	G4double fZ_sm_tf;
	G4double fP_sm_tf;
	G4double fTheta_sm_tf;
	G4double fPhi_sm_tf;


	G4double fX_sex;
	G4double fY_sex;
	G4double fZ_sex;
	G4double fP_sex;
	G4double fTheta_sex;
	G4double fPhi_sex;
	
	G4double fX_sex_tr;
	G4double fY_sex_tr;
	G4double fZ_sex_tr;
	G4double fP_sex_tr;
	G4double fTheta_sex_tr;
	G4double fPhi_sex_tr;

	G4double fX_sex_tf;
	G4double fY_sex_tf;
	G4double fZ_sex_tf;
	G4double fP_sex_tf;
	G4double fTheta_sex_tf;
	G4double fPhi_sex_tf;


	G4double fX_coil;
	G4double fY_coil;
	G4double fZ_coil;
	G4double fP_coil;
	G4double fTheta_coil;
	G4double fPhi_coil;

	G4double fX_mid;
	G4double fY_mid;
	G4double fZ_mid;
	G4double fP_mid;
	G4double fTheta_mid;
	G4double fPhi_mid;

	G4double fX_col;
	G4double fY_col;
	G4double fZ_col;
	G4double fP_col;
	G4double fTheta_col;
	G4double fPhi_col;

	G4double fX_q1en_L;
	G4double fY_q1en_L;
	G4double fZ_q1en_L;
	G4double fP_q1en_L;
	G4double fTheta_q1en_L;
	G4double fPhi_q1en_L;

	G4double fX_q1ex_L;
	G4double fY_q1ex_L;
	G4double fZ_q1ex_L;
	G4double fP_q1ex_L;
	G4double fTheta_q1ex_L;
	G4double fPhi_q1ex_L;
	
	G4double fX_q1en_R;
	G4double fY_q1en_R;
	G4double fZ_q1en_R;
	G4double fP_q1en_R;
	G4double fTheta_q1en_R;
	G4double fPhi_q1en_R;

	G4double fX_q1ex_R;
	G4double fY_q1ex_R;
	G4double fZ_q1ex_R;
	G4double fP_q1ex_R;
	G4double fTheta_q1ex_R;
	G4double fPhi_q1ex_R;
	
	G4double fX_q2en_L;
	G4double fY_q2en_L;
	G4double fZ_q2en_L;
	G4double fP_q2en_L;
	G4double fTheta_q2en_L;
	G4double fPhi_q2en_L;

	G4double fX_q2ex_L;
	G4double fY_q2ex_L;
	G4double fZ_q2ex_L;
	G4double fP_q2ex_L;
	G4double fTheta_q2ex_L;
	G4double fPhi_q2ex_L;
	
	G4double fX_q2en_R;
	G4double fY_q2en_R;
	G4double fZ_q2en_R;
	G4double fP_q2en_R;
	G4double fTheta_q2en_R;
	G4double fPhi_q2en_R;

	G4double fX_q2ex_R;
	G4double fY_q2ex_R;
	G4double fZ_q2ex_R;
	G4double fP_q2ex_R;
	G4double fTheta_q2ex_R;
	G4double fPhi_q2ex_R;

	G4double fX_den;
	G4double fY_den;
	G4double fZ_den;
	G4double fP_den;
	G4double fTheta_den;
	G4double fPhi_den;

	G4double fX_dex;
	G4double fY_dex;
	G4double fZ_dex;
	G4double fP_dex;
	G4double fTheta_dex;
	G4double fPhi_dex;

	G4double fX_q3en;
	G4double fY_q3en;
	G4double fZ_q3en;
	G4double fP_q3en;
	G4double fTheta_q3en;
	G4double fPhi_q3en;

	G4double fX_q3ex;
	G4double fY_q3ex;
	G4double fZ_q3ex;
	G4double fP_q3ex;
	G4double fTheta_q3ex;
	G4double fPhi_q3ex;

	G4double fX_vdc;
	G4double fY_vdc;
	G4double fZ_vdc;
	G4double fP_vdc;
	G4double fTheta_vdc;
	G4double fPhi_vdc;

	G4double fX_qz1;
	G4double fY_qz1;
	G4double fZ_qz1;
	G4double fP_qz1;
	G4double fTheta_qz1;
	G4double fPhi_qz1;

	G4double fX_qz2;
	G4double fY_qz2;
	G4double fZ_qz2;
	G4double fP_qz2;
	G4double fTheta_qz2;
	G4double fPhi_qz2;

	G4double fX_fp;
	G4double fY_fp;
	G4double fZ_fp;
	G4double fP_fp;
	G4double fTheta_fp;
	G4double fPhi_fp;

	G4double fX_par;
	G4double fY_par;
	G4double fZ_par;
	G4double fP_par;
	G4double fTheta_par;
	G4double fPhi_par;


	

};

#endif//__REMOLLSTEPPINGACTION_HH
