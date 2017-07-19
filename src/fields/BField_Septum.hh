// BField_Septum.h: interface for the BField_Septum class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(BFIELD_Septum_H)
#define BFIELD_Septum_H
#include <math.h>
#include <CLHEP/Vector/Rotation.h>
#include <CLHEP/Vector/ThreeVector.h>
using namespace CLHEP;

class BField_Septum
{
public:
	//Static method which returns the singleton pointer of this class.
	static BField_Septum* GetInstance();
private:
	static BField_Septum* fInstance;

public:
	BField_Septum(double pMomentum,
		const char *mapfile);
	virtual ~BField_Septum();
	bool GetBField(double Pos[],double B[]);
	bool GetBField(float fPos[],float fB[]);
	
	void SetMomentum(double pMomentum){fMomentum = pMomentum;}

private:
	bool ReadIni(const char *filename);
	bool ReadMap(const char *filename);
	bool Interpolation(double Pos[3],double B[3],int n=2);

public:

	void   Rotate_Field2Lab(const double FieldP[3],double LabP[3]);
	void   Transform_Field2Lab(const double FieldP[3],double LabP[3]);
	void   Rotate_Lab2Field(const double LabP[3],double FieldP[3]);
	void   Transform_Lab2Field(const double LabP[3],double FieldP[3]);
	bool   IsUniformField(){ return (mUseUniformB==0)?false:true;}
	void   GetUniformField(double pB[]){pB[0]=mUniformB[0];pB[1]=mUniformB[1];pB[2]=mUniformB[2];}
	void   GetOrigin(double pX[]){pX[0]=mOrigin[0];pX[1]=mOrigin[1];pX[2]=mOrigin[2];}

	CLHEP::HepRotation *GetRotation_L2F() {return mRotL2F;}
	void GetEulerAngles_L2F(double &phi, double &theta, double &psi) 
	{		
		phi=mRotL2F->getPhi();theta=mRotL2F->getTheta();psi=mRotL2F->getPsi(); 
	}

	//Create root ntuple: if verbose != 0 will print steps out 
	void CreateNtuple(const char *rootfile="septumfield.root",int verbose=0,
		double xmin=-43.0,double xmax=43.0, double xstep=1.0,
		double ymin=-12.0,double ymax=12.0, double ystep=1.0,
		double zmin=-100.0,double zmax=100.0, double zstep=1.0);

private:
	//dynalic allocated the 4 dimensional array: double mBField[indexX][indexY][indexZ][Variables];
	//mBField[indexX][indexY][indexZ][0] x
	//mBField[indexX][indexY][indexZ][1] y
	//mBField[indexX][indexY][indexZ][2] z
	//mBField[indexX][indexY][indexZ][3] Bx
	//mBField[indexX][indexY][indexZ][4] By
	//mBField[indexX][indexY][indexZ][5] Bz
	double ****mBField;

	//parameters from ini file
	int    mUseUniformB;
	double mUniformB[3];
	double mStepX,mStepY,mStepZ;
	double mXmin,mXmax;
	double mYmin,mYmax;
	double mZmin,mZmax;
	int	   mInterpolateOutOfRange;
	double mOrigin[3];	// shift from hall center to septum center
	double mTarget[3];	// shift from target center to hall center

	double mFieldUnit;
	int    mFirstDataLine;
	int    mNPara;
	int	   mRotAxis[3];
	double mRotAngle[3];

	double mDefaultMomentum;
	double fMomentum;

	bool   mDoShift, mDoRotation;
	int    mZNum,mXNum,mYNum;
	double mStepR;
	CLHEP::HepRotation *mRotL2F;
	CLHEP::HepRotation *mRotF2L;

};
typedef BField_Septum HRSSeptumField;
typedef BField_Septum G2PSeptumField;

#endif // !defined(BFIELD_Septum_H)
