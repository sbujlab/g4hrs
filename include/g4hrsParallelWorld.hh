#ifndef __G4HRSPARALLELWORLD_HH
#define __G4HRSPARALLELWORLD_HH

#include "G4VUserParallelWorld.hh"
#include <vector>

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSensitiveDetector;

class g4hrsIO;
class g4hrsEMFieldSetup;
class g4hrsEMField;
class g4hrsMaterial;
class g4hrsTune;

class g4hrsParallelWorld : public G4VUserParallelWorld
{
  public:

    g4hrsParallelWorld(G4String);
    virtual ~g4hrsParallelWorld();

  public:

    virtual void Construct();
    virtual void ConstructSD(G4LogicalVolume *);

  private:
      g4hrsIO *fIO;
	g4hrsTune* fTune;

  public:
    void SetIO(g4hrsIO *io){ fIO = io; }
	
	double fHRSAngle;
	double fSeptumAngle;
	int    fSetupHRS;

  private:
      void CreateHRS(G4LogicalVolume *);
	void CreateSeptum(G4LogicalVolume *);	

      double fTargetW;
      double fTargetH;
      double fTargetL;

      double fTargetX;
      double fTargetY;
      double fTargetZ;

      int fSeptuf_UseUniformB;
      int fUseSeptumPlusStdHRS;

      double fPivot2SieveFace;
      double fPivotXOffset;
      double fPivotYOffset;
      double fPivotZOffset;

      bool fSetupSieveSlit;
      bool fSetupCREXGeometry;
      bool fSeptum_UseUniformB;

      int fUseSeptufPlusStdHRS;
      int    fSnakeModel;

      int fSetupStdScatChamber;


      double fScatChamberRin;
      double fScatChamberRout;
      double fScatChamberL;

      double fScatChamberXOffset;
      double fScatChamberYOffset;
      double fScatChamberZOffset;


      double fFieldX;
      double fFieldY;
      double fFieldZ;
 

};

#endif//_G4HRSPARALLELWORLD_HH
