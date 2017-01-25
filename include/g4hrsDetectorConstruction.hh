#ifndef __MOLLERDETECTORCONSTRUCTION_HH
#define __MOLLERDETECTORCONSTRUCTION_HH

#include "G4GDMLParser.hh"
#include "G4VUserDetectorConstruction.hh"
#include <vector>

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSensitiveDetector;

class g4hrsIO;
class g4hrsEMFieldSetup;

class g4hrsDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    g4hrsDetectorConstruction();
    virtual ~g4hrsDetectorConstruction();

  public:

    G4VPhysicalVolume* Construct();

  private:
      g4hrsIO *fIO;
    //----------------------
    // global magnet section
    //----------------------
    //

    G4FieldManager*         fGlobalFieldManager;
    G4String                fDetFileName;

    G4VPhysicalVolume*      fWorldVolume;

  public:
    void SetIO(g4hrsIO *io){ fIO = io; }


  private:
      void CreateTarget(G4LogicalVolume *);
      void CreateTargetChamber(G4LogicalVolume *);
      void CreateSeptum(G4LogicalVolume *);
      void CreateHRS(G4LogicalVolume *);

      double fHRSAngle;
      double fSeptumAngle;

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

      int    fSetupHRS;

      double fScatChamberRin;
      double fScatChamberRout;
      double fScatChamberL;

      double fScatChamberXOffset;
      double fScatChamberYOffset;
      double fScatChamberZOffset;

      g4hrsEMFieldSetup *fEMFieldSetup;

      double fFieldX;
      double fFieldY;
      double fFieldZ;

};

#endif//__MOLLERDETECTORCONSTRUCTION_HH
