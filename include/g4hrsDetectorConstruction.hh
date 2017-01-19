#ifndef __MOLLERDETECTORCONSTRUCTION_HH
#define __MOLLERDETECTORCONSTRUCTION_HH

#include "G4GDMLParser.hh"
#include "G4VUserDetectorConstruction.hh"
#include "g4hrsGlobalField.hh"
#include <vector>

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSensitiveDetector;

class g4hrsIO;

class g4hrsDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    g4hrsDetectorConstruction();
    virtual ~g4hrsDetectorConstruction();

  public:

    G4VPhysicalVolume* Construct();

    void CreateGlobalMagneticField();

    g4hrsGlobalField* GetGlobalField(){ return fGlobalField; }

  private:
    //----------------------
    // global magnet section
    //----------------------
    //

    G4FieldManager*         fGlobalFieldManager;
    g4hrsGlobalField*       fGlobalField;
    G4String                fDetFileName;

    G4VPhysicalVolume*      fWorldVolume;

  public:


  private:
      void CreateTarget(G4LogicalVolume *);
      void CreateTargetChamber(G4LogicalVolume *);
      void CreateSeptum(G4LogicalVolume *);
      void CreateHRS(G4LogicalVolume *);


};

#endif//__MOLLERDETECTORCONSTRUCTION_HH
