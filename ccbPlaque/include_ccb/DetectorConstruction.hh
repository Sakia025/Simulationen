#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    G4LogicalVolume* GetScoringVolume_0() const
    {
        return fScoringVolume_0;
    }
    G4LogicalVolume* GetScoringVolume_1() const
    {
        return fScoringVolume_1;
    }
    G4LogicalVolume* GetScoringVolume_2() const
    {
        return fScoringVolume_2;
    }
    G4LogicalVolume* GetScoringVolume_3() const
    {
        return fScoringVolume_3;
    }
    G4LogicalVolume* GetScoringVolume_4() const
    {
        return fScoringVolume_4;
    }
    G4LogicalVolume* GetScoringVolume_5() const
    {
        return fScoringVolume_5;
    }
    G4LogicalVolume* GetScoringVolume_6() const
    {
        return fScoringVolume_6;
    }
    G4LogicalVolume* GetScoringVolume_7() const
    {
        return fScoringVolume_7;
    }
    G4LogicalVolume* GetScoringVolume_8() const
    {
        return fScoringVolume_8;
    }
    G4LogicalVolume* GetScoringVolume_9() const
    {
        return fScoringVolume_9;
    }
    G4LogicalVolume* GetScoringVolume_10() const
    {
        return fScoringVolume_10;
    }
    G4LogicalVolume* GetScoringVolume_11() const
    {
        return fScoringVolume_11;
    }
    G4LogicalVolume* GetScoringVolume_12() const
    {
        return fScoringVolume_12;
    }
    G4LogicalVolume* GetScoringVolume_13() const
    {
        return fScoringVolume_13;
    }
    G4LogicalVolume* GetScoringVolume_14() const
    {
        return fScoringVolume_14;
    }
    G4LogicalVolume* GetScoringVolume_15() const
    {
        return fScoringVolume_15;
    }
    G4LogicalVolume* GetScoringVolume_16() const
    {
        return fScoringVolume_16;
    }
    G4LogicalVolume* GetScoringVolume_17() const
    {
        return fScoringVolume_17;
    }
    G4LogicalVolume* GetScoringVolume_18() const
    {
        return fScoringVolume_18;
    }
    G4LogicalVolume* GetScoringVolume_19() const
    {
        return fScoringVolume_19;
    }
    G4LogicalVolume* GetScoringVolume_20() const
    {
        return fScoringVolume_20;
    }
    G4LogicalVolume* GetScoringVolume_21() const
    {
        return fScoringVolume_21;
    }
    G4LogicalVolume* GetScoringVolume_22() const
    {
        return fScoringVolume_22;
    }
    G4LogicalVolume* GetScoringVolume_23() const
    {
        return fScoringVolume_23;
    }
    G4LogicalVolume* GetScoringVolume_24() const
    {
        return fScoringVolume_24;
    }
    G4LogicalVolume* GetScoringVolume_25() const
    {
        return fScoringVolume_25;
    }
    G4LogicalVolume* GetScoringVolume_26() const
    {
        return fScoringVolume_26;
    }
    G4LogicalVolume* GetScoringVolume_27() const
    {
        return fScoringVolume_27;
    }
    G4LogicalVolume* GetScoringVolume_28() const
    {
        return fScoringVolume_28;
    }
    G4LogicalVolume* GetScoringVolume_29() const
    {
        return fScoringVolume_29;
    }
    G4LogicalVolume* GetScoringVolume_30() const
    {
        return fScoringVolume_30;
    }
    G4LogicalVolume* GetScoringVolume_31() const
    {
        return fScoringVolume_31;
    }
    G4LogicalVolume* GetScoringVolume_32() const
    {
        return fScoringVolume_32;
    }
    G4LogicalVolume* GetScoringVolume_33() const
    {
        return fScoringVolume_33;
    }
    G4LogicalVolume* GetScoringVolume_34() const
    {
        return fScoringVolume_34;
    }
protected:
    G4LogicalVolume*  fScoringVolume_0;
    G4LogicalVolume*  fScoringVolume_1;
    G4LogicalVolume*  fScoringVolume_2;
    G4LogicalVolume*  fScoringVolume_3;
    G4LogicalVolume*  fScoringVolume_4;
    G4LogicalVolume*  fScoringVolume_5;
    G4LogicalVolume*  fScoringVolume_6;
    G4LogicalVolume*  fScoringVolume_7;
    G4LogicalVolume*  fScoringVolume_8;
    G4LogicalVolume*  fScoringVolume_9;
    G4LogicalVolume*  fScoringVolume_10;
    G4LogicalVolume*  fScoringVolume_11;
    G4LogicalVolume*  fScoringVolume_12;
    G4LogicalVolume*  fScoringVolume_13;
    G4LogicalVolume*  fScoringVolume_14;
    G4LogicalVolume*  fScoringVolume_15;
    G4LogicalVolume*  fScoringVolume_16;
    G4LogicalVolume*  fScoringVolume_17;
    G4LogicalVolume*  fScoringVolume_18;
    G4LogicalVolume*  fScoringVolume_19;
    G4LogicalVolume*  fScoringVolume_20;
    G4LogicalVolume*  fScoringVolume_21;
    G4LogicalVolume*  fScoringVolume_22;
    G4LogicalVolume*  fScoringVolume_23;
    G4LogicalVolume*  fScoringVolume_24;
    G4LogicalVolume*  fScoringVolume_25;
    G4LogicalVolume*  fScoringVolume_26;
    G4LogicalVolume*  fScoringVolume_27;
    G4LogicalVolume*  fScoringVolume_28;
    G4LogicalVolume*  fScoringVolume_29;
    G4LogicalVolume*  fScoringVolume_30;
    G4LogicalVolume*  fScoringVolume_31;
    G4LogicalVolume*  fScoringVolume_32;
    G4LogicalVolume*  fScoringVolume_33;
    G4LogicalVolume*  fScoringVolume_34;
private:
    void DefineMaterials();
    G4Material* MatchMaterialToSTL(std::string model_name);
  //  std::vector<std::string> ReadModelNames(DIR *dr);
};

#endif
