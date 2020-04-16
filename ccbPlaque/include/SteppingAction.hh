#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction;
class G4LogicalVolume;

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(EventAction* eventAction);
    virtual ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);

private:
    EventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume_0;
    G4LogicalVolume* fScoringVolume_1;
    G4LogicalVolume* fScoringVolume_2;
    G4LogicalVolume* fScoringVolume_3;
    G4LogicalVolume* fScoringVolume_4;
    G4LogicalVolume* fScoringVolume_5;
    G4LogicalVolume* fScoringVolume_6;
    G4LogicalVolume* fScoringVolume_7;
    G4LogicalVolume* fScoringVolume_8;
    G4LogicalVolume* fScoringVolume_9;
    G4LogicalVolume* fScoringVolume_10;
    G4LogicalVolume* fScoringVolume_11;
    G4LogicalVolume* fScoringVolume_12;
    G4LogicalVolume* fScoringVolume_13;
    G4LogicalVolume* fScoringVolume_14;
    G4LogicalVolume* fScoringVolume_15;
    G4LogicalVolume* fScoringVolume_16;
    G4LogicalVolume* fScoringVolume_17;
    G4LogicalVolume* fScoringVolume_18;
    G4LogicalVolume* fScoringVolume_19;
    G4LogicalVolume* fScoringVolume_20;
    G4LogicalVolume* fScoringVolume_21;
    G4LogicalVolume* fScoringVolume_22;
    G4LogicalVolume* fScoringVolume_23;
    G4LogicalVolume* fScoringVolume_24;
    G4LogicalVolume* fScoringVolume_25;
    G4LogicalVolume* fScoringVolume_26;
    G4LogicalVolume* fScoringVolume_27;
    G4LogicalVolume* fScoringVolume_28;
    G4LogicalVolume* fScoringVolume_29;
    G4LogicalVolume* fScoringVolume_30;
    G4LogicalVolume* fScoringVolume_31;
    G4LogicalVolume* fScoringVolume_32;
    G4LogicalVolume* fScoringVolume_33;
    G4LogicalVolume* fScoringVolume_34;
};

#endif
