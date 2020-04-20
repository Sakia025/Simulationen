#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
    : G4UserSteppingAction(),
      fEventAction(eventAction),
      fScoringVolume_0(0),
      fScoringVolume_1(0),
      fScoringVolume_2(0),
      fScoringVolume_3(0),
      fScoringVolume_4(0),
      fScoringVolume_5(0),
      fScoringVolume_6(0),
      fScoringVolume_7(0),
      fScoringVolume_8(0),
      fScoringVolume_9(0),
      fScoringVolume_10(0),
      fScoringVolume_11(0),
      fScoringVolume_12(0),
      fScoringVolume_13(0),
      fScoringVolume_14(0),
      fScoringVolume_15(0),
      fScoringVolume_16(0),
      fScoringVolume_17(0),
      fScoringVolume_18(0),
      fScoringVolume_19(0),
      fScoringVolume_20(0),
      fScoringVolume_21(0),
      fScoringVolume_22(0),
      fScoringVolume_23(0),
      fScoringVolume_24(0),
      fScoringVolume_25(0),
      fScoringVolume_26(0),
      fScoringVolume_27(0),
      fScoringVolume_28(0),
      fScoringVolume_29(0),
      fScoringVolume_30(0),
      fScoringVolume_31(0),
      fScoringVolume_32(0),
      fScoringVolume_33(0),
      fScoringVolume_34(0)
{}

SteppingAction::~SteppingAction()
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    if (!fScoringVolume_0)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_0 = detectorConstruction->GetScoringVolume_0();
    }

    if (!fScoringVolume_1)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_1 = detectorConstruction->GetScoringVolume_1();
    }

    if (!fScoringVolume_2)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_2 = detectorConstruction->GetScoringVolume_2();
    }

    if (!fScoringVolume_3)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_3 = detectorConstruction->GetScoringVolume_3();
    }

    if (!fScoringVolume_4)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_4 = detectorConstruction->GetScoringVolume_4();
    }

    if (!fScoringVolume_5)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_5 = detectorConstruction->GetScoringVolume_5();
    }

    if (!fScoringVolume_6)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_6 = detectorConstruction->GetScoringVolume_6();
    }

    if (!fScoringVolume_7)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_7 = detectorConstruction->GetScoringVolume_7();
    }

    if (!fScoringVolume_8)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_8 = detectorConstruction->GetScoringVolume_8();
    }

    if (!fScoringVolume_9)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_9 = detectorConstruction->GetScoringVolume_9();
    }

    if (!fScoringVolume_10)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_10 = detectorConstruction->GetScoringVolume_10();
    }

    if (!fScoringVolume_11)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_11 = detectorConstruction->GetScoringVolume_11();
    }

    if (!fScoringVolume_12)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_12 = detectorConstruction->GetScoringVolume_12();
    }

    if (!fScoringVolume_13)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_13 = detectorConstruction->GetScoringVolume_13();
    }

    if (!fScoringVolume_14)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_14 = detectorConstruction->GetScoringVolume_14();
    }

    if (!fScoringVolume_15)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_15 = detectorConstruction->GetScoringVolume_15();
    }

    if (!fScoringVolume_16)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_16 = detectorConstruction->GetScoringVolume_16();
    }

    if (!fScoringVolume_17)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_17 = detectorConstruction->GetScoringVolume_17();
    }

    if (!fScoringVolume_18)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_18 = detectorConstruction->GetScoringVolume_18();
    }

    if (!fScoringVolume_19)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_19 = detectorConstruction->GetScoringVolume_19();
    }

    if (!fScoringVolume_20)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_20 = detectorConstruction->GetScoringVolume_20();
    }

    if (!fScoringVolume_21)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_21 = detectorConstruction->GetScoringVolume_21();
    }

    if (!fScoringVolume_22)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_22 = detectorConstruction->GetScoringVolume_22();
    }

    if (!fScoringVolume_23)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_23 = detectorConstruction->GetScoringVolume_23();
    }

    if (!fScoringVolume_24)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_24 = detectorConstruction->GetScoringVolume_24();
    }

    if (!fScoringVolume_25)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_25 = detectorConstruction->GetScoringVolume_25();
    }

    if (!fScoringVolume_26)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_26 = detectorConstruction->GetScoringVolume_26();
    }

    if (!fScoringVolume_27)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_27 = detectorConstruction->GetScoringVolume_27();
    }

    if (!fScoringVolume_28)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_28 = detectorConstruction->GetScoringVolume_28();
    }

    if (!fScoringVolume_29)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_29 = detectorConstruction->GetScoringVolume_29();
    }

    if (!fScoringVolume_30)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_30 = detectorConstruction->GetScoringVolume_30();
    }

    if (!fScoringVolume_31)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_31 = detectorConstruction->GetScoringVolume_31();
    }

    if (!fScoringVolume_32)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_32 = detectorConstruction->GetScoringVolume_32();
    }

    if (!fScoringVolume_33)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_33 = detectorConstruction->GetScoringVolume_33();
    }

    if (!fScoringVolume_34)
    {
        const DetectorConstruction* detectorConstruction
            = static_cast<const DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume_34 = detectorConstruction->GetScoringVolume_34();
    }

    G4LogicalVolume* volume
        = step->GetPreStepPoint()->GetTouchableHandle()
          ->GetVolume()->GetLogicalVolume();



    if (volume == fScoringVolume_0)
    {
        G4double edepStep_0 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_0(edepStep_0);
    }
    else if (volume == fScoringVolume_1)
    {
        G4double edepStep_1 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_1(edepStep_1);
    }
    else if (volume == fScoringVolume_2)
    {
        G4double edepStep_2 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_2(edepStep_2);
    }
    else if (volume == fScoringVolume_3)
    {
        G4double edepStep_3 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_3(edepStep_3);
    }
    else if (volume == fScoringVolume_4)
    {
        G4double edepStep_4 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_4(edepStep_4);
    }
    else if (volume == fScoringVolume_5)
    {
        G4double edepStep_5 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_5(edepStep_5);
    }
    else if (volume == fScoringVolume_6)
    {
        G4double edepStep_6 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_6(edepStep_6);
    }
    else if (volume == fScoringVolume_7)
    {
        G4double edepStep_7 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_7(edepStep_7);
    }
    else if (volume == fScoringVolume_8)
    {
        G4double edepStep_8 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_8(edepStep_8);
    }
    else if (volume == fScoringVolume_9)
    {
        G4double edepStep_9 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_9(edepStep_9);
    }
    else if (volume == fScoringVolume_10)
    {
        G4double edepStep_10 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_10(edepStep_10);
    }
    else if (volume == fScoringVolume_11)
    {
        G4double edepStep_11 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_11(edepStep_11);
    }
    else if (volume == fScoringVolume_12)
    {
        G4double edepStep_12 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_12(edepStep_12);
    }
    else if (volume == fScoringVolume_13)
    {
        G4double edepStep_13 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_13(edepStep_13);
    }
    else if (volume == fScoringVolume_14)
    {
        G4double edepStep_14 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_14(edepStep_14);
    }
    else if (volume == fScoringVolume_15)
    {
        G4double edepStep_15 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_15(edepStep_15);
    }
    else if (volume == fScoringVolume_16)
    {
        G4double edepStep_16 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_16(edepStep_16);
    }
    else if (volume == fScoringVolume_17)
    {
        G4double edepStep_17 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_17(edepStep_17);
    }
    else if (volume == fScoringVolume_18)
    {
        G4double edepStep_18 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_18(edepStep_18);
    }
    else if (volume == fScoringVolume_19)
    {
        G4double edepStep_19 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_19(edepStep_19);
    }
    else if (volume == fScoringVolume_20)
    {
        G4double edepStep_20 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_20(edepStep_20);
    }
    else if (volume == fScoringVolume_21)
    {
        G4double edepStep_21 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_21(edepStep_21);
    }
    else if (volume == fScoringVolume_22)
    {
        G4double edepStep_22 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_22(edepStep_22);
    }
    else if (volume == fScoringVolume_23)
    {
        G4double edepStep_23 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_23(edepStep_23);
    }
    else if (volume == fScoringVolume_24)
    {
        G4double edepStep_24 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_24(edepStep_24);
    }
    else if (volume == fScoringVolume_25)
    {
        G4double edepStep_25 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_25(edepStep_25);
    }
    else if (volume == fScoringVolume_26)
    {
        G4double edepStep_26 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_26(edepStep_26);
    }
    else if (volume == fScoringVolume_27)
    {
        G4double edepStep_27 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_27(edepStep_27);
    }
    else if (volume == fScoringVolume_28)
    {
        G4double edepStep_28 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_28(edepStep_28);
    }
    else if (volume == fScoringVolume_29)
    {
        G4double edepStep_29 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_29(edepStep_29);
    }
    else if (volume == fScoringVolume_30)
    {
        G4double edepStep_30 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_30(edepStep_30);
    }
    else if (volume == fScoringVolume_31)
    {
        G4double edepStep_31 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_31(edepStep_31);
    }
    else if (volume == fScoringVolume_32)
    {
        G4double edepStep_32 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_32(edepStep_32);
    }
    else if (volume == fScoringVolume_33)
    {
        G4double edepStep_33 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_33(edepStep_33);
    }
    else if (volume == fScoringVolume_34)
    {
        G4double edepStep_34 = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep_34(edepStep_34);
    }
    else return;
}
