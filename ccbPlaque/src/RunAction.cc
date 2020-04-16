#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

RunAction::RunAction()
    : G4UserRunAction(),
      fEdep_0("Edep_0", 0.),
      fEdep_1("Edep_1", 0.),
      fEdep_2("Edep_2", 0.),
      fEdep_3("Edep_3", 0.),
      fEdep_4("Edep_4", 0.),
      fEdep_5("Edep_5", 0.),
      fEdep_6("Edep_6", 0.),
      fEdep_7("Edep_7", 0.),
      fEdep_8("Edep_8", 0.),
      fEdep_9("Edep_9", 0.),
      fEdep_10("Edep_10", 0.),
      fEdep_11("Edep_11", 0.),
      fEdep_12("Edep_12", 0.),
      fEdep_13("Edep_13", 0.),
      fEdep_14("Edep_14", 0.),
      fEdep_15("Edep_15", 0.),
      fEdep_16("Edep_16", 0.),
      fEdep_17("Edep_17", 0.),
      fEdep_18("Edep_18", 0.),
      fEdep_19("Edep_19", 0.),
      fEdep_20("Edep_20", 0.),
      fEdep_21("Edep_21", 0.),
      fEdep_22("Edep_22", 0.),
      fEdep_23("Edep_23", 0.),
      fEdep_24("Edep_24", 0.),
      fEdep_25("Edep_25", 0.),
      fEdep_26("Edep_26", 0.),
      fEdep_27("Edep_27", 0.),
      fEdep_28("Edep_28", 0.),
      fEdep_29("Edep_29", 0.),
      fEdep_30("Edep_30", 0.),
      fEdep_31("Edep_31", 0.),
      fEdep_32("Edep_32", 0.),
      fEdep_33("Edep_33", 0.),
      fEdep_34("Edep_34", 0.)
{
    G4AccumulableManager* parameterManager = G4AccumulableManager::Instance();
    parameterManager->RegisterAccumulable(fEdep_0);
    parameterManager->RegisterAccumulable(fEdep_1);
    parameterManager->RegisterAccumulable(fEdep_2);
    parameterManager->RegisterAccumulable(fEdep_3);
    parameterManager->RegisterAccumulable(fEdep_4);
    parameterManager->RegisterAccumulable(fEdep_5);
    parameterManager->RegisterAccumulable(fEdep_6);
    parameterManager->RegisterAccumulable(fEdep_7);
    parameterManager->RegisterAccumulable(fEdep_8);
    parameterManager->RegisterAccumulable(fEdep_9);
    parameterManager->RegisterAccumulable(fEdep_10);
    parameterManager->RegisterAccumulable(fEdep_11);
    parameterManager->RegisterAccumulable(fEdep_12);
    parameterManager->RegisterAccumulable(fEdep_13);
    parameterManager->RegisterAccumulable(fEdep_14);
    parameterManager->RegisterAccumulable(fEdep_15);
    parameterManager->RegisterAccumulable(fEdep_16);
    parameterManager->RegisterAccumulable(fEdep_17);
    parameterManager->RegisterAccumulable(fEdep_18);
    parameterManager->RegisterAccumulable(fEdep_19);
    parameterManager->RegisterAccumulable(fEdep_20);
    parameterManager->RegisterAccumulable(fEdep_21);
    parameterManager->RegisterAccumulable(fEdep_22);
    parameterManager->RegisterAccumulable(fEdep_23);
    parameterManager->RegisterAccumulable(fEdep_24);
    parameterManager->RegisterAccumulable(fEdep_25);
    parameterManager->RegisterAccumulable(fEdep_26);
    parameterManager->RegisterAccumulable(fEdep_27);
    parameterManager->RegisterAccumulable(fEdep_28);
    parameterManager->RegisterAccumulable(fEdep_29);
    parameterManager->RegisterAccumulable(fEdep_30);
    parameterManager->RegisterAccumulable(fEdep_31);
    parameterManager->RegisterAccumulable(fEdep_32);
    parameterManager->RegisterAccumulable(fEdep_33);
    parameterManager->RegisterAccumulable(fEdep_34);
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run*)
{
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    G4AccumulableManager* parameterManager = G4AccumulableManager::Instance();
    parameterManager->Reset();
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    G4int nofEvents = run->GetNumberOfEvent();

    if (nofEvents == 0) return;

    G4AccumulableManager* parameterManager = G4AccumulableManager::Instance();
    parameterManager->Merge();

    G4double edep_0  = fEdep_0.GetValue();
    G4double edep_1  = fEdep_1.GetValue();
    G4double edep_2  = fEdep_2.GetValue();
    G4double edep_3  = fEdep_3.GetValue();
    G4double edep_4  = fEdep_4.GetValue();
    G4double edep_5  = fEdep_5.GetValue();
    G4double edep_6  = fEdep_6.GetValue();
    G4double edep_7  = fEdep_7.GetValue();
    G4double edep_8  = fEdep_8.GetValue();
    G4double edep_9  = fEdep_9.GetValue();
    G4double edep_10  = fEdep_10.GetValue();
    G4double edep_11  = fEdep_11.GetValue();
    G4double edep_12  = fEdep_12.GetValue();
    G4double edep_13  = fEdep_13.GetValue();
    G4double edep_14  = fEdep_14.GetValue();
    G4double edep_15  = fEdep_15.GetValue();
    G4double edep_16  = fEdep_16.GetValue();
    G4double edep_17  = fEdep_17.GetValue();
    G4double edep_18  = fEdep_18.GetValue();
    G4double edep_19  = fEdep_19.GetValue();
    G4double edep_20  = fEdep_20.GetValue();
    G4double edep_21  = fEdep_21.GetValue();
    G4double edep_22  = fEdep_22.GetValue();
    G4double edep_23  = fEdep_23.GetValue();
    G4double edep_24  = fEdep_24.GetValue();
    G4double edep_25  = fEdep_25.GetValue();
    G4double edep_26  = fEdep_26.GetValue();
    G4double edep_27  = fEdep_27.GetValue();
    G4double edep_28  = fEdep_28.GetValue();
    G4double edep_29  = fEdep_29.GetValue();
    G4double edep_30  = fEdep_30.GetValue();
    G4double edep_31  = fEdep_31.GetValue();
    G4double edep_32  = fEdep_32.GetValue();
    G4double edep_33  = fEdep_33.GetValue();
    G4double edep_34  = fEdep_34.GetValue();

    const DetectorConstruction* detectorConstruction
        = static_cast<const DetectorConstruction*>
          (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    G4double mass_0 = detectorConstruction->GetScoringVolume_0()->GetMass();
    G4double dose_0 = edep_0 / mass_0;
    G4double mass_1 = detectorConstruction->GetScoringVolume_1()->GetMass();
    G4double dose_1 = edep_1 / mass_1;
    G4double mass_2 = detectorConstruction->GetScoringVolume_2()->GetMass();
    G4double dose_2 = edep_2 / mass_2;
    G4double mass_3 = detectorConstruction->GetScoringVolume_3()->GetMass();
    G4double dose_3 = edep_3 / mass_3;
    G4double mass_4 = detectorConstruction->GetScoringVolume_4()->GetMass();
    G4double dose_4 = edep_4 / mass_4;
    G4double mass_5 = detectorConstruction->GetScoringVolume_5()->GetMass();
    G4double dose_5 = edep_5 / mass_5;
    G4double mass_6 = detectorConstruction->GetScoringVolume_6()->GetMass();
    G4double dose_6 = edep_6 / mass_6;
    G4double mass_7 = detectorConstruction->GetScoringVolume_7()->GetMass();
    G4double dose_7 = edep_7 / mass_7;
    G4double mass_8 = detectorConstruction->GetScoringVolume_8()->GetMass();
    G4double dose_8 = edep_8 / mass_8;
    G4double mass_9 = detectorConstruction->GetScoringVolume_9()->GetMass();
    G4double dose_9 = edep_9 / mass_9;
    G4double mass_10 = detectorConstruction->GetScoringVolume_10()->GetMass();
    G4double dose_10 = edep_10 / mass_10;
    G4double mass_11 = detectorConstruction->GetScoringVolume_11()->GetMass();
    G4double dose_11 = edep_11 / mass_11;
    G4double mass_12 = detectorConstruction->GetScoringVolume_12()->GetMass();
    G4double dose_12 = edep_12 / mass_12;
    G4double mass_13 = detectorConstruction->GetScoringVolume_13()->GetMass();
    G4double dose_13 = edep_13 / mass_13;
    G4double mass_14 = detectorConstruction->GetScoringVolume_14()->GetMass();
    G4double dose_14 = edep_14 / mass_14;
    G4double mass_15 = detectorConstruction->GetScoringVolume_15()->GetMass();
    G4double dose_15 = edep_15 / mass_15;
    G4double mass_16 = detectorConstruction->GetScoringVolume_16()->GetMass();
    G4double dose_16 = edep_16 / mass_16;
    G4double mass_17 = detectorConstruction->GetScoringVolume_17()->GetMass();
    G4double dose_17 = edep_17 / mass_17;
    G4double mass_18 = detectorConstruction->GetScoringVolume_18()->GetMass();
    G4double dose_18 = edep_18 / mass_18;
    G4double mass_19 = detectorConstruction->GetScoringVolume_19()->GetMass();
    G4double dose_19 = edep_19 / mass_19;
    G4double mass_20 = detectorConstruction->GetScoringVolume_20()->GetMass();
    G4double dose_20 = edep_20 / mass_20;
    G4double mass_21 = detectorConstruction->GetScoringVolume_21()->GetMass();
    G4double dose_21 = edep_21 / mass_21;
    G4double mass_22 = detectorConstruction->GetScoringVolume_22()->GetMass();
    G4double dose_22 = edep_22 / mass_22;
    G4double mass_23 = detectorConstruction->GetScoringVolume_23()->GetMass();
    G4double dose_23 = edep_23 / mass_23;
    G4double mass_24 = detectorConstruction->GetScoringVolume_24()->GetMass();
    G4double dose_24 = edep_24 / mass_24;
    G4double mass_25 = detectorConstruction->GetScoringVolume_25()->GetMass();
    G4double dose_25 = edep_25 / mass_25;
    G4double mass_26 = detectorConstruction->GetScoringVolume_26()->GetMass();
    G4double dose_26 = edep_26 / mass_26;
    G4double mass_27 = detectorConstruction->GetScoringVolume_27()->GetMass();
    G4double dose_27 = edep_27 / mass_27;
    G4double mass_28 = detectorConstruction->GetScoringVolume_28()->GetMass();
    G4double dose_28 = edep_28 / mass_28;
    G4double mass_29 = detectorConstruction->GetScoringVolume_29()->GetMass();
    G4double dose_29 = edep_29 / mass_29;
    G4double mass_30 = detectorConstruction->GetScoringVolume_30()->GetMass();
    G4double dose_30 = edep_30 / mass_30;
    G4double mass_31 = detectorConstruction->GetScoringVolume_31()->GetMass();
    G4double dose_31 = edep_31 / mass_31;
    G4double mass_32 = detectorConstruction->GetScoringVolume_32()->GetMass();
    G4double dose_32 = edep_32 / mass_32;
    G4double mass_33 = detectorConstruction->GetScoringVolume_33()->GetMass();
    G4double dose_33 = edep_33 / mass_33;
    G4double mass_34 = detectorConstruction->GetScoringVolume_34()->GetMass();
    G4double dose_34 = edep_34 / mass_34;

    const PrimaryGeneratorAction* generatorAction
        = static_cast<const PrimaryGeneratorAction*>
          (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    G4String runCondition;

    if (generatorAction)
    {
        const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
        runCondition += particleGun->GetParticleDefinition()->GetParticleName();
        runCondition += " of ";
        G4double particleEnergy = particleGun->GetParticleEnergy();
        runCondition += G4BestUnit(particleEnergy, "Energy");
    }

    if (IsMaster())
    {
        G4cout
                << "#**********************************************************************************" << G4endl
                << "# " << nofEvents << " events were processed in this run." << runCondition
                << G4endl
                << "# depth (mm)    dose (gray)"   << G4endl
                << "0.5 "  << dose_0 / gray  << " " << G4endl
                << "1.0 "  << dose_1 / gray  << " " << G4endl
                << "1.5 "  << dose_2 / gray  << " " << G4endl
                << "2.0 "  << dose_3 / gray  << " " << G4endl
                << "2.5 "  << dose_4 / gray  << " " << G4endl
                << "3.0 "  << dose_5 / gray  << " " << G4endl
                << "3.5 "  << dose_6 / gray  << " " << G4endl
                << "4.0 "  << dose_7 / gray  << " " << G4endl
                << "4.5 "  << dose_8 / gray  << " " << G4endl
                << "5.0 "  << dose_9 / gray  << " " << G4endl
                << "5.5 "  << dose_10 / gray << " " << G4endl
                << "6.0 "  << dose_11 / gray << " " << G4endl
                << "6.5 "  << dose_12 / gray << " " << G4endl
                << "7.0 "  << dose_13 / gray << " " << G4endl
                << "7.5 "  << dose_14 / gray << " " << G4endl
                << "8.0 "  << dose_15 / gray << " " << G4endl
                << "8.5 "  << dose_16 / gray << " " << G4endl
                << "9.0 "  << dose_17 / gray << " " << G4endl
                << "9.5 "  << dose_18 / gray << " " << G4endl
                << "10.0 " << dose_19 / gray << " " << G4endl
                << "11.0 " << dose_20 / gray << " " << G4endl
                << "12.0 " << dose_21 / gray << " " << G4endl
                << "13.0 " << dose_22 / gray << " " << G4endl
                << "14.0 " << dose_23 / gray << " " << G4endl
                << "15.0 " << dose_24 / gray << " " << G4endl
                << "16.0 " << dose_25 / gray << " " << G4endl
                << "17.0 " << dose_26 / gray << " " << G4endl
                << "18.0 " << dose_27 / gray << " " << G4endl
                << "19.0 " << dose_28 / gray << " " << G4endl
                << "20.0 " << dose_29 / gray << " " << G4endl
                << "21.0 " << dose_30 / gray << " " << G4endl
                << "22.0 " << dose_31 / gray << " " << G4endl
                << "23.0 " << dose_32 / gray << " " << G4endl
                << "24.0 " << dose_33 / gray << " " << G4endl
                << "25.0 " << dose_34 / gray << " " << G4endl;
    }

}

void RunAction::AddEdep_0(G4double edep_0)
{
    fEdep_0  += edep_0;
}
void RunAction::AddEdep_1(G4double edep_1)
{
    fEdep_1  += edep_1;
}
void RunAction::AddEdep_2(G4double edep_2)
{
    fEdep_2  += edep_2;
}
void RunAction::AddEdep_3(G4double edep_3)
{
    fEdep_3  += edep_3;
}
void RunAction::AddEdep_4(G4double edep_4)
{
    fEdep_4  += edep_4;
}
void RunAction::AddEdep_5(G4double edep_5)
{
    fEdep_5  += edep_5;
}
void RunAction::AddEdep_6(G4double edep_6)
{
    fEdep_6  += edep_6;
}
void RunAction::AddEdep_7(G4double edep_7)
{
    fEdep_7  += edep_7;
}
void RunAction::AddEdep_8(G4double edep_8)
{
    fEdep_8  += edep_8;
}
void RunAction::AddEdep_9(G4double edep_9)
{
    fEdep_9  += edep_9;
}
void RunAction::AddEdep_10(G4double edep_10)
{
    fEdep_10  += edep_10;
}
void RunAction::AddEdep_11(G4double edep_11)
{
    fEdep_11  += edep_11;
}
void RunAction::AddEdep_12(G4double edep_12)
{
    fEdep_12  += edep_12;
}
void RunAction::AddEdep_13(G4double edep_13)
{
    fEdep_13  += edep_13;
}
void RunAction::AddEdep_14(G4double edep_14)
{
    fEdep_14  += edep_14;
}
void RunAction::AddEdep_15(G4double edep_15)
{
    fEdep_15  += edep_15;
}
void RunAction::AddEdep_16(G4double edep_16)
{
    fEdep_16  += edep_16;
}
void RunAction::AddEdep_17(G4double edep_17)
{
    fEdep_17  += edep_17;
}
void RunAction::AddEdep_18(G4double edep_18)
{
    fEdep_18  += edep_18;
}
void RunAction::AddEdep_19(G4double edep_19)
{
    fEdep_19  += edep_19;
}
void RunAction::AddEdep_20(G4double edep_20)
{
    fEdep_20  += edep_20;
}
void RunAction::AddEdep_21(G4double edep_21)
{
    fEdep_21  += edep_21;
}
void RunAction::AddEdep_22(G4double edep_22)
{
    fEdep_22  += edep_22;
}
void RunAction::AddEdep_23(G4double edep_23)
{
    fEdep_23  += edep_23;
}
void RunAction::AddEdep_24(G4double edep_24)
{
    fEdep_24  += edep_24;
}
void RunAction::AddEdep_25(G4double edep_25)
{
    fEdep_25  += edep_25;
}
void RunAction::AddEdep_26(G4double edep_26)
{
    fEdep_26  += edep_26;
}
void RunAction::AddEdep_27(G4double edep_27)
{
    fEdep_27  += edep_27;
}
void RunAction::AddEdep_28(G4double edep_28)
{
    fEdep_28  += edep_28;
}
void RunAction::AddEdep_29(G4double edep_29)
{
    fEdep_29  += edep_29;
}
void RunAction::AddEdep_30(G4double edep_30)
{
    fEdep_30  += edep_30;
}
void RunAction::AddEdep_31(G4double edep_31)
{
    fEdep_31  += edep_31;
}
void RunAction::AddEdep_32(G4double edep_32)
{
    fEdep_32  += edep_32;
}
void RunAction::AddEdep_33(G4double edep_33)
{
    fEdep_33  += edep_33;
}
void RunAction::AddEdep_34(G4double edep_34)
{
    fEdep_34  += edep_34;
}
