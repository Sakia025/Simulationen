#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
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
{ }

DetectorConstruction::~DetectorConstruction()
{ }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* myWater = nist->FindOrBuildMaterial("G4_WATER");
    G4Material* mySilver = nist->FindOrBuildMaterial("G4_Ag");
    G4bool checkOverlaps = true;

    //
    // World
    //
    G4Orb* solidWorld =
        new G4Orb("World", 50.0 * mm);


    G4LogicalVolume* logicWorld =
        new G4LogicalVolume(solidWorld,
                            myWater,
                            "World");

    G4VPhysicalVolume* physWorld =
        new G4PVPlacement(0,
                          G4ThreeVector(),
                          logicWorld,
                          "World",
                          0,
                          false,
                          0,
                          checkOverlaps);

    //
    // CCB Plaque
    //
    G4Sphere*
    solidPlaque =
        new G4Sphere("silverPlaque",
                     12.0 * mm,
                     13.0 * mm,
                     0.0 * deg,
                     360.0 * deg,
                     0.0 * deg,
                     50.98 * deg );

    G4LogicalVolume*
    logicPlaque =
        new G4LogicalVolume(solidPlaque,
                            mySilver,
                            "silverPlaque");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 0.0 * mm),
                      logicPlaque,
                      "silverPlaque",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);

    //
    // Scorer
    //
    G4Box* solidScorer =
        new G4Box("solidScorer",
                  0.25 * mm, 0.25 * mm, 0.25 * mm);


    G4Box* solidScorerBig =
        new G4Box("solidScorerBig",
                  0.5 * mm, 0.5 * mm, 0.5 * mm);

    // Scorer 0.5 mm
    G4LogicalVolume* logicScorer_0 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_0");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 11.5 * mm), // the inner surcafe of the plaque has the radius 12.0 mm
                      logicScorer_0,
                      "physScorer_0",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_0 = logicScorer_0;

    // Scorer 1.0 mm
    G4LogicalVolume* logicScorer_1 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_1");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 11.0 * mm),
                      logicScorer_1,
                      "physScorer_1",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_1 = logicScorer_1;

    // Scorer 1.5 mm
    G4LogicalVolume* logicScorer_2 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_2");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 10.5 * mm),
                      logicScorer_2,
                      "physScorer_2",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_2 = logicScorer_2;

    // Scorer 2.0 mm
    G4LogicalVolume* logicScorer_3 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_3");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 10.0 * mm),
                      logicScorer_3,
                      "physScorer_3",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_3 = logicScorer_3;

    // Scorer 2.5 mm
    G4LogicalVolume* logicScorer_4 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_4");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 9.5 * mm),
                      logicScorer_4,
                      "physScorer_4",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_4 = logicScorer_4;

    // Scorer 3.0 mm
    G4LogicalVolume* logicScorer_5 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_5");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 9.0 * mm),
                      logicScorer_5,
                      "physScorer_5",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_5 = logicScorer_5;

    // Scorer 3.5 mm
    G4LogicalVolume* logicScorer_6 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_6");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 8.5 * mm),
                      logicScorer_6,
                      "physScorer_6",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_6 = logicScorer_6;

    // Scorer 4.0 mm
    G4LogicalVolume* logicScorer_7 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_7");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 8.0 * mm),
                      logicScorer_7,
                      "physScorer_7",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_7 = logicScorer_7;

    // Scorer 4.5 mm
    G4LogicalVolume* logicScorer_8 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_8");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 7.5 * mm),
                      logicScorer_8,
                      "physScorer_8",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_8 = logicScorer_8;

    // Scorer 5.0 mm
    G4LogicalVolume* logicScorer_9 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_9");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 7.0 * mm),
                      logicScorer_9,
                      "physScorer_9",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_9 = logicScorer_9;

    // Scorer 5.5 mm
    G4LogicalVolume* logicScorer_10 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_10");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 6.5 * mm),
                      logicScorer_10,
                      "physScorer_10",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_10 = logicScorer_10;

    // Scorer 6.0 mm
    G4LogicalVolume* logicScorer_11 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_11");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 6.0 * mm),
                      logicScorer_11,
                      "physScorer_11",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_11 = logicScorer_11;

    // Scorer 6.5 mm
    G4LogicalVolume* logicScorer_12 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_12");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 5.5 * mm),
                      logicScorer_12,
                      "physScorer_12",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_12 = logicScorer_12;

    // Scorer 7.0 mm
    G4LogicalVolume* logicScorer_13 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_13");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 5.0 * mm),
                      logicScorer_13,
                      "physScorer_13",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_13 = logicScorer_13;

    // Scorer 7.5 mm
    G4LogicalVolume* logicScorer_14 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_14");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 4.5 * mm),
                      logicScorer_14,
                      "physScorer_14",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_14 = logicScorer_14;

    // Scorer 8.0 mm
    G4LogicalVolume* logicScorer_15 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_15");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 4.0 * mm),
                      logicScorer_15,
                      "physScorer_15",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_15 = logicScorer_15;

    // Scorer 8.5 mm
    G4LogicalVolume* logicScorer_16 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_16");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 3.5 * mm),
                      logicScorer_16,
                      "physScorer_16",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_16 = logicScorer_16;

    // Scorer 9.0 mm
    G4LogicalVolume* logicScorer_17 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_17");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 3.0 * mm),
                      logicScorer_17,
                      "physScorer_17",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_17 = logicScorer_17;

    // Scorer 9.5 mm
    G4LogicalVolume* logicScorer_18 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_18");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 2.5 * mm),
                      logicScorer_18,
                      "physScorer_18",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_18 = logicScorer_18;

    // Scorer 10.0 mm
    G4LogicalVolume* logicScorer_19 =
        new G4LogicalVolume(solidScorer,
                            myWater,
                            "logicScorer_19");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 2.0 * mm),
                      logicScorer_19,
                      "physScorer_19",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_19 = logicScorer_19;

    // Scorer 11.0 mm
    G4LogicalVolume* logicScorer_20 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_20");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 1.0 * mm),
                      logicScorer_20,
                      "physScorer_20",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_20 = logicScorer_20;

    // Scorer 12.0 mm
    G4LogicalVolume* logicScorer_21 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_21");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, 0.0 * mm),
                      logicScorer_21,
                      "physScorer_21",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_21 = logicScorer_21;

    // Scorer 13.0 mm
    G4LogicalVolume* logicScorer_22 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_22");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -1.0 * mm),
                      logicScorer_22,
                      "physScorer_22",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_22 = logicScorer_22;

    // Scorer 14.0 mm
    G4LogicalVolume* logicScorer_23 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_23");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -2.0 * mm),
                      logicScorer_23,
                      "physScorer_23",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_23 = logicScorer_23;

    // Scorer 15.0 mm
    G4LogicalVolume* logicScorer_24 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_24");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -3.0 * mm),
                      logicScorer_24,
                      "physScorer_24",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_24 = logicScorer_24;

    // Scorer 16.0 mm
    G4LogicalVolume* logicScorer_25 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_25");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -4.0 * mm),
                      logicScorer_25,
                      "physScorer_25",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_25 = logicScorer_25;

    // Scorer 17.0 mm
    G4LogicalVolume* logicScorer_26 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_26");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -5.0 * mm),
                      logicScorer_26,
                      "physScorer_26",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_26 = logicScorer_26;

    // Scorer 18.0 mm
    G4LogicalVolume* logicScorer_27 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_27");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -6.0 * mm),
                      logicScorer_27,
                      "physScorer_27",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_27 = logicScorer_27;

    // Scorer 19.0 mm
    G4LogicalVolume* logicScorer_28 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_28");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -7.0 * mm),
                      logicScorer_28,
                      "physScorer_28",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_28 = logicScorer_28;

    // Scorer 20.0 mm
    G4LogicalVolume* logicScorer_29 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_29");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -8.0 * mm),
                      logicScorer_29,
                      "physScorer_29",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_29 = logicScorer_29;

    // Scorer 21.0 mm
    G4LogicalVolume* logicScorer_30 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_30");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -9.0 * mm),
                      logicScorer_30,
                      "physScorer_30",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_30 = logicScorer_30;

    // Scorer 22.0 mm
    G4LogicalVolume* logicScorer_31 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_31");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -10.0 * mm),
                      logicScorer_31,
                      "physScorer_31",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_31 = logicScorer_31;

    // Scorer 23.0 mm
    G4LogicalVolume* logicScorer_32 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_32");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -11.0 * mm),
                      logicScorer_32,
                      "physScorer_32",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_32 = logicScorer_32;

    // Scorer 24.0 mm
    G4LogicalVolume* logicScorer_33 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_33");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -12.0 * mm),
                      logicScorer_33,
                      "physScorer_33",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    fScoringVolume_33 = logicScorer_33;

    // Scorer 25.0 mm
    G4LogicalVolume* logicScorer_34 =
        new G4LogicalVolume(solidScorerBig,
                            myWater,
                            "logicScorer_34");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, -13.0 * mm),
                      logicScorer_34,
                      "physScorer_34",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);

    fScoringVolume_34 = logicScorer_34;

    return physWorld;
}
