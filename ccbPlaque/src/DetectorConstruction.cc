#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"

// CADMESH //
//#define CADMESH_LEXER_VERBOSE
#include "CADMesh.hh"


#include <iostream>
#include <list>
//using System.IO;
#include <string>
#include <stdio.h>
#include <dirent.h>
#include <vector>

//template <typename T>;
template<typename T>
void get_type_name(T&&);

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
{
      DefineMaterials();//dont forget to define this in detectorconstruction.hh
}

DetectorConstruction::~DetectorConstruction()
{}


  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
  {
      G4NistManager* nist = G4NistManager::Instance();

      G4bool isotopes = false;
      G4double density, a;
      G4String name, symbol;
      G4int ncomponents, natoms, n, z, nIso;
      G4State state;
      G4Element* element;
      G4Material* material;
      G4double fractionmass;

      //Eye Elements and Ruthenium106 Isotope
      G4Element*  oxygen = nist->FindOrBuildElement("O", isotopes);
      G4Element*  hydrogen = nist->FindOrBuildElement("H", isotopes);
      G4Element*  carbon = nist->FindOrBuildElement("C", isotopes);
      G4Element*  phosphorus = nist->FindOrBuildElement("P", isotopes);
      G4Element*  nitrogen = nist->FindOrBuildElement("N", isotopes);
      G4Element*  sulphur = nist->FindOrBuildElement("S", isotopes);
      G4Element*  sodium = nist->FindOrBuildElement("Na", isotopes);
      G4Element*  chlorine = nist->FindOrBuildElement("Cl", isotopes);
      G4Element*  potassium = nist->FindOrBuildElement("K", isotopes);
      G4Element*  iron = nist->FindOrBuildElement("Fe", isotopes);
      G4Element*  magnesium = nist->FindOrBuildElement("Mg", isotopes);
      G4Element*  calcium = nist->FindOrBuildElement("Ca", isotopes);


      G4Isotope* Ru106 = new G4Isotope(name = "Ru106",
                           z = 44,    // atomic number
                           n = 62,    // number of nucleons
                           a = 105.907329*g/mole );  // mass of mole

      G4Element* Ru = new G4Element(name = "Ru", symbol = "bla", nIso = 1);
        Ru->AddIsotope(Ru106, 1.0);




      //Eye Materials including Aminosäuren und Zucker und co
      G4Material* Water = nist->FindOrBuildMaterial("G4_WATER");

      G4Material* H2O = new G4Material(name = "H2O", density = 0.998*g/cm3, ncomponents = 3);
      H2O->AddElement(hydrogen,2);
      H2O->AddElement(oxygen,1);

      G4Material* bicarbonate = new G4Material(name = "bicarbonate", density = 1.67*g/cm3, ncomponents = 3);
      bicarbonate->AddElement(hydrogen,1);
      bicarbonate->AddElement(oxygen,3);
      bicarbonate->AddElement(carbon,1);

      G4Material* lactate = new G4Material(name = "lactate", density = 1.209*g/cm3, ncomponents = 3);
      lactate->AddElement(hydrogen,3);
      lactate->AddElement(oxygen,3);
      lactate->AddElement(carbon,3);

      G4Material* glucose = new G4Material(name = "glucose", density = 1.54*g/cm3, ncomponents = 3);
      glucose->AddElement(hydrogen,12);
      glucose->AddElement(oxygen,6);
      glucose->AddElement(carbon,6);

      G4Material* ascorbate = new G4Material(name = "ascorbate", density = 1.65*g/cm3, ncomponents = 3);
      ascorbate->AddElement(hydrogen,7);
      ascorbate->AddElement(oxygen,6);
      ascorbate->AddElement(carbon,6);

      G4Material* phosphate = new G4Material(name = "phosphate", density = 1.87*g/cm3, ncomponents = 2);
      phosphate->AddElement(element = oxygen, natoms = 4);
      phosphate->AddElement(phosphorus,1);


     G4Material* citrate = new G4Material(name = "citrate", density = 1.665*g/cm3, ncomponents = 3);
     citrate->AddElement(element =hydrogen, natoms = 5);
     citrate->AddElement(element =oxygen, natoms = 7);
     citrate->AddElement(element =carbon, natoms = 6);

     G4Material* alanine= new G4Material(name = "alanine", density = 1.432*g/cm3, ncomponents = 4);
      alanine->AddElement(element =hydrogen ,natoms = 7);
      alanine->AddElement(element =oxygen ,natoms = 2);
      alanine->AddElement(element =carbon ,natoms = 3);
      alanine->AddElement(element =nitrogen ,natoms = 1);


     G4Material* leucine= new G4Material(name = "leucine", density = 1.293*g/cm3, ncomponents = 4);
      leucine->AddElement(element =hydrogen, natoms = 13);
      leucine->AddElement(element =oxygen, natoms = 2);
      leucine->AddElement(element =carbon, natoms = 6);
      leucine->AddElement(element =nitrogen, natoms = 1);


    G4Material* lysine= new G4Material(name = "lysine", density = 1.1*g/cm3, ncomponents = 4);
      lysine->AddElement(element =hydrogen ,natoms = 14);
      lysine->AddElement(element =oxygen ,natoms = 2 );
      lysine->AddElement(element =carbon ,natoms = 6 );
      lysine->AddElement(element =nitrogen ,natoms = 2 );


     G4Material* threonine= new G4Material(name = "threonine", density = 1.3*g/cm3, ncomponents = 4);
      threonine->AddElement(element =hydrogen ,natoms = 9 );
      threonine->AddElement(element =oxygen ,natoms = 3 );
      threonine->AddElement(element =carbon ,natoms = 4 );
      threonine->AddElement(element =nitrogen ,natoms = 1 );


     G4Material* valine= new G4Material(name = " valine", density = 1.23*g/cm3, ncomponents = 4);
       valine->AddElement( element =hydrogen ,natoms = 11);
       valine->AddElement(element =oxygen ,natoms = 2 );
       valine->AddElement(element =carbon ,natoms = 5 );
       valine->AddElement(element =nitrogen ,natoms = 1 );


     G4Material* pyruvicacid= new G4Material(name = "pyruvicacid", density = 1.25*g/cm3, ncomponents = 3);
      pyruvicacid->AddElement(element =hydrogen, natoms = 4);
      pyruvicacid->AddElement(element =oxygen, natoms = 3);
      pyruvicacid->AddElement(element =carbon, natoms = 3);


     G4Material* arginine= new G4Material(name = "arginine", density = 1.5*g/cm3, ncomponents = 4);
      arginine->AddElement( element =hydrogen,natoms = 14);
      arginine->AddElement(element =oxygen,natoms = 4 );
      arginine->AddElement(element =carbon,natoms = 6 );
      arginine->AddElement(element =nitrogen,natoms = 2 );


     G4Material* asparticacid= new G4Material(name = "asparticacid", density = 1.66*g/cm3, ncomponents = 4);
      asparticacid->AddElement(element =carbon,natoms = 4 );
      asparticacid->AddElement(element =hydrogen,natoms = 7 );
      asparticacid->AddElement(element =nitrogen,natoms = 1 );
      asparticacid->AddElement(element =oxygen,natoms = 4 );


     G4Material* cysteine= new G4Material(name = "cysteine", density = 1.3*g/cm3, ncomponents = 5);
      cysteine->AddElement(element =carbon, natoms = 3 );
      cysteine->AddElement(element =hydrogen, natoms = 7 );
      cysteine->AddElement(element =nitrogen, natoms = 1 );
      cysteine->AddElement(element =oxygen, natoms = 2 );
      cysteine->AddElement(element =sulphur, natoms = 1 );


     G4Material* glutamicacid= new G4Material(name = "glutamicacid", density = 1.525*g/cm3, ncomponents = 4);
      glutamicacid->AddElement(element =carbon ,natoms = 5 );
      glutamicacid->AddElement(element =hydrogen ,natoms = 9 );
      glutamicacid->AddElement(element =nitrogen ,natoms = 1 );
      glutamicacid->AddElement(element =oxygen ,natoms = 4 );


     G4Material* glycine= new G4Material(name = "glycine", density = 1.595*g/cm3, ncomponents = 4);
      glycine->AddElement(element =carbon, natoms = 2 );
      glycine->AddElement(element =hydrogen, natoms = 5 );
      glycine->AddElement(element =nitrogen, natoms = 1 );
      glycine->AddElement(element =oxygen, natoms = 2 );


    G4Material* histidine= new G4Material(name = "histidine", density = 1.4*g/cm3, ncomponents = 4);
      histidine->AddElement(element =carbon ,natoms = 6 );
      histidine->AddElement(element =hydrogen ,natoms = 9 );
      histidine->AddElement(element =nitrogen ,natoms = 3 );
      histidine->AddElement(element =oxygen ,natoms = 2 );


     G4Material* hydroxylysine= new G4Material(name = "hydroxylysine", density = 1.3*g/cm3, ncomponents = 4);
      hydroxylysine->AddElement(element =carbon, natoms = 6 );
      hydroxylysine->AddElement(element =hydrogen, natoms = 14);
      hydroxylysine->AddElement(element =nitrogen, natoms = 2 );
      hydroxylysine->AddElement(element =oxygen, natoms = 3 );


     G4Material* hydroxyproline= new G4Material(name = "hydroxyproline", density = 1.4*g/cm3, ncomponents = 4);
      hydroxyproline->AddElement(element =carbon, natoms = 5 );
      hydroxyproline->AddElement(element =hydrogen, natoms = 9 );
      hydroxyproline->AddElement(element =nitrogen, natoms = 1 );
      hydroxyproline->AddElement(element =oxygen, natoms = 3 );


     G4Material* isoleucine= new G4Material(name = "isoleucine", density = 1.4*g/cm3, ncomponents = 4);
      isoleucine->AddElement(element =carbon ,natoms = 6 );
      isoleucine->AddElement( element =hydrogen ,natoms = 13);
      isoleucine->AddElement(element =nitrogen ,natoms = 1 );
      isoleucine->AddElement(element =oxygen ,natoms = 2 );


     G4Material* methionine= new G4Material(name = "methionine", density = 1.34*g/cm3, ncomponents = 5);
      methionine->AddElement(element =carbon, natoms = 5 );
      methionine->AddElement( element =hydrogen, natoms = 11);
      methionine->AddElement(element =nitrogen, natoms = 1 );
      methionine->AddElement(element =oxygen, natoms = 2 );
      methionine->AddElement(element =sulphur, natoms = 1 );


     G4Material* phenylalanine= new G4Material(name = "phenylalanine", density = 1.2*g/cm3, ncomponents = 4);
      phenylalanine->AddElement(element =carbon ,natoms = 9 );
      phenylalanine->AddElement( element =hydrogen ,natoms = 11);
      phenylalanine->AddElement(element =nitrogen ,natoms = 1 );
      phenylalanine->AddElement(element =oxygen ,natoms = 2 );


     G4Material* proline= new G4Material(name = "proline", density = 1.36*g/cm3, ncomponents = 4);
      proline->AddElement(element =carbon ,natoms = 5 );
      proline->AddElement(element =hydrogen ,natoms = 9 );
      proline->AddElement(element =nitrogen ,natoms = 1 );
      proline->AddElement(element =oxygen ,natoms = 2 );


     G4Material* serine= new G4Material(name = "serine", density = 1.537*g/cm3, ncomponents = 4);
      serine->AddElement(element =carbon ,natoms = 3);
      serine->AddElement(element =hydrogen ,natoms = 7);
      serine->AddElement(element =nitrogen ,natoms = 1);
      serine->AddElement(element =oxygen ,natoms = 3);


     G4Material* tyrosine= new G4Material(name = "tyrosine", density = 1.46*g/cm3, ncomponents = 4);
      tyrosine->AddElement(element =carbon ,natoms = 9 );
      tyrosine->AddElement(element =hydrogen ,natoms = 11);
      tyrosine->AddElement(element =nitrogen ,natoms = 1 );
      tyrosine->AddElement(element =oxygen ,natoms = 3 );



     G4Material* collagen= new G4Material(name = "collagen ", density = 1.343*g/cm3, ncomponents = 19);
      collagen->AddMaterial(material =alanine, fractionmass =0.096*100*perCent);
      collagen->AddMaterial(material =arginine, fractionmass =0.046*100*perCent);
      collagen->AddMaterial(material =asparticacid, fractionmass =0.042*100*perCent);
      collagen->AddMaterial(material =cysteine, fractionmass =0.002*100*perCent);
      collagen->AddMaterial(material =glutamicacid, fractionmass =0.071*100*perCent);
      collagen->AddMaterial(material =glycine, fractionmass =0.350*100*perCent);
      collagen->AddMaterial(material =histidine, fractionmass =0.006*100*perCent);
      collagen->AddMaterial(material =hydroxylysine, fractionmass =0.005*100*perCent);
      collagen->AddMaterial(material =hydroxyproline, fractionmass =0.125*100*perCent);
      collagen->AddMaterial(material =isoleucine, fractionmass =0.013*100*perCent);
      collagen->AddMaterial(material =leucine, fractionmass =0.022*100*perCent);
      collagen->AddMaterial(material =lysine, fractionmass =0.030*100*perCent);
      collagen->AddMaterial(material =methionine, fractionmass =0.008*100*perCent);
      collagen->AddMaterial(material =phenylalanine, fractionmass =0.008*100*perCent);
      collagen->AddMaterial(material =proline, fractionmass =0.107*100*perCent);
      collagen->AddMaterial(material =serine, fractionmass =0.039*100*perCent);
      collagen->AddMaterial(material =threonine, fractionmass =0.013*100*perCent);
      collagen->AddMaterial(material =tyrosine, fractionmass =0.003*100*perCent);
      collagen->AddMaterial(material =valine, fractionmass =0.014*100*perCent);


    G4Material* elastin= new G4Material(name = "elastin", density = 1.3*g/cm3, ncomponents = 16 );
     elastin->AddMaterial(material =alanine , fractionmass =0.17915*100*perCent);
     elastin->AddMaterial(material =arginine , fractionmass =0.00506*100*perCent);
     elastin->AddMaterial(material =asparticacid , fractionmass =0.00405*100*perCent);
     elastin->AddMaterial(material =glutamicacid , fractionmass =0.01316*100*perCent);
     elastin->AddMaterial(material =glycine , fractionmass =0.36538*100*perCent);
     elastin->AddMaterial(material =histidine , fractionmass =0.00101*100*perCent);
     elastin->AddMaterial(material =hydroxyproline , fractionmass =0.01923*100*perCent);
     elastin->AddMaterial(material =isoleucine , fractionmass =0.01923*100*perCent);
     elastin->AddMaterial(material =leucine , fractionmass =0.05061*100*perCent);
     elastin->AddMaterial(material =lysine , fractionmass =0.00304*100*perCent);
     elastin->AddMaterial(material =phenylalanine , fractionmass =0.02024*100*perCent);
     elastin->AddMaterial(material =proline , fractionmass =0.12854*100*perCent);
     elastin->AddMaterial(material =serine , fractionmass =0.00506*100*perCent);
     elastin->AddMaterial(material =threonine , fractionmass =0.00709*100*perCent);
     elastin->AddMaterial(material =tyrosine , fractionmass =0.01113*100*perCent);
     elastin->AddMaterial(material =valine , fractionmass =0.16802*100*perCent);


    G4Material* blood= new G4Material(name = "blood", density = 1.06*g/cm3, ncomponents = 10);
     blood->AddElement(element =hydrogen , fractionmass =0.102*100*perCent );
     blood->AddElement(element =carbon , fractionmass =0.110*100*perCent );
     blood->AddElement(element =nitrogen , fractionmass =0.033*100*perCent );
     blood->AddElement(element =oxygen , fractionmass =0.745*100*perCent );
     blood->AddElement(element =sulphur , fractionmass =0.002*100*perCent );
     blood->AddElement(element =sodium , fractionmass =0.001*100*perCent );
     blood->AddElement(element =chlorine , fractionmass =0.003*100*perCent );
     blood->AddElement(element =phosphorus , fractionmass =0.001*100*perCent );
     blood->AddElement(element =potassium , fractionmass =0.002*100*perCent );
     blood->AddElement(element =iron , fractionmass =0.001*100*perCent );


    G4Material* bonecranium= new G4Material(name = "bonecranium ", density = 1.61*g/cm3, ncomponents = 9);
     bonecranium->AddElement(element =hydrogen , fractionmass =0.050*100*perCent);
     bonecranium->AddElement(element =carbon , fractionmass =0.212*100*perCent);
     bonecranium->AddElement(element =nitrogen , fractionmass =0.040*100*perCent);
     bonecranium->AddElement(element =oxygen , fractionmass =0.435*100*perCent);
     bonecranium->AddElement(element =sulphur , fractionmass =0.003*100*perCent);
     bonecranium->AddElement(element =sodium , fractionmass =0.001*100*perCent);
     bonecranium->AddElement(element =phosphorus , fractionmass =0.081*100*perCent);
     bonecranium->AddElement(element =magnesium , fractionmass =0.002*100*perCent);
     bonecranium->AddElement(element =calcium , fractionmass =0.176*100*perCent);


    G4Material* eyelens= new G4Material(name = "eyelens ", density = 1.07*g/cm3, ncomponents = 8);
     eyelens->AddElement(element =hydrogen , fractionmass =0.096*100*perCent);
     eyelens->AddElement(element =carbon , fractionmass =0.195*100*perCent);
     eyelens->AddElement(element =nitrogen , fractionmass =0.057*100*perCent);
     eyelens->AddElement(element =oxygen , fractionmass =0.646*100*perCent);
     eyelens->AddElement(element =sulphur , fractionmass =0.003*100*perCent);
     eyelens->AddElement(element =sodium , fractionmass =0.001*100*perCent);
     eyelens->AddElement(element =chlorine , fractionmass =0.001*100*perCent);
     eyelens->AddElement(element =phosphorus , fractionmass =0.001*100*perCent);


    G4Material* eyenerve= new G4Material(name = "eyenerve ", density = 1.04*g/cm3, ncomponents = 9);
     eyenerve->AddElement(element =hydrogen , fractionmass =0.107*100*perCent);
     eyenerve->AddElement(element =carbon , fractionmass =0.145*100*perCent);
     eyenerve->AddElement(element =nitrogen , fractionmass =0.022*100*perCent);
     eyenerve->AddElement(element =oxygen , fractionmass =0.712*100*perCent);
     eyenerve->AddElement(element =sulphur , fractionmass =0.002*100*perCent);
     eyenerve->AddElement(element =sodium , fractionmass =0.002*100*perCent);
     eyenerve->AddElement(element =chlorine , fractionmass =0.003*100*perCent);
     eyenerve->AddElement(element =phosphorus , fractionmass =0.004*100*perCent);
     eyenerve->AddElement(element =potassium , fractionmass =0.003*100*perCent);


    G4Material* muscle= new G4Material(name = "muscle", density = 1.05*g/cm3, ncomponents = 9);
     muscle->AddElement(element =hydrogen , fractionmass =0.102*100*perCent);
     muscle->AddElement(element =carbon , fractionmass =0.143*100*perCent);
     muscle->AddElement(element =nitrogen , fractionmass =0.034*100*perCent);
     muscle->AddElement(element =oxygen , fractionmass =0.710*100*perCent);
     muscle->AddElement(element =sulphur , fractionmass =0.003*100*perCent);
     muscle->AddElement(element =sodium , fractionmass =0.001*100*perCent);
     muscle->AddElement(element =chlorine , fractionmass =0.001*100*perCent);
     muscle->AddElement(element =phosphorus , fractionmass =0.002*100*perCent);
     muscle->AddElement(element =potassium , fractionmass =0.004*100*perCent);


    G4Material* skin= new G4Material(name = "skin", density = 1.09*g/cm3, ncomponents = 9);
     skin->AddElement(element =hydrogen , fractionmass =0.100*100*perCent);
     skin->AddElement(element =carbon , fractionmass =0.204*100*perCent);
     skin->AddElement(element =nitrogen , fractionmass =0.042*100*perCent);
     skin->AddElement(element =oxygen , fractionmass =0.645*100*perCent);
     skin->AddElement(element =sulphur , fractionmass =0.002*100*perCent);
     skin->AddElement(element =sodium , fractionmass =0.002*100*perCent);
     skin->AddElement(element =chlorine , fractionmass =0.003*100*perCent);
     skin->AddElement(element =phosphorus , fractionmass =0.001*100*perCent);
     skin->AddElement(element =potassium , fractionmass =0.001*100*perCent);


    G4Material* aqueoushumor= new G4Material(name = "aqueoushumor ", density = 1.00*g/cm3, ncomponents = 17);
     aqueoushumor->AddMaterial(material =H2O , fractionmass =0.98687*100*perCent);
     aqueoushumor->AddElement(element =sodium , fractionmass =0.00619*100*perCent);
     aqueoushumor->AddElement(element =chlorine , fractionmass =0.00541*100*perCent);
     aqueoushumor->AddMaterial(material =bicarbonate , fractionmass =0.00085*100*perCent);
     aqueoushumor->AddMaterial(material =lactate , fractionmass =0.00015*100*perCent);
     aqueoushumor->AddMaterial(material =glucose , fractionmass =0.00014*100*perCent);
     aqueoushumor->AddElement(element =potassium , fractionmass =0.00013*100*perCent);
     aqueoushumor->AddElement(element =calcium , fractionmass =0.00008*100*perCent);
     aqueoushumor->AddElement(element =magnesium , fractionmass =0.00005*100*perCent);
     aqueoushumor->AddMaterial(material =ascorbate , fractionmass =0.00004*100*perCent);
     aqueoushumor->AddMaterial(material =phosphate , fractionmass =0.00003*100*perCent);
     aqueoushumor->AddMaterial(material =citrate , fractionmass =0.00001*100*perCent);
     aqueoushumor->AddMaterial(material =alanine , fractionmass =0.00001*100*perCent);
     aqueoushumor->AddMaterial(material =leucine , fractionmass =0.00001*100*perCent);
     aqueoushumor->AddMaterial(material =lysine , fractionmass =0.00001*100*perCent);
     aqueoushumor->AddMaterial(material =threonine , fractionmass =0.00001*100*perCent);
     aqueoushumor->AddMaterial(material =valine , fractionmass =0.00001*100*perCent);


    G4Material* cornea= new G4Material(name = "cornea", density = 1.024*g/cm3, ncomponents = 13);
     cornea->AddMaterial(material =H2O , fractionmass =0.78000*100*perCent);
     cornea->AddMaterial(material =collagen , fractionmass =0.15000*100*perCent);
     cornea->AddElement(element =chlorine , fractionmass =0.01977*100*perCent);
     cornea->AddElement(element =sodium , fractionmass =0.01725*100*perCent);
     cornea->AddElement(element =oxygen , fractionmass =0.00985*100*perCent);
     cornea->AddElement(element =potassium , fractionmass =0.00719*100*perCent);
     cornea->AddElement(element =carbon , fractionmass =0.00660*100*perCent);
     cornea->AddElement(element =sulphur , fractionmass =0.00659*100*perCent);
     cornea->AddElement(element =hydrogen , fractionmass =0.00090*100*perCent);
     cornea->AddElement(element =nitrogen , fractionmass =0.00055*100*perCent);
     cornea->AddElement(element =phosphorus , fractionmass =0.00076*100*perCent);
     cornea->AddElement(element =calcium , fractionmass =0.00030*100*perCent);
     cornea->AddElement(element =magnesium , fractionmass =0.00024*100*perCent);


    G4Material* sclera= new G4Material(name = "sclera", density = 1.049*g/cm3, ncomponents = 14);
     sclera->AddMaterial(material =H2O , fractionmass =0.68000*100*perCent);
     sclera->AddMaterial(material =collagen , fractionmass =0.27000*100*perCent);
     sclera->AddMaterial(material =elastin , fractionmass =0.01000*100*perCent);
     sclera->AddElement(element =oxygen , fractionmass =0.00985*100*perCent);
     sclera->AddElement(element =chlorine , fractionmass =0.00920*100*perCent);
     sclera->AddElement(element =carbon , fractionmass =0.00660*100*perCent);
     sclera->AddElement(element =sodium , fractionmass =0.00583*100*perCent);
     sclera->AddElement(element =sulphur , fractionmass =0.00333*100*perCent);
     sclera->AddElement(element =potassium , fractionmass =0.00297*100*perCent);
     sclera->AddElement(element =hydrogen , fractionmass =0.00090*100*perCent);
     sclera->AddElement(element =nitrogen , fractionmass =0.00055*100*perCent);
     sclera->AddElement(element =phosphorus , fractionmass =0.00043*100*perCent);
     sclera->AddElement(element =calcium , fractionmass =0.00027*100*perCent);
     sclera->AddElement(element =magnesium , fractionmass =0.00007*100*perCent);


    G4Material* tumor= new G4Material(name = "tumor", density = 1.049*g/cm3, ncomponents = 9);
     tumor->AddElement(element =oxygen , fractionmass =0.61500*100*perCent);
     tumor->AddElement(element =chlorine , fractionmass =0.21200*100*perCent);
     tumor->AddElement(element =hydrogen , fractionmass =0.09400*100*perCent);
     tumor->AddElement(element =nitrogen , fractionmass =0.05600*100*perCent);
     tumor->AddElement(element =sulphur , fractionmass =0.00644*100*perCent);
     tumor->AddElement(element =potassium , fractionmass =0.00506*100*perCent);
     tumor->AddElement(element =phosphorus , fractionmass =0.00506*100*perCent);
     tumor->AddElement(element =chlorine , fractionmass =0.00391*100*perCent);
     tumor->AddElement(element =sodium , fractionmass =0.00253*100*perCent);


    G4Material* vitreoushumor= new G4Material(name = "vitreoushumor", density = 1.0065*g/cm3, ncomponents = 10);
     vitreoushumor->AddMaterial(material =H2O , fractionmass =0.99064*100*perCent);
     vitreoushumor->AddElement(element =chlorine , fractionmass =0.00431*100*perCent);
     vitreoushumor->AddElement(element =sodium , fractionmass =0.00337*100*perCent);
     vitreoushumor->AddMaterial(material =glucose , fractionmass =0.00054*100*perCent);
     vitreoushumor->AddMaterial(material =collagen , fractionmass =0.00040*100*perCent);
     vitreoushumor->AddMaterial(material =lactate , fractionmass =0.00036*100*perCent);
     vitreoushumor->AddElement(element =potassium , fractionmass =0.00022*100*perCent);
     vitreoushumor->AddElement(element =calcium , fractionmass =0.00007*100*perCent);
     vitreoushumor->AddMaterial(material =pyruvicacid , fractionmass =0.00007*100*perCent);
     vitreoushumor->AddElement(element =magnesium , fractionmass =0.00002*100*perCent);


    G4Material* choroidea= new G4Material(name = "choroidea", density = 1.002*g/cm3, ncomponents = 14);
     choroidea->AddMaterial(material =H2O , fractionmass =0.68000*100*perCent);
     choroidea->AddMaterial(material =collagen , fractionmass =0.27000*100*perCent);
     choroidea->AddMaterial(material =elastin , fractionmass =0.01000*100*perCent);
     choroidea->AddElement(element =oxygen , fractionmass =0.00985*100*perCent);
     choroidea->AddElement(element =chlorine , fractionmass =0.00920*100*perCent);
     choroidea->AddElement(element =carbon , fractionmass =0.00660*100*perCent);
     choroidea->AddElement(element =sodium , fractionmass =0.00583*100*perCent);
     choroidea->AddElement(element =sulphur , fractionmass =0.00333*100*perCent);
     choroidea->AddElement(element =potassium , fractionmass =0.00297*100*perCent);
     choroidea->AddElement(element =hydrogen , fractionmass =0.00090*100*perCent);
     choroidea->AddElement(element =nitrogen , fractionmass =0.00055*100*perCent);
     choroidea->AddElement(element =phosphorus , fractionmass =0.00043*100*perCent);
     choroidea->AddElement(element =calcium , fractionmass =0.00027*100*perCent);
     choroidea->AddElement(element =magnesium , fractionmass =0.00007*100*perCent);


    G4Material* netzhaut= new G4Material(name = "netzhaut", density = 1.008*g/cm3, ncomponents = 10);
     netzhaut->AddMaterial(material =H2O , fractionmass =0.99064*100*perCent);
     netzhaut->AddElement(element =chlorine , fractionmass =0.00431*100*perCent);
     netzhaut->AddElement(element =sodium , fractionmass =0.00337*100*perCent);
     netzhaut->AddMaterial(material =glucose , fractionmass =0.00054*100*perCent);
     netzhaut->AddMaterial(material =collagen , fractionmass =0.00040*100*perCent);
     netzhaut->AddMaterial(material =lactate , fractionmass =0.00036*100*perCent);
     netzhaut->AddElement(element =potassium , fractionmass =0.00022*100*perCent);
     netzhaut->AddElement(element =calcium , fractionmass =0.00007*100*perCent);
     netzhaut->AddMaterial(material =pyruvicacid , fractionmass =0.00007*100*perCent);
     netzhaut->AddElement(element =magnesium , fractionmass =0.00002*100*perCent);

     //hier stand in der .py Datei state = solid, weiß nicht, ob das wichtig ist
     G4Material* ruthenium = new G4Material(name = "ruthenium", density = 12.36*g/cm3, ncomponents = 1, state = kStateSolid);
      ruthenium->AddElement(element =Ru, fractionmass =1.*100*perCent);

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




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
        new G4Orb("World", 5000.0 * mm);


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
    /*G4Box* solidScorer =
        new G4Box("solidScorer",
                  0.25 * mm, 0.25 * mm, 0.25 * mm);


    G4Box* solidScorerBig =
        new G4Box("solidScorerBig",
                  0.5 * mm, 0.5 * mm, 0.5 * mm);
    */
    // Scorer 0.5 mm
    /*
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


    */
    ////////////////////
    // CADMesh :: STL //
    ////////////////////

    //auto bunny_mesh = CADMesh::TessellatedMesh::FromSTL("./bunny.stl");

    //auto bunny_logical = new G4LogicalVolume( bunny_mesh->GetSolid()
    //                                         , myWater
    //                                         , "logical"
    //                                         , 0, 0, 0
    //);

    //new G4PVPlacement( 0
    //                 , G4ThreeVector()
    //                 , bunny_logical
    //                 , "physical"
    //                 , logicWorld
    //                 , false, 0
    //);




    struct dirent *de;  // Pointer for directory entry

    // opendir() returns a pointer of DIR type.
    DIR *dr = opendir("/Users/smller/Simulationen/ccbPlaque/models/");

    if (dr == NULL)  // opendir returns NULL if couldn't open directory
    {
      G4cout << ("Could not open current directory") << G4endl;
      return 0;
    }

    // Refer http://pubs.opengroup.org/onlinepubs/7990989775/xsh/readdir.html
    // for readdir()
    std::vector<std::string> model_names;
    while ((de = readdir(dr)) != NULL){
      if(strstr(de->d_name,".stl") != NULL)
      {
          model_names.push_back(de->d_name);
      }
    }
    closedir(dr);



    //realsise and place all .stl objects in geant4
    if (model_names.size() == 0)
    {
        G4cout << ("No models found, empty physworld") << G4endl;
    }
    else
    {
        //Material bestimmen aus dem Namen der .stl File, später mit MatchMaterialtoSTL
        G4Material* logical_material;


        //std::vector<std::shared_ptr<CADMesh::TessellatedMesh>> model_mesh_vector; //works but decltype is better since ..FromSTL has a template that can be saved this way instead of saying it is always CADMEsh
        std::vector<decltype(CADMesh::TessellatedMesh::FromSTL("../ccbPlaque/models/"+ model_names[0]))> model_mesh_vector;
        std::vector<G4LogicalVolume*> model_logical_vector;

        for(unsigned int i = 0; i < model_names.size(); i++)
        {
            logical_material = DetectorConstruction::MatchMaterialToSTL(model_names[i]);
            model_mesh_vector.push_back(CADMesh::TessellatedMesh::FromSTL("../ccbPlaque/models/"+ model_names[i]));
            model_logical_vector.push_back(new G4LogicalVolume(model_mesh_vector[i]->GetSolid()
                                                , myWater
                                                , "logical"
                                                , 0, 0, 0));
            new G4PVPlacement( 0
              , G4ThreeVector()
              , model_logical_vector[i]
              , "physical"
              , logicWorld
              , false, 0
            );
        }
    }


    return physWorld;
}
//{}

//------------------------------------------------------------------------------
G4Material* DetectorConstruction::MatchMaterialToSTL(std::string model_name)
{
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* model_material;
    if (model_name.find("Sclera")){model_material = nist->FindOrBuildMaterial("sclera");}
    else if (model_name.find("Choroidea")){model_material = nist->FindOrBuildMaterial("choroidea");}
    else if (model_name.find("Netzhaut")){model_material = nist->FindOrBuildMaterial("netzhaut");}
    else if (model_name.find("Glaskörper")){model_material = nist->FindOrBuildMaterial("vitreoushumor");}
    else if (model_name.find("Cornea")){model_material = nist->FindOrBuildMaterial("cornea");}
    else if (model_name.find("Kammerwasser")){model_material = nist->FindOrBuildMaterial("aqueoushumor");}
    else if (model_name.find("Linse")){model_material = nist->FindOrBuildMaterial("eyelens");}
    else if (model_name.find("fovea")){model_material = nist->FindOrBuildMaterial("netzhaut");}
    else if (model_name.find("Fovea")){model_material = nist->FindOrBuildMaterial("netzhaut");}
    else if (model_name.find("Sehnerv")){model_material = nist->FindOrBuildMaterial("eyenerve");}
    else if (model_name.find("Applikator")){model_material = nist->FindOrBuildMaterial("Ru");}
    else if (model_name.find("Target")){model_material = nist->FindOrBuildMaterial("Ru");}
    else
    {
        G4cout<< "no material found for this eye component, instead water was used" << G4endl;
        model_material = nist->FindOrBuildMaterial("H20");
    }
    return model_material;
}
