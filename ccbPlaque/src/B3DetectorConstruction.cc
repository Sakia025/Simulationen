//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

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




int idontknow = -1;
//must be here, its global!
std::vector<std::string> model_names;
//G4cout << idontknow << G4endl;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::B3DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true)
{
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();

  G4bool isotopes_old = false;

  G4Element*  O = man->FindOrBuildElement("O" , isotopes_old);
  G4Element* Si = man->FindOrBuildElement("Si", isotopes_old);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes_old);

  G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
  LSO->AddElement(Lu, 2);
  LSO->AddElement(Si, 1);
  LSO->AddElement(O , 5);




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
  //G4Material* Water = nist->FindOrBuildMaterial("G4_WATER");

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

G4VPhysicalVolume* B3DetectorConstruction::Construct()
{
  // Gamma detector Parameters
  //
  G4double cryst_dX = 6*cm, cryst_dY = 6*cm, cryst_dZ = 3*cm;
  G4int nb_cryst = 32;
  G4int nb_rings = 9;
  //
  G4double dPhi = twopi/nb_cryst, half_dPhi = 0.5*dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  //
  G4double ring_R1 = 0.5*cryst_dY/tandPhi;
  G4double ring_R2 = (ring_R1+cryst_dZ)/cosdPhi;
  //
  G4double detector_dZ = nb_rings*cryst_dX;
  //
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2SiO5");


//meins
  //G4NistManager* nist = G4NistManager::Instance();
  G4Material* myWater = nist->FindOrBuildMaterial("G4_WATER");
  G4Material* mySilver = nist->FindOrBuildMaterial("G4_Ag");
  G4bool checkOverlaps = true;

  //
  // World
  //
/*
  G4double world_sizeXY = 2.4*ring_R2;
  G4double world_sizeZ  = 1.2*detector_dZ;

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        default_mat,         //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps
*/

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
  // ring
  //
  G4Tubs* solidRing =
    new G4Tubs("Ring", ring_R1, ring_R2, 0.5*cryst_dX, 0., twopi);

  G4LogicalVolume* logicRing =
    new G4LogicalVolume(solidRing,           //its solid
                        default_mat,         //its material
                        "Ring");             //its name

  //
  // define crystal
  //
  G4double gap = 0.5*mm;        //a gap for wrapping
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
  G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);

  G4LogicalVolume* logicCryst =
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name

  // place crystals within a ring
  //
  for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {
    G4double phi = icrys*dPhi;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg);
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
    G4ThreeVector position = (ring_R1+0.5*cryst_dZ)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);

    new G4PVPlacement(transform,             //rotation,position
                      logicCryst,            //its logical volume
                      "crystal",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      fCheckOverlaps);       // checking overlaps
  }

  //
  // full detector
  //
  G4Tubs* solidDetector =
    new G4Tubs("Detector", ring_R1, ring_R2, 0.5*detector_dZ, 0., twopi);

  G4LogicalVolume* logicDetector =
    new G4LogicalVolume(solidDetector,       //its solid
                        default_mat,         //its material
                        "Detector");         //its name

  //
  // place rings within detector
  //
  G4double OG = -0.5*(detector_dZ + cryst_dX);
  for (G4int iring = 0; iring < nb_rings ; iring++) {
    OG += cryst_dX;
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,OG), //position
                      logicRing,             //its logical volume
                      "ring",                //its name
                      logicDetector,         //its mother  volume
                      false,                 //no boolean operation
                      iring,                 //copy number
                      fCheckOverlaps);       // checking overlaps
  }

  //
  // place detector in world
  //
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicDetector,           //its logical volume
                    "Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps

  //
  // patient
  //
  G4double patient_radius = 8*cm;
  G4double patient_dZ = 10*cm;
  G4Material* patient_mat = nist->FindOrBuildMaterial("G4_BRAIN_ICRP");

  G4Tubs* solidPatient =
    new G4Tubs("Patient", 0., patient_radius, 0.5*patient_dZ, 0., twopi);

  G4LogicalVolume* logicPatient =
    new G4LogicalVolume(solidPatient,        //its solid
                        patient_mat,         //its material
                        "PatientLV");        //its name

  //
  // place patient in world
  //
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicPatient,            //its logical volume
                    "Patient",               //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps


//meins
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Funktion draus machen

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
    //model_names defined and saved in global_variables.hh, so it can be used in EventAction.cc
    extern int wastuichhier;
    std::cout << wastuichhier << std::endl;
    //std::vector<std::string> model_names;
    while ((de = readdir(dr)) != NULL){
      if(strstr(de->d_name,".stl") != NULL)
      {
          model_names.push_back(de->d_name);
      }
    }
    std::vector<int> model_names_test(5, 0);
    std::vector<int>* model_names_ptr = &model_names_test;
    G4cout << model_names_ptr->size() << G4endl;
    closedir(dr);

  //extern std::vector<std::string> model_names;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

    //realsise and place all .stl objects in geant4
    if (model_names.empty())
    {
        G4cout << ("No models found, empty physworld") << G4endl;
        return 0;
    }
    else
    {
        //Material bestimmen aus dem Namen der .stl File, später mit MatchMaterialtoSTL
        G4Material* logical_material;


        //std::vector<std::shared_ptr<CADMesh::TessellatedMesh>> model_mesh_vector; //works but decltype is better since ..FromSTL has a template that can be saved this way instead of saying it is always CADMEsh
        std::vector<decltype(CADMesh::TessellatedMesh::FromSTL("../ccbPlaque/models/"+ model_names[0]))> model_mesh_vector;
        std::vector<G4LogicalVolume*> model_logical_vector;
        std::vector<G4MultiFunctionalDetector*> my_Scorers_vector;

        for(unsigned int i = 0; i < model_names.size(); i++)
        {
           //construct Logical Volumes from .stl files via CADMeesh 2.0
            logical_material = B3DetectorConstruction::MatchMaterialToSTL(model_names[i]);
            model_mesh_vector.push_back(CADMesh::TessellatedMesh::FromSTL("../ccbPlaque/models/"+ model_names[i]));
            model_logical_vector.push_back(new G4LogicalVolume(model_mesh_vector[i]->GetSolid()
                                                , logical_material
                                                , "logical"
                                                , 0, 0, 0));



           //MultiFunctional Detector(link: https://agenda.infn.it/event/5981/sessions/10089/attachments/43814/52014/SensitiveDetector_alghero2013.pdfg)
           //!!!Do NOT change the MultiDetector name, if you do you need to change it in EventAction.cc as well!!!!!
            my_Scorers_vector.push_back(new G4MultiFunctionalDetector("my_" + std::string(model_names[i]) + "_Scorer"));
            G4SDManager::GetSDMpointer() -> AddNewDetector(my_Scorers_vector[i]);
            //kann hier oder auch unten stehen, eine version im pdf eine in examples
            model_logical_vector[i]->SetSensitiveDetector(my_Scorers_vector[i]);
            //in dem oben genannten Link steht nur -> Register und PrimitiveSensitivity,
            //die sind jewils nicht bekannt , PrimitiveScorer und RegisterPrimitive
            //jedoch schon, zu sehen in example B3a, daher auch SetSensitiveDetector,
            //nicht sicher, ob das muss und was es tut
            G4VPrimitiveScorer* totalDose = new G4PSDoseDeposit("TotalDose");
            my_Scorers_vector[i]->RegisterPrimitive(totalDose);




            //Place Volume in the World
            new G4PVPlacement( 0
              , G4ThreeVector()
              , model_logical_vector[i]
              , "physical"
              , logicWorld
              , false, 0
            );
        }
    }




  // Visualization attributes
  //
  logicRing->SetVisAttributes (G4VisAttributes::GetInvisible());
  logicDetector->SetVisAttributes (G4VisAttributes::GetInvisible());

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare crystal as a MultiFunctionalDetector scorer
  //
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV",cryst);

  // declare patient as a MultiFunctionalDetector scorer
  //
  G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
  G4SDManager::GetSDMpointer()->AddNewDetector(patient);
  G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
  patient->RegisterPrimitive(primitiv2);
  SetSensitiveDetector("PatientLV",patient);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4Material* B3DetectorConstruction::MatchMaterialToSTL(std::string model_name)
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
