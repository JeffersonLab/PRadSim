//===================================================//
//    Event Generator for ep and ee scattering       //
//    Target: Hydrogen gas                           //
//    Incident beam: electron beam                   //
//                                                   //
//    Developer: Chao Peng                           //
//    Date: 22-Jan-2013                              //
//===================================================//

#include <iostream>
#include <fstream>
#include <string>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TFile.h>
#include <TImage.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TString.h>
#include <TLatex.h>
#include <TRandom.h>

#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom.h>

using namespace std;

void BGanalysis() 
{
// gROOT->Reset(); 
gRandom->SetSeed();

ifstream datafile1("./Recon_1G_150umKap_35cm.dat");

hfile = new TFile("noCol_150umKap.root","RECREATE","acsii");
//Tree = new TTree("T","ascii",0);
//Tree->Branch("E",    &E,     "E/D");

C1_E = new TH1F("Energy_NoTrigger", "No Collimator",  110, 0., 1100.);
C1_A = new TH1F("Angle_NoTrigger", "No Collimator",  16, 0., 4.);
C2_E = new TH1F("Energy_Trigger", "No Collimator",  110, 0., 1100.);
C2_A = new TH1F("Angle_Trigger", "No Collimator",  16, 0., 4.);

 double Ecut = 1000., Anglecut = 0.;
 double E, angle, deg;
 double temp1, temp2, temp3, temp4;

 deg = 180./3.14159265358979;

 while (true) {
   datafile1 >> temp1 >> temp2;
   if( datafile1.eof() ) break;
   C1_E->Fill(temp2);
   C1_A->Fill(temp1);
   if(temp2 >= 500) {
     C2_E->Fill(temp2);
     C2_A->Fill(temp1);
   }
 }

hfile->Write();
hfile->Close();
}
