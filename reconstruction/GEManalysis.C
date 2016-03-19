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
#include <math.h>

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

using namespace std;

void GEManalysis() 
{
  gROOT->Reset(); 

  double tmp1, tmp2, tmp3, tmp4, tmp5;
  double x, y, rcx, rcy, E;
  double x_tr, y_tr, r_tr, r_last, theta;
  double HyCal_z, GEM_z, Target_z;
  double GEM_x[100], GEM_y[100];
  int xindex, yindex, n, GEM_n, n_tr, i_tr, Nhits;
  double DEG = 180./3.14159265358979;

  Target_z = -3000.;
  HyCal_z = 2410.;
  GEM_z = 2220.;
  
  ifstream GEMData;
  GEMData.open("../GEM_32s_FandT.dat");
  hfile = new TFile("GEMData.root", "RECREATE", "acsii");
  Tree = new TTree("T_all", "ascii", 0);

  Tree->Branch("Nhits",   &Nhits,    "Nhits/I");

  Angle = new TH1F("Theta", "N", 1000, 0.,  10.);
  Energy = new TH1F("Energy", "N", 2400, 0., 1200.);
  AngEne = new TH2F("Reconstruction", "N", 1000, 0., 10., 1500, 0., 1500.);

  n_tr = 0;

  while(true) {
    GEMData >> tmp1 >> tmp2 >> tmp3 >> tmp4;
    if(GEMData.eof()) break;
    if(tmp1 != 0) Energy->Fill(tmp4);
    else {Nhits = (int)tmp4; n_tr += Nhits; Tree->Fill(); }
  } 

  
//  EnergyDeposit.close();
//  TruePosition.close();
  cout << n_tr << endl;
  GEMData.close();
  hfile->Write();  
  hfile->Close();
}
