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

void ReconstructionAnalysis() 
{
  gROOT->Reset(); 

  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  double x, y, rcx, rcy, E;
  double x_tr, y_tr, r_tr, r_last, theta, Final_x, Final_y;
  double HyCal_z, Target_z;
  double HyCal_x[100], HyCal_y[100], HyCal_E[100], GEM_x[100], GEM_y[100], GEM_z[100];
  int xindex, yindex, n, GEM_n, HyCal_n, n_tr, i_tr;
  double DEG = 180./3.14159265358979;

  Target_z = -3000.;
  HyCal_z = 2650.;

  ifstream HyCalData;
  ifstream GEMData;
  ifstream Events;

  HyCalData.open("./Recon.dat");
  GEMData.open("../GEMPosition.dat");
  hfile = new TFile("ReconData.root", "RECREATE", "acsii");

  TTree* Tree = new TTree("T_HyCal", "Hits reconstructed on HyCal");
  TTree* Tree2 = new TTree("T_Coincidence", "Hits by HyCal and GEM");

  Tree->Branch("rcx",       &rcx,           "rcx/D");
  Tree->Branch("rcy",       &rcy,           "rcy/D");
  Tree->Branch("E",             &E,             "E/D");
  Tree->Branch("theta",     &theta,         "theta/D");

  Tree2->Branch("theta", &theta, "theta/D");
  Tree2->Branch("E", &E, "E/D");
  Tree2->Branch("x", &Final_x, "x/D");
  Tree2->Branch("y", &Final_y, "y/D");

  Angle = new TH1F("Theta", "N", 1000, 0.5,  1.3);
  Energy = new TH1F("Energy", "N", 2000, 0., 2000.);
  AngEne = new TH2F("Reconstruction", "Reconstructed Hits", 1000, 0., 10., 1300, 0., 1300.);
  ReconXY = new TH2F("PositionProjection", "Reconstructed Hits on HyCal", 1400, -700., 700., 1400, -700., 700.);

  Original_AngEne = new TH2F("Generated_Events", "Generated Events", 1000, 0., 10., 1300, 0., 1300.);
  Original_XY = new TH2F("GeneratedEvents_Position on HyCal", "Generated Events (Projection on HyCal)", 1200, -600., 600., 1200, -600., 600.);

  while(true) {
    HyCal_n = 0;
    do{
      HyCalData >> tmp1 >> tmp2 >> tmp3;
      if(HyCalData.eof()) break;
      if(tmp1 != -9999) {
        HyCal_x[HyCal_n] = tmp1;
        HyCal_y[HyCal_n] = tmp2;
        HyCal_E[HyCal_n] = tmp3;
        HyCal_n += 1;
      }
    }while(tmp1 != -9999);
//    cout << tmp2 << endl;
    if(HyCalData.eof()) break;


    n = 0;
    do {
      GEMData >> tmp4 >> tmp5 >> tmp6;
      if(GEMData.eof()) {cout << "Warning!" << endl; break;}
      if(tmp4 != 0) {
        GEM_x[n] = tmp4;
        GEM_y[n] = tmp5;
        GEM_z[n] = tmp6;
        n += 1;
      }
    } while(tmp4 != 0);



  for(int l = 0; l < HyCal_n; ++l) {
    rcx = 10.*HyCal_x[l];
    rcy = 10.*HyCal_y[l];
    E = HyCal_E[l];
    theta = atan(sqrt(rcx*rcx + rcy*rcy)/(HyCal_z - Target_z))*DEG;
//    cout << rcx << endl;
    Tree->Fill();

    n_tr = 0;  
    r_last = 1000.; //350 for central part, 900 for lead glass part

/*
    if(n > 0) {
      for(int i = 0; i < n; ++i) {
         x_tr = GEM_x[i] + (-GEM_x[i] + rcx)/(-HyCal_z + GEM_z)*(GEM_z - Target_z);
         y_tr = GEM_y[i] + (-GEM_y[i] + rcy)/(-HyCal_z + GEM_z)*(GEM_z - Target_z);
         r_tr = sqrt(x_tr*x_tr + y_tr*y_tr);
//         cout << rcx << "  " << rcy << "  " << GEM_x[i] << "  " << GEM_y[i] << "  " << r_tr << endl;
         if(r_tr < 1000.) { // 350 for central part, 900 for lead glass part
           n_tr += 1;
           if(r_tr < r_last) {
             i_tr = i; 
             r_last = r_tr;
             theta = atan(sqrt(GEM_x[i]*GEM_x[i] + GEM_y[i]*GEM_y[i])/(GEM_z - Target_z))*DEG;
           }
         }     
      }
    }
*/

  if(n > 0) {
    for(int i = 0; i < n; ++i) {
      x_tr = GEM_x[i]/(GEM_z[i] - Target_z)*(HyCal_z - Target_z);
      y_tr = GEM_y[i]/(GEM_z[i] - Target_z)*(HyCal_z - Target_z);
//      ReconXY->Fill(x_tr - rcx, y_tr - rcy);
      if(fabs(x_tr - rcx) < 40. && fabs(y_tr - rcy) < 40.) {
        n_tr += 1;
        r_tr = sqrt((x_tr - rcx)*(x_tr - rcx) + (y_tr - rcy)*(y_tr - rcy));
        if(r_tr < r_last) {
          i_tr = i;
          r_last = r_tr;
          theta = atan(sqrt(GEM_x[i]*GEM_x[i] + GEM_y[i]*GEM_y[i])/(GEM_z[i] - Target_z))*DEG;
          Final_x = x_tr;
          Final_y = y_tr;
        }
      }
    }
  }

  if(n_tr != 0){  if(theta >= 0.) Angle->Fill(theta);  Energy->Fill(E);  AngEne->Fill(theta,E); ReconXY->Fill(Final_x, Final_y); Tree2->Fill();}
  }   
  }
  
//  EnergyDeposit.close();
//  TruePosition.close();
  HyCalData.close();
  GEMData.close();

  /*
  Events.open("../RCEP.dat");
  double ev_theta[3], ev_phi[3], ev_energy[3];
  double ev_x, ev_y;
//  double DEG = 180./3.141592653589793;
  while (true) {
    Events >> ev_energy[0] >> ev_theta[0] >> ev_phi[0] >> ev_energy[1] >> ev_theta[1] >> ev_phi[1] >> ev_energy[2] >> ev_theta[2] >> ev_phi[2];
    if( Events.eof() ) break;
    for(int i = 0 ; i < 3; ++i) {
      if(ev_energy[i] > 0) {
        Original_AngEne->Fill(ev_theta[i]*DEG, ev_energy[i]);
        ev_x = tan(ev_theta[i])*(HyCal_z - Target_z)*cos(ev_phi[i]);
        ev_y = tan(ev_theta[i])*(HyCal_z - Target_z)*sin(ev_phi[i]);
        Original_XY->Fill(ev_x, ev_y);
      }
    }
  }
*/

  gStyle->SetPadLeftMargin(0.08);  // Margin left axis 
  gStyle->SetPadRightMargin(0.16); // Margin for palettes in 2D histos
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPalette(1);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptStat(0);

  TCanvas *C = new TCanvas("C","",1200,1200);
  C->Divide(2,2);
  C->cd(1);
  gPad->SetLogz();
  AngEne->GetXaxis()->SetTitle("Scat. Angle (degree)");
  AngEne->GetYaxis()->SetTitle("Energy (MeV)");
  AngEne->GetYaxis()->SetTitleOffset(1.15);
  AngEne->Draw("COLZ");
  C->cd(2);
  gPad->SetLogz();
  ReconXY->GetXaxis()->SetTitle("x (mm)");
  ReconXY->GetYaxis()->SetTitle("y (mm)");
  ReconXY->GetYaxis()->SetTitleOffset(1.15);
  ReconXY->Draw("COLZ");
  /*
  C->cd(3);
  gPad->SetLogz();
  Original_AngEne->GetXaxis()->SetTitle("Scat. Angle (degree)");
  Original_AngEne->GetYaxis()->SetTitle("Energy (MeV)");
  Original_AngEne->GetYaxis()->SetTitleOffset(1.15);
  Original_AngEne->Draw("COLZ");
  C->cd(4);
  gPad->SetLogz();
  Original_XY->GetXaxis()->SetTitle("x (mm)");
  Original_XY->GetYaxis()->SetTitle("y (mm)");
  Original_XY->GetYaxis()->SetTitleOffset(1.15);
  Original_XY->Draw("COLZ");
  */
  C->Print("ReconPlots.pdf");

  hfile->Delete("T_HyCal;1");
  hfile->Write();  
 // hfile->Close();
}
