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
#include <sstream>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <TObject.h>
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
#include <TPaletteAxis.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom.h>
#include "posrecon.h"

using namespace std;

int Eventviewer(Int_t Number) 
{
// gROOT->Reset(); 
gRandom->SetSeed();
ifstream HyCaldata("../EnergyDeposit.dat");

hfile = new TFile("EView.root","RECREATE","acsii");

hpxpy = new TH2F("hpxpy", "Events Viewer", 2400, -600., 600., 2400 ,-600., 600.);
TLegend* leg = new TLegend(0.80, 0.90, 0.99, 0.99);


 double tmp1, tmp2;
 int BlockNo, Nfill, EventNb, Show;
 double Energy, TotalE, ShownE;
 double x, y;
 double xmin, ymin, xmax, ymax;
 double step = 1200./2400.; 
 int Nx = (int)(697./step);
 int Ny = (int)(697./step);

 Show = 1;
 if(Number) Show = Number;

// draw frames
 for(int i = 0; i < Nx; ++i) {
   for(int j = 0; j < 35; ++j) {
     x = -348.5 + step*i;
     y = -348.5 + 20.5*(double)j;
     if(fabs(x) >= 20.5 || fabs(y) >= 20.5)hpxpy->Fill(x,y);
   }
 }
 for(int i = 0; i < Nx; ++i) {
   for(int j = 0; j < 35; ++j) {
     y = -348.5 + step*i;
     x = -348.5 + 20.5*(double)j;
     if(fabs(x) >= 20.5 || fabs(y) >= 20.5)hpxpy->Fill(x,y);
   }
 }
 
 Nx = (int)(38.2*24./step);
 for(int i = 0; i < Nx; ++ i) {
   for(int j = 0; j < 7; ++ j) {
     x = -(348.5 - 38.2*24. + (double)step*i);
     y = -(348.5 + 38.2*(double)j);
     if(fabs(x) > 348.5 || fabs(y) > 348.5)hpxpy->Fill(x,y);
   }
 }
 Ny = (int)(38.2*6./step);
 for(int i = 0; i < 25; ++ i) {
   for(int j = 0; j < Ny; ++ j) {
     x = -(348.5 - 38.2*(double)(24-i));
     y = -(348.5 + (double)j*step);
     hpxpy->Fill(x,y);
   }
 }

 Nx = (int)(38.2*24./step);
 for(int i = 0; i < Nx; ++ i) {
   for(int j = 0; j < 7; ++ j) {
     x = 348.5 - 38.2*24. + (double)step*i;
     y = 348.5 + 38.2*(double)j;
     if(fabs(x) > 348.5 || fabs(y) > 348.5)hpxpy->Fill(x,y);
   }
 }
 Ny = (int)(38.2*6./step);
 for(int i = 0; i < 25; ++ i) {
   for(int j = 0; j < Ny; ++ j) {
     x = 348.5 - 38.2*(double)(24-i);
     y = 348.5 + (double)j*step;
     hpxpy->Fill(x,y);
   }
 }

  Nx = (int)(38.2*24./step);
 for(int i = 0; i < Nx; ++ i) {
   for(int j = 1; j < 7; ++ j) {
     y = 348.5 - 38.2*24. + (double)step*i;
     x = -(348.5 + 38.2*(double)j);
     if(fabs(x) > 348.5 || fabs(y) > 348.5)hpxpy->Fill(x,y);
   }
 }
 Ny = (int)(38.2*6./step);
 for(int i = 0; i < 25; ++ i) {
   for(int j = 0; j < Ny; ++ j) {
     y = 348.5 - 38.2*(double)(24-i);
     x = -(348.5 + (double)j*step);
     if(i < 24)hpxpy->Fill(x,y);
     else if(fabs(x) > 38.2*24.-348.5) hpxpy->Fill(x,y);
   }
 }

 Nx = (int)(38.2*24./step);
 for(int i = 0; i < Nx; ++ i) {
   for(int j = 1; j < 7; ++ j) {
     y = -(348.5 - 38.2*24. + (double)step*i);
     x = (348.5 + 38.2*(double)j);
     if(fabs(x) > 348.5 || fabs(y) > 348.5)hpxpy->Fill(x,y);
   }
 }
 Ny = (int)(38.2*6./step);
 for(int i = 0; i < 25; ++ i) {
   for(int j = 0; j < Ny; ++ j) {
     y = -(348.5 - 38.2*(double)(24-i));
     x = (348.5 + (double)j*step);
     if(i < 24)hpxpy->Fill(x,y);
     else if(fabs(x) > 38.2*24.-348.5) hpxpy->Fill(x,y);
   }
 }


 
 EventNb = 1;
 while(true) {
    TotalE = 0.;
    do {
      HyCaldata >> tmp1 >> tmp2;
      if(HyCaldata.eof()) break;
      BlockNo = tmp1 - 1;
      if(EventNb == Show) {
       if( tmp1 != 0 && tmp2 > 0.) {
        Nfill = (int)(tmp2 + 0.5);
        if(BlockNo < 1152) {
          xmin = -GetXCoord(BlockNo)*10. - 10.25;
          ymin = GetYCoord(BlockNo)*10. - 10.25;
          xmax = -GetXCoord(BlockNo)*10. + 10.25;
          ymax = GetYCoord(BlockNo)*10. + 10.25;
          Nx = (int)(20.5/step);
          Ny = (int)(20.5/step);
        }
        else {
          xmin = -GetXCoord(BlockNo)*10. - 19.10;
          ymin = GetYCoord(BlockNo)*10. - 19.10;
          xmax = -GetXCoord(BlockNo)*10. + 19.10;
          ymax = GetYCoord(BlockNo)*10. + 19.10;
          Nx = (int)(38.2/step);
          Ny = (int)(38.2/step);
        }
        for(int i = 0; i < Nfill; ++i) {
          for(int j = 1; j < Nx; ++j) {
            for(int k = 1; k < Ny; ++k) {
              x = xmin + j*step;
              y = ymin + k*step;
              hpxpy->Fill(x, y); 
            }
          }         
        }
//        cout << BlockNo <<"  " << Nfill << endl;
        TotalE += tmp2;
       }
      }
    } while (tmp1 != 0);

    if(EventNb == Show) {ShownE = TotalE; cout << "Event " << Show << " is shown, its total energy is " << TotalE << " MeV." << endl;}
    if(EventNb > Show) break;
    if(HyCaldata.eof()) {break;}
    EventNb += 1;
  }
  if(Show > EventNb - 1) cout << "###Error### : Cannot find this event, the total events number is " << EventNb - 1 << endl;
// Tree->Fill();

  //Draw, set style
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

  //Canvas and axis setting
  TCanvas *C = new TCanvas("C","",600,600);
  C->cd();
  C->SetLogz();

  hpxpy->GetXaxis()->SetTitle("X (mm)");
  hpxpy->GetXaxis()->SetTitleSize(0.024);
  hpxpy->GetXaxis()->SetLabelSize(0.024);
  hpxpy->GetYaxis()->SetTitle("Y (mm)");
  hpxpy->GetYaxis()->SetTitleSize(0.024);
  hpxpy->GetYaxis()->SetTitleOffset(1.45);
  hpxpy->GetYaxis()->SetLabelSize(0.024);
  hpxpy->GetZaxis()->SetTitle("Energy (MeV)");
  hpxpy->GetZaxis()->SetTitleSize(0.024);
  hpxpy->GetZaxis()->SetLabelSize(0.024);
  hpxpy->GetZaxis()->SetTitleOffset(1.45);
  hpxpy->GetZaxis()->SetRangeUser(0.5, 300.);
 
  //Legend
  stringstream ss;
  ss << "Events  " << Number;
  string s = ss.str();
  const char*  c = s.c_str();
  leg->AddEntry((TObject*)0, c, "1");  
  ss.str("");
  ss << ShownE << " MeV ";
  s = ss.str();
  c = s.c_str();
  leg->AddEntry((TObject*)0, c, "1");

  hpxpy->Draw("COLZ");
  leg->Draw();
/*
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hpxpy->GetListOfFunctions()->FindObject("palette");
  palette->GetAxis()->SetLabelSize(0.024);
*/
  C->Print("h2D.pdf");
  return 0;
  hfile->Write();
  hfile->Close();
}
