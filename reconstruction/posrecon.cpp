#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream> 
#include "posrecon.h"

using namespace std;

int main(){

  int ntotal, tmp1, BlockNo, xcheck, nid, Ci;
  int ClusterCenter, Group, xid, yid, NewCluster, ExistCluster[10];
  int ClusterN, Cluster_BlockN[50];
  int NbEvents = 0;
  double ClusterEnergy, CenterEnergy;
  double TotEne, reconx[10], recony[10], reconE[10], weightx, weighty, weight, Totalweight;
  double tmp2, tmp3, tmp4;
  double ClusterRadius, Distance;
  double RadL_Crystal = 0.8903;
  double RadL_Leadglass = 2.577;
  double Moliere_Crystal = 2.05;
  double Moliere_Leadglass = 3.82;
  double MoliereRatio = Moliere_Leadglass/Moliere_Crystal;
  double baseR = 6.0;

//  ifstream enefile("../100M_Col.dat");
  ifstream enefile("../EnergyDeposit.dat");
//  ifstream trueposition("../position2.dat");
  ofstream output("Recon.dat");

  double averagex = 0., averagey = 0.;
  for(int j = 0; j < 1728; ++j) {
      Position.x[j] = GetXCoord(j);
      Position.y[j] = GetYCoord(j);
   }
  
  while(true){
    for(int j = 0; j < 1728; ++j) {
      Energy[j] = 0.;
    }
    
    TotEne = 0.;
    ClusterCenter = 0;
    CenterEnergy = 0.;
    NewCluster = 1;
    Ci = 0;

    do {
      enefile >> tmp1 >> tmp2;
      if(enefile.eof()) break;
      BlockNo = tmp1 - 1;
       if( tmp1 != 0 ) { 
        Energy[BlockNo] = tmp2;
        if(tmp2 > CenterEnergy) {CenterEnergy = tmp2; ClusterCenter = BlockNo;}
        TotEne += tmp2;
      }
    } while (tmp1 != 0); 

    if(enefile.eof()) break;
//ClusterCenter = 1;
//cout << GetXCoord(ClusterCenter) <<"  " << GetYCoord(ClusterCenter) <<endl;
  for(int t = 0; t < 10; t ++) {

    ClusterEnergy = 0.;

    if(Ci > 0) {
      NewCluster = 0;
      tmp1 = 0;
      tmp2 = 0.;
      for(int j = 0; j < 1728; ++j) {
        if(Energy[j] > 0.) {
        ClusterRadius = baseR;
        if(j > 1152) ClusterRadius = MoliereRatio*baseR;
        Distance = 120.;
          for(int k = 0; k < Ci; ++k) {
            tmp1 = ExistCluster[k];
            Distance = fmin(Distance, sqrt((Position.x[j] - Position.x[tmp1])*(Position.x[j] - Position.x[tmp1]) + (Position.y[j] - Position.y[tmp1])*(Position.y[j] - Position.y[tmp1])));
          }
//        cout << Ci << "  " << Distance << endl;
          if(Energy[j] > tmp2 && Distance >= 2.*ClusterRadius) {tmp2 = Energy[j]; ClusterCenter = j; CenterEnergy = Energy[j]; NewCluster += 1;}
        }
      }
    }                                    // search if there are another cluster

/*
    cout << Ci << "  " << ClusterCenter << "  " << Energy[ClusterCenter] << endl;
    for(int k = 0; k < Ci; k ++) {
    cout << ExistCluster[k] << "  "  << GetDistance(ExistCluster[k], ClusterCenter) << "  " << Distance << endl;
    }
*/

    if(NewCluster == 0 || CenterEnergy < 10.) break; // No more clusters

    // Find which group the cluster center is
    if(ClusterCenter <= 1152) Group = 0;  // 0, crystal part
    if(ClusterCenter > 1152) Group = (ClusterCenter - 1153)/144; // 1~4, 4 groups of Pb-glass part

    ClusterN = 0; // Reset number of blocks in this cluster
    // Crystal part, find the blocks in this cluster
    if(Group == 0) {
      ClusterRadius = baseR; // slightly larger than sqrt(2)*2*2.05, but less than 3*2.05
      for(int j = 0; j < 1152; j++) {
        if(Energy[j]>0.&&(((Position.x[j] - Position.x[ClusterCenter])*(Position.x[j] - Position.x[ClusterCenter]) + 
        (Position.y[j] - Position.y[ClusterCenter])*(Position.y[j] - Position.y[ClusterCenter])) <= ClusterRadius*ClusterRadius)) {
        Cluster_BlockN[ClusterN] = j;
        ClusterN += 1;
        ClusterEnergy += Energy[j];
//        cout << j << "  " ; 
//        cout << GetXCoord(j) <<"  " << GetYCoord(j) <<endl;}
        }
      }
/*
      if(Position.x[ClusterCenter] - baseR < -34.85) ClusterRadius = fmax(ClusterRadius, MoliereRatio*baseR - (MoliereRatio - 1.)*(Position.x[ClusterCenter] + 34.85));
      if(Position.x[ClusterCenter] + baseR > 34.85) ClusterRadius = fmax(ClusterRadius, MoliereRatio*baseR - (MoliereRatio - 1.)*(Position.x[ClusterCenter] - 34.85));
      if(Position.y[ClusterCenter] - baseR < -34.85) ClusterRadius = fmax(ClusterRadius, MoliereRatio*baseR - (MoliereRatio - 1.)*(Position.y[ClusterCenter] + 34.85));
      if(Position.y[ClusterCenter] + baseR > 34.85) ClusterRadius = fmax(ClusterRadius, MoliereRatio*baseR - (MoliereRatio - 1.)*(Position.y[ClusterCenter] - 34.85));   
*/    
      ClusterRadius = baseR*MoliereRatio;
      for(int j = 1152; j < 1728; j++) {
        if(Energy[j]>0.&&(((Position.x[j] - Position.x[ClusterCenter])*(Position.x[j] - Position.x[ClusterCenter]) +
        (Position.y[j] - Position.y[ClusterCenter])*(Position.y[j] - Position.y[ClusterCenter])) <= ClusterRadius*ClusterRadius)) {
          Cluster_BlockN[ClusterN] = j;
          ClusterN += 1;
          ClusterEnergy += Energy[j];
        } 
      }

    }

    if(Group != 0) {
      ClusterRadius = MoliereRatio*baseR; // slightly larger than sqrt(2)*2*3.82, but less than 3*3.82
      for(int j = 1152; j < 1728; j++) {
        if(Energy[j]>0.&&(((Position.x[j] - Position.x[ClusterCenter])*(Position.x[j] - Position.x[ClusterCenter]) +
        (Position.y[j] - Position.y[ClusterCenter])*(Position.y[j] - Position.y[ClusterCenter])) <= ClusterRadius*ClusterRadius)) {
          Cluster_BlockN[ClusterN] = j;
          ClusterN += 1;
          ClusterEnergy += Energy[j];
        }
      }
/*
      if(fabs(Position.x[ClusterCenter]) - MoliereRatio*baseR < 34.85) ClusterRadius = MoliereRatio*baseR + (1./MoliereRatio - 1.)*(34.85 + MoliereRatio*baseR - fabs(Position.x[ClusterCenter]));
      if(fabs(Position.y[ClusterCenter]) - MoliereRatio*baseR < 34.85) ClusterRadius = fmin(ClusterRadius, MoliereRatio*baseR + (1./MoliereRatio - 1.)*(34.85 + MoliereRatio*baseR - fabs(Position.y[ClusterCenter])));
*/
      for(int j = 0; j < 1152; j++) {
        if(Energy[j]>0.&&(((Position.x[j] - Position.x[ClusterCenter])*(Position.x[j] - Position.x[ClusterCenter]) +
        (Position.y[j] - Position.y[ClusterCenter])*(Position.y[j] - Position.y[ClusterCenter])) <= ClusterRadius*ClusterRadius)) {
          Cluster_BlockN[ClusterN] = j;
          ClusterN += 1;
          ClusterEnergy += Energy[j];
        }
      }

    }

//    cout << ClusterEnergy << endl;
     // Cluster Reconstruction
     Totalweight = 0.;
     weightx = 0.;
     weighty = 0.;
     double CorrectEnergy = 0.;
     double CorrectFactor = 1;

     if(ClusterEnergy >= 50.) {
       int LogMethod = 0;
       if(fabs(Position.x[ClusterCenter]) > (34.85 - 4.) && fabs(Position.x[ClusterCenter]) < (34.85 + 4.*MoliereRatio)) LogMethod = 1;
       if(fabs(Position.y[ClusterCenter]) > (34.85 - 4.) && fabs(Position.y[ClusterCenter]) < (34.85 + 4.*MoliereRatio)) LogMethod = 1;
       for(int cn = 0; cn < ClusterN; ++cn) {
         nid = Cluster_BlockN[cn];
         CorrectEnergy += Energy[nid]/CorrectFactor;
         weight = Energy[nid]/ClusterEnergy;
         if(LogMethod == 1) weight = 10. + log(Energy[nid]/ClusterEnergy); //logarithmic method for transition area
         if(weight > 0.) {         
           weightx += weight*Position.x[nid];
           weighty += weight*Position.y[nid];
           Totalweight += weight;
         }
       }
       reconx[Ci] = weightx/Totalweight;
       recony[Ci] = weighty/Totalweight;
       reconE[Ci] = (CorrectEnergy - 2.35909)/0.91439;
       ExistCluster[Ci] = ClusterCenter;
       Ci += 1;
     }
     else {reconE[Ci] = 0.; ExistCluster[Ci] = ClusterCenter; Ci += 1;}
   }

    //Place holder, for position correction from center of shower to HyCal surface
    //
    //
    //Place holder, for energy compensation at the edge
    //
    //
    //

    // Output
    double recontheta;
    int Noutput = 0;
    for(int k = 0; k < Ci; ++ k) {
      if(reconE[k] > 0.) {
        Noutput += 1;
        recontheta = atan(sqrt(reconx[k]*reconx[k] + recony[k]*recony[k])/541.)*180./3.14159265358979;
//        output << recontheta << "  " << tmp2 << endl;// "  " << tmp3 << "  " << tmp4 << endl;
        output << reconx[k] << "  " << recony[k] <<"  "<< reconE[k] << endl;
//      output << recontheta << "  " << reconE << endl;
      }
    }
    output << -9999 << "  " << Noutput << "  " << Ci << endl;
    NbEvents += 1;
    cout << "\r---N = " << NbEvents  << flush;

  }
  cout << endl;
  enefile.close();
  output.close();
  return 0;
}
