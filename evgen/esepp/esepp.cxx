//----------------------------------------------------------------------------------------------
// The ESEPP event generator (Elastic Scattering of Electrons and Positrons on Protons)
// Version 1.4 (April 27, 2014)
// Detailed description: A.V. Gramolin, et al., arXiv:1401.2959
//
// ESEPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ESEPP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ESEPP.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (c) Alexander Gramolin, 2014. E-mail: gramolin (at) inp.nsk.su 
// http://gramolin.com/esepp/
//==============================================================================================


// To support large output files:
#define _FILE_OFFSET_BITS 64

#include "const.h"  // Header file with some constants
#include "esepp.h"  // Main header file
#include "dialog.h" // Header file with the initial dialog

// Accurate QED calculation of first-order bremsstrahlung
// (V.S. Fadin, A.L. Feldman and R.E. Gerasimov):
#include "lepton.h"       // Lepton term
#include "interference.h" // Interference term
#include "proton.h"       // Proton term


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR0 (1D, Rosenbluth events):
class TFDISTR0: public TFoamIntegrand {
public:
  TFDISTR0(){};
  Double_t Density(int nDim, Double_t *arg) {

  theta_l = theta_min + arg[0]*(theta_max - theta_min); // Theta angle for the lepton

  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

  E_g = 0.; theta_g = 0.; phi_g = 0.; // Elastic scattering
  SetFinalFourMomenta(); // Set four-momenta for the final particles
  
  return (phi_max - phi_min)*(theta_max - theta_min)*inter_ros_sin.Eval(theta_l)*mkb;
}
};
// The end of TFDISTR0
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR1e (1D, elastic scattering, e-/mu-):
class TFDISTR1e: public TFoamIntegrand {
public:
  TFDISTR1e(){};
  Double_t Density(int nDim, Double_t *arg) {

  theta_l = theta_min + arg[0]*(theta_max - theta_min); // Theta angle for the lepton
  
  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

  E_g = 0.; theta_g = 0.; phi_g = 0.; // Elastic scattering
  SetFinalFourMomenta(); // Set four-momenta for the final particles
  
  // Mode of calculation for bremsstrahlung:
  switch (flag_mode)
    {
    case 1: case 4: case 7: case 10: // Only lepton bremsstrahlung
      {
      delta_sum = inter_brem_ee.Eval(theta_l) + inter_virt.Eval(theta_l) + inter_prime.Eval(theta_l);
      break;
      }
      
    case 2: case 5: case 8: case 11: // Only proton bremsstrahlung
      {
      delta_sum = inter_brem_pp.Eval(theta_l) + inter_virt.Eval(theta_l) + inter_prime.Eval(theta_l);
      break;
      }
      
    case 3: case 6: case 9: case 12: // All the bremsstrahlung terms
      {
      delta_sum = inter_brem_ee.Eval(theta_l) + inter_brem_pp.Eval(theta_l) + inter_brem_ep.Eval(theta_l) + inter_virt.Eval(theta_l) + inter_prime.Eval(theta_l);
      break;
      }
    }
  
  return (phi_max - phi_min)*(theta_max - theta_min)*(1. + delta_sum)*inter_ros_sin.Eval(theta_l)*mkb;
}
};
// The end of TFDISTR1e
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR1p (1D, elastic scatterring, e+/mu+):
class TFDISTR1p: public TFoamIntegrand {
public:
  TFDISTR1p(){};
  Double_t Density(int nDim, Double_t *arg) {

  theta_l = theta_min + arg[0]*(theta_max - theta_min); // Theta angle for the lepton

  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

  E_g = 0.; theta_g = 0.; phi_g = 0.; // Elastic scattering
  SetFinalFourMomenta(); // Set four-momenta for the final particles

  // Mode of calculation for bremsstrahlung:
  switch (flag_mode)
    {
    case 1: case 4: case 7: case 10: // Only lepton bremsstrahlung
      {
      delta_sum = inter_brem_ee.Eval(theta_l) + inter_virt.Eval(theta_l) - inter_prime.Eval(theta_l);
      break;
      }
      
    case 2: case 5: case 8: case 11: // Only proton bremsstrahlung
      {
      delta_sum = inter_brem_pp.Eval(theta_l) + inter_virt.Eval(theta_l) - inter_prime.Eval(theta_l);
      break;
      }
      
    case 3: case 6: case 9: case 12: // All the bremsstrahlung terms
      {
      delta_sum = inter_brem_ee.Eval(theta_l) + inter_brem_pp.Eval(theta_l) - inter_brem_ep.Eval(theta_l) + inter_virt.Eval(theta_l) - inter_prime.Eval(theta_l);
      break;
      }
    }
  
  return (phi_max - phi_min)*(theta_max - theta_min)*(1. + delta_sum)*inter_ros_sin.Eval(theta_l)*mkb;
}
};
// The end of TFDISTR1p
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR2e (4D, bremsstrahlung, e-/mu-):
class TFDISTR2e: public TFoamIntegrand {
public:
  TFDISTR2e(){};
  Double_t Density(int nDim, Double_t *arg) {

  // Four arguments (basic kinematic variables):
  theta_g = arg[0]*Pi; // Theta angle for the photon
  phi_g = arg[1]*2.*Pi; // Phi angle for the photon
  E_g = E_g_cut + arg[2]*(E_g_max - E_g_cut); // Energy for the photon
  theta_l = theta_min + arg[3]*(theta_max - theta_min); // Theta angle for the lepton
  
  if (flag_mode > 3)
    {
    // Checking the photon energy:
    if (E_g > M*(E_li - m)/(M + E_li - Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_g))) return 0.;
      
    // Evaluation of the final lepton energy:
    if (EvalEnergy() == 0) return 0.;

    EvalKinematicParams(E_g); // Evaluation of some kinematic parameters
    SetFinalFourMomenta(); // Set four-momenta for the final particles
    }
  else // Primary soft-photon approximation
    {
    E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton
      
    // Four-momenta of the final lepton and proton (elastic kinematics):
    v_lf.SetPxPyPzE(Sqrt(Pow2(E_lf) - Pow2(m))*Sin(theta_l), 0., Sqrt(Pow2(E_lf) - Pow2(m))*Cos(theta_l), E_lf);
    v_pf = v_li + v_pi - v_lf;
    E_p = v_pf.E(); theta_p = v_pf.Theta(); phi_p = v_pf.Phi();
      
    // Four-momentum of the photon:
    v_kf.SetPxPyPzE(E_g*Sin(theta_g)*Cos(phi_g), E_g*Sin(theta_g)*Sin(phi_g), E_g*Cos(theta_g), E_g);
    }

  // Mode of calculation for bremsstrahlung:
  switch (flag_mode)
    {
    case 1: // Primary soft-photon approximation, only lepton bremsstrahlung
      {
      ret =  2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
      // Back to the elastic kinematics:
      E_g = 0.; theta_g = 0.; phi_g = 0.;
      v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
      return ret;
      }

    case 2: // Primary soft-photon approximation, only proton bremsstrahlung
      {
      ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
      // Back to the elastic kinematics:
      E_g = 0.; theta_g = 0.; phi_g = 0.;
      v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
      return ret;
      }
      
    case 3: // Primary soft-photon approximation, all the terms
      {
      ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() + BremSoftInterference())*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
      // Back to the elastic kinematics:
      E_g = 0.; theta_g = 0.; phi_g = 0.;
      v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
      return ret;
      }

    case 4: // Modified soft-photon approximation, only lepton bremsstrahlung
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      }

    case 5: // Modified soft-photon approximation, only proton bremsstrahlung
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      }

    case 6: // Modified soft-photon approximation, all the terms
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() + BremSoftInterference())*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      }

    case 7: // Improved soft-photon approximation, only lepton bremsstrahlung
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
      }

    case 8: // Improved soft-photon approximation, only proton bremsstrahlung
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
      }

    case 9: // Improved soft-photon approximation, all the terms
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() + BremSoftInterference())*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
      }

    case 10: // Accurate QED calculation, only lepton bremsstrahlung
      {
      EvalAllProducts(); // Evaluation of the four-momentum products
	
      // Square of the total amplitude:
      M_sum = lterm();

      // Checking the the square of the amplitude:
      if (M_sum < 0.)
        {
	cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
        flag_warn = true;
        return 0.;
	}

      return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m))*M))*E_g*((Pow2(E_lf) - Pow2(m))/Abs((Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
      }

    case 11: // Accurate QED calculation, only proton bremsstrahlung
      {
      EvalAllProducts(); // Evaluation of the four-momentum products
	
      // Square of the total amplitude:
      M_sum = pterm();

      // Checking the the square of the amplitude:
      if (M_sum < 0.)
        {
        cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
        flag_warn = true;
        return 0.;
	}

      return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m))*M))*E_g*((Pow2(E_lf) - Pow2(m))/Abs((Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
      }
      
    case 12: // Accurate QED calculation, all the terms
      {
      EvalAllProducts(); // Evaluation of the four-momentum products
	
      // Square of the total amplitude:
      M_sum = lterm() + pterm() + lpterm();

      // Checking the the square of the amplitude:
      if (M_sum < 0.)
        {
	cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
        flag_warn = true;
        return 0.;
	}
      
      return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m))*M))*E_g*((Pow2(E_lf) - Pow2(m))/Abs((Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
      }
    }
  
  return 0.;
}
};
// The end of TFDISTR2e
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR2p (4D, bremsstrahlung, e+/mu+):
class TFDISTR2p: public TFoamIntegrand {
public:
  TFDISTR2p(){};
  Double_t Density(int nDim, Double_t *arg) {

  // Four arguments (basic kinematic variables):
  theta_g = arg[0]*Pi; // Theta angle for the photon
  phi_g = arg[1]*2.*Pi; // Phi angle for the photon
  E_g = E_g_cut + arg[2]*(E_g_max - E_g_cut); // Energy for the photon
  theta_l = theta_min + arg[3]*(theta_max - theta_min); // Theta angle for the lepton
  
  if (flag_mode > 3)
    {
    // Checking the photon energy:
    if (E_g > M*(E_li - m)/(M + E_li - Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_g))) return 0.;

    // Evaluation of the final lepton energy:
    if (EvalEnergy() == 0) return 0.;

    EvalKinematicParams(E_g); // Evaluation of some kinematic parameters
    SetFinalFourMomenta(); // Set four-momenta for the final particles
    }
  else // Primary soft-photon approximation
    {
    E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton
      
    // Four-momenta of the final lepton and proton (elastic kinematics):
    v_lf.SetPxPyPzE(Sqrt(Pow2(E_lf) - Pow2(m))*Sin(theta_l), 0., Sqrt(Pow2(E_lf) - Pow2(m))*Cos(theta_l), E_lf);
    v_pf = v_li + v_pi - v_lf;
    E_p = v_pf.E(); theta_p = v_pf.Theta(); phi_p = v_pf.Phi();
      
    // Four-momentum of the photon:
    v_kf.SetPxPyPzE(E_g*Sin(theta_g)*Cos(phi_g), E_g*Sin(theta_g)*Sin(phi_g), E_g*Cos(theta_g), E_g);
    }

  // Mode of calculation for bremsstrahlung:
  switch (flag_mode)
    {
    case 1: // Primary soft-photon approximation, only lepton bremsstrahlung
      {
      ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
      // Back to the elastic kinematics:
      E_g = 0.; theta_g = 0.; phi_g = 0.;
      v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
      return ret;
      }

    case 2: // Primary soft-photon approximation, only proton bremsstrahlung
      {
      ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
      // Back to the elastic kinematics:
      E_g = 0.; theta_g = 0.; phi_g = 0.;
      v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
      return ret;
      }
      
    case 3: // Primary soft-photon approximation, all the terms
      {
      ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() - BremSoftInterference())*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
      // Back to the elastic kinematics:
      E_g = 0.; theta_g = 0.; phi_g = 0.;
      v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
      return ret;
      }

    case 4: // Modified soft-photon approximation, only lepton bremsstrahlung
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      }
      
    case 5: // Modified soft-photon approximation, only proton bremsstrahlung
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      }
      
    case 6: // Modified soft-photon approximation, all the terms
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() - BremSoftInterference())*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      }
      
    case 7: // Improved soft-photon approximation, only lepton bremsstrahlung
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
      }
      
    case 8: // Improved soft-photon approximation, only proton bremsstrahlung
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
      }
      
    case 9: // Improved soft-photon approximation, all the terms
      {
      return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() - BremSoftInterference())*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
      }
      
    case 10: // Accurate QED calculation, only lepton bremsstrahlung
      {
      EvalAllProducts(); // Evaluation of the four-momentum products
	
      // Square of the total amplitude:
      M_sum = lterm();

      // Checking the the square of the amplitude:
      if (M_sum < 0.)
        {
	cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
        flag_warn = true;
        return 0.;
	}

      return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m))*M))*E_g*((Pow2(E_lf) - Pow2(m))/Abs((Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
      }

    case 11: // Accurate QED calculation, only proton bremsstrahlung
      {
      EvalAllProducts(); // Evaluation of the four-momentum products
	
      // Square of the total amplitude:
      M_sum = pterm();

      // Checking the the square of the amplitude:
      if (M_sum < 0.)
        {
	cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
        flag_warn = true;
        return 0.;
	}

      return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m))*M))*E_g*((Pow2(E_lf) - Pow2(m))/Abs((Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
      }
      
    case 12: // Accurate QED calculation, all the terms
      {
      EvalAllProducts(); // Evaluation of the four-momentum products
	
      // Square of the total amplitude:
      M_sum = lterm() + pterm() - lpterm();

      // Checking the the square of the amplitude:
      if (M_sum < 0.)
        {
	cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
        flag_warn = true;
        return 0.;
	}

      return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m))*M))*E_g*((Pow2(E_lf) - Pow2(m))/Abs((Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
      }
    }
  
  return 0.;
}
};
// The end of TFDISTR2p
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Main function:
int main(int argc, char **argv)
{

// Initial values for the flags:
flag_quick = false;   // Flag to speed up the calculations
flag_init = false;    // Flag to set the number of events for initialization
flag_info = false;    // Flag to print additional information
flag_vepp = false;    // Flag to generate events for the VEPP-3 experiment
flag_target = false;  // Flag to work with a storage cell
flag_warn = false;    // Flag indicating the presence/absence of warnings


//----------------------------------------------------------------------------------------------
// Processing of the command line options:
if (argc > 1)
  {
  TString str;

  for (int j = 1; j < argc; j++)
    {
    str = argv[j];
  
    if (str == "-quick") { flag_quick = true; if (flag_init == false) { flag_init = true; ninit = 1000000; } }
    if (str == "-init") { flag_init = true; ninit = atoi(argv[j+1]); }
    if (str == "-info") { flag_info = true; }
    if (str == "-vepp" || str == "-VEPP") { flag_vepp = true; flag_target = true; }
    if (str == "-target") { flag_target = true; }
    
    if (str == "-c" || str == "-C") // Key "-c"
      {
      cout << endl << "This program is free software: you can redistribute it and/or modify" << endl;
      cout << "it under the terms of the GNU General Public License as published by" << endl;
      cout << "the Free Software Foundation, either version 3 of the License, or" << endl;
      cout << "(at your option) any later version." << endl;

      cout << endl << "You should have received a copy of the GNU General Public License" << endl;
      cout << "along with this program. If not, see http://www.gnu.org/licenses/." << endl;
      }
    
    if (str == "-w" || str == "-W") // Key "-w"
      {
      cout << endl << "This program is distributed in the hope that it will be useful," << endl;
      cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
      cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the" << endl;
      cout << "GNU General Public License for more details." << endl;

      cout << endl << "You should have received a copy of the GNU General Public License" << endl;
      cout << "along with this program. If not, see http://www.gnu.org/licenses/." << endl;
      }
      
    if (str == "-h" || str == "-help" || str == "--help") // Help (keys "-h", "-help" or "--help")
      {
      cout << endl << "You can run esepp with the following options:" << endl << endl;
      cout << "-help       to print this help;" << endl;
      cout << "-info       to print additional information;" << endl;
      cout << "-init N     to set the number N of events for initialization;" << endl;
      cout << "-quick      to speed up the calculations (only for debugging!);" << endl;
      cout << "-target     to work with a storage cell;" << endl;
      cout << "-vepp       to generate events for the VEPP-3 experiment." << endl << endl;
      return EXIT_SUCCESS;
      }
    }
  }
//==============================================================================================


if (Dialog() == EXIT_FAILURE) return EXIT_FAILURE; // Dialog for entering the initial data

time(&starttime); // Starting the timer

// Set the initial lepton four-momentum (in the laboratory frame):
v_li.SetPxPyPzE(0., 0., Sqrt(Pow2(E_li) - Pow2(m)), E_li);

// Set the initial proton four-momentum (in the laboratory frame):
v_pi.SetPxPyPzE(0., 0., 0., M);

omega = (phi_max - phi_min)*(Cos(theta_min) - Cos(theta_max)); // Solid angle (steradian)

// Initialization of the random number generator:
PseRan->SetSeed(0);


//----------------------------------------------------------------------------------------------
// Setting the functions and tolerance for numerical integration:
i_li_lf.SetFunction(func_li_lf); // To calculate B(v_li, v_lf, E_g_cut)
i_li_lf.SetRelTolerance(IntTol);

i_li_pi.SetFunction(func_li_pi); // To calculate B(v_li, v_pi, E_g_cut)
i_li_pi.SetRelTolerance(IntTol);

i_li_pf.SetFunction(func_li_pf); // To calculate B(v_li, v_pf, E_g_cut)
i_li_pf.SetRelTolerance(IntTol);

i_lf_pi.SetFunction(func_lf_pi); // To calculate B(v_lf, v_pi, E_g_cut)
i_lf_pi.SetRelTolerance(IntTol);

i_lf_pf.SetFunction(func_lf_pf); // To calculate B(v_lf, v_pf, E_g_cut)
i_lf_pf.SetRelTolerance(IntTol);

i_pi_pf.SetFunction(func_pi_pf); // To calculate B(v_pi, v_pf, E_g_cut)
i_pi_pf.SetRelTolerance(IntTol);

// To integrate the "elastic" part of cross section (with E_g < E_g_cut):
i_el_e.SetFunction(func_el_e); // e-/mu-
i_el_e.SetRelTolerance(IntTol);

i_el_p.SetFunction(func_el_p); // e+/mu+
i_el_p.SetRelTolerance(IntTol);

// To integrate the Rosenbluth differential cross section:
i_ros.SetFunction(func_ros);
i_ros.SetRelTolerance(IntTol);
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Interpolation of bremsstrahlung and virtual-photon corrections:
if (flag_vpol == 3) // In the case of accurate calculation of the vacuum polarization
  {
  int npoints = 0;
  char str[128];

  while (!feof(fvpol)) // Reading from the file "vpol.dat"
    {
    str[0] = 0;
    fgets(str, 128, fvpol);
    if (feof(fvpol) || strlen(str) == 0) break; // The end or empty string
    
    if (str[0] != '/')
      {
      sscanf(str, "%lf %lf", &s[npoints], &rep[npoints]);
      rep[npoints] = 2.*rep[npoints];
      npoints++;
      }
    }
    
  fclose(fvpol); // Closing the file "vpol.dat"
  
  inter_vpol.SetData(npoints, s, rep); // Interpolating
  }

E_g = 0.; theta_g = 0.; phi_g = 0.; // Elastic scattering
for (i = 0; i < InterpolPoints; i++)
  {
  // Theta angle for the lepton:
  theta_l = theta_min + i*(theta_max - theta_min)/(InterpolPoints - 1);

  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

  EvalKinematicParams(0.); // Evaluation of some kinematic parameters
  SetFinalFourMomenta();   // Set four-momenta for the final particles

  xx[i] = theta_l;               // Theta angle for the lepton
  y_ros_sin[i] = f_ros(theta_l); // RosenbluthCS()*Sin(theta_l)
  y_brem_ee[i] = d_brem_ee();    // Bremsstrahlung correction from the lepton term
  y_brem_pp[i] = d_brem_pp();    // Bremsstrahlung correction from the proton term
  y_brem_ep[i] = d_brem_ep();    // Bremsstrahlung correction from the interference term
  y_prime[i] = d_prime();        // TPE contribution by Maximon & Tjon
  
  // Virtual-photon correction:
  if (flag_vpol == 1) y_virt[i] = d_vac_e() + d_vertex();
  if (flag_vpol == 2) y_virt[i] = d_vac_e() + d_vac_mu() + d_vac_tau() + d_vertex();
  if (flag_vpol == 3) y_virt[i] = inter_vpol.Eval(-qq) + d_vertex();
  }

// Interpolating functions:
inter_ros_sin.SetData(InterpolPoints, xx, y_ros_sin); // RosenbluthCS()*Sin(theta_l)
inter_brem_ee.SetData(InterpolPoints, xx, y_brem_ee); // Lepton bremsstrahlung term
inter_brem_pp.SetData(InterpolPoints, xx, y_brem_pp); // Proton bremsstrahlung term
inter_brem_ep.SetData(InterpolPoints, xx, y_brem_ep); // Interference bremsstrahlung term
inter_virt.SetData(InterpolPoints, xx, y_virt);       // Virtual-photon correction
inter_prime.SetData(InterpolPoints, xx, y_prime);     // TPE contribution by Maximon & Tjon
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Some additional information:
if (flag_info == true)
  {
  theta_l = (theta_min + theta_max)/2.; // The average polar angle
  
  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

  EvalKinematicParams(0.); // Evaluation of some kinematic parameters
  SetFinalFourMomenta();   // Set four-momenta for the final particles

  cout << endl << "theta_l = " << theta_l/degrad << " degree, Q^2 = " << -qq << " GeV^2, epsilon = " << eps << endl;
  cout << "Mott differential cross section: " << Pow2(alpha*Cos(theta_l/2.)/(2.*E_li*Pow2(Sin(theta_l/2.))))*(E_lf/E_li)*mkb << " microbarn / steradian" << endl;
  cout << "d_virt = " << inter_virt.Eval(theta_l) << ", d_brem = " << d_brem_ee() + d_brem_ep() + d_brem_pp() << ", d_prime = " << d_prime() << endl;
  cout << "d_brem: d_brem_ee = " << d_brem_ee() << ", d_brem_ep = " << d_brem_ep() << ", d_brem_pp = " << d_brem_pp() << endl;
  cout << "d_virt: d_vertex = " << d_vertex() << ", d_vac_e = " << d_vac_e() << ", d_vac_mu = " << d_vac_mu() << ", d_vac_tau = " << d_vac_tau();
  if (flag_vpol == 3) cout << ", d_vac = " << inter_vpol.Eval(-qq);
  cout << endl;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// mFOAM simulators (S. Jadach, P. Sawicki, arXiv:physics/0506084):

// Some options for the simulators:
if (flag_quick == false) // Default options
  {
  nCells_1D = 1000;  // Number of cells for 1D case
  nSampl_1D = 200;   // Number of samples for 1Dcase
  nCells_4D = 30000; // Number of cells for 4D case
  nSampl_4D = 1500;  // Number of samples for 4D case
  }
else // Options for the case of quick calculation
  {
  nCells_1D = 500;   // Number of cells for 1D case
  nSampl_1D = 100;   // Number of samples for 1Dcase
  nCells_4D = 10000; // Number of cells for 4D case
  nSampl_4D = 600;   // Number of samples for 4D case
  }
  
// Creating Foam simulators:
TFoam *Foam_Ros = new TFoam("Foam_Ros"); // Rosenbluth events
TFoam *Foam_ElE = new TFoam("Foam_ElE"); // Elastic scattering, e-/mu-
TFoam *Foam_BrE1 = new TFoam("Foam_BrE1"); // First-order bremsstrahlung, e-/mu-, en_sign = -1
TFoam *Foam_BrE2 = new TFoam("Foam_BrE2"); // First-order bremsstrahlung, e-/mu-, en_sign = +1
TFoam *Foam_ElP = new TFoam("Foam_ElP"); // Elastic scattering, e+/mu+
TFoam *Foam_BrP1 = new TFoam("Foam_BrP1"); // First-order bremsstrahlung, e+/mu+, en_sign = -1
TFoam *Foam_BrP2 = new TFoam("Foam_BrP2"); // First-order bremsstrahlung, e+/mu+, en_sign = +1

if (flag_rosen == true) // Rosenbluth events
  {
  // Distribution function:
  TFoamIntegrand *Rho_Ros = new TFDISTR0();

  // Initialization of the Foam_Ros simulator (Rosenbluth events):
  cout << endl << "Initialization (Rosenbluth scattering):" << endl;
  Foam_Ros->SetkDim(1);           // Set number of dimensions
  Foam_Ros->SetnCells(nCells_1D); // Set number of cells
  Foam_Ros->SetnSampl(nSampl_1D); // Set number os samples
  Foam_Ros->SetOptRej(1);         // Unweighted events in MC generation
  Foam_Ros->SetRho(Rho_Ros);      // Set distribution function
  Foam_Ros->SetPseRan(PseRan);    // Set random number generator
  Foam_Ros->SetChat(1);           // Set "chat level" in the standard output
  Foam_Ros->Initialize();         // Initialization
  }

if (flag_lepton != 2 && flag_lepton != 5)  // For negatively charged leptons
  {
  // Distribution functions:
  TFoamIntegrand *Rho_ElE = new TFDISTR1e(); // Elastic scattering, e-/mu-
  TFoamIntegrand *Rho_BrE = new TFDISTR2e(); // First-order bremsstrahlung, e-/mu-
  
  cout << endl;
  
  // Initialization of the Foam_ElE1 simulator (elastic scattering, e-/mu-):
  if (flag_lepton == 1 || flag_lepton == 3) cout << "Initialization (elastic scattering, e-):" << endl;
  if (flag_lepton == 4 || flag_lepton == 6) cout << "Initialization (elastic scattering, mu-):" << endl;
  Foam_ElE->SetkDim(1);           // Set number of dimensions
  Foam_ElE->SetnCells(nCells_1D); // Set number of cells
  Foam_ElE->SetnSampl(nSampl_1D); // Set number os samples
  Foam_ElE->SetOptRej(1);         // Unweighted events in MC generation
  Foam_ElE->SetRho(Rho_ElE);      // Set distribution function
  Foam_ElE->SetPseRan(PseRan);    // Set random number generator
  Foam_ElE->SetChat(1);           // Set "chat level" in the standard output
  Foam_ElE->Initialize();         // Initialization
  
  cout << endl;
  
  en_sign = -1.; // The "-" root of the equation for the scattered lepton energy
  
  // Initialization of the Foam_BrE1 simulator (bremsstrahlung, e-/mu-, root "-"):
  if (flag_lepton == 1 || flag_lepton == 3) cout << "Initialization (bremsstrahlung, e-, root \"-\"):" << endl;
  if (flag_lepton == 4 || flag_lepton == 6) cout << "Initialization (bremsstrahlung, mu-, root \"-\"):" << endl;
  Foam_BrE1->SetkDim(4);           // Set number of dimensions
  Foam_BrE1->SetnCells(nCells_4D); // Set number of cells
  Foam_BrE1->SetnSampl(nSampl_4D); // Set number os samples
  Foam_BrE1->SetOptRej(1);         // Unweighted events in MC generation
  Foam_BrE1->SetRho(Rho_BrE);      // Set distribution function
  Foam_BrE1->SetPseRan(PseRan);    // Set random number generator
  Foam_BrE1->SetChat(1);           // Set "chat level" in the standard output
  Foam_BrE1->Initialize();         // Initialization
  
  if (flag_mode > 3)
    {
    cout << endl;
      
    en_sign = +1.; // The "+" root of the equation for the scattered lepton energy  

    // Initialization of the Foam_BrE2 simulator (bremsstrahlung, e-/mu-, root "+"):
    if (flag_lepton == 1 || flag_lepton == 3) cout << "Initialization (bremsstrahlung, e-, root \"+\"):" << endl;
    if (flag_lepton == 4 || flag_lepton == 6) cout << "Initialization (bremsstrahlung, mu-, root \"+\"):" << endl;
    Foam_BrE2->SetkDim(4);           // Set number of dimensions
    Foam_BrE2->SetnCells(nCells_4D); // Set number of cells
    Foam_BrE2->SetnSampl(nSampl_4D); // Set number os samples
    Foam_BrE2->SetOptRej(1);         // Unweighted events in MC generation
    Foam_BrE2->SetRho(Rho_BrE);      // Set distribution function
    Foam_BrE2->SetPseRan(PseRan);    // Set random number generator
    Foam_BrE2->SetChat(1);           // Set "chat level" in the standard output
    Foam_BrE2->Initialize();         // Initialization
  
    bre2 = Foam_BrE2->GetPrimary();
    if (bre2 == 0) cout << "Foam_BrE2 (bremsstrahlung, e-/mu-, root \"+\"): no need to use." << endl;
    }
  else bre2 = 0;
  }
  
if (flag_lepton != 1 && flag_lepton != 4)  // For positively charged leptons
  {
  // Distribution functions:
  TFoamIntegrand *Rho_ElP = new TFDISTR1p(); // Elastic scattering, e+/mu+
  TFoamIntegrand *Rho_BrP = new TFDISTR2p(); // First-order bremsstrahlung, e+/mu+

  cout << endl;
  
  // Initialization of the Foam_ElP simulator (elastic scattering, e+/mu+):
  if (flag_lepton == 2 || flag_lepton == 3) cout << "Initialization (elastic scattering, e+):" << endl;
  if (flag_lepton == 5 || flag_lepton == 6) cout << "Initialization (elastic scattering, mu+):" << endl;
  Foam_ElP->SetkDim(1);           // Set number of dimensions
  Foam_ElP->SetnCells(nCells_1D); // Set number of cells
  Foam_ElP->SetnSampl(nSampl_1D); // Set number os samples
  Foam_ElP->SetOptRej(1);         // Unweighted events in MC generation
  Foam_ElP->SetRho(Rho_ElP);      // Set distribution function
  Foam_ElP->SetPseRan(PseRan);    // Set random number generator
  Foam_ElP->SetChat(1);           // Set "chat level" in the standard output
  Foam_ElP->Initialize();         // Initialization

  cout << endl;
  
  en_sign = -1.; // The "-" root of the equation for the scattered lepton energy
  
  // Initialization of the Foam_BrP simulator (bremsstrahlung, e+/mu+, root "-"):
  if (flag_lepton == 2 || flag_lepton == 3) cout << "Initialization (bremsstrahlung, e+, root \"-\"):" << endl;
  if (flag_lepton == 5 || flag_lepton == 6) cout << "Initialization (bremsstrahlung, mu+, root \"-\"):" << endl;
  Foam_BrP1->SetkDim(4);           // Set number of dimensions
  Foam_BrP1->SetnCells(nCells_4D); // Set number of cells
  Foam_BrP1->SetnSampl(nSampl_4D); // Set number os samples
  Foam_BrP1->SetOptRej(1);         // Unweighted events in MC generation
  Foam_BrP1->SetRho(Rho_BrP);      // Set distribution function
  Foam_BrP1->SetPseRan(PseRan);    // Set random number generator
  Foam_BrP1->SetChat(1);           // Set "chat level" in the standard output
  Foam_BrP1->Initialize();         // Initialization
  
  if (flag_mode > 3)
    {
    cout << endl;

    en_sign = +1.; // The "+" root of the equation for the scattered lepton energy
  
    // Initialization of the Foam_BrP simulator (bremsstrahlung, e+/mu+, root "+"):
    if (flag_lepton == 2 || flag_lepton == 3) cout << "Initialization (bremsstrahlung, e+, root \"+\"):" << endl;
    if (flag_lepton == 5 || flag_lepton == 6) cout << "Initialization (bremsstrahlung, mu+, root \"+\"):" << endl;
    Foam_BrP2->SetkDim(4);           // Set number of dimensions
    Foam_BrP2->SetnCells(nCells_4D); // Set number of cells
    Foam_BrP2->SetnSampl(nSampl_4D); // Set number os samples
    Foam_BrP2->SetOptRej(1);         // Unweighted events in MC generation
    Foam_BrP2->SetRho(Rho_BrP);      // Set distribution function
    Foam_BrP2->SetPseRan(PseRan);    // Set random number generator
    Foam_BrP2->SetChat(1);           // Set "chat level" in the standard output
    Foam_BrP2->Initialize();         // Initialization
  
    brp2 = Foam_BrP2->GetPrimary();
    if (brp2 == 0) cout << "Foam_BrP2 (bremsstrahlung, e+/mu+, root \"+\"): no need to use." << endl;
    }
  else brp2 = 0;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Initialization of the generator (calculation of the bremsstrahlung cross sections):

// The number of events for initialization (default: 10 million):
if (flag_init == false)
  {
  ninit = long (1e7);
  if (nevents >= 4e7) ninit = long (2e7);
  if (nevents >= 8e7) ninit = long (4e7);
  }

i = 0;
cout << setprecision(3) << endl;
en_sign = -1.; // The "-" root of the equation for the scattered lepton energy

// Negatively charged leptons:------------------------------------------------------------------
if (flag_lepton == 1 || flag_lepton == 4)
  {
  for (loop = 0; loop < ninit; loop++)
    {
    Foam_BrE1->MakeEvent(); // Generating a bremsstrahlung event for e-/mu-, root "-"
    
    // The progress bar:
    if (loop == long (0.01*i*(ninit - 1)))
      {
      cout << "Initialization 1 (" << ninit << " events): ";
      cout << i << "%" << '\r' << flush;
      i += 1;
      }
    }
    
  cout << endl;  
  Foam_BrE1->GetIntegMC(bre1, bre_error); // Cross section for e-/mu-, root "-"
  }

// Positively charged leptons:------------------------------------------------------------------
if (flag_lepton == 2 || flag_lepton == 5)
  {
  for (loop = 0; loop < ninit; loop++)
    {
    Foam_BrP1->MakeEvent(); // Generating a bremsstrahlung event for e+/mu+, root "-"
    
    // The progress bar:
    if (loop == long (0.01*i*(ninit - 1)))
      {
      cout << "Initialization 1 (" << ninit << " events): ";
      cout << i << "%" << '\r' << flush;
      i += 1;
      }
    }
    
  cout << endl;
  Foam_BrP1->GetIntegMC(brp1, brp_error); // Cross section for e+/mu+, root "-"
  }

// Both negatively and positively charged leptons:----------------------------------------------
if (flag_lepton == 3 || flag_lepton == 6)
  {
  for (loop = 0; loop < ninit; loop++)
    {
    Foam_BrE1->MakeEvent(); // Generating a bremsstrahlung event for e-/mu-, root "-"
    Foam_BrP1->MakeEvent(); // Generating a bremsstrahlung event for e+/mu+, root "-"
    
    // The progress status:
    if (loop == long (0.01*i*(ninit - 1)))
      {
      cout << "Initialization 1 (" << ninit << " events): ";
      cout << i << "%" << '\r' << flush;
      i += 1;
      }
    }

  cout << endl;
  Foam_BrE1->GetIntegMC(bre1, bre_error); // Cross section for e-/mu-, root "-"
  Foam_BrP1->GetIntegMC(brp1, brp_error); // Cross section for e+/mu+, root "-"
  }

// Initialization 2 (for "+" roots):
if (bre2 == 0 && brp2 == 0) cout << "Initialization 2: no need to perform";
else
  {
  i = 0;
  en_sign = +1.; // The "+" root of the equation for the scattered lepton energy

  for (loop = 0; loop < long (0.05*ninit); loop++)
    {
    // Generating a bremsstrahlung event for e-/mu-, root "+":
    if (bre2 > 0) Foam_BrE2->MakeEvent();

    // Generating a bremsstrahlung event for e+/mu+, root "+":
    if (brp2 > 0) Foam_BrP2->MakeEvent();

    // The progress status:
    if (loop == long (0.01*i*(long (0.05*ninit) - 1)))
      {
      cout << "Initialization 2 (" << long (0.05*ninit) << " events): ";
      cout << i << "%" << '\r' << flush;
      i += 1;
      }
    }
  }

cout << endl;
if (bre2 > 0) Foam_BrE2->GetIntegMC(bre2, bre_error); // Cross section for e-/mu-, root "+"
if (brp2 > 0) Foam_BrP2->GetIntegMC(brp2, brp_error); // Cross section for e+/mu+, root "+"
//==============================================================================================


// Elastic scattering cross sections:
E_g = 0.; theta_g = 0.; phi_g = 0.;
double elast_e = (phi_max - phi_min)*mkb*i_el_e.Integral(theta_min, theta_max); // e-/mu-
double elast_p = (phi_max - phi_min)*mkb*i_el_p.Integral(theta_min, theta_max); // e+/mu+

// Rosenbluth cross section:
double rosen = (phi_max - phi_min)*mkb*i_ros.Integral(theta_min, theta_max);


//----------------------------------------------------------------------------------------------
// The numbers of events to generate:
long nev_ElE = 0; // Elastic scattering, e-/mu-
long nev_ElP = 0; // Elastic scattering, e+/mu+

long nev_BrE1 = 0; // Bremsstrahlung, e-/mu-, root "-"
long nev_BrE2 = 0; // Bremsstrahlung, e-/mu-, root "+"

long nev_BrP1 = 0; // Bremsstrahlung, e+/mu+, root "-"
long nev_BrP2 = 0; // Bremsstrahlung, e+/mu+, root "+"

long nev_0 = 0; // Number of Rosenbluth events
long nev_e = 0; // Total number of events for e-/mu-
long nev_p = 0; // Total number of events for e+/mu+

if (flag_lepton == 1 || flag_lepton == 4) // Only e-/mu-
  {
  sum_e = elast_e + bre1 + bre2;
  nev_0 = long (nevents*rosen/sum_e); // Rosenbluth scattering

  nev_ElE = long (nevents*elast_e/sum_e); // Elastic scattering

  nev_BrE2 = long (nevents*bre2/sum_e); // Bremsstrahlung, root "+"
  nev_BrE1 = nevents - (nev_ElE + nev_BrE2); // Bremsstrahlung, root "-"
  }
  
if (flag_lepton == 2 || flag_lepton == 5) // Only e+/mu+
  {
  sum_p = elast_p + brp1 + brp2;
  nev_0 = long (nevents*rosen/sum_p); // Rosenbluth scattering

  nev_ElP = long (nevents*elast_p/sum_p); // Elastic scattering

  nev_BrP2 = long (nevents*brp2/sum_p); // Bremsstrahlung, root "+"
  nev_BrP1 = nevents - (nev_ElP + nev_BrP2); // Bremsstrahlung, root "-"
  }

if ((flag_lepton == 3 || flag_lepton == 6) && flag_bint == false) // Both e-/mu- and e+/mu+
  {
  sum_e = elast_e + bre1 + bre2;
  sum_p = elast_p + brp1 + brp2;

  nev_0 = long (nevents*rosen/(sum_e + sum_p)); // Rosenbluth scattering

  nev_ElE = long (nevents*elast_e/(sum_e + sum_p)); // Elastic scattering, e-/mu-
  nev_ElP = long (nevents*elast_p/(sum_e + sum_p)); // Elastic scattering, e+/mu+

  nev_BrE1 = long (nevents*bre1/(sum_e + sum_p)); // Bremsstrahlung, e-/mu-, root "-"
  nev_BrE2 = long (nevents*bre2/(sum_e + sum_p)); // Bremsstrahlung, e-/mu-, root "+"
  
  nev_BrP2 = long (nevents*brp2/(sum_e + sum_p)); // Bremsstrahlung, e+/mu+, root "+"
  nev_BrP1 = nevents - (nev_ElE + nev_ElP + nev_BrE1 + nev_BrE2 + nev_BrP2); // Bremsstrahlung, e+/mu+, root "-"
  }
  
if (flag_bint == true) // If the beam current integral is specified (for VEPP-3)
  {
  nev_0 = long (kiloc*rosen/1.602e-7); // Rosenbluth scattering
    
  nev_ElE = long (kiloc*elast_e/1.602e-7); // Elastic scattering, e-
  nev_ElP = long (kiloc*elast_p/1.602e-7); // Elastic scattering, e+

  nev_BrE1 = long (kiloc*bre1/1.602e-7); // Bremsstrahlung, e-, root "-"
  nev_BrE2 = long (kiloc*bre2/1.602e-7); // Bremsstrahlung, e-, root "+"
  nev_BrP1 = long (kiloc*brp1/1.602e-7); // Bremsstrahlung, e+, root "-"
  nev_BrP2 = long (kiloc*brp2/1.602e-7); // Bremsstrahlung, e+, root "+"
  
  nevents = nev_ElE + nev_ElP + nev_BrE1 + nev_BrE2 + nev_BrP1 + nev_BrP2;
  }
  
if (flag_rosen == false) nev_0 = 0;
nev_e = nev_ElE + nev_BrE1 + nev_BrE2; // Total number of events for e-/mu-
nev_p = nev_ElP + nev_BrP1 + nev_BrP2; // Total number of events for e+/mu+
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Some additional information (to compare analytical and numerical integrations):
if (flag_info == true)
  {
  double E_tmp = E_g_cut;
  E_g_cut = E_g_max;

  double check_e = i_el_e.Integral(theta_min, theta_max); // e-/mu-
  double check_p = i_el_p.Integral(theta_min, theta_max); // e+/mu+
  
  E_g_cut = E_tmp;

  check_e = (phi_max-phi_min)*mkb*(check_e - i_el_e.Integral(theta_min, theta_max))/omega; // e-/mu-
  check_p = (phi_max-phi_min)*mkb*(check_p - i_el_p.Integral(theta_min, theta_max))/omega; //e+/mu+
  
  cout << endl << "Bremsstrahlung cross section, microbarn / steradian" << endl;
  cout << "(E_g from " << 1000.*E_g_cut << " to " << 1000.*E_g_max << " MeV):" << endl;
  cout << "analytical integration, e-/mu-: " << scientific << check_e << ", e+/mu+: " << check_p << endl;
  cout << " numerical integration, e-/mu-: " << (bre1 + bre2)/omega << ", e+/mu+: " << (brp1 + brp2)/omega << endl << endl;
  cout.unsetf(ios_base::floatfield); // The default format for cout
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Generating events for the Rosenbluth scattering:

i = 0;
for (loop = 0; loop < nev_0; loop++)
  {
  Foam_Ros->MakeEvent(); // Generating a Rosenbluth event
  
  // Generating azimuthal angle (phi) for the lepton:
  phi_l = phi_min + (phi_max - phi_min)*(PseRan->Rndm());

  // Azimuthal angle (phi) for the proton:
  if (phi_l < 0.) phi_p = phi_l + Pi;
  else phi_p = phi_l - Pi;

  // If we want to have (0 < phi < 2*Pi) instead of (-Pi < phi < Pi):
  if (flag_phi == 1) { if (phi_p < 0.) phi_p += 2.*Pi; }

  // Generating z-coordinate of the event (if we use a storage cell):
  if (flag_target == true) zcoord = scell();

  // Writing to the output files: --------------------------------------------------------------
  if (flag_root == true) // Writing the event to the *.root file
    {
    if (flag_target == false) ntp_0->Fill(1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g);
    else ntp_0->Fill(1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g, zcoord);
    }

  if (flag_dat == true)  // Writing the event to the *.dat file
    {
    if (flag_target == false) fprintf(f0, "%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf\n", 1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g);
    else fprintf(f0, "%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %6.1lf\n", 1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g, zcoord);
    }

  // The progress status for Rosenbluth events:
  if (loop == long (0.01*i*(nev_0 - 1)))
    {
    cout << "Generation of " << nev_0 << " Rosenbluth events: " << i << "%" << '\r' << flush;
    i += 1;
    }
  }
  
if (nev_0 > 0) cout << endl;
  
// The end of generating events for the Rosenbluth scattering
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Generating events for the case of e-/mu-:

i = 0;
for (loop = 0; loop < nev_e; loop++)
  {
  // "Elastic" scattering events: --------------------------------------------------------------
  if ((nev_ElE - count_ElE) > 0 && PseRan->Rndm() < 1.*(nev_ElE - count_ElE)/(nev_e - count_ElE - count_BrE1 - count_BrE2))
    {
    Foam_ElE->MakeEvent(); // Generating an elastic event, e-/mu-
    count_ElE++; // Counter for elastic scattering events, e-/mu-

    // Generating azimuthal angle (phi) for the lepton:
    phi_l = phi_min + (phi_max - phi_min)*(PseRan->Rndm());

    // Azimuthal angle (phi) for the proton:
    if (phi_l < 0.) phi_p = phi_l + Pi;
    else phi_p = phi_l - Pi;

    // If we want to have (0 < phi < 2*Pi) instead of (-Pi < phi < Pi):
    if (flag_phi == 1) { if (phi_p < 0.) phi_p += 2.*Pi; }
    }
  // Bremsstrahlung events: --------------------------------------------------------------------
  else
    {
    if (nev_BrE2 - count_BrE2 > 0 && count_BrE2/nev_BrE2 < count_BrE1/nev_BrE1)
      {
      en_sign = +1.; // The "+" root of the equation for the scattered lepton energy
      Foam_BrE2->MakeEvent(); // Generating a bremsstrahlung event, e-/mu-, root "+"
      count_BrE2++; // Counter for bremsstrahlung events, e-/mu-, root "+"
      }
    else
      {
      en_sign = -1.; // The "-" root of the equation for the scattered lepton energy
      Foam_BrE1->MakeEvent(); // Generating a bremsstrahlung event, e-/mu-, root "-"
      count_BrE1++; // Counter for bremsstrahlung events, e-/mu-, root "-"
      }
      
    // Generating azimuthal angle (phi) for the lepton:
    phi_l = phi_min + (phi_max - phi_min)*(PseRan->Rndm());

    // Rotation the final particle four-momenta around the z-axis:
    v_lf.RotateZ(phi_l); // Four-momentum of the final lepton
    v_pf.RotateZ(phi_l); // Four-momentum of the final proton
    v_kf.RotateZ(phi_l); // Four-momentum of the photon

    // Kinematic parameters of the final particles after rotation:
    E_lf = v_lf.E(); theta_l = v_lf.Theta(); phi_l = v_lf.Phi(); // Final lepton
    E_p = v_pf.E(); theta_p = v_pf.Theta(); phi_p = v_pf.Phi();  // Final proton
    E_g = v_kf.E(); theta_g = v_kf.Theta(); phi_g = v_kf.Phi();  // Photon

    // If we want to have (0 < phi < 2*Pi) instead of (-Pi < phi < Pi):
    if (flag_phi == 1)
      {
      if (phi_l < 0.) phi_l += 2.*Pi; // Phi angle for the final lepton
      if (phi_p < 0.) phi_p += 2.*Pi; // Phi angle for the final proton
      if (phi_g < 0.) phi_g += 2.*Pi; // Phi angle for the photon
      }
    }

  // Generating z-coordinate of the event (if we use a storage cell):
  if (flag_target == true) zcoord = scell();

  // Writing to the output files: --------------------------------------------------------------
  if (flag_root == true) // Writing the event to the *.root file
    {
    if (flag_target == false) ntp_e->Fill(1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g);
    else ntp_e->Fill(1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g, zcoord);
    }

  if (flag_dat == true)  // Writing the event to the *.dat file
    {
    if (flag_target == false) fprintf(fe, "%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf\n", 1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g);
    else fprintf(fe, "%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %6.1lf\n", 1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g, zcoord);
    }

  // The progress status for e-/mu- events:
  if (loop == long (0.01*i*(nev_e - 1)))
    {
    if (flag_lepton == 1 || flag_lepton == 3) cout << "Generation of " << nev_e << " events for e-: ";
    if (flag_lepton == 4 || flag_lepton == 6) cout << "Generation of " << nev_e << " events for mu-: ";
    cout << i << "%" << '\r' << flush;
    i += 1;
    }
    
  }

if (nev_e > 0) cout << endl;
  
// The end of generating events for the case of e-/mu-
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Generating events for the case of e+/mu+:

i = 0;
for (loop = 0; loop < nev_p; loop++)
  {
  // "Elastic" scattering events: --------------------------------------------------------------
  if ((nev_ElP - count_ElP) > 0 && PseRan->Rndm() < 1.*(nev_ElP - count_ElP)/(nev_p - count_ElP - count_BrP1 - count_BrP2))
    {
    Foam_ElP->MakeEvent(); // Generating an elastic event, e+/mu+
    count_ElP++; // Counter for elastic scattering events, e+/mu+

    // Generating azimuthal angle (phi) for the lepton:
    phi_l = phi_min + (phi_max - phi_min)*(PseRan->Rndm());

    // Azimuthal angle (phi) for the proton:
    if (phi_l < 0.) phi_p = phi_l + Pi;
    else phi_p = phi_l - Pi;

    // If we want to have (0 < phi < 2*Pi) instead of (-Pi < phi < Pi):
    if (flag_phi == 1) { if (phi_p < 0.) phi_p += 2.*Pi; }
    }
  // Bremsstrahlung events: --------------------------------------------------------------------
  else
    { 
    if (nev_BrP2 - count_BrP2 > 0 && count_BrP2/nev_BrP2 < count_BrP1/nev_BrP1)
      {
      en_sign = +1.; // The "+" root of the equation for the scattered lepton energy
      Foam_BrP2->MakeEvent(); // Generating a bremsstrahlung event, e+/mu+, root "+"
      count_BrP2++; // Counter for bremsstrahlung events, e+/mu+, root "+"
      }
    else
      {
      en_sign = -1.; // The "-" root of the equation for the scattered lepton energy
      Foam_BrP1->MakeEvent(); // Generating a bremsstrahlung event, e+/mu+, root "-"
      count_BrP1++; // Counter for bremsstrahlung events, e+/mu+, root "-"
      }

    // Generating azimuthal angle (phi) for the lepton:
    phi_l = phi_min + (phi_max - phi_min)*(PseRan->Rndm());

    // Rotation the final particle four-momenta around the z-axis:
    v_lf.RotateZ(phi_l); // Four-momentum of the final lepton
    v_pf.RotateZ(phi_l); // Four-momentum of the final proton
    v_kf.RotateZ(phi_l); // Four-momentum of the photon

    // Kinematic parameters of the final particles after rotation:
    E_lf = v_lf.E(); theta_l = v_lf.Theta(); phi_l = v_lf.Phi(); // Final lepton
    E_p = v_pf.E(); theta_p = v_pf.Theta(); phi_p = v_pf.Phi();  // Final proton
    E_g = v_kf.E(); theta_g = v_kf.Theta(); phi_g = v_kf.Phi();  // Photon

    // If we want to have (0 < phi < 2*Pi) instead of (-Pi < phi < Pi):
    if (flag_phi == 1)
      {
      if (phi_l < 0.) phi_l += 2.*Pi; // Phi angle for the final lepton
      if (phi_p < 0.) phi_p += 2.*Pi; // Phi angle for the final proton
      if (phi_g < 0.) phi_g += 2.*Pi; // Phi angle for the photon
      }
    }

  // Generating z-coordinate of the event (if we use a storage cell):
  if (flag_target == true) zcoord = scell();

  // Writing to the output files: --------------------------------------------------------------
  if (flag_root == true) // Writing the event to the *.root file
    {
    if (flag_target == false) ntp_p->Fill(1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g);
    else ntp_p->Fill(1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g, zcoord);
    }

  if (flag_dat == true)  // Writing the event to the *.dat file
    {
    if (flag_target == false) fprintf(fp, "%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf\n", 1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g);
    else fprintf(fp, "%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %6.1lf\n", 1000.*E_lf, theta_l, phi_l, 1000.*E_p, theta_p, phi_p, 1000.*E_g, theta_g, phi_g, zcoord);
    }

  // The progress status for e+/mu+ events:
  if (loop == long (0.01*i*(nev_p - 1)))
    {
    if (flag_lepton == 2 || flag_lepton == 3) cout << "Generation of " << nev_p << " events for e+: ";
    if (flag_lepton == 5 || flag_lepton == 6) cout << "Generation of " << nev_p << " events for mu+: ";
    cout << i << "%" << '\r' << flush;
    i += 1;
    }
    
  }
  
if (nev_p > 0) cout << endl;
  
// The end of generating events for the case of e+/mu+
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Closing the output files:
if (flag_dat == true) // Files *.dat
  {
  if (flag_rosen == true) fclose(f0); // Rosenbluth events
  if (flag_lepton != 2 && flag_lepton != 5) fclose(fe); // e-/mu-
  if (flag_lepton != 1 && flag_lepton != 4) fclose(fp); // e+/mu+
  }

if (flag_root == true) // Files *.root
  {
  if (flag_rosen == true) // Rosenbluth events
    {
    froot_0->cd();
    cout << endl;
    ntp_0->Print();
    froot_0->Write();
    froot_0->Close();
    }
    
  if (flag_lepton != 2 && flag_lepton != 5) // e-/mu-
    {
    froot_e->cd();
    cout << endl;
    ntp_e->Print();
    froot_e->Write();
    froot_e->Close();
    }

  if (flag_lepton != 1 && flag_lepton != 4) // e+/mu+
    {
    froot_p->cd();
    cout << endl;
    ntp_p->Print();
    froot_p->Write();
    froot_p->Close();
    }    
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Finalization of Foam simulators:
double x1, x2;
if (flag_rosen == true) // For Rosenbluth events
  {
  cout << endl << "Finalization (Rosenbluth scattering):" << endl;
  Foam_Ros->Finalize(x1, x2); // Foam_Ros: Rosenbluth scattering
  }

if (flag_lepton != 2 && flag_lepton != 5) // For negatively charged leptons
  {
  if (flag_lepton == 1 || flag_lepton == 3) cout << endl << "Finalization (elastic scattering, e-):" << endl;
  if (flag_lepton == 4 || flag_lepton == 6) cout << endl << "Finalization (elastic scattering, mu-):" << endl;
  Foam_ElE->Finalize(x1, x2); // Foam_ElE: elastic scattering, e-/mu-
  
  if (flag_lepton == 1 || flag_lepton == 3) cout << endl << "Finalization (bremsstrahlung, e-, root \"-\"):" << endl;
  if (flag_lepton == 4 || flag_lepton == 6) cout << endl << "Finalization (bremsstrahlung, mu-, root \"-\"):" << endl;
  Foam_BrE1->Finalize(x1, x2); // Foam_BrE1: bremsstrahlung, e-/mu-, root "-"
  
  if (bre2 > 0)
    {
    if (flag_lepton == 1 || flag_lepton == 3) cout << endl << "Finalization (bremsstrahlung, e-, root \"+\"):" << endl;
    if (flag_lepton == 4 || flag_lepton == 6) cout << endl << "Finalization (bremsstrahlung, mu-, root \"+\"):" << endl;
    Foam_BrE2->Finalize(x1, x2); // Foam_BrE2: bremsstrahlung, e-/mu-, root "+"
    }
  }

if (flag_lepton != 1 && flag_lepton != 4) // For positively charged leptons
  {
  if (flag_lepton == 2 || flag_lepton == 3) cout << endl << "Finalization (elastic scattering, e+):" << endl;
  if (flag_lepton == 5 || flag_lepton == 6) cout << endl << "Finalization (elastic scattering, mu+):" << endl;
  Foam_ElP->Finalize(x1, x2); // Foam_ElP: elastic scattering, e+/mu+
  
  if (flag_lepton == 2 || flag_lepton == 3) cout << endl << "Finalization (bremsstrahlung, e+, root \"-\"):" << endl;
  if (flag_lepton == 5 || flag_lepton == 6) cout << endl << "Finalization (bremsstrahlung, mu+, root \"-\"):" << endl;
  Foam_BrP1->Finalize(x1, x2); // Foam_BrP1: bremsstrahlung, e+/mu+, root "-"
  
  if (brp2 > 0)
    {
    if (flag_lepton == 2 || flag_lepton == 3) cout << endl << "Finalization (bremsstrahlung, e+, root \"+\"):" << endl;
    if (flag_lepton == 5 || flag_lepton == 6) cout << endl << "Finalization (bremsstrahlung, mu+, root \"+\"):" << endl;
    Foam_BrP2->Finalize(x1, x2); // Foam_BrP2: bremsstrahlung, e+/mu+, root "+"
    }
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Writing to the standard output and the *.info file:
ofstream finfo;
if (flag_lepton < 4) finfo.open (filename + "_e.info");  // Opening the *.info file for e-/e+
if (flag_lepton > 3) finfo.open (filename + "_mu.info"); // Opening the *.info file for mu-/mu+

if (finfo.is_open())
  {    
  cout << endl << "INPUT PARAMETERS AND SOME INFORMATION ON THE EVENTS:" << endl << endl;
  finfo << "INPUT PARAMETERS AND SOME INFORMATION ON THE EVENTS:" << endl << endl;

  // Type of the scattered lepton:
  switch (flag_lepton)
    {
    case 1: // e-
      {
      cout << "* Elastic scattering of electrons (e-) on protons;" << endl;
      finfo << "* Elastic scattering of electrons (e-) on protons;" << endl;
      break;
      }
    case 2: // e+
      {
      cout << "* Elastic scattering of positrons (e+) on protons;" << endl;
      finfo << "* Elastic scattering of positrons (e+) on protons;" << endl;
      break;
      }
    case 3: // e- and e+
      {
      cout << "* Elastic scattering of electrons and positrons (e- and e+) on protons;" << endl;
      finfo << "* Elastic scattering of electrons and positrons (e- and e+) on protons;" << endl;
      break;
      }
    case 4: // mu-
      {
      cout << "* Elastic scattering of muons (mu-) on protons;" << endl;
      finfo << "* Elastic scattering of muons (mu-) on protons;" << endl;
      break;
      }
    case 5: // mu+
      {
      cout << "* Elastic scattering of antimuons (mu+) on protons;" << endl;
      finfo << "* Elastic scattering of antimuons (mu+) on protons;" << endl;
      break;
      }
    case 6: // mu- and mu+
      {
      cout << "* Elastic scattering of muons and antimuons (mu- and mu+) on protons;" << endl;
      finfo << "* Elastic scattering of muons and antimuons (mu- and mu+) on protons;" << endl;
      break;
      }
    }

  // Internal structure of the proton:
  switch (flag_struct)
    {
    case 1: // A point-like proton
      {
      cout << "* A point-like proton;" << endl;
      finfo << "* A point-like proton;" << endl;
      break;
      }
    case 2: // A proton with the dipole form factors
      {
      cout << "* A proton with the dipole form factors;" << endl;
      finfo << "* A proton with the dipole form factors;" << endl;
      break;
      }
    case 3: // A proton with the Kelly form factors
      {
      cout << "* A proton with the Kelly form factors;" << endl;
      finfo << "* A proton with the Kelly form factors;" << endl;
      break;
      }
    case 4: // A proton with the Puckett form factors
      {
      cout << "* A proton with the Puckett form factors;" << endl;
      finfo << "* A proton with the Puckett form factors;" << endl;
      break;
      }
    case 5: // A proton with the form factors specified in the file "const.h"
      {
      cout << "* A proton with the form factors specified in \"const.h\";" << endl;
      finfo << "* A proton with the form factors specified in \"const.h\";" << endl;
      break;
      }
    }
    
  // Mode of calculation for bremsstrahlung:
  switch (flag_mode)
    {
    case 1: // Primary soft-photon approximation, only lepton bremsstrahlung
      {
      cout << "* Primary soft-photon approximation, only lepton bremsstrahlung;" << endl;
      finfo << "* Primary soft-photon approximation, only lepton bremsstrahlung;" << endl;
      break;
      }
    case 2: // Primary soft-photon approximation, only proton bremsstrahlung
      {
      cout << "* Primary soft-photon approximation, only proton bremsstrahlung;" << endl;
      finfo << "* Primary soft-photon approximation, only proton bremsstrahlung;" << endl;
      break;
      }
    case 3: // Primary soft-photon approximation, all the terms
      {
      cout << "* Primary soft-photon approximation, all the terms;" << endl;
      finfo << "* Primary soft-photon approximation, all the terms;" << endl;
      break;
      }
    case 4: // Modified soft-photon approximation, only lepton bremsstrahlung
      {
      cout << "* Modified soft-photon approximation, only lepton bremsstrahlung;" << endl;
      finfo << "* Modified soft-photon approximation, only lepton bremsstrahlung;" << endl;
      break;
      }
    case 5: // Modified soft-photon approximation, only proton bremsstrahlung
      {
      cout << "* Modified soft-photon approximation, only proton bremsstrahlung;" << endl;
      finfo << "* Modified soft-photon approximation, only proton bremsstrahlung;" << endl;
      break;
      }
    case 6: // Modified soft-photon approximation, all the terms
      {
      cout << "* Modified soft-photon approximation, all the terms;" << endl;
      finfo << "* Modified soft-photon approximation, all the terms;" << endl;
      break;
      }
    case 7: // Improved soft-photon approximation, only lepton bremsstrahlung
      {
      cout << "* Improved soft-photon approximation, only lepton bremsstrahlung;" << endl;
      finfo << "* Improved soft-photon approximation, only lepton bremsstrahlung;" << endl;
      break;
      }
    case 8: // Improved soft-photon approximation, only proton bremsstrahlung
      {
      cout << "* Improved soft-photon approximation, only proton bremsstrahlung;" << endl;
      finfo << "* Improved soft-photon approximation, only proton bremsstrahlung;" << endl;
      break;
      }
    case 9: // Improved soft-photon approximation, all the terms
      {
      cout << "* Improved soft-photon approximation, all the terms;" << endl;
      finfo << "* Improved soft-photon approximation, all the terms;" << endl;
      break;
      }
    case 10: // Accurate QED calculation, only lepton bremsstrahlung
      {
      cout << "* Accurate QED calculation, only lepton bremsstrahlung;" << endl;
      finfo << "* Accurate QED calculation, only lepton bremsstrahlung;" << endl;
      break;
      }
    case 11: // Accurate QED calculation, only proton bremsstrahlung
      {
      cout << "* Accurate QED calculation, only proton bremsstrahlung;" << endl;
      finfo << "* Accurate QED calculation, only proton bremsstrahlung;" << endl;
      break;
      }
    case 12: // Accurate QED calculation, all the terms
      {
      cout << "* Accurate QED calculation, all the terms;" << endl;
      finfo << "* Accurate QED calculation, all the terms;" << endl;
      break;
      }
    }

  // Vacuum polarization correction:
  switch (flag_vpol)
    {
    case 1: // Only electron-positron loops
      {
      cout << "* Vacuum polarization: only electron-positron loops;" << endl;
      finfo << "* Vacuum polarization: only electron-positron loops;" << endl;
      break;
      }
    case 2: // Full leptonic contribution
      {
      cout << "* Vacuum polarization: full leptonic contribution;" << endl;
      finfo << "* Vacuum polarization: full leptonic contribution;" << endl;
      break;
      }
    case 3: // More accurate data (see http://cmd.inp.nsk.su/~ignatov/vpl/)
      {
      cout << "* Vacuum polarization: full vacuum polarization correction;" << endl;
      finfo << "* Vacuum polarization: full vacuum polarization correction;" << endl;
      break;
      }
    }
    
  // Two-photon exchange amplitudes:
  switch (flag_tpe)
    {
    case 1: // Approach of Mo & Tsai
      {
      cout << "* Approach of Mo & Tsai for the TPE amplitudes;" << endl;
      finfo << "* Approach of Mo & Tsai for the TPE amplitudes;" << endl;
      break;
      }
    case 2: // Approach of Maximon & Tjon
      {
      cout << "* Approach of Maximon & Tjon for the TPE amplitudes;" << endl;
      finfo << "* Approach of Maximon & Tjon for the TPE amplitudes;" << endl;
      break;
      }
    }
    
  // Full energy of incident leptons (MeV):
  cout << "* Full energy of incident leptons: " << E_li*1000. << " MeV;" << endl;
  finfo << "* Full energy of incident leptons: " << E_li*1000. << " MeV;" << endl;

  // Cut-off and maximum energies of bremsstrahlung photons (MeV):
  cout << "* Energies of bremsstrahlung photons: from " << E_g_cut*1000. << " to " << E_g_max*1000. << " MeV;" << endl;
  finfo << "* Energies of bremsstrahlung photons: from " << E_g_cut*1000. << " to " << E_g_max*1000. << " MeV;" << endl;

  // Minimum and maximum polar angles (theta) of scattered leptons (degree):
  cout << "* Polar angles of scattered leptons (theta): from " << theta_min/degrad << " to " << theta_max/degrad << " degree;" << endl;
  finfo << "* Polar angles of scattered leptons (theta): from " << theta_min/degrad << " to " << theta_max/degrad << " degree;" << endl;

  // Minimum and maximum azimuthal angles (phi) of scattered leptons (degree):
  cout << "* Azimuthal angles of scattered leptons (phi): from " << phi_min/degrad << " to " << phi_max/degrad << " degree;" << endl;
  finfo << "* Azimuthal angles of scattered leptons (phi): from " << phi_min/degrad << " to " << phi_max/degrad << " degree;" << endl;
  
  // Storage cell parameters (if it's used):
  if (flag_target == true)
    {
    if (flag_cell == 1) // In the case of uniform distribution of gas pressure
      {
      cout << "* Uniform distribution of gas pressure; Z from " << cell_zmin << " to " << cell_zmax << " mm;" << endl;
      finfo << "* Uniform distribution of gas pressure; Z from " << cell_zmin << " to " << cell_zmax << " mm;" << endl;
      }
  
    if (flag_cell == 2) // In the case of triangle distribution of gas pressure
      {
      cout << "* Triangle distribution of gas pressure; Z from " << cell_zmin << " to " << cell_zmax << " mm, injection at Z = " << cell_zinj << " mm;" << endl;
      finfo << "* Triangle distribution of gas pressure; Z from " << cell_zmin << " to " << cell_zmax << " mm, injection at Z = " << cell_zinj << " mm;" << endl;
      }
    }
  
  // Precison for numbers to output:
  cout.precision(3);
  finfo.precision(3);
  
  // Solid angle:
  cout << "* Solid angle: " << scientific << omega << " steradian;" << endl;
  finfo << "* Solid angle: " << scientific << omega << " steradian;" << endl;

  // Evaluation of Q^2 and epsilon values for the minimum scattering angle:
  E_g = 0.; theta_g = 0.; phi_g = 0.; // Elastic scattering
  theta_l = theta_min; // Minimum scattering angle
  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton
  EvalKinematicParams(0.); // Evaluation of some kinematic parameters
  double QQ_min = -qq; double eps_max = eps;

  // Evaluation of Q^2 and epsilon values for the maximum scattering angle:
  theta_l = theta_max; // Maximum scattering angle
  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton
  EvalKinematicParams(0.); // Evaluation of some kinematic parameters
  double QQ_max = -qq; double eps_min = eps;
  
  // Writing the minimum and maximum values for Q^2 and epsilon:
  cout << "* Values of Q^2: from " << QQ_min << " to " << QQ_max << " GeV^2;" << endl; // Writing to the standard output
  cout << "* Values of epsilon: from " << fixed << eps_min << " to " << fixed << eps_max << ";" << endl;
  
  finfo << "* Values of Q^2: from " << QQ_min << " to " << QQ_max << " GeV^2;" << endl; // Writing to the *.info file
  finfo << "* Values of epsilon: from " << fixed << eps_min << " to " << fixed << eps_max << ";" << endl;
  
  // Are there any warnings?
  if (flag_warn == false) // There are no warnings
    {
    cout << "* Warnings: no." << endl;
    finfo << "* Warnings: no." << endl;
    }
  else // There are some warnings
    {
    cout << "* Warnings: yes. Please be aware of them!" << endl;
    finfo << "* Warnings: yes. Please be aware of them!" << endl;
    }

  // Writing the table containing the number of events:
  cout << endl << "|---------------------------------------|" << endl; // Writing to the standard output
  cout << "|           NUMBERS OF EVENTS           |" << endl;
  cout << "|---------------------------------------|" << endl;
  if (flag_lepton < 4) cout << "|             |     e-     |     e+     |" << endl;
  if (flag_lepton > 3) cout << "|             |     mu-    |     mu+    |" << endl;
  cout << "|---------------------------------------|" << endl;
  cout << "|    Elastic  |  ";
  if (flag_lepton != 2 && flag_lepton != 5) cout << setw(9) << count_ElE << " |  ";
  else cout << "       -- |  ";
  if (flag_lepton != 1 && flag_lepton != 4) cout << setw(9) << count_ElP << " |" << endl;
  else cout << "       -- |" << endl;
  cout << "|  Inelastic  |  ";
  if (flag_lepton != 2 && flag_lepton != 5) cout << setw(9) << count_BrE1 + count_BrE2 << " |  ";
  else cout << "       -- |  ";
  if (flag_lepton != 1 && flag_lepton != 4) cout << setw(9) << count_BrP1 + count_BrP2 << " |" << endl;
  else cout << "       -- |" << endl;
  cout << "|---------------------------------------|" << endl;
  cout << "|      Total  |  ";
  if (flag_lepton != 2 && flag_lepton != 5) cout << setw(9) << count_ElE + count_BrE1 + count_BrE2 << " |  ";
  else cout << "       -- |  ";
  if (flag_lepton != 1 && flag_lepton != 4) cout << setw(9) << count_ElP + count_BrP1 + count_BrP2 << " |" << endl;
  else cout << "       -- |" << endl;
  cout << "|---------------------------------------|" << endl;
  if (flag_rosen == true) { cout << "| Rosenbluth  |               " << setw(9) << nev_0 << " |" << endl;
  cout << "|---------------------------------------|" << endl; }

  finfo << endl << "|---------------------------------------|" << endl; // Writing to the *.info file
  finfo << "|            NUMBER OF EVENTS           |" << endl;
  finfo << "|---------------------------------------|" << endl;
  if (flag_lepton < 4) finfo << "|             |     e-     |     e+     |" << endl;
  if (flag_lepton > 3) finfo << "|             |     mu-    |     mu+    |" << endl;
  finfo << "|---------------------------------------|" << endl;
  finfo << "|    Elastic  |  ";
  if (flag_lepton != 2 && flag_lepton != 5) finfo << setw(9) << count_ElE << " |  ";
  else finfo << "       -- |  ";
  if (flag_lepton != 1 && flag_lepton != 4) finfo << setw(9) << count_ElP << " |" << endl;
  else finfo << "       -- |" << endl;
  finfo << "|  Inelastic  |  ";
  if (flag_lepton != 2 && flag_lepton != 5) finfo << setw(9) << count_BrE1 + count_BrE2 << " |  ";
  else finfo << "       -- |  ";
  if (flag_lepton != 1 && flag_lepton != 4) finfo << setw(9) << count_BrP1 + count_BrP2 << " |" << endl;
  else finfo << "       -- |" << endl;
  finfo << "|---------------------------------------|" << endl;
  finfo << "|      Total  |  ";
  if (flag_lepton != 2 && flag_lepton != 5) finfo << setw(9) << count_ElE + count_BrE1 + count_BrE2 << " |  ";
  else finfo << "       -- |  ";
  if (flag_lepton != 1 && flag_lepton != 4) finfo << setw(9) << count_ElP + count_BrP1 + count_BrP2 << " |" << endl;
  else finfo << "       -- |" << endl;
  finfo << "|---------------------------------------|" << endl;
  if (flag_rosen == true) { finfo << "| Rosenbluth  |               " << setw(9) << nev_0 << " |" << endl;
  finfo << "|---------------------------------------|" << endl; }

  // Precison for numbers to output:
  cout.precision(4);
  finfo.precision(4);
  
  // Writing the table containing the values for integrated cross sections:
  cout << endl << "|---------------------------------------|" << endl; // Writing to the standard output
  cout << "| INTEGRATED CROSS SECTIONS (microbarn) |" << endl;
  cout << "|---------------------------------------|" << endl;
  if (flag_lepton < 4) cout << "|             |     e-     |     e+     |" << endl;
  if (flag_lepton > 3) cout << "|             |     mu-    |     mu+    |" << endl;
  cout << "|---------------------------------------|" << endl;
  cout << "|    Elastic  | ";
  if (flag_lepton != 2 && flag_lepton != 5) cout << scientific << elast_e << " | ";
  else cout << "        -- | ";
  if (flag_lepton != 1 && flag_lepton != 4) cout << scientific << elast_p << " |" << endl;
  else cout << "        -- |" << endl;
  cout << "|  Inelastic  | ";
  if (flag_lepton != 2 && flag_lepton != 5) cout << bre1 + bre2 << " | ";
  else cout << "        -- | ";
  if (flag_lepton != 1 && flag_lepton != 4) cout << brp1 + brp2 << " |" << endl;
  else cout << "        -- |" << endl;
  cout << "|---------------------------------------|" << endl;
  cout << "|      Total  | ";
  if (flag_lepton != 2 && flag_lepton != 5) cout << elast_e + bre1 + bre2 << " | ";
  else cout << "        -- | ";
  if (flag_lepton != 1 && flag_lepton != 4) cout << elast_p + brp1 + brp2 << " |" << endl;
  else cout << "        -- |" << endl;
  cout << "|---------------------------------------|" << endl;
  if (flag_rosen == true) { cout << "| Rosenbluth  |              " << rosen << " |" << endl;
  cout << "|---------------------------------------|" << endl; }
  
  finfo << endl << "|---------------------------------------|" << endl; // Writing to the *.info file
  finfo << "| INTEGRATED CROSS SECTIONS (microbarn) |" << endl;
  finfo << "|---------------------------------------|" << endl;
  if (flag_lepton < 4) finfo << "|             |     e-     |     e+     |" << endl;
  if (flag_lepton > 3) finfo << "|             |     mu-    |     mu+    |" << endl;
  finfo << "|---------------------------------------|" << endl;
  finfo << "|    Elastic  | ";
  if (flag_lepton != 2 && flag_lepton != 5) finfo << scientific << elast_e << " | ";
  else finfo << "        -- | ";
  if (flag_lepton != 1 && flag_lepton != 4) finfo << scientific << elast_p << " |" << endl;
  else finfo << "        -- |" << endl;
  finfo << "|  Inelastic  | ";
  if (flag_lepton != 2 && flag_lepton != 5) finfo << bre1 + bre2 << " | ";
  else finfo << "        -- | ";
  if (flag_lepton != 1 && flag_lepton != 4) finfo << brp1 + brp2 << " |" << endl;
  else finfo << "        -- |" << endl;
  finfo << "|---------------------------------------|" << endl;
  finfo << "|      Total  | ";
  if (flag_lepton != 2 && flag_lepton != 5) finfo << elast_e + bre1 + bre2 << " | ";
  else finfo << "        -- | ";
  if (flag_lepton != 1 && flag_lepton != 4) finfo << elast_p + brp1 + brp2 << " |" << endl;
  else finfo << "        -- |" << endl;
  finfo << "|---------------------------------------|" << endl;
  if (flag_rosen == true) { finfo << "| Rosenbluth  |              " << rosen << " |" << endl;
  finfo << "|---------------------------------------|" << endl; }

  // Writing the Rosenbluth differential cross section (averaged over the solid angle):
  cout << endl << "ROSENBLUTH DIFFERENTIAL CROSS SECTION (averaged over the solid angle):" << endl;
  cout << rosen/omega << " microbarn / steradian" << endl;

  finfo << endl << "ROSENBLUTH DIFFERENTIAL CROSS SECTION (averaged over the solid angle):" << endl;
  finfo << rosen/omega << " microbarn / steradian" << endl;
  
  // Writing the real differential cross sections averaged over the solid angle:
  cout << endl << "ACTUAL DIFFERENTIAL CROSS SECTION (averaged over the solid angle):" << endl; // Writing to the standard output
  if (flag_lepton == 1 || flag_lepton == 3) cout << "e-: ";
  if (flag_lepton == 4 || flag_lepton == 6) cout << "mu-: ";
  if (flag_lepton != 2 && flag_lepton != 5) cout << (elast_e + bre1 + bre2)/omega << " microbarn / steradian" << fixed << " (delta = " << (elast_e + bre1 + bre2)/rosen - 1. << ")" << endl;
  if (flag_lepton == 2 || flag_lepton == 3) cout << "e+: ";
  if (flag_lepton == 5 || flag_lepton == 6) cout << "mu+: ";
  if (flag_lepton != 1 && flag_lepton != 4) cout << scientific << (elast_p + brp1 + brp2)/omega << " microbarn / steradian" << fixed << " (delta = " << (elast_p + brp1 + brp2)/rosen - 1. << ")" << endl;
  
  finfo << endl << "ACTUAL DIFFERENTIAL CROSS SECTION (averaged over the solid angle):" << endl; // Writing to the *.info file
  if (flag_lepton == 1 || flag_lepton == 3) finfo << "e-: ";
  if (flag_lepton == 4 || flag_lepton == 6) finfo << "mu-: ";
  if (flag_lepton != 2 && flag_lepton != 5) finfo << (elast_e + bre1 + bre2)/omega << " microbarn / steradian" << fixed << " (delta = " << (elast_e + bre1 + bre2)/rosen - 1. << ")" << endl;
  if (flag_lepton == 2 || flag_lepton == 3) finfo << "e+: ";
  if (flag_lepton == 5 || flag_lepton == 6) finfo << "mu+: ";
  if (flag_lepton != 1 && flag_lepton != 4) finfo << scientific << (elast_p + brp1 + brp2)/omega << " microbarn / steradian" << fixed << " (delta = " << (elast_p + brp1 + brp2)/rosen - 1. << ")" << endl;

  // Writing the integrated luminosity (inverse picobarn):
  cout << endl << "INTEGRATED LUMINOSITY:" << endl; // Writing to the standard output
  if (flag_rosen == true && flag_lepton < 4) cout << "e0: ";
  if (flag_rosen == true && flag_lepton > 3) cout << "mu0: ";
  if (flag_rosen == true) cout << scientific << 1.e-6*nev_0/rosen << " inverse picobarn" << endl;
  if (flag_lepton == 1 || flag_lepton == 3) cout << "e-: ";
  if (flag_lepton == 4 || flag_lepton == 6) cout << "mu-: ";
  if (flag_lepton != 2 && flag_lepton != 5) cout << scientific << 1.e-6*nev_e/(elast_e + bre1 + bre2) << " inverse picobarn" << endl;
  if (flag_lepton == 2 || flag_lepton == 3) cout << "e+: ";
  if (flag_lepton == 5 || flag_lepton == 6) cout << "mu+: ";
  if (flag_lepton != 1 && flag_lepton != 4) cout << scientific << 1.e-6*nev_p/(elast_p + brp1 + brp2) << " inverse picobarn" << endl;

  finfo << endl << "INTEGRATED LUMINOSITY:" << endl; // Writing to the *.info file
  if (flag_rosen == true && flag_lepton < 4) finfo << "e0: ";
  if (flag_rosen == true && flag_lepton > 3) finfo << "mu0: ";
  if (flag_rosen == true) finfo << scientific << 1.e-6*nev_0/rosen << " inverse picobarn" << endl;
  if (flag_lepton == 1 || flag_lepton == 3) finfo << "e-: ";
  if (flag_lepton == 4 || flag_lepton == 6) finfo << "mu-: ";
  if (flag_lepton != 2 && flag_lepton != 5) finfo << scientific << 1.e-6*nev_e/(elast_e + bre1 + bre2) << " inverse picobarn" << endl;
  if (flag_lepton == 2 || flag_lepton == 3) finfo << "e+: ";
  if (flag_lepton == 5 || flag_lepton == 6) finfo << "mu+: ";
  if (flag_lepton != 1 && flag_lepton != 4) finfo << scientific << 1.e-6*nev_p/(elast_p + brp1 + brp2) << " inverse picobarn" << endl;
  
  if (flag_vepp == true) // In the case of the experiment at VEPP-3
    {
    cout << endl << "BEAM CURRENT INTEGRAL" << endl; // Writing to the standard output
    cout << "(at the target thickness of 10^15 atoms / cm^2):" << endl;
    if (flag_rosen == true) cout << "e0: " << 1.602e-7*nev_0/rosen << " kC" << endl;
    cout << "e-: " << 1.602e-7*nev_e/(elast_e + bre1 + bre2) << " kC" << endl;
    cout << "e+: " << 1.602e-7*nev_p/(elast_p + brp1 + brp2) << " kC" << endl;

    finfo << endl << "BEAM CURRENT INTEGRAL" << endl; // Writing to the *.info file
    finfo << "(at the target thickness of 10^15 atoms / cm^2):" << endl;
    if (flag_rosen == true) finfo << "e0: " << 1.602e-7*nev_0/rosen << " kC" << endl;
    finfo << "e-: " << 1.602e-7*nev_e/(elast_e + bre1 + bre2) << " kC" << endl;
    finfo << "e+: " << 1.602e-7*nev_p/(elast_p + brp1 + brp2) << " kC" << endl;
    }
 
  time (&stoptime); // Stopping the timer
  
  // How many seconds it took to generate events:
  cout << endl << "It took " << stoptime - starttime << " seconds to generate " << nevents << " events" << endl << endl;
  finfo << endl << "It took " << stoptime - starttime << " seconds to generate " << nevents << " events" << endl << endl;
  
  finfo.close(); // Closing the *.info file
  }
else { cout << endl << "Error opening *.info file!"; return EXIT_FAILURE; }
// The end of writing
//==============================================================================================

return EXIT_SUCCESS;
}

