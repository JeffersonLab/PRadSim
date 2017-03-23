//----------------------------------------------------------------------------------------------
// This file is part of ESEPP (version 1.4).
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


#include "TMath.h"


//----------------------------------------------------------------------------------------------
// Some mathematical and physical constants:
const double Pi = TMath::Pi();
const double m_e = 0.51099893e-3;          // Mass of the electron/positron (in GeV)
const double m_mu = 105.658372e-3;         // Mass of the muon/antimuon (in GeV)
const double m_tau = 1.77682;              // Mass of the tau/antitau (in GeV)
const double M = 938.272046e-3;            // Mass of the proton (in GeV)
const double M2 = TMath::Power(M,2);
const double M4 = TMath::Power(M,4);
const double mu = 2.79284736;              // Magnetic moment of the proton
const double alpha = 1./137.036;           // Fine-structure constant
const double e = TMath::Sqrt(4.*Pi*alpha); // Electron charge magnitude

const double degrad = TMath::Pi()/180.; // Degree to radian conversion

const double mb = 0.389379304; // GeV^{-2} to mbarn conversion (GeV^{-2} = 0.389379304 millibarn)
const double mkb = 389.379404; // GeV^{-2} to mkbarn conversion


//----------------------------------------------------------------------------------------------
// Kelly parametrization of the proton form factors (J.J. Kelly, PRC 70 (2004) 068202):
const double a11_K = -0.24; // Electric form factor
const double b11_K = 10.98;
const double b12_K = 12.82;
const double b13_K = 21.97;

const double a21_K = 0.12; // Magnetic form factor
const double b21_K = 10.97;
const double b22_K = 18.86;
const double b23_K = 6.55;


//----------------------------------------------------------------------------------------------
// Puckett parametrization of the proton form factors (A.J.R. Puckett, arXiv:1008.0855):
const double a11_P = -0.299; // Electric form factor
const double b11_P = 11.11;
const double b12_P = 14.11;
const double b13_P = 15.7;

const double a21_P = 0.081; // Magnetic form factor
const double b21_P = 11.15;
const double b22_P = 18.45;
const double b23_P = 5.31;


//----------------------------------------------------------------------------------------------
// You can use your own parametrization of the proton form factors:
const double a11 = 2.90966; // Electric form factor
const double a12 = -1.11542229;
const double a13 = 3.866171e-2;
const double b11 = 14.5187212;
const double b12 = 40.88333;
const double b13 = 99.999998;
const double b14 = 4.579e-5;
const double b15 = 10.3580447;

const double a21 = -1.43573; // Magnetic form factor
const double a22 = 1.19052066;
const double a23 = 2.5455841e-1;
const double b21 = 9.70703681;
const double b22 = 3.7357e-4;
const double b23 = 6.0e-8;
const double b24 = 9.9527277;
const double b25 = 12.7977739;

