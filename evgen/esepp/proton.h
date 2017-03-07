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


//----------------------------------------------------------------------------------------------
// First-order bremsstrahlung, proton term:
double pterm ()
{
res = -((pow(F22,2)*pow(kfpi,4)*(-lilf + m2)*M2*
       ((-2*pow(F10,2) + pow(F20,2))*kfpf + 2*pow(F10,2)*M2) + 2*pow(F10,2)*pow(kfpf,2)*M4*
       (pow(F22,2)*pow(kfpf,2)*(-lilf + m2) - 
       4*pow(F12,2)*kflf*lipf*M2 - 
       4*F12*F22*kflf*lipf*M2 - 
       pow(F22,2)*kflf*lipf*M2 - 
       4*F12*F22*lfpf*lipf*M2 - 
       3*pow(F22,2)*lfpf*lipf*M2 + 
       4*pow(F12,2)*lfpi*lipf*M2 + 
       4*F12*F22*lfpi*lipf*M2 + 
       pow(F22,2)*lfpi*lipf*M2 + 
       4*F12*F22*kflf*lipi*M2 + 
       3*pow(F22,2)*kflf*lipi*M2 + 
       4*pow(F12,2)*lfpf*lipi*M2 + 
       4*F12*F22*lfpf*lipi*M2 + 
       pow(F22,2)*lfpf*lipi*M2 - 
       4*F12*F22*lfpi*lipi*M2 - 
       3*pow(F22,2)*lfpi*lipi*M2 - 
       4*pow(F12,2)*lilf*M4 - 4*F12*F22*lilf*M4 - 
       pow(F22,2)*lilf*M4 + 
       8*pow(F12,2)*m2*M4 + 
       12*F12*F22*m2*M4 + 
       5*pow(F22,2)*m2*M4 - 
       pow(F22,2)*kflf*lipf*pipf + pow(F22,2)*lfpf*lipf*pipf + 
       pow(F22,2)*lfpi*lipf*pipf - pow(F22,2)*kflf*lipi*pipf + 
       pow(F22,2)*lfpf*lipi*pipf + pow(F22,2)*lfpi*lipi*pipf + 
       4*F12*F22*lilf*M2*pipf + 
       2*pow(F22,2)*lilf*M2*pipf - 
       4*pow(F12,2)*m2*M2*pipf - 
       12*F12*F22*m2*M2*pipf - 
       6*pow(F22,2)*m2*M2*pipf - 
       pow(F22,2)*lilf*pow(pipf,2) + 
       pow(F22,2)*m2*pow(pipf,2) + 
       kfpf*(-(pow(F22,2)*lfpf*lipf) - pow(F22,2)*lfpi*lipf - 
       pow(F22,2)*lfpf*lipi - pow(F22,2)*lfpi*lipi + 
       pow(F22,2)*kflf*(lipf + lipi) + 
       4*pow(F12,2)*m2*M2 + 
       12*F12*F22*m2*M2 + 
       6*pow(F22,2)*m2*M2 - 
       2*pow(F22,2)*m2*pipf - 
       2*F22*lilf*((2*F12 + F22)*M2 - F22*pipf)) - 
       kfli*(F22*kflf*(F22*kfpf + (4*F12 + 3*F22)*M2 - 
       F22*pipf) - F22*lfpi*(F22*kfpf + (4*F12 + 3*F22)*M2 - F22*pipf) + 
       lfpf*(-(pow(F22,2)*kfpf) + pow(2*F12 + F22,2)*M2 + 
       pow(F22,2)*pipf))) - 
       pow(kfpi,3)*(-4*pow(F20,2)*pow(F22,2)*pow(kfpf,3)*
       m2 + F22*kfli*(-(F22*kflf*M2*((-2*pow(F10,2) + pow(F20,2))*kfpf + 
       2*pow(F10,2)*M2)) - 
       F22*lfpi*M2*((-2*pow(F10,2) + pow(F20,2))*kfpf + 
       2*pow(F10,2)*M2) + 
       lfpf*(2*pow(F20,2)*F22*pow(kfpf,2) + 
       (4*pow(F10,2)*F22 + 2*F10*F20*(F12 + F22) + 
       pow(F20,2)*(2*F12 + F22))*kfpf*M2 - 
       2*pow(F10,2)*F22*M4)) + 
       F22*pow(kfpf,2)*(2*pow(F20,2)*F22*kflf*lipf + 
       4*pow(F20,2)*F22*lfpf*lipf + 2*pow(F20,2)*F22*lfpi*lipf + 
       2*pow(F20,2)*F22*lfpf*lipi - 
       8*F10*F12*F20*m2*M2 - 
       8*F12*pow(F20,2)*m2*M2 + 
       2*pow(F10,2)*F22*m2*M2 - 
       8*F10*F20*F22*m2*M2 - 
       7*pow(F20,2)*F22*m2*M2 + 
       2*pow(F20,2)*F22*m2*pipf - 
       F22*lilf*((2*pow(F10,2) + pow(F20,2))*M2 + 
       2*pow(F20,2)*pipf)) + 
       2*pow(F10,2)*M4*(-(pow(F22,2)*lfpf*lipf) - pow(F22,2)*lfpi*lipf - 
       pow(F22,2)*lfpf*lipi - pow(F22,2)*lfpi*lipi - 
       pow(F22,2)*kflf*(lipf + lipi) + 
       4*pow(F12,2)*m2*M2 + 
       12*F12*F22*m2*M2 + 
       6*pow(F22,2)*m2*M2 - 
       2*pow(F22,2)*m2*pipf - 
       2*F22*lilf*((2*F12 + F22)*M2 - F22*pipf)) + 
       kfpf*M2*(-8*pow(F12,2)*pow(F20,2)*lfpf*lipf - 
       12*F10*F12*F20*F22*lfpf*lipf - 
       12*F12*pow(F20,2)*F22*lfpf*lipf + 
       6*pow(F10,2)*pow(F22,2)*lfpf*lipf - 
       12*F10*F20*pow(F22,2)*lfpf*lipf - 
       5*pow(F20,2)*pow(F22,2)*lfpf*lipf + 
       2*F10*F12*F20*F22*lfpi*lipf + 
       2*F12*pow(F20,2)*F22*lfpi*lipf + 
       4*pow(F10,2)*pow(F22,2)*lfpi*lipf + 
       2*F10*F20*pow(F22,2)*lfpi*lipf + 
       pow(F20,2)*pow(F22,2)*lfpi*lipf + 
       2*F10*F12*F20*F22*lfpf*lipi + 
       2*F12*pow(F20,2)*F22*lfpf*lipi + 
       4*pow(F10,2)*pow(F22,2)*lfpf*lipi + 
       2*F10*F20*pow(F22,2)*lfpf*lipi + 
       pow(F20,2)*pow(F22,2)*lfpf*lipi + 
       2*pow(F10,2)*pow(F22,2)*lfpi*lipi - 
       pow(F20,2)*pow(F22,2)*lfpi*lipi + 
       F22*kflf*((4*pow(F10,2)*F22 + 2*F10*F20*(F12 + F22) + 
       pow(F20,2)*(2*F12 + F22))*lipf + 
       (2*pow(F10,2) - pow(F20,2))*F22*lipi) - 
       8*pow(F10,2)*pow(F12,2)*m2*M2 - 
       8*F10*pow(F12,2)*F20*m2*M2 - 
       8*pow(F10,2)*F12*F22*m2*M2 - 
       4*F10*F12*F20*F22*m2*M2 + 
       4*F12*pow(F20,2)*F22*m2*M2 + 
       2*pow(F10,2)*pow(F22,2)*m2*M2 + 
       4*F10*F20*pow(F22,2)*m2*M2 + 
       3*pow(F20,2)*pow(F22,2)*m2*M2 + 
       6*pow(F10,2)*pow(F22,2)*m2*pipf - 
       pow(F20,2)*pow(F22,2)*m2*pipf + 
       lilf*((4*pow(F12,2)*pow(F20,2) + 
       4*F12*F20*(F10 + F20)*F22 + 
       (-2*pow(F10,2) + 4*F10*F20 + pow(F20,2))*pow(F22,2))*M2 + 
       (-6*pow(F10,2) + pow(F20,2))*pow(F22,2)*pipf))) + 
       kfpf*kfpi*M2*((2*pow(F10,2) - pow(F20,2))*pow(F22,2)*
       pow(kfpf,3)*(-lilf + m2) - 
       4*pow(F10,2)*pipf*(-4*F12*F22*lfpf*lipf*M2 - 
       3*pow(F22,2)*lfpf*lipf*M2 + 
       4*pow(F12,2)*lfpi*lipf*M2 + 
       4*F12*F22*lfpi*lipf*M2 + 
       pow(F22,2)*lfpi*lipf*M2 + 
       4*pow(F12,2)*lfpf*lipi*M2 + 
       4*F12*F22*lfpf*lipi*M2 + 
       pow(F22,2)*lfpf*lipi*M2 - 
       4*F12*F22*lfpi*lipi*M2 - 
       3*pow(F22,2)*lfpi*lipi*M2 + 
       2*pow(F12 + F22,2)*kflf*(-lipf + lipi)*M2 + 
       8*pow(F12,2)*m2*M4 + 
       12*F12*F22*m2*M4 + 
       5*pow(F22,2)*m2*M4 + 
       pow(F22,2)*lfpf*lipf*pipf + pow(F22,2)*lfpi*lipf*pipf + 
       pow(F22,2)*lfpf*lipi*pipf + pow(F22,2)*lfpi*lipi*pipf - 
       4*pow(F12,2)*m2*M2*pipf - 
       12*F12*F22*m2*M2*pipf - 
       6*pow(F22,2)*m2*M2*pipf + 
       pow(F22,2)*m2*pow(pipf,2) - 
       lilf*pow((2*F12 + F22)*M2 - F22*pipf,2)) - 
       pow(kfpf,2)*(2*pow(F10,2)*pow(F22,2)*lfpf*lipf - 
       pow(F20,2)*pow(F22,2)*lfpf*lipf + 
       2*F10*F12*F20*F22*lfpi*lipf + 
       2*F12*pow(F20,2)*F22*lfpi*lipf + 
       4*pow(F10,2)*pow(F22,2)*lfpi*lipf + 
       2*F10*F20*pow(F22,2)*lfpi*lipf + 
       pow(F20,2)*pow(F22,2)*lfpi*lipf + 
       2*F10*F12*F20*F22*lfpf*lipi + 
       2*F12*pow(F20,2)*F22*lfpf*lipi + 
       4*pow(F10,2)*pow(F22,2)*lfpf*lipi + 
       2*F10*F20*pow(F22,2)*lfpf*lipi + 
       pow(F20,2)*pow(F22,2)*lfpf*lipi - 
       8*pow(F12,2)*pow(F20,2)*lfpi*lipi - 
       12*F10*F12*F20*F22*lfpi*lipi - 
       12*F12*pow(F20,2)*F22*lfpi*lipi + 
       6*pow(F10,2)*pow(F22,2)*lfpi*lipi - 
       12*F10*F20*pow(F22,2)*lfpi*lipi - 
       5*pow(F20,2)*pow(F22,2)*lfpi*lipi - 
       F22*kflf*((2*pow(F10,2) - pow(F20,2))*F22*lipf + 
       (4*pow(F10,2)*F22 + 2*F10*F20*(F12 + F22) + 
       pow(F20,2)*(2*F12 + F22))*lipi) - 
       8*pow(F10,2)*pow(F12,2)*m2*M2 - 
       8*F10*pow(F12,2)*F20*m2*M2 - 
       8*pow(F10,2)*F12*F22*m2*M2 - 
       4*F10*F12*F20*F22*m2*M2 + 
       4*F12*pow(F20,2)*F22*m2*M2 + 
       2*pow(F10,2)*pow(F22,2)*m2*M2 + 
       4*F10*F20*pow(F22,2)*m2*M2 + 
       3*pow(F20,2)*pow(F22,2)*m2*M2 + 
       6*pow(F10,2)*pow(F22,2)*m2*pipf - 
       pow(F20,2)*pow(F22,2)*m2*pipf + 
       lilf*((4*pow(F12,2)*pow(F20,2) + 
       4*F12*F20*(F10 + F20)*F22 + 
       (-2*pow(F10,2) + 4*F10*F20 + pow(F20,2))*
       pow(F22,2))*M2 + 
       (-6*pow(F10,2) + pow(F20,2))*pow(F22,2)*pipf)) + 
       kfli*(F22*pow(kfpf,2)*
       ((2*pow(F10,2) - pow(F20,2))*F22*lfpf + 
       (4*pow(F10,2)*F22 + 2*F10*F20*(F12 + F22) + 
       pow(F20,2)*(2*F12 + F22))*lfpi) - 
       8*pow(F10,2)*pow(F12 + F22,2)*(-lfpf + lfpi)*M2*
       pipf + 2*kfpf*(lfpi*((2*pow(F10,2)*F22*(2*F12 + F22) + 
       pow(F20,2)*(2*pow(F12,2) + 3*F12*F22 + pow(F22,2)) + 
       F10*F20*(4*pow(F12,2) + 7*F12*F22 + 3*pow(F22,2)))*M2 - 
       (4*pow(F12,2)*pow(F20,2) + 
       7*F12*F20*(F10 + F20)*F22 + 
       (3*pow(F10,2) + 7*F10*F20 + 3*pow(F20,2))*
       pow(F22,2))*pipf) + 
       lfpf*(-((F12 + F22)*(4*pow(F10,2)*F12 + F10*F20*(4*F12 - 3*F22) - 
       pow(F20,2)*(2*F12 + F22))*M2) + 
       F22*(-(pow(F10,2)*F22) + F10*F20*(F12 + F22) + 
       pow(F20,2)*(F12 + F22))*pipf)) + 
       kflf*((-2*pow(F10,2) + pow(F20,2))*pow(F22,2)*pow(kfpf,2) - 
       2*F22*kfpf*((pow(F10,2)*F22 - 2*F10*F20*(F12 + F22) - 
       2*pow(F20,2)*(F12 + F22))*M2 + 
       2*F20*(F10 + F20)*(F12 + F22)*pipf) + 
       4*(4*F10*F12*(F10 + F20)*(F12 + F22)*M4 + 
       (pow(F10,2)*pow(F22,2) - 
       4*F10*F20*pow(F12 + F22,2) - 
       2*pow(F20,2)*pow(F12 + F22,2))*M2*pipf + 
       (2*pow(F12,2)*pow(F20,2) + 
       4*F12*F20*(F10 + F20)*F22 + 
       (pow(F10,2) + 4*F10*F20 + 2*pow(F20,2))*
       pow(F22,2))*pow(pipf,2)))) - 
       2*kfpf*(kflf*(-(lipi*((2*pow(F10,2)*F22*(2*F12 + F22) + 
       pow(F20,2)*(2*pow(F12,2) + 3*F12*F22 + pow(F22,2)) + 
       F10*F20*(4*pow(F12,2) + 7*F12*F22 + 3*pow(F22,2)))*M2 - 
       (4*pow(F12,2)*pow(F20,2) + 
       7*F12*F20*(F10 + F20)*F22 + 
       (3*pow(F10,2) + 7*F10*F20 + 3*pow(F20,2))*
       pow(F22,2))*pipf)) + 
       lipf*((F12 + F22)*(4*pow(F10,2)*F12 + F10*F20*(4*F12 - 3*F22) - 
       pow(F20,2)*(2*F12 + F22))*M2 + 
       F22*(pow(F10,2)*F22 - F10*F20*(F12 + F22) - 
       pow(F20,2)*(F12 + F22))*pipf)) + 
       pow(F10,2)*(-(pow(F22,2)*lfpf*lipf*M2) - 
       4*pow(F12,2)*lfpi*lipf*M2 + 
       pow(F22,2)*lfpi*lipf*M2 + 
       4*pow(F12,2)*m2*M4 + 
       12*F12*F22*m2*M4 + 
       6*pow(F22,2)*m2*M4 - 
       2*pow(F22,2)*lfpf*lipf*pipf - 
       4*pow(F22,2)*lfpi*lipf*pipf + 
       8*pow(F12,2)*m2*M2*pipf + 
       24*F12*F22*m2*M2*pipf + 
       10*pow(F22,2)*m2*M2*pipf - 
       4*pow(F22,2)*m2*pow(pipf,2) - 
       2*F22*lilf*(M2 + 2*pipf)*
       ((2*F12 + F22)*M2 - F22*pipf) + 
       lipi*(lfpi*((-8*pow(F12,2) + 3*pow(F22,2))*
       M2 - 6*pow(F22,2)*pipf) + 
       lfpf*((-4*pow(F12,2) + pow(F22,2))*M2 - 
       4*pow(F22,2)*pipf))))) + 
       pow(kfpi,2)*(F22*pow(kfpf,3)*
       (2*pow(F20,2)*F22*lfpi*lipf - 2*pow(F20,2)*F22*kflf*lipi + 
       2*pow(F20,2)*F22*lfpf*lipi + 4*pow(F20,2)*F22*lfpi*lipi - 
       8*F10*F12*F20*m2*M2 - 
       8*F12*pow(F20,2)*m2*M2 + 
       2*pow(F10,2)*F22*m2*M2 - 
       8*F10*F20*F22*m2*M2 - 
       7*pow(F20,2)*F22*m2*M2 + 
       2*pow(F20,2)*F22*m2*pipf - 
       F22*lilf*((2*pow(F10,2) + pow(F20,2))*M2 + 
       2*pow(F20,2)*pipf)) - 
       2*pow(kfpf,2)*(-2*F10*F12*F20*F22*lfpf*lipf*M2 - 
       2*F12*pow(F20,2)*F22*lfpf*lipf*M2 - 
       4*pow(F10,2)*pow(F22,2)*lfpf*lipf*M2 - 
       2*F10*F20*pow(F22,2)*lfpf*lipf*M2 - 
       2*pow(F20,2)*pow(F22,2)*lfpf*lipf*M2 + 
       4*pow(F12,2)*pow(F20,2)*lfpi*lipf*M2 + 
       6*F10*F12*F20*F22*lfpi*lipf*M2 + 
       6*F12*pow(F20,2)*F22*lfpi*lipf*M2 - 
       4*pow(F10,2)*pow(F22,2)*lfpi*lipf*M2 + 
       6*F10*F20*pow(F22,2)*lfpi*lipf*M2 + 
       2*pow(F20,2)*pow(F22,2)*lfpi*lipf*M2 + 
       4*pow(F12,2)*pow(F20,2)*lfpf*lipi*M2 + 
       6*F10*F12*F20*F22*lfpf*lipi*M2 + 
       6*F12*pow(F20,2)*F22*lfpf*lipi*M2 - 
       4*pow(F10,2)*pow(F22,2)*lfpf*lipi*M2 + 
       6*F10*F20*pow(F22,2)*lfpf*lipi*M2 + 
       2*pow(F20,2)*pow(F22,2)*lfpf*lipi*M2 - 
       2*F10*F12*F20*F22*lfpi*lipi*M2 - 
       2*F12*pow(F20,2)*F22*lfpi*lipi*M2 - 
       4*pow(F10,2)*pow(F22,2)*lfpi*lipi*M2 - 
       2*F10*F20*pow(F22,2)*lfpi*lipi*M2 - 
       2*pow(F20,2)*pow(F22,2)*lfpi*lipi*M2 - 
       F22*(pow(F10,2)*F22 - 3*F10*F20*(F12 + F22) - 
       3*pow(F20,2)*(F12 + F22))*kflf*(-lipf + lipi)*M2 +
       8*F10*pow(F12,2)*F20*m2*M4 + 
       4*pow(F12,2)*pow(F20,2)*m2*M4 + 
       8*pow(F10,2)*F12*F22*m2*M4 + 
       16*F10*F12*F20*F22*m2*M4 + 
       8*F12*pow(F20,2)*F22*m2*M4 + 
       8*F10*F20*pow(F22,2)*m2*M4 + 
       4*pow(F20,2)*pow(F22,2)*m2*M4 + 
       pow(F20,2)*pow(F22,2)*lfpf*lipf*pipf + 
       pow(F20,2)*pow(F22,2)*lfpi*lipf*pipf + 
       pow(F20,2)*pow(F22,2)*lfpf*lipi*pipf + 
       pow(F20,2)*pow(F22,2)*lfpi*lipi*pipf - 
       4*pow(F12,2)*pow(F20,2)*m2*M2*pipf - 
       12*F10*F12*F20*F22*m2*M2*pipf - 
       12*F12*pow(F20,2)*F22*m2*M2*pipf - 
       6*pow(F10,2)*pow(F22,2)*m2*M2*pipf - 
       12*F10*F20*pow(F22,2)*m2*M2*pipf - 
       7*pow(F20,2)*pow(F22,2)*m2*M2*pipf + 
       pow(F20,2)*pow(F22,2)*m2*pow(pipf,2) + 
       lilf*(4*pow(F10,2)*pow(F22,2)*M4 - 
       (4*pow(F12,2)*pow(F20,2) + 
       4*F12*F20*(F10 + F20)*F22 + 
       (-2*pow(F10,2) + 4*F10*F20 + pow(F20,2))*
       pow(F22,2))*M2*pipf - 
       pow(F20,2)*pow(F22,2)*pow(pipf,2))) + 
       2*pow(F10,2)*M4*
       (-4*F12*F22*lfpf*lipf*M2 - 
       3*pow(F22,2)*lfpf*lipf*M2 + 
       4*pow(F12,2)*lfpi*lipf*M2 + 
       4*F12*F22*lfpi*lipf*M2 + 
       pow(F22,2)*lfpi*lipf*M2 + 
       4*pow(F12,2)*lfpf*lipi*M2 + 
       4*F12*F22*lfpf*lipi*M2 + 
       pow(F22,2)*lfpf*lipi*M2 - 
       4*F12*F22*lfpi*lipi*M2 - 
       3*pow(F22,2)*lfpi*lipi*M2 + 
       8*pow(F12,2)*m2*M4 + 
       12*F12*F22*m2*M4 + 
       5*pow(F22,2)*m2*M4 + 
       pow(F22,2)*lfpf*lipf*pipf + pow(F22,2)*lfpi*lipf*pipf + 
       pow(F22,2)*lfpf*lipi*pipf + pow(F22,2)*lfpi*lipi*pipf - 
       4*pow(F12,2)*m2*M2*pipf - 
       12*F12*F22*m2*M2*pipf - 
       6*pow(F22,2)*m2*M2*pipf + 
       pow(F22,2)*m2*pow(pipf,2) - 
       lilf*pow((2*F12 + F22)*M2 - F22*pipf,2) + 
       kflf*(F22*lipf*(-((4*F12 + 3*F22)*M2) + F22*pipf) + 
       lipi*(pow(2*F12 + F22,2)*M2 + pow(F22,2)*pipf))) +
       2*kfli*(F22*kflf*(pow(F10,2)*M4*
       (-((4*F12 + 3*F22)*M2) + F22*pipf) + 
       pow(kfpf,2)*(-2*pow(F10,2)*F22*M2 + 
       pow(F20,2)*F22*pipf) + 
       kfpf*M2*((pow(F10,2)*F22 - 2*F10*F20*(F12 + F22) - 
       2*pow(F20,2)*(F12 + F22))*M2 + 
       2*F20*(F10 + F20)*(F12 + F22)*pipf)) + 
       lfpf*M2*(F22*(-(pow(F10,2)*F22) + 3*F10*F20*(F12 + F22) + 
       3*pow(F20,2)*(F12 + F22))*pow(kfpf,2) + 
       pow(F10,2)*F22*M2*
       (-((4*F12 + 3*F22)*M2) + F22*pipf) + 
       kfpf*((2*pow(F10,2)*F22*(2*F12 + F22) + 
       pow(F20,2)*(2*pow(F12,2) + 3*F12*F22 + pow(F22,2)) + 
       F10*F20*(4*pow(F12,2) + 7*F12*F22 + 3*pow(F22,2)))*M2 - 
       (4*pow(F12,2)*pow(F20,2) + 
       7*F12*F20*(F10 + F20)*F22 + 
       (3*pow(F10,2) + 7*F10*F20 + 3*pow(F20,2))*pow(F22,2))*pipf)) + 
       lfpi*(-(pow(F20,2)*pow(F22,2)*pow(kfpf,3)) + 
       F22*(pow(F10,2)*F22 - 3*F10*F20*(F12 + F22) - 
       3*pow(F20,2)*(F12 + F22))*pow(kfpf,2)*M2 + 
       pow(F10,2)*M4*
       (pow(2*F12 + F22,2)*M2 + pow(F22,2)*pipf) + 
       kfpf*M2*(-((F12 + F22)*(4*pow(F10,2)*F12 + F10*F20*(4*F12 - 3*F22) - 
       pow(F20,2)*(2*F12 + F22))*M2) + 
       F22*(-(pow(F10,2)*F22) + F10*F20*(F12 + F22) + 
       pow(F20,2)*(F12 + F22))*pipf))) + 
       2*kfpf*M2*(kflf*
       (lipf*((2*pow(F10,2)*F22*(2*F12 + F22) + 
       pow(F20,2)*(2*pow(F12,2) + 3*F12*F22 + pow(F22,2)) + 
       F10*F20*(4*pow(F12,2) + 7*F12*F22 + 3*pow(F22,2)))*
       M2 - (4*pow(F12,2)*pow(F20,2) + 
       7*F12*F20*(F10 + F20)*F22 + 
       (3*pow(F10,2) + 7*F10*F20 + 3*pow(F20,2))*pow(F22,2))*pipf) + 
       lipi*(-((F12 + F22)*(4*pow(F10,2)*F12 + F10*F20*(4*F12 - 3*F22) - 
       pow(F20,2)*(2*F12 + F22))*M2) + 
       F22*(-(pow(F10,2)*F22) + F10*F20*(F12 + F22) + 
       pow(F20,2)*(F12 + F22))*pipf)) + 
       pow(F10,2)*(-8*pow(F12,2)*lfpf*lipf*M2 + 
       3*pow(F22,2)*lfpf*lipf*M2 - 
       4*pow(F12,2)*lfpi*lipf*M2 + 
       pow(F22,2)*lfpi*lipf*M2 + 
       4*pow(F12,2)*m2*M4 + 
       12*F12*F22*m2*M4 + 
       6*pow(F22,2)*m2*M4 - 
       6*pow(F22,2)*lfpf*lipf*pipf - 
       4*pow(F22,2)*lfpi*lipf*pipf + 
       8*pow(F12,2)*m2*M2*pipf + 
       24*F12*F22*m2*M2*pipf + 
       10*pow(F22,2)*m2*M2*pipf - 
       4*pow(F22,2)*m2*pow(pipf,2) - 
       2*F22*lilf*(M2 + 2*pipf)*
       ((2*F12 + F22)*M2 - F22*pipf) - 
       lipi*(pow(F22,2)*lfpi*(M2 + 2*pipf) + 
       lfpf*((4*pow(F12,2) - pow(F22,2))*M2 + 
       4*pow(F22,2)*pipf))))))/(pow(kfpf,2)*pow(kfpi,2)*M4));

return (Pow6(e)/Pow2(q_22))*res;
}

