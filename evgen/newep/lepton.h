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
// First-order bremsstrahlung, lepton term:
double lterm ()
{
res = (2*(pow(kfli,3)*(-kflf + m2)*pow((2*F11 + F21)*M2 - F21*pipf,2) - 
       pow(kflf,2)*m2*(-4*pow(F11,2)*kfpf*lfpi*M2 - 
       4*F11*F21*kfpf*lfpi*M2 - 
       pow(F21,2)*kfpf*lfpi*M2 + 
       4*F11*F21*kfpi*lfpi*M2 + 
       3*pow(F21,2)*kfpi*lfpi*M2 + 
       4*pow(F11,2)*lfpi*lipf*M2 + 
       4*F11*F21*lfpi*lipf*M2 + 
       pow(F21,2)*lfpi*lipf*M2 - 
       4*F11*F21*lfpi*lipi*M2 - 
       3*pow(F21,2)*lfpi*lipi*M2 + 
       (8*pow(F11,2) + 12*F11*F21 + 5*pow(F21,2))*m2*
       M4 - pow(F21,2)*kfpf*lfpi*pipf - 
       pow(F21,2)*kfpi*lfpi*pipf + pow(F21,2)*lfpi*lipf*pipf + 
       pow(F21,2)*lfpi*lipi*pipf - 
       4*pow(F11,2)*m2*M2*pipf - 
       12*F11*F21*m2*M2*pipf - 
       6*pow(F21,2)*m2*M2*pipf + 
       pow(F21,2)*m2*pow(pipf,2) + 
       kflf*pow((2*F11 + F21)*M2 - F21*pipf,2) - 
       lilf*pow((2*F11 + F21)*M2 - F21*pipf,2) + 
       lfpf*((F21*(4*F11 + 3*F21)*kfpf - pow(2*F11 + F21,2)*kfpi - 
       F21*(4*F11 + 3*F21)*lipf + pow(2*F11 + F21,2)*lipi)*
       M2 + pow(F21,2)*(-kfpf - kfpi + lipf + lipi)*pipf))
       + pow(kfli,2)*(-(m2*
       (4*pow(F11,2)*lfpi*lipf*M2 + 
       4*F11*F21*lfpi*lipf*M2 + 
       pow(F21,2)*lfpi*lipf*M2 - 
       4*F11*F21*lfpi*lipi*M2 - 
       3*pow(F21,2)*lfpi*lipi*M2 - 
       4*pow(F11,2)*lilf*M4 - 
       4*F11*F21*lilf*M4 - pow(F21,2)*lilf*M4 + 
       (8*pow(F11,2) + 12*F11*F21 + 5*pow(F21,2))*m2*
       M4 + pow(F21,2)*lfpi*lipf*pipf + 
       pow(F21,2)*lfpi*lipi*pipf + 
       4*F11*F21*lilf*M2*pipf + 
       2*pow(F21,2)*lilf*M2*pipf - 
       4*pow(F11,2)*m2*M2*pipf - 
       12*F11*F21*m2*M2*pipf - 
       6*pow(F21,2)*m2*M2*pipf - 
       pow(F21,2)*lilf*pow(pipf,2) + 
       pow(F21,2)*m2*pow(pipf,2) + 
       kfpf*((-(F21*(4*F11 + 3*F21)*lipf) + 
       pow(2*F11 + F21,2)*lipi)*M2 + 
       pow(F21,2)*(lipf + lipi)*pipf) + 
       lfpf*((-(F21*(4*F11 + 3*F21)*lipf) + 
       pow(2*F11 + F21,2)*lipi)*M2 + 
       pow(F21,2)*(lipf + lipi)*pipf) + 
       kfpi*((pow(2*F11 + F21,2)*lipf - 
       F21*(4*F11 + 3*F21)*lipi)*M2 + 
       pow(F21,2)*(lipf + lipi)*pipf))) + 
       kflf*(-4*F11*F21*pow(lfpi,2)*M2 - 
       3*pow(F21,2)*pow(lfpi,2)*M2 + 
       4*pow(F11,2)*lfpi*lipf*M2 + 
       4*F11*F21*lfpi*lipf*M2 + 
       pow(F21,2)*lfpi*lipf*M2 - 
       4*F11*F21*lfpi*lipi*M2 - 
       3*pow(F21,2)*lfpi*lipi*M2 - 
       8*pow(F11,2)*lilf*M4 - 8*F11*F21*lilf*M4 - 
       2*pow(F21,2)*lilf*M4 - 
       pow(2*F11 + F21,2)*m2*M4 + 
       pow(F21,2)*pow(lfpi,2)*pipf + 
       pow(F21,2)*lfpi*lipf*pipf + pow(F21,2)*lfpi*lipi*pipf + 
       8*F11*F21*lilf*M2*pipf + 
       4*pow(F21,2)*lilf*M2*pipf + 
       4*F11*F21*m2*M2*pipf + 
       2*pow(F21,2)*m2*M2*pipf - 
       2*pow(F21,2)*lilf*pow(pipf,2) - 
       pow(F21,2)*m2*pow(pipf,2) + 
       F21*pow(lfpf,2)*(-((4*F11 + 3*F21)*M2) + F21*pipf) + 
       kfpf*((-(F21*(4*F11 + 3*F21)*lipf) + 
       pow(2*F11 + F21,2)*lipi)*M2 + 
       pow(F21,2)*(lipf + lipi)*pipf) + 
       kfpi*((pow(2*F11 + F21,2)*lipf - F21*(4*F11 + 3*F21)*lipi)*
       M2 + pow(F21,2)*(lipf + lipi)*pipf) + 
       lfpf*((2*pow(2*F11 + F21,2)*lfpi - 
       F21*(4*F11 + 3*F21)*lipf + pow(2*F11 + F21,2)*lipi)*
       M2 + pow(F21,2)*(2*lfpi + lipf + lipi)*pipf))) + 
       kflf*kfli*(-4*pow(F11,2)*kfpf*lfpi*lilf*M2 - 
       4*F11*F21*kfpf*lfpi*lilf*M2 - 
       pow(F21,2)*kfpf*lfpi*lilf*M2 - 
       4*F11*F21*kfpf*lilf*lipf*M2 - 
       3*pow(F21,2)*kfpf*lilf*lipf*M2 + 
       8*pow(F11,2)*lfpi*lilf*lipf*M2 + 
       8*F11*F21*lfpi*lilf*lipf*M2 + 
       2*pow(F21,2)*lfpi*lilf*lipf*M2 + 
       4*pow(F11,2)*kfpf*lilf*lipi*M2 + 
       4*F11*F21*kfpf*lilf*lipi*M2 + 
       pow(F21,2)*kfpf*lilf*lipi*M2 - 
       8*F11*F21*lfpi*lilf*lipi*M2 - 
       6*pow(F21,2)*lfpi*lilf*lipi*M2 + 
       4*F11*F21*pow(kfpf,2)*m2*M2 + 
       3*pow(F21,2)*pow(kfpf,2)*m2*M2 - 
       8*pow(F11,2)*pow(lilf,2)*M4 - 
       8*F11*F21*pow(lilf,2)*M4 - 
       2*pow(F21,2)*pow(lilf,2)*M4 + 
       16*pow(F11,2)*lilf*m2*M4 + 
       24*F11*F21*lilf*m2*M4 + 
       10*pow(F21,2)*lilf*m2*M4 - 
       pow(F21,2)*kfpf*lfpi*lilf*pipf + 
       pow(F21,2)*kfpf*lilf*lipf*pipf + 
       2*pow(F21,2)*lfpi*lilf*lipf*pipf + 
       pow(F21,2)*kfpf*lilf*lipi*pipf + 
       2*pow(F21,2)*lfpi*lilf*lipi*pipf - 
       pow(F21,2)*pow(kfpf,2)*m2*pipf + 
       8*F11*F21*pow(lilf,2)*M2*pipf + 
       4*pow(F21,2)*pow(lilf,2)*M2*pipf - 
       8*pow(F11,2)*lilf*m2*M2*pipf - 
       24*F11*F21*lilf*m2*M2*pipf - 
       12*pow(F21,2)*lilf*m2*M2*pipf - 
       2*pow(F21,2)*pow(lilf,2)*pow(pipf,2) + 
       2*pow(F21,2)*lilf*m2*pow(pipf,2) - 
       pow(kflf,2)*pow((2*F11 + F21)*M2 - F21*pipf,2) + 
       F21*pow(kfpi,2)*m2*((4*F11 + 3*F21)*M2 - F21*pipf) + 
       lfpf*lilf*((F21*(4*F11 + 3*F21)*kfpf - 
       2*F21*(4*F11 + 3*F21)*lipf + 2*pow(2*F11 + F21,2)*lipi)*
       M2 + pow(F21,2)*(-kfpf + 2*(lipf + lipi))*pipf) + 
       kfpi*(-4*F11*F21*kflf*lfpi*M2 - 
       3*pow(F21,2)*kflf*lfpi*M2 + 
       4*F11*F21*lfpi*lilf*M2 + 
       3*pow(F21,2)*lfpi*lilf*M2 + 
       4*pow(F11,2)*lilf*lipf*M2 + 
       4*F11*F21*lilf*lipf*M2 + 
       pow(F21,2)*lilf*lipf*M2 + 
       pow(F21,2)*kflf*lfpi*pipf - pow(F21,2)*lfpi*lilf*pipf + 
       pow(F21,2)*lilf*lipf*pipf + 
       F21*lilf*lipi*(-((4*F11 + 3*F21)*M2) + F21*pipf) + 
       lfpf*(kflf - lilf)*(pow(2*F11 + F21,2)*M2 + pow(F21,2)*pipf) - 
       2*kfpf*m2*(pow(2*F11 + F21,2)*M2 + 
       pow(F21,2)*pipf)) + 
       kflf*(4*pow(F11,2)*kfpf*lfpi*M2 + 
       4*F11*F21*kfpf*lfpi*M2 + 
       pow(F21,2)*kfpf*lfpi*M2 - 
       4*pow(F11,2)*lfpi*lipf*M2 - 
       4*F11*F21*lfpi*lipf*M2 - 
       pow(F21,2)*lfpi*lipf*M2 + 
       4*F11*F21*pow(lipf,2)*M2 + 
       3*pow(F21,2)*pow(lipf,2)*M2 + 
       pow(2*F11 + F21,2)*m2*M4 + 
       pow(F21,2)*kfpf*lfpi*pipf - pow(F21,2)*lfpi*lipf*pipf - 
       pow(F21,2)*pow(lipf,2)*pipf - 
       4*F11*F21*m2*M2*pipf - 
       2*pow(F21,2)*m2*M2*pipf + 
       pow(F21,2)*m2*pow(pipf,2) + 
       2*lilf*pow((2*F11 + F21)*M2 - F21*pipf,2) + 
       F21*lfpf*(-kfpf + lipf)*((4*F11 + 3*F21)*M2 - F21*pipf) + 
       F21*pow(lipi,2)*((4*F11 + 3*F21)*M2 - F21*pipf) + 
       lipi*(F21*lfpi*((4*F11 + 3*F21)*M2 - F21*pipf) + 
       2*lipf*(-(pow(2*F11 + F21,2)*M2) - 
       pow(F21,2)*pipf) - 
       lfpf*(pow(2*F11 + F21,2)*M2 + pow(F21,2)*pipf))))))/(pow(kflf,2)*pow(kfli,2)*M2);

return (Pow6(e)/Pow2(q_12))*res;
}

