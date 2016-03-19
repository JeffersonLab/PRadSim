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
// First-order bremsstrahlung, interference term:
double lpterm ()
{
res1 = (F20*F21*F22*kflf*pow(kfpf,2)*lilf + 
       2*F20*F21*F22*pow(kfpf,2)*lfpf*lilf - 
       4*F11*F20*F22*pow(kfpf,2)*lfpi*lilf + 
       2*F10*F21*F22*pow(kfpf,2)*lfpi*lilf - 
       2*F20*F21*F22*pow(kfpf,2)*lfpi*lilf - 
       F20*F21*F22*pow(kflf,2)*kfpf*lipf - 
       4*F20*F21*F22*kflf*pow(kfpf,2)*lipf - 
       4*F20*F21*F22*kflf*kfpf*lfpf*lipf + 2*F12*F20*F21*kflf*kfpf*lfpi*lipf - 
       2*F11*F20*F22*kflf*kfpf*lfpi*lipf - 2*F10*F21*F22*kflf*kfpf*lfpi*lipf - 
       4*F20*F21*F22*kflf*kfpf*lfpi*lipf + 
       4*F12*F20*F21*pow(kfpf,2)*lfpi*lipf - 
       4*F11*F20*F22*pow(kfpf,2)*lfpi*lipf - 
       2*F10*F21*F22*pow(kfpf,2)*lfpi*lipf + 
       8*F12*F20*F21*kfpf*lfpf*lfpi*lipf - 8*F11*F20*F22*kfpf*lfpf*lfpi*lipf + 
       4*F12*F20*F21*kfpf*pow(lfpi,2)*lipf - 
       4*F11*F20*F22*kfpf*pow(lfpi,2)*lipf + 
       2*F11*F20*F22*pow(kflf,2)*kfpf*lipi + 
       F20*F21*F22*pow(kflf,2)*kfpf*lipi + 
       12*F11*F20*F22*kflf*pow(kfpf,2)*lipi - 
       2*F10*F21*F22*kflf*pow(kfpf,2)*lipi + 
       8*F20*F21*F22*kflf*pow(kfpf,2)*lipi + 
       2*F12*F20*F21*kflf*kfpf*lfpf*lipi + 
       10*F11*F20*F22*kflf*kfpf*lfpf*lipi + 
       8*F20*F21*F22*kflf*kfpf*lfpf*lipi - 
       4*F10*F21*F22*pow(kfpf,2)*lfpf*lipi - 
       4*F10*F21*F22*kfpf*pow(lfpf,2)*lipi + 
       4*F12*F20*F21*kflf*kfpf*lfpi*lipi + 4*F11*F20*F22*kflf*kfpf*lfpi*lipi - 
       2*F10*F21*F22*kflf*kfpf*lfpi*lipi + 4*F20*F21*F22*kflf*kfpf*lfpi*lipi - 
       4*F11*F20*F22*pow(kfpf,2)*lfpi*lipi - 
       6*F10*F21*F22*pow(kfpf,2)*lfpi*lipi - 
       4*F20*F21*F22*pow(kfpf,2)*lfpi*lipi + 
       4*F12*F20*F21*kfpf*lfpf*lfpi*lipi - 4*F11*F20*F22*kfpf*lfpf*lfpi*lipi - 
       8*F10*F21*F22*kfpf*lfpf*lfpi*lipi - 
       4*F10*F21*F22*kfpf*pow(lfpi,2)*lipi - 
       F20*F21*F22*kflf*pow(kfpf,2)*m2 + 
       2*F20*F21*F22*pow(kfpf,3)*m2 + 
       4*F11*F20*F22*pow(kfpf,2)*lfpi*m2 - 
       2*F10*F21*F22*pow(kfpf,2)*lfpi*m2 + 
       4*F20*F21*F22*pow(kfpf,2)*lfpi*m2 + 
       2*F20*F21*F22*pow(kfpf,2)*lipf*m2 - 
       8*F11*F20*F22*pow(kfpf,2)*lipi*m2 - 
       6*F20*F21*F22*pow(kfpf,2)*lipi*m2 + 
       8*F11*F12*F20*kflf*kfpf*lilf*M2 + 
       4*F12*F20*F21*kflf*kfpf*lilf*M2 + 
       4*F10*F11*F22*kflf*kfpf*lilf*M2 + 
       4*F11*F20*F22*kflf*kfpf*lilf*M2 + 
       2*F10*F21*F22*kflf*kfpf*lilf*M2 + 
       2*F20*F21*F22*kflf*kfpf*lilf*M2 - 
       8*F11*F12*F20*pow(kfpf,2)*lilf*M2 - 
       6*F12*F20*F21*pow(kfpf,2)*lilf*M2 - 
       12*F10*F11*F22*pow(kfpf,2)*lilf*M2 - 
       6*F11*F20*F22*pow(kfpf,2)*lilf*M2 - 
       8*F10*F21*F22*pow(kfpf,2)*lilf*M2 - 
       4*F20*F21*F22*pow(kfpf,2)*lilf*M2 + 
       4*F10*F12*F21*kfpf*lfpf*lilf*M2 - 
       4*F10*F11*F22*kfpf*lfpf*lilf*M2 - 
       16*F10*F11*F12*kfpf*lfpi*lilf*M2 - 
       12*F10*F12*F21*kfpf*lfpi*lilf*M2 - 
       4*F10*F11*F22*kfpf*lfpi*lilf*M2 + 
       2*F12*F20*F21*pow(kflf,2)*lipf*M2 - 
       4*F10*F11*F22*pow(kflf,2)*lipf*M2 + 
       2*F11*F20*F22*pow(kflf,2)*lipf*M2 - 
       2*F10*F21*F22*pow(kflf,2)*lipf*M2 + 
       4*F20*F21*F22*pow(kflf,2)*lipf*M2 - 
       8*F11*F12*F20*kflf*kfpf*lipf*M2 - 
       8*F10*F12*F21*kflf*kfpf*lipf*M2 - 
       10*F12*F20*F21*kflf*kfpf*lipf*M2 - 
       4*F10*F11*F22*kflf*kfpf*lipf*M2 - 
       6*F11*F20*F22*kflf*kfpf*lipf*M2 - 
       6*F10*F21*F22*kflf*kfpf*lipf*M2 - 
       8*F20*F21*F22*kflf*kfpf*lipf*M2 - 
       12*F10*F12*F21*kflf*lfpf*lipf*M2 - 
       20*F10*F11*F22*kflf*lfpf*lipf*M2 - 
       24*F10*F21*F22*kflf*lfpf*lipf*M2 - 
       16*F10*F12*F21*kfpf*lfpf*lipf*M2 - 
       16*F10*F11*F22*kfpf*lfpf*lipf*M2 - 
       24*F10*F21*F22*kfpf*lfpf*lipf*M2 - 
       16*F10*F12*F21*pow(lfpf,2)*lipf*M2 - 
       16*F10*F11*F22*pow(lfpf,2)*lipf*M2 - 
       24*F10*F21*F22*pow(lfpf,2)*lipf*M2 + 
       16*F10*F11*F12*kflf*lfpi*lipf*M2 + 
       4*F10*F12*F21*kflf*lfpi*lipf*M2 - 
       4*F12*F20*F21*kflf*lfpi*lipf*M2 + 
       12*F10*F11*F22*kflf*lfpi*lipf*M2 + 
       4*F11*F20*F22*kflf*lfpi*lipf*M2 + 
       4*F10*F21*F22*kflf*lfpi*lipf*M2 + 
       16*F10*F11*F12*kfpf*lfpi*lipf*M2 + 
       8*F10*F12*F21*kfpf*lfpi*lipf*M2 + 
       8*F10*F11*F22*kfpf*lfpi*lipf*M2 + 
       4*F10*F21*F22*kfpf*lfpi*lipf*M2 + 
       32*F10*F11*F12*lfpf*lfpi*lipf*M2 + 
       16*F10*F12*F21*lfpf*lfpi*lipf*M2 + 
       16*F10*F11*F22*lfpf*lfpi*lipf*M2 + 
       8*F10*F21*F22*lfpf*lfpi*lipf*M2 - 
       8*F11*F12*F20*pow(kflf,2)*lipi*M2 - 
       6*F12*F20*F21*pow(kflf,2)*lipi*M2 - 
       4*F10*F11*F22*pow(kflf,2)*lipi*M2 - 
       6*F11*F20*F22*pow(kflf,2)*lipi*M2 - 
       2*F10*F21*F22*pow(kflf,2)*lipi*M2 - 
       4*F20*F21*F22*pow(kflf,2)*lipi*M2 + 
       32*F10*F11*F12*kflf*kfpf*lipi*M2 + 
       8*F11*F12*F20*kflf*kfpf*lipi*M2 + 
       24*F10*F12*F21*kflf*kfpf*lipi*M2 + 
       6*F12*F20*F21*kflf*kfpf*lipi*M2 + 
       4*F10*F11*F22*kflf*kfpf*lipi*M2 + 
       6*F11*F20*F22*kflf*kfpf*lipi*M2 + 
       4*F10*F21*F22*kflf*kfpf*lipi*M2 + 
       4*F20*F21*F22*kflf*kfpf*lipi*M2 + 
       32*F10*F11*F12*kflf*lfpf*lipi*M2 + 
       20*F10*F12*F21*kflf*lfpf*lipi*M2 + 
       12*F10*F11*F22*kflf*lfpf*lipi*M2 + 
       8*F10*F21*F22*kflf*lfpf*lipi*M2 + 
       32*F10*F11*F12*kfpf*lfpf*lipi*M2 + 
       16*F10*F12*F21*kfpf*lfpf*lipi*M2 + 
       16*F10*F11*F22*kfpf*lfpf*lipi*M2 + 
       8*F10*F21*F22*kfpf*lfpf*lipi*M2 + 
       32*F10*F11*F12*pow(lfpf,2)*lipi*M2 + 
       16*F10*F12*F21*pow(lfpf,2)*lipi*M2 + 
       16*F10*F11*F22*pow(lfpf,2)*lipi*M2 + 
       8*F10*F21*F22*pow(lfpf,2)*lipi*M2 - 
       12*F10*F12*F21*kflf*lfpi*lipi*M2 - 
       4*F12*F20*F21*kflf*lfpi*lipi*M2 - 
       4*F10*F11*F22*kflf*lfpi*lipi*M2 + 
       4*F11*F20*F22*kflf*lfpi*lipi*M2 - 
       12*F10*F21*F22*kflf*lfpi*lipi*M2 - 
       8*F10*F12*F21*kfpf*lfpi*lipi*M2 - 
       8*F10*F11*F22*kfpf*lfpi*lipi*M2 - 
       12*F10*F21*F22*kfpf*lfpi*lipi*M2 - 
       16*F10*F12*F21*lfpf*lfpi*lipi*M2 - 
       16*F10*F11*F22*lfpf*lfpi*lipi*M2 - 
       24*F10*F21*F22*lfpf*lfpi*lipi*M2 - 
       16*F11*F12*F20*kflf*kfpf*m2*M2 - 
       12*F12*F20*F21*kflf*kfpf*m2*M2 - 
       12*F11*F20*F22*kflf*kfpf*m2*M2 + 
       2*F10*F21*F22*kflf*kfpf*m2*M2 - 
       10*F20*F21*F22*kflf*kfpf*m2*M2 + 
       16*F11*F12*F20*pow(kfpf,2)*m2*M2 + 
       4*F10*F12*F21*pow(kfpf,2)*m2*M2 + 
       14*F12*F20*F21*pow(kfpf,2)*m2*M2 + 
       28*F10*F11*F22*pow(kfpf,2)*m2*M2 + 
       14*F11*F20*F22*pow(kfpf,2)*m2*M2 + 
       24*F10*F21*F22*pow(kfpf,2)*m2*M2 + 
       12*F20*F21*F22*pow(kfpf,2)*m2*M2 + 
       24*F10*F11*F22*kfpf*lfpf*m2*M2 + 
       20*F10*F21*F22*kfpf*lfpf*m2*M2 + 
       16*F10*F11*F12*kfpf*lfpi*m2*M2 + 
       16*F10*F12*F21*kfpf*lfpi*m2*M2 + 
       8*F10*F11*F22*kfpf*lfpi*m2*M2 + 
       4*F10*F21*F22*kfpf*lfpi*m2*M2 + 
       4*F10*F12*F21*kfpf*lipf*m2*M2 + 
       4*F10*F11*F22*kfpf*lipf*m2*M2 + 
       4*F10*F21*F22*kfpf*lipf*m2*M2 - 
       16*F10*F11*F12*kfpf*lipi*m2*M2 - 
       12*F10*F12*F21*kfpf*lipi*m2*M2 + 
       4*F10*F11*F22*kfpf*lipi*m2*M2 + 
       4*F10*F21*F22*kfpf*lipi*m2*M2 - 
       16*F10*F11*F12*kflf*lilf*M4 - 
       8*F10*F12*F21*kflf*lilf*M4 - 
       8*F10*F11*F22*kflf*lilf*M4 - 
       4*F10*F21*F22*kflf*lilf*M4 - 
       16*F10*F11*F12*kfpf*lilf*M4 - 
       8*F10*F12*F21*kfpf*lilf*M4 - 
       8*F10*F11*F22*kfpf*lilf*M4 - 
       4*F10*F21*F22*kfpf*lilf*M4 - 
       32*F10*F11*F12*lfpf*lilf*M4 - 
       16*F10*F12*F21*lfpf*lilf*M4 - 
       16*F10*F11*F22*lfpf*lilf*M4 - 
       8*F10*F21*F22*lfpf*lilf*M4 + 
       16*F10*F11*F12*kflf*lipf*M4 + 
       16*F10*F12*F21*kflf*lipf*M4 + 
       16*F10*F11*F22*kflf*lipf*M4 + 
       16*F10*F21*F22*kflf*lipf*M4 - 
       16*F10*F11*F12*kflf*lipi*M4 - 
       8*F10*F12*F21*kflf*lipi*M4 - 
       8*F10*F11*F22*kflf*lipi*M4 - 
       4*F10*F21*F22*kflf*lipi*M4 + 
       32*F10*F11*F12*kflf*m2*M4 + 
       24*F10*F12*F21*kflf*m2*M4 + 
       24*F10*F11*F22*kflf*m2*M4 + 
       20*F10*F21*F22*kflf*m2*M4 + 
       32*F10*F11*F12*kfpf*m2*M4 + 
       24*F10*F12*F21*kfpf*m2*M4 + 
       24*F10*F11*F22*kfpf*m2*M4 + 
       20*F10*F21*F22*kfpf*m2*M4 + 
       64*F10*F11*F12*lfpf*m2*M4 + 
       48*F10*F12*F21*lfpf*m2*M4 + 
       48*F10*F11*F22*lfpf*m2*M4 + 
       40*F10*F21*F22*lfpf*m2*M4 + 
       pow(kfpi,2)*(F21*F22*kflf*(F20*lilf + 2*F10*lipf - F20*m2) + 
       2*(-2*F12*F20*F21*lfpf*lipf + F10*F21*F22*lfpf*lipf - 
       2*F20*F21*F22*lfpf*lipf + F10*F21*F22*lfpf*lipi + 
       F20*F21*F22*kfpf*m2 + F10*F21*F22*lfpf*m2 + 
       2*F10*F12*F21*m2*M2 + 
       F12*F20*F21*m2*M2 + 
       2*F10*F11*F22*m2*M2 + 
       F11*F20*F22*m2*M2 + 
       4*F10*F21*F22*m2*M2 + 
       2*F20*F21*F22*m2*M2 + 
       lilf*(-(F10*F21*F22*lfpf) + 
       (F12*F20*F21 - (F11*F20 + 2*F10*(F11 + F21))*F22)*M2))) -
       4*F12*F20*F21*kflf*kfpf*lilf*pipf - 2*F10*F21*F22*kflf*kfpf*lilf*pipf - 
       2*F20*F21*F22*kflf*kfpf*lilf*pipf + 
       4*F11*F20*F22*pow(kfpf,2)*lilf*pipf + 
       4*F10*F21*F22*pow(kfpf,2)*lilf*pipf + 
       4*F20*F21*F22*pow(kfpf,2)*lilf*pipf - 
       4*F12*F20*F21*kfpf*lfpf*lilf*pipf + 4*F11*F20*F22*kfpf*lfpf*lilf*pipf - 
       4*F12*F20*F21*kfpf*lfpi*lilf*pipf + 4*F11*F20*F22*kfpf*lfpi*lilf*pipf + 
       2*F12*F20*F21*pow(kflf,2)*lipf*pipf - 
       2*F11*F20*F22*pow(kflf,2)*lipf*pipf + 
       2*F10*F21*F22*pow(kflf,2)*lipf*pipf + 
       4*F12*F20*F21*kflf*kfpf*lipf*pipf + 6*F10*F21*F22*kflf*kfpf*lipf*pipf + 
       4*F20*F21*F22*kflf*kfpf*lipf*pipf + 4*F12*F20*F21*kflf*lfpf*lipf*pipf - 
       4*F11*F20*F22*kflf*lfpf*lipf*pipf + 8*F10*F21*F22*kflf*lfpf*lipf*pipf + 
       8*F10*F21*F22*kfpf*lfpf*lipf*pipf + 
       8*F10*F21*F22*pow(lfpf,2)*lipf*pipf + 
       4*F10*F21*F22*kflf*lfpi*lipf*pipf + 4*F10*F21*F22*kfpf*lfpi*lipf*pipf + 
       8*F10*F21*F22*lfpf*lfpi*lipf*pipf + 
       2*F12*F20*F21*pow(kflf,2)*lipi*pipf - 
       2*F11*F20*F22*pow(kflf,2)*lipi*pipf + 
       2*F10*F21*F22*pow(kflf,2)*lipi*pipf + 
       4*F12*F20*F21*kflf*kfpf*lipi*pipf - 4*F11*F20*F22*kflf*kfpf*lipi*pipf + 
       8*F10*F21*F22*kflf*kfpf*lipi*pipf + 4*F12*F20*F21*kflf*lfpf*lipi*pipf - 
       4*F11*F20*F22*kflf*lfpf*lipi*pipf + 8*F10*F21*F22*kflf*lfpf*lipi*pipf + 
       8*F10*F21*F22*kfpf*lfpf*lipi*pipf + 
       8*F10*F21*F22*pow(lfpf,2)*lipi*pipf + 
       4*F10*F21*F22*kflf*lfpi*lipi*pipf + 4*F10*F21*F22*kfpf*lfpi*lipi*pipf + 
       8*F10*F21*F22*lfpf*lfpi*lipi*pipf + 
       4*F12*F20*F21*kflf*kfpf*m2*pipf + 
       4*F11*F20*F22*kflf*kfpf*m2*pipf + 
       2*F10*F21*F22*kflf*kfpf*m2*pipf + 
       6*F20*F21*F22*kflf*kfpf*m2*pipf - 
       4*F12*F20*F21*pow(kfpf,2)*m2*pipf - 
       4*F11*F20*F22*pow(kfpf,2)*m2*pipf - 
       8*F10*F21*F22*pow(kfpf,2)*m2*pipf - 
       8*F20*F21*F22*pow(kfpf,2)*m2*pipf - 
       4*F10*F21*F22*kfpf*lfpf*m2*pipf - 
       4*F10*F21*F22*kfpf*lfpi*m2*pipf - 
       4*F12*F20*F21*kfpf*lipf*m2*pipf + 
       4*F11*F20*F22*kfpf*lipf*m2*pipf - 
       4*F10*F21*F22*kfpf*lipf*m2*pipf - 
       4*F12*F20*F21*kfpf*lipi*m2*pipf + 
       4*F11*F20*F22*kfpf*lipi*m2*pipf - 
       4*F10*F21*F22*kfpf*lipi*m2*pipf + 
       8*F10*F12*F21*kflf*lilf*M2*pipf + 
       8*F10*F11*F22*kflf*lilf*M2*pipf + 
       8*F10*F21*F22*kflf*lilf*M2*pipf + 
       8*F10*F12*F21*kfpf*lilf*M2*pipf + 
       8*F10*F11*F22*kfpf*lilf*M2*pipf + 
       8*F10*F21*F22*kfpf*lilf*M2*pipf + 
       16*F10*F12*F21*lfpf*lilf*M2*pipf + 
       16*F10*F11*F22*lfpf*lilf*M2*pipf + 
       16*F10*F21*F22*lfpf*lilf*M2*pipf - 
       16*F10*F11*F12*kflf*lipf*M2*pipf - 
       16*F10*F12*F21*kflf*lipf*M2*pipf - 
       16*F10*F11*F22*kflf*lipf*M2*pipf - 
       16*F10*F21*F22*kflf*lipf*M2*pipf + 
       8*F10*F12*F21*kflf*lipi*M2*pipf + 
       8*F10*F11*F22*kflf*lipi*M2*pipf + 
       8*F10*F21*F22*kflf*lipi*M2*pipf - 
       16*F10*F11*F12*kflf*m2*M2*pipf - 
       24*F10*F12*F21*kflf*m2*M2*pipf - 
       24*F10*F11*F22*kflf*m2*M2*pipf - 
       24*F10*F21*F22*kflf*m2*M2*pipf - 
       16*F10*F11*F12*kfpf*m2*M2*pipf - 
       24*F10*F12*F21*kfpf*m2*M2*pipf - 
       24*F10*F11*F22*kfpf*m2*M2*pipf - 
       24*F10*F21*F22*kfpf*m2*M2*pipf - 
       32*F10*F11*F12*lfpf*m2*M2*pipf - 
       48*F10*F12*F21*lfpf*m2*M2*pipf - 
       48*F10*F11*F22*lfpf*m2*M2*pipf - 
       48*F10*F21*F22*lfpf*m2*M2*pipf - 
       4*F10*F21*F22*kflf*lilf*pow(pipf,2) - 
       4*F10*F21*F22*kfpf*lilf*pow(pipf,2) - 
       8*F10*F21*F22*lfpf*lilf*pow(pipf,2) - 
       4*F10*F21*F22*kflf*lipi*pow(pipf,2) + 
       4*F10*F21*F22*kflf*m2*pow(pipf,2) + 
       4*F10*F21*F22*kfpf*m2*pow(pipf,2) + 
       8*F10*F21*F22*lfpf*m2*pow(pipf,2) + 
       kfli*(4*F20*(F11 + F21)*F22*pow(kflf,2)*(M2 - pipf) + 
       kfpi*(F22*kflf*(2*F11*F20*lfpf + F20*F21*lfpf - F20*F21*lfpi + 
       8*F10*F11*M2 + 4*F11*F20*M2 + 
       6*F10*F21*M2 + 4*F20*F21*M2 - 2*F10*F21*pipf) +
       2*((2*F11*F20 + (F10 + 2*F20)*F21)*F22*pow(lfpf,2) - 
       2*(F10 + F20)*(F11 + F21)*F22*m2*M2 + 
       lfpi*(F10*F21*F22*lfpf + 
       (-(F12*F20*F21) + (F11*F20 + 2*F10*(F11 + F21))*F22)*M2) + 
       lfpf*((F11*(-4*F12*F20 + 2*F10*F22 - 3*F20*F22) + 
       F21*(-3*F12*F20 + F10*F22 - 2*F20*F22))*M2 + 
       F21*(2*F12*F20 - F10*F22 + 2*F20*F22)*pipf))) - 
       kflf*(F22*kfpf*(F20*F21*lfpf - F20*(2*F11 + F21)*lfpi + 
       8*F10*F11*M2 + 8*F11*F20*M2 + 
       6*F10*F21*M2 + 4*F20*F21*M2 + 
       (4*F11*F20 - 2*F10*F21)*pipf) + 
       2*(M2*(-((F12*F20*F21 + 
       (-6*F10*F11 + F11*F20 - 5*F10*F21 + 2*F20*F21)*F22)*lfpf) +
       2*(2*F11 + F21)*(F20*(F12 + F22) + F10*(2*F12 + F22))*M2) - 
       ((F12*F20*F21 + (-5*F11*F20 + F10*F21 - 4*F20*F21)*F22)*
       lfpf + 4*((F10 + F20)*F21*(F12 + F22) + 
       F11*(F12*F20 + (F10 + F20)*F22))*M2)*pipf + 
       2*F21*(F12*F20 + (F10 + F20)*F22)*pow(pipf,2) + 
       lfpi*((F21*(3*F12*F20 - 3*F10*F22 - 2*F20*F22) + 
       F11*(4*F12*F20 - (2*F10 + F20)*F22))*M2 - 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf))) - 
       2*(kfpf*(F20*F21*F22*pow(lfpf,2) + 
       (F10 + F20)*F21*F22*pow(lfpi,2) - 
       4*F11*F12*F20*lfpf*M2 - 
       3*F12*F20*F21*lfpf*M2 + 
       2*F10*F11*F22*lfpf*M2 - 
       3*F11*F20*F22*lfpf*M2 + 
       2*F10*F21*F22*lfpf*M2 - 
       2*F20*F21*F22*lfpf*M2 - 
       2*F10*F11*F22*m2*M2 - 
       2*F10*F21*F22*m2*M2 - 
       2*F20*(F11 + F21)*F22*(-lfpf + m2)*pipf + 
       lfpi*((-2*F11*F20*F22 + F10*F21*F22)*lfpf + 
       (F12*F20*F21 + 2*F10*F11*F22 - F11*F20*F22 + 
       F10*F21*F22)*M2 - F10*F21*F22*pipf)) + 
       2*((F10 + F20)*(F12*F21 - F11*F22)*pow(lfpi,2)*M2 - 
       lfpf*lfpi*((-(F12*F20*F21) + F11*F20*F22 + 
       F10*(2*F11 + F21)*(2*F12 + F22))*M2 + 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) + 
       lfpf*(F10*M2*((F12*F21 + 3*(F11 + F21)*F22)*lfpf + 
       (2*F11 + F21)*(2*F12 + F22)*M2) - 
       ((F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*lfpf + 
       2*F10*(F12*F21 + (F11 + F21)*F22)*M2)*pipf + 
       F10*F21*F22*pow(pipf,2))))) + 
       kfpi*(F20*F22*pow(kflf,2)*((2*F11 + F21)*lipf - F21*lipi) + 
       2*kflf*(-2*F12*F20*F21*lfpf*lipf + 4*F11*F20*F22*lfpf*lipf + 
       F10*F21*F22*lfpf*lipf + F20*F21*F22*lfpf*lipf - 
       F12*F20*F21*lfpi*lipf + F11*F20*F22*lfpi*lipf - 
       F20*F21*F22*lfpi*lipf + 
       F22*kfpf*(-(F20*(2*F11 + F21)*lilf) + 2*F11*F20*lipf + 
       F10*F21*lipf - (F10 + 2*F20)*F21*lipi + 
       2*F11*F20*m2 + F20*F21*m2) + 
       2*F10*F11*F22*lilf*M2 - 2*F11*F20*F22*lilf*M2 + 
       F10*F21*F22*lilf*M2 - F20*F21*F22*lilf*M2 - 
       4*F10*F12*F21*lipf*M2 - F12*F20*F21*lipf*M2 - 
       2*F10*F11*F22*lipf*M2 + F11*F20*F22*lipf*M2 - 
       4*F10*F21*F22*lipf*M2 + 
       4*F11*F12*F20*m2*M2 + 
       4*F12*F20*F21*m2*M2 - 
       4*F10*F11*F22*m2*M2 + 
       4*F11*F20*F22*m2*M2 - 
       3*F10*F21*F22*m2*M2 + 
       3*F20*F21*F22*m2*M2 - F10*F21*F22*lilf*pipf + 
       F20*F21*F22*lilf*pipf + 2*F12*F20*F21*lipf*pipf + 
       2*F10*F21*F22*lipf*pipf + 2*F20*F21*F22*lipf*pipf + 
       F10*F21*F22*m2*pipf - F20*F21*F22*m2*pipf - 
       lipi*(F12*F20*F21*lfpf - F11*F20*F22*lfpf - F10*F21*F22*lfpf + 
       F20*F21*F22*lfpf + F20*F21*F22*lfpi + 
       4*F10*F12*F21*M2 + 3*F12*F20*F21*M2 - 
       2*F10*F11*F22*M2 - F11*F20*F22*M2 + 
       3*F10*F21*F22*M2 + 2*F20*F21*F22*M2 + 
       F10*F21*F22*pipf)) - 
       2*(2*F20*(2*F11 + F21)*F22*pow(kfpf,2)*m2 + 
       kfpf*(2*F12*F20*F21*lfpf*lipf - 2*F11*F20*F22*lfpf*lipf - 
       3*F10*F21*F22*lfpf*lipf - 2*F12*F20*F21*lfpi*lipf - 
       F10*F21*F22*lfpi*lipf - 2*F20*F21*F22*lfpi*lipf + 
       2*F11*F20*F22*lfpf*m2 - F10*F21*F22*lfpf*m2 + 
       2*F20*F21*F22*lfpf*m2 + F10*F21*F22*lfpi*m2 - 
       F20*F21*F22*lipf*m2 - 
       F22*lipi*((2*F11*F20 + F10*F21 + 2*F20*F21)*lfpf - 
       F10*F21*lfpi + F20*F21*m2) + 
       8*F10*F11*F12*m2*M2 + 
       8*F11*F12*F20*m2*M2 + 
       4*F10*F12*F21*m2*M2 + 
       6*F12*F20*F21*m2*M2 + 
       16*F10*F11*F22*m2*M2 + 
       6*F11*F20*F22*m2*M2 + 
       12*F10*F21*F22*m2*M2 + 
       4*F20*F21*F22*m2*M2 - 
       lilf*(-2*F11*F20*F22*lfpf - F10*F21*F22*lfpf - 
       F20*F21*F22*lfpf + (F10 + F20)*F21*F22*lfpi + 
       4*F11*F12*F20*M2 + 4*F12*F20*F21*M2 + 
       8*F10*F11*F22*M2 + 2*F11*F20*F22*M2 + 
       6*F10*F21*F22*M2 + 2*F20*F21*F22*M2 - 
       2*F21*(F12*F20 + (F10 + F20)*F22)*pipf)) - 
       2*(F10*F12*F21*lfpf*lilf*M2 + 
       F12*F20*F21*lfpf*lilf*M2 + 
       3*F10*F11*F22*lfpf*lilf*M2 - 
       F11*F20*F22*lfpf*lilf*M2 + 
       2*F10*F21*F22*lfpf*lilf*M2 + 
       F10*F12*F21*lfpi*lilf*M2 + 
       F12*F20*F21*lfpi*lilf*M2 - 
       F10*F11*F22*lfpi*lilf*M2 - 
       F11*F20*F22*lfpi*lilf*M2 - 
       4*F10*F11*F12*lfpf*m2*M2 - 
       4*F10*F12*F21*lfpf*m2*M2 - 
       8*F10*F11*F22*lfpf*m2*M2 - 
       6*F10*F21*F22*lfpf*m2*M2 - 
       2*F10*F21*F22*lfpf*lilf*pipf + 
       2*F10*F21*F22*lfpf*m2*pipf + 
       lipi*((-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*
       pow(lfpf,2) + (F10 + F20)*(F12*F21 - F11*F22)*m2*M2 + 
       F10*lfpf*(F21*F22*lfpi - 
       (2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2 + 
       F21*F22*pipf)) + 
       lipf*((-2*F12*F20*F21 + 2*F11*F20*F22 + F10*F21*F22)*pow(lfpf,2) + 
       (F10 + F20)*(F12*F21 - F11*F22)*m2*M2 + 
       lfpf*((-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*lfpi + 
       F10*((2*F11 + F21)*(2*F12 + F22)*M2 + 
       F21*F22*pipf)))))))/(2.*kfpf*M2);

res2 = (2*F20*F21*F22*pow(kfpi,3)*m2 + 
       pow(kfpi,2)*(F22*kflf*(-(F20*F21*lilf) + 12*F11*F20*lipf - 
       2*F10*F21*lipf + 8*F20*F21*lipf - 4*F20*F21*lipi + 
       F20*F21*m2) - 
       2*(2*F11*F20*F22*lfpf*lilf - F10*F21*F22*lfpf*lilf + 
       F20*F21*F22*lfpf*lilf - F20*F21*F22*lfpi*lilf - 
       2*F11*F20*F22*lfpf*lipf - 3*F10*F21*F22*lfpf*lipf - 
       2*F20*F21*F22*lfpf*lipf - 2*F10*F21*F22*lfpi*lipf + 
       2*F20*(2*F11 + F21)*F22*kfpf*m2 - 
       2*F11*F20*F22*lfpf*m2 + F10*F21*F22*lfpf*m2 - 
       2*F20*F21*F22*lfpf*m2 + 4*F11*F20*F22*lipf*m2 + 
       3*F20*F21*F22*lipf*m2 - 
       lipi*((-2*F12*F20*F21 + 2*F11*F20*F22 + F10*F21*F22)*lfpf + 
       F20*F21*F22*m2) - 4*F11*F12*F20*lilf*M2 - 
       3*F12*F20*F21*lilf*M2 - 6*F10*F11*F22*lilf*M2 - 
       3*F11*F20*F22*lilf*M2 - 4*F10*F21*F22*lilf*M2 - 
       2*F20*F21*F22*lilf*M2 + 
       8*F11*F12*F20*m2*M2 + 
       2*F10*F12*F21*m2*M2 + 
       7*F12*F20*F21*m2*M2 + 
       14*F10*F11*F22*m2*M2 + 
       7*F11*F20*F22*m2*M2 + 
       12*F10*F21*F22*m2*M2 + 
       6*F20*F21*F22*m2*M2 + 2*F11*F20*F22*lilf*pipf + 
       2*F10*F21*F22*lilf*pipf + 2*F20*F21*F22*lilf*pipf - 
       2*F12*F20*F21*m2*pipf - 2*F11*F20*F22*m2*pipf - 
       4*F10*F21*F22*m2*pipf - 4*F20*F21*F22*m2*pipf)) - 
       pow(kflf,2)*(4*F20*(F11 + F21)*F22*kfli*(M2 - pipf) + 
       lipf*(-(F20*F21*F22*kfpf) + 
       2*((4*F11*F12*F20 + 3*F12*F20*F21 + 2*F10*F11*F22 + 
       3*F11*F20*F22 + F10*F21*F22 + 2*F20*F21*F22)*M2 - 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf)) + 
       lipi*(F20*(2*F11 + F21)*F22*kfpf - 
       2*((F12*F20*F21 + (-(F10*(2*F11 + F21)) + F20*(F11 + 2*F21))*F22)*M2 + 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf))) + 
       kflf*(F21*F22*pow(kfpf,2)*(2*F10*lipi + F20*(-lilf + m2)) + 
       2*kfpf*(-(F20*F21*F22*lfpf*lipf) - F12*F20*F21*lfpi*lipf + 
       F11*F20*F22*lfpi*lipf + F10*F21*F22*lfpi*lipf - 
       F20*F21*F22*lfpi*lipf + 4*F10*F12*F21*lipf*M2 + 
       3*F12*F20*F21*lipf*M2 - 2*F10*F11*F22*lipf*M2 - 
       F11*F20*F22*lipf*M2 + 3*F10*F21*F22*lipf*M2 + 
       2*F20*F21*F22*lipf*M2 + 
       4*F11*F12*F20*m2*M2 + 
       4*F12*F20*F21*m2*M2 - 
       4*F10*F11*F22*m2*M2 + 
       4*F11*F20*F22*m2*M2 - 
       3*F10*F21*F22*m2*M2 + 
       3*F20*F21*F22*m2*M2 + F10*F21*F22*lipf*pipf + 
       F10*F21*F22*m2*pipf - F20*F21*F22*m2*pipf + 
       (F10 - F20)*F22*lilf*((2*F11 + F21)*M2 - F21*pipf) + 
       lipi*(-(F12*F20*F21*lfpf) + F11*F20*F22*lfpf - 
       F20*F21*F22*lfpf + 
       (-2*F12*F20*F21 + (4*F11*F20 + (F10 + F20)*F21)*F22)*lfpi + 
       4*F10*F12*F21*M2 + F12*F20*F21*M2 + 
       2*F10*F11*F22*M2 - F11*F20*F22*M2 + 
       4*F10*F21*F22*M2 - 
       2*F21*(F12*F20 + (F10 + F20)*F22)*pipf)) - 
       4*(-3*F10*F12*F21*lfpf*lipf*M2 - 
       F12*F20*F21*lfpf*lipf*M2 - 
       F10*F11*F22*lfpf*lipf*M2 + 
       F11*F20*F22*lfpf*lipf*M2 - 
       3*F10*F21*F22*lfpf*lipf*M2 + 
       8*F10*F11*F12*lfpi*lipf*M2 + 
       5*F10*F12*F21*lfpi*lipf*M2 + 
       3*F10*F11*F22*lfpi*lipf*M2 + 
       2*F10*F21*F22*lfpi*lipf*M2 - 
       4*F10*F11*F12*lilf*M4 - 2*F10*F12*F21*lilf*M4 - 
       2*F10*F11*F22*lilf*M4 - F10*F21*F22*lilf*M4 + 
       4*F10*F11*F12*lipf*M4 + 2*F10*F12*F21*lipf*M4 + 
       2*F10*F11*F22*lipf*M4 + F10*F21*F22*lipf*M4 + 
       8*F10*F11*F12*m2*M4 + 
       6*F10*F12*F21*m2*M4 + 
       6*F10*F11*F22*m2*M4 + 
       5*F10*F21*F22*m2*M4 + 
       F10*F21*F22*lfpf*lipf*pipf + F12*F20*F21*lfpi*lipf*pipf - 
       F11*F20*F22*lfpi*lipf*pipf + 2*F10*F21*F22*lfpi*lipf*pipf + 
       2*F10*F12*F21*lilf*M2*pipf + 
       2*F10*F11*F22*lilf*M2*pipf + 
       2*F10*F21*F22*lilf*M2*pipf - 
       2*F10*F12*F21*lipf*M2*pipf - 
       2*F10*F11*F22*lipf*M2*pipf - 
       2*F10*F21*F22*lipf*M2*pipf - 
       4*F10*F11*F12*m2*M2*pipf - 
       6*F10*F12*F21*m2*M2*pipf - 
       6*F10*F11*F22*m2*M2*pipf - 
       6*F10*F21*F22*m2*M2*pipf - 
       F10*F21*F22*lilf*pow(pipf,2) + 
       F10*F21*F22*lipf*pow(pipf,2) + 
       F10*F21*F22*m2*pow(pipf,2) + 
       lipi*(M2*((-(F12*F20*F21) + F11*F20*F22 + 
       F10*(F21*(F12 + F22) + F11*(4*F12 + 3*F22)))*lfpf - 
       4*F10*(F11 + F21)*(F12 + F22)*M2) + 
       F10*(F21*F22*lfpf + 4*(F11 + F21)*(F12 + F22)*M2)*
       pipf - lfpi*(F10*(3*F12*F21 + 5*F11*F22 + 6*F21*F22)*
       M2 + (-(F12*F20*F21) + F11*F20*F22 - 2*F10*F21*F22)*pipf))) + 
       kfli*(F22*kfpf*(F20*F21*lfpf - F20*(2*F11 + F21)*lfpi + 
       8*F10*F11*M2 + 4*F11*F20*M2 + 
       6*F10*F21*M2 + 4*F20*F21*M2 - 2*F10*F21*pipf) +
       2*(M2*((F21*(-3*F12*F20 + 3*F10*F22 + 2*F20*F22) + 
       F11*(-4*F12*F20 + (2*F10 + F20)*F22))*lfpf + 
       2*(2*F11 + F21)*(F20*(F12 + F22) + F10*(2*F12 + F22))*M2) - 
       (-((F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*lfpf) + 
       4*((F10 + F20)*F21*(F12 + F22) + 
       F11*(F12*F20 + (F10 + F20)*F22))*M2)*pipf + 
       2*F21*(F12*F20 + (F10 + F20)*F22)*pow(pipf,2) + 
       lfpi*((F12*F20*F21 + 
       (-6*F10*F11 + F11*F20 - 5*F10*F21 + 2*F20*F21)*F22)*M2 + (F12*F20*F21 + 
       (-5*F11*F20 + F10*F21 - 4*F20*F21)*F22)*pipf)))) + 
       2*(-(pow(kfpf,2)*(F21*lfpi*((-2*F12*F20 + (F10 - 2*F20)*F22)*lipi + 
       F10*F22*(lipf - m2)) + 
       (2*F10 + F20)*(F12*F21 + (F11 + 2*F21)*F22)*m2*
       M2 + lilf*(F10*F21*F22*lfpi + 
       (F12*F20*F21 - (F11*F20 + 2*F10*(F11 + F21))*F22)*M2))) + 
       4*F10*lfpi*(-2*F12*F21*lfpf*lipf*M2 - 
       2*F11*F22*lfpf*lipf*M2 - 
       3*F21*F22*lfpf*lipf*M2 + 
       4*F11*F12*lfpi*lipf*M2 + 
       2*F12*F21*lfpi*lipf*M2 + 
       2*F11*F22*lfpi*lipf*M2 + F21*F22*lfpi*lipf*M2 + 
       8*F11*F12*m2*M4 + 
       6*F12*F21*m2*M4 + 
       6*F11*F22*m2*M4 + 
       5*F21*F22*m2*M4 + F21*F22*lfpf*lipf*pipf + 
       F21*F22*lfpi*lipf*pipf - 4*F11*F12*m2*M2*pipf - 
       6*F12*F21*m2*M2*pipf - 
       6*F11*F22*m2*M2*pipf - 
       6*F21*F22*m2*M2*pipf + 
       F21*F22*m2*pow(pipf,2) - 
       lilf*((2*F11 + F21)*M2 - F21*pipf)*
       ((2*F12 + F22)*M2 - F22*pipf) + 
       lipi*(lfpf*((2*F11 + F21)*(2*F12 + F22)*M2 + 
       F21*F22*pipf) + lfpi*(-((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2) + 
       F21*F22*pipf))) - 2*kfpf*(F10*F12*F21*lfpf*lilf*M2 + 
       F12*F20*F21*lfpf*lilf*M2 - 
       F10*F11*F22*lfpf*lilf*M2 - 
       F11*F20*F22*lfpf*lilf*M2 + 
       F10*F12*F21*lfpi*lilf*M2 + 
       F12*F20*F21*lfpi*lilf*M2 + 
       3*F10*F11*F22*lfpi*lilf*M2 - 
       F11*F20*F22*lfpi*lilf*M2 + 
       2*F10*F21*F22*lfpi*lilf*M2 - 
       4*F10*F11*F12*lfpi*m2*M2 - 
       4*F10*F12*F21*lfpi*m2*M2 - 
       8*F10*F11*F22*lfpi*m2*M2 - 
       6*F10*F21*F22*lfpi*m2*M2 - 
       2*F10*F21*F22*lfpi*lilf*pipf + 
       2*F10*F21*F22*lfpi*m2*pipf + 
       lipf*((-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*pow(lfpi,2) + 
       (F10 + F20)*(F12*F21 - F11*F22)*m2*M2 + 
       F10*lfpi*(F21*F22*lfpf + 
       (2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2 - F21*F22*pipf)) + 
       lipi*((-2*F12*F20*F21 + 2*F11*F20*F22 + F10*F21*F22)*pow(lfpi,2) + 
       (F10 + F20)*(F12*F21 - F11*F22)*m2*M2 - 
       lfpi*((F12*F20*F21 - (F11*F20 + F10*F21)*F22)*lfpf + 
       F10*(2*F11 + F21)*(2*F12 + F22)*M2 + 
       F10*F21*F22*pipf))) - 
       kfli*(kfpf*(-((2*F11*F20 + (F10 + 2*F20)*F21)*F22*pow(lfpi,2)) + 
       ((-(F12*F20*F21) + (F11*F20 + 2*F10*(F11 + F21))*F22)*
       lfpf + 2*(F10 + F20)*(F11 + F21)*F22*m2)*M2 -
       lfpi*(F10*F21*F22*lfpf + 
       (4*F11*F12*F20 + 3*F12*F20*F21 - 2*F10*F11*F22 + 
       3*F11*F20*F22 - F10*F21*F22 + 2*F20*F21*F22)*M2 +
       F21*(-2*F12*F20 + (F10 - 2*F20)*F22)*pipf)) +
       2*(-((F10 + F20)*(F12*F21 - F11*F22)*pow(lfpf,2)*M2) - 
       pow(lfpi,2)*(F10*(F12*F21 + 3*(F11 + F21)*F22)*M2 - 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) + 
       lfpi*(M2*((-(F12*F20*F21) + F11*F20*F22 + 
       F10*(2*F11 + F21)*(2*F12 + F22))*lfpf + 
       F10*(2*F11 + F21)*(2*F12 + F22)*M2) + 
       ((F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*lfpf - 
       2*F10*(F12*F21 + (F11 + F21)*F22)*M2)*pipf + 
       F10*F21*F22*pow(pipf,2))))) + 
       kfpi*(F20*F22*pow(kflf,2)*(-((2*F11 + F21)*lipf) + F21*lipi) - 
       2*kflf*(-2*F12*F20*F21*lfpf*lipf - 2*F11*F20*F22*lfpf*lipf + 
       F10*F21*F22*lfpf*lipf - 2*F20*F21*F22*lfpf*lipf - 
       F12*F20*F21*lfpi*lipf - 5*F11*F20*F22*lfpi*lipf - 
       4*F20*F21*F22*lfpi*lipf + 
       F22*kfpf*(-(F20*(2*F11 + F21)*lilf) + F10*F21*lipf + 
       2*F20*F21*lipf - (2*F11*F20 + F10*F21)*lipi + 
       2*F11*F20*m2 + F20*F21*m2) - 
       4*F11*F12*F20*lilf*M2 - 2*F12*F20*F21*lilf*M2 - 
       2*F10*F11*F22*lilf*M2 - 2*F11*F20*F22*lilf*M2 - 
       F10*F21*F22*lilf*M2 - F20*F21*F22*lilf*M2 + 
       16*F10*F11*F12*lipf*M2 + 4*F11*F12*F20*lipf*M2 + 
       12*F10*F12*F21*lipf*M2 + 3*F12*F20*F21*lipf*M2 + 
       2*F10*F11*F22*lipf*M2 + 3*F11*F20*F22*lipf*M2 + 
       2*F10*F21*F22*lipf*M2 + 2*F20*F21*F22*lipf*M2 + 
       8*F11*F12*F20*m2*M2 + 
       6*F12*F20*F21*m2*M2 + 
       6*F11*F20*F22*m2*M2 - 
       F10*F21*F22*m2*M2 + 
       5*F20*F21*F22*m2*M2 + 2*F12*F20*F21*lilf*pipf + 
       F10*F21*F22*lilf*pipf + F20*F21*F22*lilf*pipf + 
       2*F12*F20*F21*lipf*pipf - 2*F11*F20*F22*lipf*pipf + 
       4*F10*F21*F22*lipf*pipf - 2*F12*F20*F21*m2*pipf - 
       2*F11*F20*F22*m2*pipf - F10*F21*F22*m2*pipf - 
       3*F20*F21*F22*m2*pipf + 
       lipi*(-(F12*F20*F21*lfpf) + F11*F20*F22*lfpf + 
       F10*F21*F22*lfpf + 2*F20*F21*F22*lfpf + 2*F20*F21*F22*lfpi - 
       4*F11*F12*F20*M2 - 4*F10*F12*F21*M2 - 
       5*F12*F20*F21*M2 - 2*F10*F11*F22*M2 - 
       3*F11*F20*F22*M2 - 3*F10*F21*F22*M2 - 
       4*F20*F21*F22*M2 + 
       F21*(2*F12*F20 + 3*F10*F22 + 2*F20*F22)*pipf)) + 
       kfli*(-(F22*kflf*(2*F11*F20*lfpf + F20*F21*lfpf - F20*F21*lfpi + 
       8*F10*F11*M2 + 8*F11*F20*M2 + 
       6*F10*F21*M2 + 4*F20*F21*M2 + 
       (4*F11*F20 - 2*F10*F21)*pipf)) + 
       2*(-(F10*F21*F22*pow(lfpf,2)) - F20*F21*F22*pow(lfpf,2) - 
       F20*F21*F22*pow(lfpi,2) + F12*F20*F21*lfpf*M2 + 
       2*F10*F11*F22*lfpf*M2 - F11*F20*F22*lfpf*M2 + 
       F10*F21*F22*lfpf*M2 + 
       2*F10*F11*F22*m2*M2 + 
       2*F10*F21*F22*m2*M2 - F10*F21*F22*lfpf*pipf + 
       2*F11*F20*F22*m2*pipf + 
       2*F20*F21*F22*m2*pipf + 
       lfpi*((2*F11*F20 - F10*F21)*F22*lfpf - 
       (4*F11*F12*F20 + 3*F12*F20*F21 - 2*F10*F11*F22 + 
       3*F11*F20*F22 - 2*F10*F21*F22 + 2*F20*F21*F22)*
       M2 + 2*F20*(F11 + F21)*F22*pipf))) + 
       2*(F20*F21*F22*pow(kfpf,2)*m2 + 
       kfpf*(F10*F21*F22*lfpf*lipf - 2*F11*F20*F22*lfpi*lipf - 
       F10*F21*F22*lfpi*lipf - 2*F20*F21*F22*lfpi*lipf - 
       F10*F21*F22*lfpf*m2 - 2*F11*F20*F22*lfpi*m2 + 
       F10*F21*F22*lfpi*m2 - 2*F20*F21*F22*lfpi*m2 + 
       F20*F21*F22*lipf*m2 + 
       lipi*((2*F12*F20*F21 - 2*F11*F20*F22 - 3*F10*F21*F22)*lfpi + 
       F21*(-((2*F12*F20 + F10*F22 + 2*F20*F22)*lfpf) + 
       F20*F22*m2)) + 
       8*F10*F11*F12*m2*M2 + 
       8*F11*F12*F20*m2*M2 + 
       4*F10*F12*F21*m2*M2 + 
       6*F12*F20*F21*m2*M2 + 
       16*F10*F11*F22*m2*M2 + 
       6*F11*F20*F22*m2*M2 + 
       12*F10*F21*F22*m2*M2 + 
       4*F20*F21*F22*m2*M2 - 
       lilf*(-(F10*F21*F22*lfpf) - F20*F21*F22*lfpf + 
       (2*F11*F20 + (F10 + F20)*F21)*F22*lfpi + 
       4*F11*F12*F20*M2 + 4*F12*F20*F21*M2 + 
       8*F10*F11*F22*M2 + 2*F11*F20*F22*M2 + 
       6*F10*F21*F22*M2 + 2*F20*F21*F22*M2 - 
       2*F21*(F12*F20 + (F10 + F20)*F22)*pipf)) + 
       2*(F10*F21*F22*pow(lfpf,2)*lipf - F12*F20*F21*lfpf*lfpi*lipf + 
       F11*F20*F22*lfpf*lfpi*lipf + 2*F10*F21*F22*lfpf*lfpi*lipf + 
       F10*F21*F22*pow(lfpi,2)*lipf - 
       2*F10*F12*F21*lfpf*lipf*M2 - 
       2*F10*F11*F22*lfpf*lipf*M2 - 
       3*F10*F21*F22*lfpf*lipf*M2 + 
       8*F10*F11*F12*lfpi*lipf*M2 + 
       4*F10*F12*F21*lfpi*lipf*M2 + 
       4*F10*F11*F22*lfpi*lipf*M2 + 
       2*F10*F21*F22*lfpi*lipf*M2 - 
       4*F10*F11*F12*lfpf*m2*M2 - 
       4*F10*F12*F21*lfpf*m2*M2 - 
       2*F10*F11*F22*lfpf*m2*M2 - 
       F10*F21*F22*lfpf*m2*M2 - 
       6*F10*F11*F22*lfpi*m2*M2 - 
       5*F10*F21*F22*lfpi*m2*M2 + 
       4*F10*F11*F12*lipf*m2*M2 + 
       3*F10*F12*F21*lipf*m2*M2 - 
       F10*F11*F22*lipf*m2*M2 - 
       F10*F21*F22*lipf*m2*M2 + 
       8*F10*F11*F12*m2*M4 + 
       6*F10*F12*F21*m2*M4 + 
       6*F10*F11*F22*m2*M4 + 
       5*F10*F21*F22*m2*M4 + 
       F10*F21*F22*lfpf*lipf*pipf + 2*F10*F21*F22*lfpi*lipf*pipf + 
       F10*F21*F22*lfpf*m2*pipf + 
       F10*F21*F22*lfpi*m2*pipf + 
       F12*F20*F21*lipf*m2*pipf - 
       F11*F20*F22*lipf*m2*pipf + 
       F10*F21*F22*lipf*m2*pipf - 
       4*F10*F11*F12*m2*M2*pipf - 
       6*F10*F12*F21*m2*M2*pipf - 
       6*F10*F11*F22*m2*M2*pipf - 
       6*F10*F21*F22*m2*M2*pipf + 
       F10*F21*F22*m2*pow(pipf,2) + 
       lilf*(F10*M2*((4*F11*F12 + 3*F12*F21 + F11*F22)*lfpf - 
       (2*F11 + F21)*(2*F12 + F22)*M2) + 
       (F20*(F12*F21 - F11*F22)*lfpf + 
       2*F10*(F12*F21 + (F11 + F21)*F22)*M2)*pipf - 
       F10*F21*F22*pow(pipf,2) - 
       (F12*F21 - F11*F22)*lfpi*(F10*M2 - F20*pipf)) + 
       lipi*(-(F12*F20*F21*pow(lfpf,2)) + 
       F11*F20*F22*pow(lfpf,2) + 4*F10*F11*F12*lfpf*M2 + 
       2*F10*F12*F21*lfpf*M2 + 
       2*F10*F11*F22*lfpf*M2 + 
       F10*F21*F22*lfpf*M2 - 
       F10*F12*F21*m2*M2 - 
       F10*F11*F22*m2*M2 - 
       F10*F21*F22*m2*M2 + 
       (F10*F21*F22*lfpf + 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*m2)*
       pipf - 2*lfpi*(F20*(F12*F21 - F11*F22)*lfpf + 
       F10*(2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2 - 
       F10*F21*F22*pipf))))))/(2.*kfpi*M2);

res3 = -(pow(kfli,2)*(4*F20*(F11 + F21)*F22*kflf*(M2 - pipf) - 
       lfpi*(-(F20*(2*F11 + F21)*F22*kfpf) + F20*F21*F22*kfpi + 
       2*((4*F11*F12*F20 + 3*F12*F20*F21 + 2*F10*F11*F22 + 
       3*F11*F20*F22 + F10*F21*F22 + 2*F20*F21*F22)*M2 - 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf)) + 
       lfpf*(-(F20*F21*F22*kfpf) + F20*(2*F11 + F21)*F22*kfpi + 
       2*((F12*F20*F21 + (-(F10*(2*F11 + F21)) + F20*(F11 + 2*F21))*F22)*M2 + 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf))) + 
       kfli*(4*F20*F21*F22*pow(kfpf,2)*lfpf - 
       12*F11*F20*F22*pow(kfpf,2)*lfpi + 
       2*F10*F21*F22*pow(kfpf,2)*lfpi - 
       8*F20*F21*F22*pow(kfpf,2)*lfpi + 
       F20*F21*F22*pow(kfpf,2)*lilf - 4*F20*F21*F22*kfpf*lfpf*lipf + 
       2*F12*F20*F21*kfpf*lfpi*lipf + 10*F11*F20*F22*kfpf*lfpi*lipf + 
       8*F20*F21*F22*kfpf*lfpi*lipf + 2*F12*F20*F21*kfpf*lfpf*lipi - 
       2*F11*F20*F22*kfpf*lfpf*lipi - 2*F10*F21*F22*kfpf*lfpf*lipi - 
       4*F20*F21*F22*kfpf*lfpf*lipi + 4*F12*F20*F21*kfpf*lfpi*lipi + 
       4*F11*F20*F22*kfpf*lfpi*lipi - 2*F10*F21*F22*kfpf*lfpi*lipi + 
       4*F20*F21*F22*kfpf*lfpi*lipi - 
       F20*F21*F22*pow(kfpf,2)*m2 + 
       F21*F22*pow(kfpi,2)*(-2*F10*lfpf + F20*lilf - F20*m2) + 
       8*F11*F12*F20*kfpf*lfpf*M2 + 
       8*F10*F12*F21*kfpf*lfpf*M2 + 
       10*F12*F20*F21*kfpf*lfpf*M2 + 
       4*F10*F11*F22*kfpf*lfpf*M2 + 
       6*F11*F20*F22*kfpf*lfpf*M2 + 
       6*F10*F21*F22*kfpf*lfpf*M2 + 
       8*F20*F21*F22*kfpf*lfpf*M2 - 
       32*F10*F11*F12*kfpf*lfpi*M2 - 
       8*F11*F12*F20*kfpf*lfpi*M2 - 
       24*F10*F12*F21*kfpf*lfpi*M2 - 
       6*F12*F20*F21*kfpf*lfpi*M2 - 
       4*F10*F11*F22*kfpf*lfpi*M2 - 
       6*F11*F20*F22*kfpf*lfpi*M2 - 
       4*F10*F21*F22*kfpf*lfpi*M2 - 
       4*F20*F21*F22*kfpf*lfpi*M2 + 
       8*F11*F12*F20*kfpf*lilf*M2 + 
       4*F12*F20*F21*kfpf*lilf*M2 + 
       4*F10*F11*F22*kfpf*lilf*M2 + 
       4*F11*F20*F22*kfpf*lilf*M2 + 
       2*F10*F21*F22*kfpf*lilf*M2 + 
       2*F20*F21*F22*kfpf*lilf*M2 - 
       12*F10*F12*F21*lfpf*lipf*M2 - 
       20*F10*F11*F22*lfpf*lipf*M2 - 
       24*F10*F21*F22*lfpf*lipf*M2 + 
       32*F10*F11*F12*lfpi*lipf*M2 + 
       20*F10*F12*F21*lfpi*lipf*M2 + 
       12*F10*F11*F22*lfpi*lipf*M2 + 
       8*F10*F21*F22*lfpi*lipf*M2 + 
       16*F10*F11*F12*lfpf*lipi*M2 + 
       4*F10*F12*F21*lfpf*lipi*M2 - 
       4*F12*F20*F21*lfpf*lipi*M2 + 
       12*F10*F11*F22*lfpf*lipi*M2 + 
       4*F11*F20*F22*lfpf*lipi*M2 + 
       4*F10*F21*F22*lfpf*lipi*M2 - 
       12*F10*F12*F21*lfpi*lipi*M2 - 
       4*F12*F20*F21*lfpi*lipi*M2 - 
       4*F10*F11*F22*lfpi*lipi*M2 + 
       4*F11*F20*F22*lfpi*lipi*M2 - 
       12*F10*F21*F22*lfpi*lipi*M2 - 
       16*F11*F12*F20*kfpf*m2*M2 - 
       12*F12*F20*F21*kfpf*m2*M2 - 
       12*F11*F20*F22*kfpf*m2*M2 + 
       2*F10*F21*F22*kfpf*m2*M2 - 
       10*F20*F21*F22*kfpf*m2*M2 - 
       16*F10*F11*F12*lfpf*M4 - 16*F10*F12*F21*lfpf*M4 - 
       16*F10*F11*F22*lfpf*M4 - 16*F10*F21*F22*lfpf*M4 + 
       16*F10*F11*F12*lfpi*M4 + 8*F10*F12*F21*lfpi*M4 + 
       8*F10*F11*F22*lfpi*M4 + 4*F10*F21*F22*lfpi*M4 - 
       16*F10*F11*F12*lilf*M4 - 8*F10*F12*F21*lilf*M4 - 
       8*F10*F11*F22*lilf*M4 - 4*F10*F21*F22*lilf*M4 + 
       32*F10*F11*F12*m2*M4 + 
       24*F10*F12*F21*m2*M4 + 
       24*F10*F11*F22*m2*M4 + 
       20*F10*F21*F22*m2*M4 - 
       4*F12*F20*F21*kfpf*lfpf*pipf - 6*F10*F21*F22*kfpf*lfpf*pipf - 
       4*F20*F21*F22*kfpf*lfpf*pipf - 4*F12*F20*F21*kfpf*lfpi*pipf + 
       4*F11*F20*F22*kfpf*lfpi*pipf - 8*F10*F21*F22*kfpf*lfpi*pipf - 
       4*F12*F20*F21*kfpf*lilf*pipf - 2*F10*F21*F22*kfpf*lilf*pipf - 
       2*F20*F21*F22*kfpf*lilf*pipf + 4*F12*F20*F21*lfpf*lipf*pipf - 
       4*F11*F20*F22*lfpf*lipf*pipf + 8*F10*F21*F22*lfpf*lipf*pipf + 
       4*F12*F20*F21*lfpi*lipf*pipf - 4*F11*F20*F22*lfpi*lipf*pipf + 
       8*F10*F21*F22*lfpi*lipf*pipf + 4*F10*F21*F22*lfpf*lipi*pipf + 
       4*F10*F21*F22*lfpi*lipi*pipf + 
       4*F12*F20*F21*kfpf*m2*pipf + 
       4*F11*F20*F22*kfpf*m2*pipf + 
       2*F10*F21*F22*kfpf*m2*pipf + 
       6*F20*F21*F22*kfpf*m2*pipf + 
       16*F10*F11*F12*lfpf*M2*pipf + 
       16*F10*F12*F21*lfpf*M2*pipf + 
       16*F10*F11*F22*lfpf*M2*pipf + 
       16*F10*F21*F22*lfpf*M2*pipf - 
       8*F10*F12*F21*lfpi*M2*pipf - 
       8*F10*F11*F22*lfpi*M2*pipf - 
       8*F10*F21*F22*lfpi*M2*pipf + 
       8*F10*F12*F21*lilf*M2*pipf + 
       8*F10*F11*F22*lilf*M2*pipf + 
       8*F10*F21*F22*lilf*M2*pipf - 
       16*F10*F11*F12*m2*M2*pipf - 
       24*F10*F12*F21*m2*M2*pipf - 
       24*F10*F11*F22*m2*M2*pipf - 
       24*F10*F21*F22*m2*M2*pipf + 
       4*F10*F21*F22*lfpi*pow(pipf,2) - 
       4*F10*F21*F22*lilf*pow(pipf,2) + 
       4*F10*F21*F22*m2*pow(pipf,2) + 
       kfpi*(-(F22*kflf*(-(F20*(2*F11 + F21)*lipf) + F20*F21*lipi + 
       8*F10*F11*M2 + 4*F11*F20*M2 + 
       6*F10*F21*M2 + 4*F20*F21*M2 - 
       2*F10*F21*pipf)) + 
       2*(-2*F12*F20*F21*lfpf*lipf + 4*F11*F20*F22*lfpf*lipf + 
       F10*F21*F22*lfpf*lipf + F20*F21*F22*lfpf*lipf - 
       F12*F20*F21*lfpi*lipf + F11*F20*F22*lfpi*lipf + 
       F10*F21*F22*lfpi*lipf - F20*F21*F22*lfpi*lipf - 
       F12*F20*F21*lfpf*lipi + F11*F20*F22*lfpf*lipi - 
       F20*F21*F22*lfpf*lipi - F20*F21*F22*lfpi*lipi + 
       F22*kfpf*(-2*F11*F20*lfpf - F10*F21*lfpf + 
       (F10 + 2*F20)*F21*lfpi - F20*(2*F11 + F21)*lilf + 
       2*F11*F20*m2 + F20*F21*m2) + 
       4*F10*F12*F21*lfpf*M2 + 
       F12*F20*F21*lfpf*M2 + 
       2*F10*F11*F22*lfpf*M2 - 
       F11*F20*F22*lfpf*M2 + 
       4*F10*F21*F22*lfpf*M2 + 
       4*F10*F12*F21*lfpi*M2 + 
       3*F12*F20*F21*lfpi*M2 - 
       2*F10*F11*F22*lfpi*M2 - 
       F11*F20*F22*lfpi*M2 + 
       3*F10*F21*F22*lfpi*M2 + 
       2*F20*F21*F22*lfpi*M2 + 
       4*F11*F12*F20*m2*M2 + 
       4*F12*F20*F21*m2*M2 - 
       4*F10*F11*F22*m2*M2 + 
       4*F11*F20*F22*m2*M2 - 
       3*F10*F21*F22*m2*M2 + 
       3*F20*F21*F22*m2*M2 - 
       2*F12*F20*F21*lfpf*pipf - 2*F10*F21*F22*lfpf*pipf - 
       2*F20*F21*F22*lfpf*pipf + F10*F21*F22*lfpi*pipf + 
       F10*F21*F22*m2*pipf - F20*F21*F22*m2*pipf + 
       (F10 - F20)*F22*lilf*((2*F11 + F21)*M2 - F21*pipf))) + 
       kflf*(F22*kfpf*(-(F20*F21*lipf) + F20*(2*F11 + F21)*lipi + 
       8*F10*F11*M2 + 8*F11*F20*M2 + 
       6*F10*F21*M2 + 4*F20*F21*M2 + 
       4*F11*F20*pipf - 2*F10*F21*pipf) + 
       2*(2*((2*F11 + F21)*M2 - F21*pipf)*
       ((F20*(F12 + F22) + F10*(2*F12 + F22))*M2 - 
       (F12*F20 + (F10 + F20)*F22)*pipf) + 
       lipi*((F11*(-4*F12*F20 + 2*F10*F22 + F20*F22) + 
       F21*(-3*F12*F20 + 3*F10*F22 + 2*F20*F22))*M2 + 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) + 
       lipf*((F12*F20*F21 + 
       (-6*F10*F11 + F11*F20 - 5*F10*F21 + 2*F20*F21)*F22)*M2 + 
       (F12*F20*F21 + (-5*F11*F20 + F10*F21 - 4*F20*F21)*F22)*pipf)))) - 
       2*(F20*F21*F22*pow(kfpf,3)*m2 + 
       pow(kfpi,2)*(-2*F12*F20*F21*lfpf*lipf + F10*F21*F22*lfpf*lipf - 
       2*F20*F21*F22*lfpf*lipf + F10*F21*F22*lfpi*lipf + 
       F20*F21*F22*kfpf*m2 - F10*F21*F22*lipf*m2 + 
       2*F10*F12*F21*m2*M2 + 
       F12*F20*F21*m2*M2 + 
       2*F10*F11*F22*m2*M2 + 
       F11*F20*F22*m2*M2 + 
       4*F10*F21*F22*m2*M2 + 
       2*F20*F21*F22*m2*M2 + 
       lilf*(F10*F21*F22*lipf + 
       (F12*F20*F21 - (F11*F20 + 2*F10*(F11 + F21))*F22)*M2)) + 
       pow(kfpf,2)*(-2*F10*F21*F22*lfpi*lipf - 
       F20*F21*F22*lfpf*m2 + 4*F11*F20*F22*lfpi*m2 + 
       3*F20*F21*F22*lfpi*m2 - 
       lipi*(-2*F12*F20*F21*lfpf + 2*F11*F20*F22*lfpf + 
       F10*F21*F22*lfpf + 
       (2*F11*F20 + 3*F10*F21 + 2*F20*F21)*F22*lfpi + 
       (-2*F11*F20 + (F10 - F20)*F21)*F22*lilf + 
       2*F11*F20*F22*m2 - F10*F21*F22*m2 + 
       2*F20*F21*F22*m2) + 
       8*F11*F12*F20*m2*M2 + 
       2*F10*F12*F21*m2*M2 + 
       7*F12*F20*F21*m2*M2 + 
       14*F10*F11*F22*m2*M2 + 
       7*F11*F20*F22*m2*M2 + 
       12*F10*F21*F22*m2*M2 + 
       6*F20*F21*F22*m2*M2 - 
       2*F12*F20*F21*m2*pipf - 2*F11*F20*F22*m2*pipf - 
       4*F10*F21*F22*m2*pipf - 4*F20*F21*F22*m2*pipf - 
       lilf*(F20*F21*F22*lipf + 
       (F21*(3*F12*F20 + 4*F10*F22 + 2*F20*F22) + 
       F11*(4*F12*F20 + 6*F10*F22 + 3*F20*F22))*M2 - 
       2*(F11*F20 + (F10 + F20)*F21)*F22*pipf)) - 
       2*(2*F10*lipf*(-2*F12*F21*lfpf*lipf*M2 - 
       2*F11*F22*lfpf*lipf*M2 - 
       3*F21*F22*lfpf*lipf*M2 + 
       4*F11*F12*lfpi*lipf*M2 + 
       2*F12*F21*lfpi*lipf*M2 + 
       2*F11*F22*lfpi*lipf*M2 + 
       F21*F22*lfpi*lipf*M2 + 
       8*F11*F12*m2*M4 + 
       6*F12*F21*m2*M4 + 
       6*F11*F22*m2*M4 + 
       5*F21*F22*m2*M4 + F21*F22*lfpf*lipf*pipf + 
       F21*F22*lfpi*lipf*pipf - 
       4*F11*F12*m2*M2*pipf - 
       6*F12*F21*m2*M2*pipf - 
       6*F11*F22*m2*M2*pipf - 
       6*F21*F22*m2*M2*pipf + 
       F21*F22*m2*pow(pipf,2) - 
       lilf*((2*F11 + F21)*M2 - F21*pipf)*
       ((2*F12 + F22)*M2 - F22*pipf) + 
       lipi*(lfpf*((2*F11 + F21)*(2*F12 + F22)*M2 + 
       F21*F22*pipf) + 
       lfpi*(-((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*
       M2) + F21*F22*pipf))) + 
       kflf*(-((F10 + F20)*(F12*F21 - F11*F22)*pow(lipi,2)*M2) + 
       lipf*lipi*((-(F12*F20*F21) + F11*F20*F22 + 
       F10*(2*F11 + F21)*(2*F12 + F22))*M2 + 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) + 
       lipf*(F10*((2*F11 + F21)*M2 - F21*pipf)*
       ((2*F12 + F22)*M2 - F22*pipf) - 
       lipf*(F10*(F12*F21 + 3*(F11 + F21)*F22)*M2 - 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf)))) + 
       kfpf*(kflf*(F20*F21*F22*pow(lipf,2) + 
       (F10 + F20)*F21*F22*pow(lipi,2) - 
       2*(F11 + F21)*F22*m2*(F10*M2 + F20*pipf) + 
       lipi*((-2*F11*F20*F22 + F10*F21*F22)*lipf - 
       (F12*F20*F21 + 2*F10*F11*F22 - F11*F20*F22 + 
       F10*F21*F22)*M2 + F10*F21*F22*pipf) + 
       lipf*((F21*(3*F12*F20 - 2*F10*F22 + 2*F20*F22) + 
       F11*(4*F12*F20 - 2*F10*F22 + 3*F20*F22))*M2 - 
       2*F20*(F11 + F21)*F22*pipf)) + 
       2*(F10*F21*F22*lfpi*pow(lipf,2) + 
       (F20*(-(F12*F21) + F11*F22)*lfpf + F10*F21*F22*lfpi)*
       pow(lipi,2) - 4*F10*F12*F21*lfpf*lipf*M2 - 
       4*F10*F11*F22*lfpf*lipf*M2 - 
       6*F10*F21*F22*lfpf*lipf*M2 + 
       8*F10*F11*F12*lfpi*lipf*M2 + 
       4*F10*F12*F21*lfpi*lipf*M2 + 
       4*F10*F11*F22*lfpi*lipf*M2 + 
       2*F10*F21*F22*lfpi*lipf*M2 - 
       F10*F12*F21*lfpf*m2*M2 - 
       F10*F11*F22*lfpf*m2*M2 - 
       F10*F21*F22*lfpf*m2*M2 + 
       4*F10*F11*F12*lfpi*m2*M2 + 
       3*F10*F12*F21*lfpi*m2*M2 - 
       F10*F11*F22*lfpi*m2*M2 - 
       F10*F21*F22*lfpi*m2*M2 - 
       6*F10*F11*F22*lipf*m2*M2 - 
       5*F10*F21*F22*lipf*m2*M2 + 
       8*F10*F11*F12*m2*M4 + 
       6*F10*F12*F21*m2*M4 + 
       6*F10*F11*F22*m2*M4 + 
       5*F10*F21*F22*m2*M4 + 
       2*F10*F21*F22*lfpf*lipf*pipf + 
       2*F10*F21*F22*lfpi*lipf*pipf + 
       F12*F20*F21*lfpf*m2*pipf - 
       F11*F20*F22*lfpf*m2*pipf + 
       F10*F21*F22*lfpf*m2*pipf + 
       F12*F20*F21*lfpi*m2*pipf - 
       F11*F20*F22*lfpi*m2*pipf + 
       F10*F21*F22*lfpi*m2*pipf + 
       F10*F21*F22*lipf*m2*pipf - 
       4*F10*F11*F12*m2*M2*pipf - 
       6*F10*F12*F21*m2*M2*pipf - 
       6*F10*F11*F22*m2*M2*pipf - 
       6*F10*F21*F22*m2*M2*pipf + 
       F10*F21*F22*m2*pow(pipf,2) - 
       lilf*((F12*F21 - F11*F22)*lipf*(F10*M2 - F20*pipf) + 
       F10*((2*F11 + F21)*M2 - F21*pipf)*
       ((2*F12 + F22)*M2 - F22*pipf)) + 
       lipi*(-2*F12*F20*F21*lfpf*lipf + 2*F11*F20*F22*lfpf*lipf + 
       4*F10*F11*F12*lfpf*M2 + 
       2*F10*F12*F21*lfpf*M2 + 
       2*F10*F11*F22*lfpf*M2 + 
       F10*F21*F22*lfpf*M2 - 
       4*F10*F11*F12*m2*M2 - 
       4*F10*F12*F21*m2*M2 - 
       2*F10*F11*F22*m2*M2 - 
       F10*F21*F22*m2*M2 + 
       F10*F21*F22*lfpf*pipf + F10*F21*F22*m2*pipf + 
       lilf*(F10*(4*F11*F12 + 3*F12*F21 + F11*F22)*M2 + 
       F20*(F12*F21 - F11*F22)*pipf) + 
       lfpi*((-(F12*F20*F21) + F11*F20*F22 + 2*F10*F21*F22)*lipf + 
       F10*(-((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2) + 
       F21*F22*pipf))))) + 
       kfpi*(-2*F20*(2*F11 + F21)*F22*pow(kfpf,2)*m2 + 
       kflf*(-((2*F11*F20 + (F10 + 2*F20)*F21)*F22*pow(lipf,2)) + 
       2*(F10 + F20)*(F11 + F21)*F22*m2*M2 + 
       lipi*(-(F10*F21*F22*lipf) + 
       (-(F12*F20*F21) + (F11*F20 + 2*F10*(F11 + F21))*F22)*M2) + 
       lipf*((F11*(-4*F12*F20 + 2*F10*F22 - 3*F20*F22) + 
       F21*(-3*F12*F20 + F10*F22 - 2*F20*F22))*M2 + 
       F21*(2*F12*F20 - F10*F22 + 2*F20*F22)*pipf)) + 
       kfpf*(-2*F12*F20*F21*lfpf*lipf + 2*F11*F20*F22*lfpf*lipf + 
       3*F10*F21*F22*lfpf*lipf + 2*F11*F20*F22*lfpi*lipf + 
       F10*F21*F22*lfpi*lipf + 2*F20*F21*F22*lfpi*lipf - 
       F20*F21*F22*lfpf*m2 - F20*F21*F22*lfpi*m2 + 
       2*F11*F20*F22*lipf*m2 - 
       F10*F21*F22*lipf*m2 + 
       2*F20*F21*F22*lipf*m2 + 
       F21*lipi*(2*F12*F20*lfpf + F10*F22*lfpf + 2*F20*F22*lfpf - 
       F10*F22*lfpi - (F10 + F20)*F22*lilf + F10*F22*m2) - 
       8*F10*F11*F12*m2*M2 - 
       8*F11*F12*F20*m2*M2 - 
       4*F10*F12*F21*m2*M2 - 
       6*F12*F20*F21*m2*M2 - 
       16*F10*F11*F22*m2*M2 - 
       6*F11*F20*F22*m2*M2 - 
       12*F10*F21*F22*m2*M2 - 
       4*F20*F21*F22*m2*M2 + 
       lilf*((2*F11*F20 + (F10 + F20)*F21)*F22*lipf + 
       2*((F21*(2*F12*F20 + 3*F10*F22 + F20*F22) + 
       F11*(2*F12*F20 + 4*F10*F22 + F20*F22))*M2 - 
       F21*(F12*F20 + (F10 + F20)*F22)*pipf))) - 
       2*(((-2*F12*F20*F21 + 2*F11*F20*F22 + F10*F21*F22)*lfpf + 
       (-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*lfpi)*pow(lipf,2) + 
       (F10 + F20)*(F12*F21 - F11*F22)*(lfpf + lfpi)*m2*M2 + 
       lipi*(((-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*lfpf + 
       F10*F21*F22*lfpi)*lipf + 
       (F10 + F20)*(F12*F21 - F11*F22)*lilf*M2) - 
       lipf*(-(lilf*((F10*F12*F21 + F12*F20*F21 + 3*F10*F11*F22 - 
       F11*F20*F22 + 2*F10*F21*F22)*M2 - 
       2*F10*F21*F22*pipf)) + 
       F10*(((2*F11 + F21)*(2*F12 + F22)*lfpf + 
       2*(2*F11*F12 + 2*F12*F21 + 4*F11*F22 + 3*F21*F22)*m2)*M2 + 
       F21*F22*(lfpf - 2*m2)*pipf + 
       lfpi*(-((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2) + 
       F21*F22*pipf)))))))/(2.*kfpf*M2);

res4 = (pow(kfli,2)*(4*F20*(F11 + F21)*F22*kflf*(M2 - pipf) + 
       lfpf*(-(F20*F21*F22*kfpf) + F20*(2*F11 + F21)*F22*kfpi + 
       8*F11*F12*F20*M2 + 6*F12*F20*F21*M2 + 
       4*F10*F11*F22*M2 + 6*F11*F20*F22*M2 + 
       2*F10*F21*F22*M2 + 4*F20*F21*F22*M2 - 
       2*F12*F20*F21*pipf + 2*F11*F20*F22*pipf - 2*F10*F21*F22*pipf) + 
       lfpi*(F20*(2*F11 + F21)*F22*kfpf - F20*F21*F22*kfpi - 
       2*((F12*F20*F21 + (-(F10*(2*F11 + F21)) + F20*(F11 + 2*F21))*F22)*M2 + 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf))) + 
       kfli*(2*F10*F21*F22*pow(kfpf,2)*lfpi + 
       F20*F21*F22*pow(kfpf,2)*lilf + 2*F20*F21*F22*kfpf*lfpf*lipf + 
       2*F12*F20*F21*kfpf*lfpi*lipf - 2*F11*F20*F22*kfpf*lfpi*lipf + 
       2*F20*F21*F22*kfpf*lfpi*lipf + 2*F12*F20*F21*kfpf*lfpf*lipi - 
       2*F11*F20*F22*kfpf*lfpf*lipi - 2*F10*F21*F22*kfpf*lfpf*lipi + 
       2*F20*F21*F22*kfpf*lfpf*lipi + 4*F12*F20*F21*kfpf*lfpi*lipi - 
       8*F11*F20*F22*kfpf*lfpi*lipi - 2*F10*F21*F22*kfpf*lfpi*lipi - 
       2*F20*F21*F22*kfpf*lfpi*lipi - 
       F20*F21*F22*pow(kfpf,2)*m2 + 
       F22*pow(kfpi,2)*(12*F11*F20*lfpf - 2*F10*F21*lfpf + 
       8*F20*F21*lfpf - 4*F20*F21*lfpi + F20*F21*lilf - 
       F20*F21*m2) + 8*F10*F12*F21*kfpf*lfpf*M2 + 
       6*F12*F20*F21*kfpf*lfpf*M2 - 
       4*F10*F11*F22*kfpf*lfpf*M2 - 
       2*F11*F20*F22*kfpf*lfpf*M2 + 
       6*F10*F21*F22*kfpf*lfpf*M2 + 
       4*F20*F21*F22*kfpf*lfpf*M2 + 
       8*F10*F12*F21*kfpf*lfpi*M2 + 
       2*F12*F20*F21*kfpf*lfpi*M2 + 
       4*F10*F11*F22*kfpf*lfpi*M2 - 
       2*F11*F20*F22*kfpf*lfpi*M2 + 
       8*F10*F21*F22*kfpf*lfpi*M2 - 
       4*F10*F11*F22*kfpf*lilf*M2 + 
       4*F11*F20*F22*kfpf*lilf*M2 - 
       2*F10*F21*F22*kfpf*lilf*M2 + 
       2*F20*F21*F22*kfpf*lilf*M2 - 
       12*F10*F12*F21*lfpf*lipf*M2 - 
       4*F12*F20*F21*lfpf*lipf*M2 - 
       4*F10*F11*F22*lfpf*lipf*M2 + 
       4*F11*F20*F22*lfpf*lipf*M2 - 
       12*F10*F21*F22*lfpf*lipf*M2 + 
       16*F10*F11*F12*lfpi*lipf*M2 + 
       4*F10*F12*F21*lfpi*lipf*M2 - 
       4*F12*F20*F21*lfpi*lipf*M2 + 
       12*F10*F11*F22*lfpi*lipf*M2 + 
       4*F11*F20*F22*lfpi*lipf*M2 + 
       4*F10*F21*F22*lfpi*lipf*M2 + 
       32*F10*F11*F12*lfpf*lipi*M2 + 
       20*F10*F12*F21*lfpf*lipi*M2 + 
       12*F10*F11*F22*lfpf*lipi*M2 + 
       8*F10*F21*F22*lfpf*lipi*M2 - 
       12*F10*F12*F21*lfpi*lipi*M2 - 
       20*F10*F11*F22*lfpi*lipi*M2 - 
       24*F10*F21*F22*lfpi*lipi*M2 - 
       8*F11*F12*F20*kfpf*m2*M2 - 
       8*F12*F20*F21*kfpf*m2*M2 + 
       8*F10*F11*F22*kfpf*m2*M2 - 
       8*F11*F20*F22*kfpf*m2*M2 + 
       6*F10*F21*F22*kfpf*m2*M2 - 
       6*F20*F21*F22*kfpf*m2*M2 - 
       16*F10*F11*F12*lfpf*M4 - 8*F10*F12*F21*lfpf*M4 - 
       8*F10*F11*F22*lfpf*M4 - 4*F10*F21*F22*lfpf*M4 + 
       16*F10*F11*F12*lfpi*M4 + 16*F10*F12*F21*lfpi*M4 + 
       16*F10*F11*F22*lfpi*M4 + 16*F10*F21*F22*lfpi*M4 - 
       16*F10*F11*F12*lilf*M4 - 8*F10*F12*F21*lilf*M4 - 
       8*F10*F11*F22*lilf*M4 - 4*F10*F21*F22*lilf*M4 + 
       32*F10*F11*F12*m2*M4 + 
       24*F10*F12*F21*m2*M4 + 
       24*F10*F11*F22*m2*M4 + 
       20*F10*F21*F22*m2*M4 + 
       2*F10*F21*F22*kfpf*lfpf*pipf - 4*F12*F20*F21*kfpf*lfpi*pipf - 
       4*F10*F21*F22*kfpf*lfpi*pipf - 4*F20*F21*F22*kfpf*lfpi*pipf + 
       2*F10*F21*F22*kfpf*lilf*pipf - 2*F20*F21*F22*kfpf*lilf*pipf + 
       4*F10*F21*F22*lfpf*lipf*pipf + 4*F10*F21*F22*lfpi*lipf*pipf + 
       4*F12*F20*F21*lfpf*lipi*pipf - 4*F11*F20*F22*lfpf*lipi*pipf + 
       8*F10*F21*F22*lfpf*lipi*pipf + 4*F12*F20*F21*lfpi*lipi*pipf - 
       4*F11*F20*F22*lfpi*lipi*pipf + 8*F10*F21*F22*lfpi*lipi*pipf - 
       2*F10*F21*F22*kfpf*m2*pipf + 
       2*F20*F21*F22*kfpf*m2*pipf + 
       8*F10*F12*F21*lfpf*M2*pipf + 
       8*F10*F11*F22*lfpf*M2*pipf + 
       8*F10*F21*F22*lfpf*M2*pipf - 
       16*F10*F11*F12*lfpi*M2*pipf - 
       16*F10*F12*F21*lfpi*M2*pipf - 
       16*F10*F11*F22*lfpi*M2*pipf - 
       16*F10*F21*F22*lfpi*M2*pipf + 
       8*F10*F12*F21*lilf*M2*pipf + 
       8*F10*F11*F22*lilf*M2*pipf + 
       8*F10*F21*F22*lilf*M2*pipf - 
       16*F10*F11*F12*m2*M2*pipf - 
       24*F10*F12*F21*m2*M2*pipf - 
       24*F10*F11*F22*m2*M2*pipf - 
       24*F10*F21*F22*m2*M2*pipf - 
       4*F10*F21*F22*lfpf*pow(pipf,2) - 
       4*F10*F21*F22*lilf*pow(pipf,2) + 
       4*F10*F21*F22*m2*pow(pipf,2) + 
       kfpi*(-(F22*kflf*(-(F20*(2*F11 + F21)*lipf) + F20*F21*lipi + 
       8*F10*F11*M2 + 8*F11*F20*M2 + 
       6*F10*F21*M2 + 4*F20*F21*M2 + 
       4*F11*F20*pipf - 2*F10*F21*pipf)) - 
       2*(2*F12*F20*F21*lfpf*lipf + 2*F11*F20*F22*lfpf*lipf - 
       F10*F21*F22*lfpf*lipf + 2*F20*F21*F22*lfpf*lipf + 
       F12*F20*F21*lfpi*lipf - F11*F20*F22*lfpi*lipf - 
       F10*F21*F22*lfpi*lipf - 2*F20*F21*F22*lfpi*lipf + 
       F12*F20*F21*lfpf*lipi + 5*F11*F20*F22*lfpf*lipi + 
       4*F20*F21*F22*lfpf*lipi - 2*F20*F21*F22*lfpi*lipi + 
       F22*kfpf*(F10*F21*lfpf + 2*F20*F21*lfpf - 
       (2*F11*F20 + F10*F21)*lfpi + F20*(2*F11 + F21)*lilf - 
       2*F11*F20*m2 - F20*F21*m2) + 
       16*F10*F11*F12*lfpf*M2 + 
       4*F11*F12*F20*lfpf*M2 + 
       12*F10*F12*F21*lfpf*M2 + 
       3*F12*F20*F21*lfpf*M2 + 
       2*F10*F11*F22*lfpf*M2 + 
       3*F11*F20*F22*lfpf*M2 + 
       2*F10*F21*F22*lfpf*M2 + 
       2*F20*F21*F22*lfpf*M2 - 
       4*F11*F12*F20*lfpi*M2 - 
       4*F10*F12*F21*lfpi*M2 - 
       5*F12*F20*F21*lfpi*M2 - 
       2*F10*F11*F22*lfpi*M2 - 
       3*F11*F20*F22*lfpi*M2 - 
       3*F10*F21*F22*lfpi*M2 - 
       4*F20*F21*F22*lfpi*M2 - 
       8*F11*F12*F20*m2*M2 - 
       6*F12*F20*F21*m2*M2 - 
       6*F11*F20*F22*m2*M2 + 
       F10*F21*F22*m2*M2 - 
       5*F20*F21*F22*m2*M2 + 
       2*F12*F20*F21*lfpf*pipf - 2*F11*F20*F22*lfpf*pipf + 
       4*F10*F21*F22*lfpf*pipf + 2*F12*F20*F21*lfpi*pipf + 
       3*F10*F21*F22*lfpi*pipf + 2*F20*F21*F22*lfpi*pipf + 
       2*F12*F20*F21*m2*pipf + 
       2*F11*F20*F22*m2*pipf + 
       F10*F21*F22*m2*pipf + 
       3*F20*F21*F22*m2*pipf + 
       (2*F12*F20 + (F10 + F20)*F22)*lilf*((2*F11 + F21)*M2 - F21*pipf))) + 
       kflf*(F22*kfpf*(-(F20*F21*lipf) + F20*(2*F11 + F21)*lipi + 
       8*F10*F11*M2 + 4*F11*F20*M2 + 
       6*F10*F21*M2 + 4*F20*F21*M2 - 2*F10*F21*pipf) + 
       2*(2*((2*F11 + F21)*M2 - F21*pipf)*
       ((F20*(F12 + F22) + F10*(2*F12 + F22))*M2 - 
       (F12*F20 + (F10 + F20)*F22)*pipf) + 
       lipf*((F21*(3*F12*F20 - 3*F10*F22 - 2*F20*F22) + 
       F11*(4*F12*F20 - (2*F10 + F20)*F22))*M2 - 
       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) - lipi*((F12*F20*F21 + 
       (-6*F10*F11 + F11*F20 - 5*F10*F21 + 2*F20*F21)*F22)*M2 + (F12*F20*F21 + 
       (-5*F11*F20 + F10*F21 - 4*F20*F21)*F22)*pipf)))) + 
       2*(-(F10*F21*F22*pow(kfpf,2)*lfpf*lipi) + 
       2*F12*F20*F21*pow(kfpf,2)*lfpi*lipi - 
       F10*F21*F22*pow(kfpf,2)*lfpi*lipi + 
       2*F20*F21*F22*pow(kfpf,2)*lfpi*lipi + 
       F10*F21*F22*pow(kfpf,2)*lilf*lipi - 
       F10*F21*F22*kflf*kfpf*lipf*lipi + 
       2*F10*F21*F22*kfpf*lfpf*lipf*lipi - 
       2*F12*F20*F21*kfpf*lfpi*lipf*lipi + 
       2*F11*F20*F22*kfpf*lfpi*lipf*lipi + 
       2*F10*F21*F22*kfpf*lfpi*lipf*lipi - 
       2*F11*F20*F22*kflf*kfpf*pow(lipi,2) - 
       F10*F21*F22*kflf*kfpf*pow(lipi,2) - 
       2*F20*F21*F22*kflf*kfpf*pow(lipi,2) - 
       2*F12*F20*F21*kfpf*lfpf*pow(lipi,2) + 
       2*F11*F20*F22*kfpf*lfpf*pow(lipi,2) + 
       2*F10*F21*F22*kfpf*lfpf*pow(lipi,2) - 
       4*F12*F20*F21*kfpf*lfpi*pow(lipi,2) + 
       4*F11*F20*F22*kfpf*lfpi*pow(lipi,2) + 
       2*F10*F21*F22*kfpf*lfpi*pow(lipi,2) + 
       F20*F21*F22*pow(kfpi,3)*m2 - 
       F10*F21*F22*pow(kfpf,2)*lipi*m2 - 
       F12*F20*F21*pow(kfpf,2)*lilf*M2 + 
       2*F10*F11*F22*pow(kfpf,2)*lilf*M2 + 
       F11*F20*F22*pow(kfpf,2)*lilf*M2 + 
       2*F10*F21*F22*pow(kfpf,2)*lilf*M2 + 
       F12*F20*F21*kflf*kfpf*lipf*M2 - 
       2*F10*F11*F22*kflf*kfpf*lipf*M2 - 
       F11*F20*F22*kflf*kfpf*lipf*M2 - 
       2*F10*F21*F22*kflf*kfpf*lipf*M2 + 
       2*F10*F12*F21*kfpf*lilf*lipf*M2 + 
       2*F12*F20*F21*kfpf*lilf*lipf*M2 - 
       2*F10*F11*F22*kfpf*lilf*lipf*M2 - 
       2*F11*F20*F22*kfpf*lilf*lipf*M2 - 
       2*F10*F12*F21*kflf*pow(lipf,2)*M2 - 
       2*F12*F20*F21*kflf*pow(lipf,2)*M2 + 
       2*F10*F11*F22*kflf*pow(lipf,2)*M2 + 
       2*F11*F20*F22*kflf*pow(lipf,2)*M2 + 
       4*F11*F12*F20*kflf*kfpf*lipi*M2 + 
       3*F12*F20*F21*kflf*kfpf*lipi*M2 - 
       2*F10*F11*F22*kflf*kfpf*lipi*M2 + 
       3*F11*F20*F22*kflf*kfpf*lipi*M2 - 
       F10*F21*F22*kflf*kfpf*lipi*M2 + 
       2*F20*F21*F22*kflf*kfpf*lipi*M2 - 
       4*F10*F12*F21*kfpf*lfpf*lipi*M2 - 
       4*F10*F11*F22*kfpf*lfpf*lipi*M2 - 
       6*F10*F21*F22*kfpf*lfpf*lipi*M2 + 
       8*F10*F11*F12*kfpf*lfpi*lipi*M2 + 
       4*F10*F12*F21*kfpf*lfpi*lipi*M2 + 
       4*F10*F11*F22*kfpf*lfpi*lipi*M2 + 
       2*F10*F21*F22*kfpf*lfpi*lipi*M2 + 
       2*F10*F12*F21*kfpf*lilf*lipi*M2 + 
       2*F12*F20*F21*kfpf*lilf*lipi*M2 + 
       6*F10*F11*F22*kfpf*lilf*lipi*M2 - 
       2*F11*F20*F22*kfpf*lilf*lipi*M2 + 
       4*F10*F21*F22*kfpf*lilf*lipi*M2 + 
       8*F10*F11*F12*kflf*lipf*lipi*M2 + 
       4*F10*F12*F21*kflf*lipf*lipi*M2 - 
       2*F12*F20*F21*kflf*lipf*lipi*M2 + 
       4*F10*F11*F22*kflf*lipf*lipi*M2 + 
       2*F11*F20*F22*kflf*lipf*lipi*M2 + 
       2*F10*F21*F22*kflf*lipf*lipi*M2 + 
       8*F10*F12*F21*lfpf*lipf*lipi*M2 + 
       8*F10*F11*F22*lfpf*lipf*lipi*M2 + 
       12*F10*F21*F22*lfpf*lipf*lipi*M2 - 
       16*F10*F11*F12*lfpi*lipf*lipi*M2 - 
       8*F10*F12*F21*lfpi*lipf*lipi*M2 - 
       8*F10*F11*F22*lfpi*lipf*lipi*M2 - 
       4*F10*F21*F22*lfpi*lipf*lipi*M2 - 
       2*F10*F12*F21*kflf*pow(lipi,2)*M2 - 
       6*F10*F11*F22*kflf*pow(lipi,2)*M2 - 
       6*F10*F21*F22*kflf*pow(lipi,2)*M2 - 
       16*F10*F11*F12*lfpf*pow(lipi,2)*M2 - 
       8*F10*F12*F21*lfpf*pow(lipi,2)*M2 - 
       8*F10*F11*F22*lfpf*pow(lipi,2)*M2 - 
       4*F10*F21*F22*lfpf*pow(lipi,2)*M2 + 
       8*F10*F12*F21*lfpi*pow(lipi,2)*M2 + 
       8*F10*F11*F22*lfpi*pow(lipi,2)*M2 + 
       12*F10*F21*F22*lfpi*pow(lipi,2)*M2 + 
       2*F10*F11*F22*kflf*kfpf*m2*M2 + 
       2*F11*F20*F22*kflf*kfpf*m2*M2 + 
       2*F10*F21*F22*kflf*kfpf*m2*M2 + 
       2*F20*F21*F22*kflf*kfpf*m2*M2 - 
       2*F10*F12*F21*pow(kfpf,2)*m2*M2 - 
       F12*F20*F21*pow(kfpf,2)*m2*M2 - 
       2*F10*F11*F22*pow(kfpf,2)*m2*M2 - 
       F11*F20*F22*pow(kfpf,2)*m2*M2 - 
       4*F10*F21*F22*pow(kfpf,2)*m2*M2 - 
       2*F20*F21*F22*pow(kfpf,2)*m2*M2 + 
       2*F10*F12*F21*kfpf*lfpf*m2*M2 + 
       2*F12*F20*F21*kfpf*lfpf*m2*M2 - 
       2*F10*F11*F22*kfpf*lfpf*m2*M2 - 
       2*F11*F20*F22*kfpf*lfpf*m2*M2 + 
       2*F10*F12*F21*kfpf*lfpi*m2*M2 + 
       2*F12*F20*F21*kfpf*lfpi*m2*M2 - 
       2*F10*F11*F22*kfpf*lfpi*m2*M2 - 
       2*F11*F20*F22*kfpf*lfpi*m2*M2 - 
       8*F10*F11*F12*kfpf*lipi*m2*M2 - 
       8*F10*F12*F21*kfpf*lipi*m2*M2 - 
       16*F10*F11*F22*kfpf*lipi*m2*M2 - 
       12*F10*F21*F22*kfpf*lipi*m2*M2 - 
       8*F10*F11*F12*kflf*lipi*M4 - 
       4*F10*F12*F21*kflf*lipi*M4 - 
       4*F10*F11*F22*kflf*lipi*M4 - 
       2*F10*F21*F22*kflf*lipi*M4 + 
       16*F10*F11*F12*lilf*lipi*M4 + 
       8*F10*F12*F21*lilf*lipi*M4 + 
       8*F10*F11*F22*lilf*lipi*M4 + 
       4*F10*F21*F22*lilf*lipi*M4 - 
       32*F10*F11*F12*lipi*m2*M4 - 
       24*F10*F12*F21*lipi*m2*M4 - 
       24*F10*F11*F22*lipi*m2*M4 - 
       20*F10*F21*F22*lipi*m2*M4 - 
       2*F12*F20*F21*kflf*kfpf*lipi*pipf + 
       F10*F21*F22*kflf*kfpf*lipi*pipf - 
       2*F20*F21*F22*kflf*kfpf*lipi*pipf + 
       2*F10*F21*F22*kfpf*lfpf*lipi*pipf + 
       2*F10*F21*F22*kfpf*lfpi*lipi*pipf - 
       4*F10*F21*F22*kfpf*lilf*lipi*pipf + 
       2*F12*F20*F21*kflf*lipf*lipi*pipf - 
       2*F11*F20*F22*kflf*lipf*lipi*pipf + 
       2*F10*F21*F22*kflf*lipf*lipi*pipf - 
       4*F10*F21*F22*lfpf*lipf*lipi*pipf - 
       4*F10*F21*F22*lfpi*lipf*lipi*pipf + 
       2*F12*F20*F21*kflf*pow(lipi,2)*pipf - 
       2*F11*F20*F22*kflf*pow(lipi,2)*pipf + 
       2*F10*F21*F22*kflf*pow(lipi,2)*pipf - 
       4*F10*F21*F22*lfpf*pow(lipi,2)*pipf - 
       4*F10*F21*F22*lfpi*pow(lipi,2)*pipf + 
       4*F10*F21*F22*kfpf*lipi*m2*pipf + 
       4*F10*F12*F21*kflf*lipi*M2*pipf + 
       4*F10*F11*F22*kflf*lipi*M2*pipf + 
       4*F10*F21*F22*kflf*lipi*M2*pipf - 
       8*F10*F12*F21*lilf*lipi*M2*pipf - 
       8*F10*F11*F22*lilf*lipi*M2*pipf - 
       8*F10*F21*F22*lilf*lipi*M2*pipf + 
       16*F10*F11*F12*lipi*m2*M2*pipf + 
       24*F10*F12*F21*lipi*m2*M2*pipf + 
       24*F10*F11*F22*lipi*m2*M2*pipf + 
       24*F10*F21*F22*lipi*m2*M2*pipf - 
       2*F10*F21*F22*kflf*lipi*pow(pipf,2) + 
       4*F10*F21*F22*lilf*lipi*pow(pipf,2) - 
       4*F10*F21*F22*lipi*m2*pow(pipf,2) + 
       pow(kfpi,2)*(2*F11*F20*F22*lfpf*lipf + 3*F10*F21*F22*lfpf*lipf + 
       2*F20*F21*F22*lfpf*lipf - 2*F12*F20*F21*lfpi*lipf + 
       2*F11*F20*F22*lfpi*lipf + F10*F21*F22*lfpi*lipf + 
       2*F10*F21*F22*lfpf*lipi - 
       2*F20*(2*F11 + F21)*F22*kfpf*m2 + 
       4*F11*F20*F22*lfpf*m2 + 3*F20*F21*F22*lfpf*m2 - 
       F20*F21*F22*lfpi*m2 - 2*F11*F20*F22*lipf*m2 + 
       F10*F21*F22*lipf*m2 - 2*F20*F21*F22*lipf*m2 - 
       8*F11*F12*F20*m2*M2 - 
       2*F10*F12*F21*m2*M2 - 
       7*F12*F20*F21*m2*M2 - 
       14*F10*F11*F22*m2*M2 - 
       7*F11*F20*F22*m2*M2 - 
       12*F10*F21*F22*m2*M2 - 
       6*F20*F21*F22*m2*M2 + 
       2*F12*F20*F21*m2*pipf + 2*F11*F20*F22*m2*pipf + 
       4*F10*F21*F22*m2*pipf + 4*F20*F21*F22*m2*pipf + 
       lilf*((2*F11*F20 - F10*F21 + F20*F21)*F22*lipf - 
       F20*F21*F22*lipi + 4*F11*F12*F20*M2 + 
       3*F12*F20*F21*M2 + 6*F10*F11*F22*M2 + 
       3*F11*F20*F22*M2 + 4*F10*F21*F22*M2 + 
       2*F20*F21*F22*M2 - 2*F11*F20*F22*pipf - 
       2*F10*F21*F22*pipf - 2*F20*F21*F22*pipf)) + 
       kfpi*(F20*F21*F22*pow(kfpf,2)*m2 + 
       kflf*((F10 + F20)*F21*F22*pow(lipf,2) + 
       F20*F21*F22*pow(lipi,2) - 
       2*(F11 + F21)*F22*m2*(F10*M2 + F20*pipf) + 
       lipf*((F12*F20*F21 + (2*F10*F11 - F11*F20 + F10*F21)*F22)*
       M2 - F10*F21*F22*pipf) + 
       lipi*((-2*F11*F20*F22 + F10*F21*F22)*lipf - 
       (4*F11*F12*F20 + 3*F12*F20*F21 - 2*F10*F11*F22 + 
       3*F11*F20*F22 - 2*F10*F21*F22 + 2*F20*F21*F22)*M2 + 
       2*F20*(F11 + F21)*F22*pipf)) + 
       kfpf*(F10*F21*F22*lfpf*lipf - 2*F12*F20*F21*lfpi*lipf - 
       F10*F21*F22*lfpi*lipf - 2*F20*F21*F22*lfpi*lipf - 
       F20*F21*F22*lfpf*m2 - F20*F21*F22*lfpi*m2 + 
       F10*F21*F22*lipf*m2 + 
       lipi*((2*F12*F20*F21 - 2*F11*F20*F22 - 3*F10*F21*F22)*lfpi + 
       (2*F11*F20 + (F10 + F20)*F21)*F22*lilf + 
       F22*(-((2*F11*F20 + F10*F21 + 2*F20*F21)*lfpf) + 
       (2*F11*F20 - F10*F21 + 2*F20*F21)*m2)) + 
       8*F10*F11*F12*m2*M2 + 
       8*F11*F12*F20*m2*M2 + 
       4*F10*F12*F21*m2*M2 + 
       6*F12*F20*F21*m2*M2 + 
       16*F10*F11*F22*m2*M2 + 
       6*F11*F20*F22*m2*M2 + 
       12*F10*F21*F22*m2*M2 + 
       4*F20*F21*F22*m2*M2 - 
       lilf*((F10 + F20)*F21*F22*lipf + 
       2*((F21*(2*F12*F20 + 3*F10*F22 + F20*F22) + 
       F11*(2*F12*F20 + 4*F10*F22 + F20*F22))*M2 - 
       F21*(F12*F20 + (F10 + F20)*F22)*pipf))) + 
       2*(-(F10*F21*F22*lfpf*pow(lipf,2)) + 
       F12*F20*F21*lfpi*pow(lipf,2) - 
       F11*F20*F22*lfpi*pow(lipf,2) - 
       F10*F21*F22*lfpf*pow(lipi,2) - 
       2*F10*F12*F21*lfpf*lipf*M2 - 
       2*F10*F11*F22*lfpf*lipf*M2 - 
       3*F10*F21*F22*lfpf*lipf*M2 + 
       4*F10*F11*F12*lfpi*lipf*M2 + 
       2*F10*F12*F21*lfpi*lipf*M2 + 
       2*F10*F11*F22*lfpi*lipf*M2 + 
       F10*F21*F22*lfpi*lipf*M2 - 
       4*F10*F11*F12*lfpf*m2*M2 - 
       3*F10*F12*F21*lfpf*m2*M2 + 
       F10*F11*F22*lfpf*m2*M2 + 
       F10*F21*F22*lfpf*m2*M2 + 
       F10*F12*F21*lfpi*m2*M2 + 
       F10*F11*F22*lfpi*m2*M2 + 
       F10*F21*F22*lfpi*m2*M2 + 
       4*F10*F11*F12*lipf*m2*M2 + 
       4*F10*F12*F21*lipf*m2*M2 + 
       2*F10*F11*F22*lipf*m2*M2 + 
       F10*F21*F22*lipf*m2*M2 + 
       8*F10*F11*F12*m2*M4 + 
       6*F10*F12*F21*m2*M4 + 
       6*F10*F11*F22*m2*M4 + 
       5*F10*F21*F22*m2*M4 + 
       F10*F21*F22*lfpf*lipf*pipf + F10*F21*F22*lfpi*lipf*pipf - 
       F12*F20*F21*lfpf*m2*pipf + 
       F11*F20*F22*lfpf*m2*pipf - 
       F10*F21*F22*lfpf*m2*pipf - 
       F12*F20*F21*lfpi*m2*pipf + 
       F11*F20*F22*lfpi*m2*pipf - 
       F10*F21*F22*lfpi*m2*pipf - 
       F10*F21*F22*lipf*m2*pipf - 
       4*F10*F11*F12*m2*M2*pipf - 
       6*F10*F12*F21*m2*M2*pipf - 
       6*F10*F11*F22*m2*M2*pipf - 
       6*F10*F21*F22*m2*M2*pipf + 
       F10*F21*F22*m2*pow(pipf,2) - 
       lilf*(F10*((2*F11 + F21)*M2 - F21*pipf)*
       ((2*F12 + F22)*M2 - F22*pipf) + 
       lipf*(F10*(4*F11*F12 + 3*F12*F21 + F11*F22)*M2 + 
       F20*(F12*F21 - F11*F22)*pipf)) + 
       lipi*(F12*F20*F21*lfpf*lipf - F11*F20*F22*lfpf*lipf - 
       2*F10*F21*F22*lfpf*lipf + 8*F10*F11*F12*lfpf*M2 + 
       4*F10*F12*F21*lfpf*M2 + 
       4*F10*F11*F22*lfpf*M2 + 
       2*F10*F21*F22*lfpf*M2 + 
       6*F10*F11*F22*m2*M2 + 
       5*F10*F21*F22*m2*M2 + 
       2*F10*F21*F22*lfpf*pipf - F10*F21*F22*m2*pipf + 
       (F12*F21 - F11*F22)*lilf*(F10*M2 - F20*pipf) - 
       2*lfpi*((-(F12*F20*F21) + F11*F20*F22)*lipf + 
       F10*((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2 - 
       F21*F22*pipf)))))))/(2.*kfpi*M2);

return (Pow6(e)/(q_12*q_22))*((res1 - res2)/kflf - (res4 - res3)/kfli);
}

