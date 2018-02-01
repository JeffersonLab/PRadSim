#include <complex>

std::complex<double> IxA(std::complex<double>, std::complex<double>, std::complex<double>);
std::complex<double> IIxA(std::complex<double>, std::complex<double>, std::complex<double>, double);
std::complex<double> IIIxA(std::complex<double>, std::complex<double>, std::complex<double>);
std::complex<double> IVxA(std::complex<double>, std::complex<double>, std::complex<double>);

std::complex<double> IxB(std::complex<double>, std::complex<double>, std::complex<double>);
std::complex<double> IIxB(std::complex<double>, std::complex<double>, std::complex<double>, double);
std::complex<double> IIIxB(std::complex<double>, std::complex<double>, std::complex<double>);
std::complex<double> IVxB(std::complex<double>, std::complex<double>, std::complex<double>);

std::complex<double> ComplexPower(double x, double p)
{
  std::complex<double> a(x, 0);
  return std::polar(pow(x, p), std::arg(a) * p);
};

std::complex<double> ComplexPower(double x, int p)
{
  double pd = (double)p;
  return ComplexPower(x, pd);
};

std::complex<double> ComplexPower(std::complex<double> x, double p)
{
  x *= ((std::arg(x) == -M_PI) ? std::polar(1.0, 2 * M_PI) : 1.0);
  return std::polar(pow(std::abs(x), p), std::arg(x) * p);
};

std::complex<double> ComplexPower(std::complex<double> x, int p)
{
  double pd = (double)p;
  return ComplexPower(x, pd);
};

std::complex<double> ComplexLog(std::complex<double> x)
{
  x *= ((std::arg(x) == -M_PI) ? std::polar(1.0, 2 * M_PI) : 1.0);
  return log(x);
};

double PolyLog(std::complex<double> arg)
{
  return DiLog(std::real(arg));
}

double Born_Terms(double SD, double TD, double UD){
  return pow(TD,-2)*(2*TD*UD + pow(TD,2) + 2*pow(-2 + UD,2)) + 
   pow(UD,-2)*(2*TD*UD + 2*pow(-2 + TD,2) + pow(UD,2)) + 
   2*pow(TD,-1)*pow(UD,-1)*(-4 + pow(TD + UD,2));
}

std::complex<double> IxA(std::complex<double> SD, std::complex<double> TD, std::complex<double> UD){
  return (0.3183098861837907*alpha*ComplexPower(TD,-1)*
   (ComplexLog(0.5*(ComplexPower(4. - 1.*TD,0.5) + ComplexPower(-1.*TD,0.5)))*ComplexPower(4. - 1.*TD,-0.5)*
      ComplexPower(-1.*TD,-0.5)*(2.*(8. - 3.*TD)*UD*(-4. + TD + UD) + 
        (4. - 1.*TD)*(16. + 3.*ComplexPower(TD,2))) + 
     (2.*TD*UD + ComplexPower(TD,2) + 2.*ComplexPower(-2. + UD,2))*
      (-2. + (-2. + TD)*ComplexPower(4. - 1.*TD,-0.5)*ComplexPower(-1.*TD,-0.5)*
         (-1.*PolyLog(1. + 0.5*(-1. + ComplexPower(-1.*TD*ComplexPower(4. - 1.*TD,-1),0.5))) + 
           PolyLog(0.5*(1. - 1.*ComplexPower(-1.*TD*ComplexPower(4. - 1.*TD,-1),0.5))) - 
           0.5*ComplexPower(ComplexLog(1. + 0.5*(-1. + ComplexPower(-1.*TD*ComplexPower(4. - 1.*TD,-1),0.5))),
             2) + 0.5*ComplexPower(ComplexLog(0.5*(1. - 1.*ComplexPower(-1.*TD*ComplexPower(4. - 1.*TD,-1),0.5))),
             2)))));
}

std::complex<double> IIxA(std::complex<double> SD, std::complex<double> TD, std::complex<double> UD, double r){
  return (0.1061032953945969*alpha*ComplexPower(TD,-2)*
   (-1.6666666666666667*TD - 4.*ComplexPower(r,-1) + 
     ComplexLog(0.5*(ComplexPower(-1.*r*TD,0.5) + ComplexPower(4. - 1.*r*TD,0.5)))*
      (2.*TD + 4.*ComplexPower(r,-1))*ComplexPower(-1.*r*TD,-0.5)*ComplexPower(4. - 1.*r*TD,0.5))*
   (2.*TD*UD + ComplexPower(TD,2) + 2.*ComplexPower(-2. + UD,2)));
}

std::complex<double> IIIxA(std::complex<double> SD, std::complex<double> TD, std::complex<double> UD){
  return (0.3183098861837907*alpha*(ComplexLog(0.5*(ComplexPower(4. - 1.*SD,0.5) + ComplexPower(-1.*SD,0.5)))*
      ComplexPower(4. - 1.*SD,-0.5)*ComplexPower(-1.*SD,-0.5)*(8. - 6.*SD + SD*TD + ComplexPower(SD,2)) + 
     ComplexLog(-1.*TD)*(-1. - 0.5*UD + 4.*UD*ComplexPower(4. - 1.*TD,-1)) + 
     (2. - 1.*SD)*ComplexLog(-1.*TD)*ComplexLog(0.5*(ComplexPower(4. - 1.*SD,0.5) + ComplexPower(-1.*SD,0.5)))*
      ComplexPower(4. - 1.*SD,-0.5)*ComplexPower(-1.*SD,-0.5)*
      (-8. + 3.*TD + 6.*UD + 4.*ComplexPower(TD,-1)*ComplexPower(-2. + UD,2)) + 
     ComplexPower(4. - 1.*TD,-0.5)*ComplexPower(-1.*TD,-0.5)*
      (24. - 14.*TD - 24.*UD + 6.*TD*UD + 16.*UD*ComplexPower(4. - 1.*TD,-1) + 
        3.*ComplexPower(TD,2) + 4.*ComplexPower(UD,2))*
      (3.2898681336964524 + PolyLog(
         0.5*(2. - 1.*TD - 1.*ComplexPower(4. - 1.*TD,0.5)*ComplexPower(-1.*TD,0.5))) + 
        ComplexPower(ComplexLog(0.5*(ComplexPower(4. - 1.*TD,0.5) + ComplexPower(-1.*TD,0.5))),2)) + 
     0.5*(2. - 1.*SD)*(2.*SD + TD)*ComplexPower(4. - 1.*SD,-0.5)*ComplexPower(-1.*SD,-0.5)*
      (-1.*PolyLog(1. + 0.5*(-1. + ComplexPower(-1.*SD*ComplexPower(4. - 1.*SD,-1),0.5))) + 
        PolyLog(0.5*(1. - 1.*ComplexPower(-1.*SD*ComplexPower(4. - 1.*SD,-1),0.5))) - 
        0.5*ComplexPower(ComplexLog(1. + 0.5*(-1. + ComplexPower(-1.*SD*ComplexPower(4. - 1.*SD,-1),0.5))),2) + 
        0.5*ComplexPower(ComplexLog(0.5*(1. - 1.*ComplexPower(-1.*SD*ComplexPower(4. - 1.*SD,-1),0.5))),2))));
}

std::complex<double> IVxA(std::complex<double> SD, std::complex<double> TD, std::complex<double> UD){
  return (0.3183098861837907*alpha*(ComplexLog(-1.*TD)*
      (-3. + 0.5*SD + 4.*UD*ComplexPower(4. - 1.*TD,-1)) + 
     (-2. + UD)*ComplexLog(-1.*TD)*ComplexLog(0.5*(ComplexPower(4. - 1.*UD,0.5) + ComplexPower(-1.*UD,0.5)))*
      ComplexPower(4. - 1.*UD,-0.5)*(TD + 2.*UD + 4.*ComplexPower(TD,-1)*ComplexPower(-2. + UD,2))*
      ComplexPower(-1.*UD,-0.5) + ComplexLog(0.5*(ComplexPower(4. - 1.*UD,0.5) + ComplexPower(-1.*UD,0.5)))*
      ComplexPower(4. - 1.*UD,-0.5)*ComplexPower(-1.*UD,-0.5)*
      (-8. + 6.*UD - 1.*TD*UD - 1.*ComplexPower(UD,2)) + 
     ComplexPower(4. - 1.*TD,-0.5)*ComplexPower(-1.*TD,-0.5)*
      (-8. - 2.*TD + 8.*UD - 2.*TD*UD + 16.*UD*ComplexPower(4. - 1.*TD,-1) - 
        1.*ComplexPower(TD,2) - 4.*ComplexPower(UD,2))*
      (3.2898681336964524 + PolyLog(
         0.5*(2. - 1.*TD - 1.*ComplexPower(4. - 1.*TD,0.5)*ComplexPower(-1.*TD,0.5))) + 
        ComplexPower(ComplexLog(0.5*(ComplexPower(4. - 1.*TD,0.5) + ComplexPower(-1.*TD,0.5))),2)) + 
     0.5*(-2. + UD)*(TD + 2.*UD)*ComplexPower(4. - 1.*UD,-0.5)*ComplexPower(-1.*UD,-0.5)*
      (-1.*PolyLog(1. + 0.5*(-1. + ComplexPower(-1.*UD*ComplexPower(4. - 1.*UD,-1),0.5))) + 
        PolyLog(0.5*(1. - 1.*ComplexPower(-1.*UD*ComplexPower(4. - 1.*UD,-1),0.5))) - 
        0.5*ComplexPower(ComplexLog(1. + 0.5*(-1. + ComplexPower(-1.*UD*ComplexPower(4. - 1.*UD,-1),0.5))),2) + 
        0.5*ComplexPower(ComplexLog(0.5*(1. - 1.*ComplexPower(-1.*UD*ComplexPower(4. - 1.*UD,-1),0.5))),2))));
}

std::complex<double> IxB(std::complex<double> SD, std::complex<double> TD, std::complex<double> UD){
  return (0.3183098861837907*alpha*ComplexPower(TD,-1)*
   (ComplexLog(0.5*(ComplexPower(4. - 1.*TD,0.5) + ComplexPower(-1.*TD,0.5)))*ComplexPower(4. - 1.*TD,-0.5)*
      ComplexPower(-1.*TD,-0.5)*(-32. + 2.*ComplexPower(TD,2) - 4.*ComplexPower(UD,2) + 
        3.*(4. - 1.*TD)*ComplexPower(TD + UD,2)) + 
     (-4. + ComplexPower(TD + UD,2))*(-2. + 
        (-2. + TD)*ComplexPower(4. - 1.*TD,-0.5)*ComplexPower(-1.*TD,-0.5)*
         (-1.*PolyLog(1. + 0.5*(-1. + ComplexPower(-1.*TD*ComplexPower(4. - 1.*TD,-1),0.5))) + 
           PolyLog(0.5*(1. - 1.*ComplexPower(-1.*TD*ComplexPower(4. - 1.*TD,-1),0.5))) - 
           0.5*ComplexPower(ComplexLog(1. + 0.5*(-1. + ComplexPower(-1.*TD*ComplexPower(4. - 1.*TD,-1),0.5))),
             2) + 0.5*ComplexPower(ComplexLog(0.5*(1. - 1.*ComplexPower(-1.*TD*ComplexPower(4. - 1.*TD,-1),0.5))),
             2)))));
}

std::complex<double> IIxB(std::complex<double> SD, std::complex<double> TD, std::complex<double> UD, double r){
  return (0.1061032953945969*alpha*ComplexPower(TD,-2)*
   (-1.6666666666666667*TD - 4.*ComplexPower(r,-1) + 
     ComplexLog(0.5*(ComplexPower(-1.*r*TD,0.5) + ComplexPower(4. - 1.*r*TD,0.5)))*
      (2.*TD + 4.*ComplexPower(r,-1))*ComplexPower(-1.*r*TD,-0.5)*ComplexPower(4. - 1.*r*TD,0.5))*
   (-4. + ComplexPower(TD + UD,2)));
}

std::complex<double> IIIxB(std::complex<double> SD, std::complex<double> TD, std::complex<double> UD){
  return (0.3183098861837907*alpha*(2.*UD*ComplexLog(-1.*TD)*ComplexPower(4. - 1.*TD,-1) + 
     2.*ComplexLog(0.5*(ComplexPower(4. - 1.*SD,0.5) + ComplexPower(-1.*SD,0.5)))*ComplexPower(4. - 1.*SD,-0.5)*
      ComplexPower(-1.*SD,-0.5)*(2.*UD + ComplexLog(-1.*TD)*
         (-2.*SD - 1.*UD + (6. - 1.*SD)*ComplexPower(-2. + SD,2)*ComplexPower(TD,-1))) + 
     2.*ComplexPower(4. - 1.*TD,-0.5)*ComplexPower(-1.*TD,-0.5)*
      (-6. - 1.*TD - 2.*UD + 4.*UD*ComplexPower(4. - 1.*TD,-1) + ComplexPower(TD + UD,2))*
      (3.2898681336964524 + PolyLog(
         0.5*(2. - 1.*TD - 1.*ComplexPower(4. - 1.*TD,0.5)*ComplexPower(-1.*TD,0.5))) + 
        ComplexPower(ComplexLog(0.5*(ComplexPower(4. - 1.*TD,0.5) + ComplexPower(-1.*TD,0.5))),2)) + 
     (2.*SD + UD)*ComplexPower(4. - 1.*SD,-0.5)*ComplexPower(-1.*SD,-0.5)*
      (-1.*PolyLog(1. + 0.5*(-1. + ComplexPower(-1.*SD*ComplexPower(4. - 1.*SD,-1),0.5))) + 
        PolyLog(0.5*(1. - 1.*ComplexPower(-1.*SD*ComplexPower(4. - 1.*SD,-1),0.5))) - 
        0.5*ComplexPower(ComplexLog(1. + 0.5*(-1. + ComplexPower(-1.*SD*ComplexPower(4. - 1.*SD,-1),0.5))),2) + 
        0.5*ComplexPower(ComplexLog(0.5*(1. - 1.*ComplexPower(-1.*SD*ComplexPower(4. - 1.*SD,-1),0.5))),2))));
}

std::complex<double> IVxB(std::complex<double> SD, std::complex<double> TD, std::complex<double> UD){
  return (0.3183098861837907*alpha*(ComplexLog(-1.*TD)*
      (-1. + 0.5*SD + 2.*UD*ComplexPower(4. - 1.*TD,-1)) + 
     ComplexLog(0.5*(ComplexPower(4. - 1.*UD,0.5) + ComplexPower(-1.*UD,0.5)))*ComplexPower(4. - 1.*UD,-0.5)*
      ComplexPower(-1.*UD,-0.5)*((2. - 1.*TD)*(4. + UD) - 1.*ComplexPower(UD,2) + 
        ComplexLog(-1.*TD)*(-8. - 6.*UD + TD*UD + 
           2.*(2. + UD)*ComplexPower(TD,-1)*ComplexPower(-2. + UD,2) + 2.*ComplexPower(UD,2))) + 
     ComplexPower(4. - 1.*TD,-0.5)*ComplexPower(-1.*TD,-0.5)*
      (4. - 4.*UD - 2.*TD*UD + 8.*UD*ComplexPower(4. - 1.*TD,-1) - 1.*ComplexPower(TD,2) - 
        2.*ComplexPower(UD,2))*(3.2898681336964524 + 
        PolyLog(0.5*(2. - 1.*TD - 1.*ComplexPower(4. - 1.*TD,0.5)*ComplexPower(-1.*TD,0.5))) + 
        ComplexPower(ComplexLog(0.5*(ComplexPower(4. - 1.*TD,0.5) + ComplexPower(-1.*TD,0.5))),2)) + 
     0.5*ComplexPower(4. - 1.*UD,-0.5)*ComplexPower(-1.*UD,-0.5)*
      ((2. - 1.*TD)*(4. - 1.*UD) + 2.*ComplexPower(UD,2))*
      (-1.*PolyLog(1. + 0.5*(-1. + ComplexPower(-1.*UD*ComplexPower(4. - 1.*UD,-1),0.5))) + 
        PolyLog(0.5*(1. - 1.*ComplexPower(-1.*UD*ComplexPower(4. - 1.*UD,-1),0.5))) - 
        0.5*ComplexPower(ComplexLog(1. + 0.5*(-1. + ComplexPower(-1.*UD*ComplexPower(4. - 1.*UD,-1),0.5))),2) + 
        0.5*ComplexPower(ComplexLog(0.5*(1. - 1.*ComplexPower(-1.*UD*ComplexPower(4. - 1.*UD,-1),0.5))),2))));
}

double MollerLoop(double SD, double TD, double UD){
  double rVal = pow(m/mmu,2);
  double Line2 = 2.*std::real(\
    (IxA(SD,TD,UD) + IIxA(SD,TD,UD,1.) + IIxA(SD,TD,UD,rVal) + IIIxA(SD,TD,UD) + IVxA(SD,TD,UD))/TD+
    (IxB(SD,TD,UD) + IIxB(SD,TD,UD,1.) + IIxB(SD,TD,UD,rVal) + IIIxB(SD,TD,UD) + IVxB(SD,TD,UD))/UD+
    (IxA(SD,UD,TD) + IIxA(SD,UD,TD,1.) + IIxA(SD,UD,TD,rVal) + IIIxA(SD,UD,TD) + IVxA(SD,UD,TD))/UD+
    (IxB(SD,UD,TD) + IIxB(SD,UD,TD,1.) + IIxB(SD,UD,TD,rVal) + IIIxB(SD,UD,TD) + IVxB(SD,UD,TD))/TD);
  double Line1 = Born_Terms(SD,TD,UD);
  return Line2/Line1;
}
