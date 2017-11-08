#include <Rcpp.h>
#include <math.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]

double vert_int(double Ea, double tC){
  //DESCRIPTION:
  //The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
  // tC: in situ temperature
  // Tr: reference temperature
  //
  //INPUT PARAMETERS:
  // boltzman constant constant [ eV /K ]
  double kb = 8.62E-5;
  double Tr = 15.;
  double Y;
  Y = exp(-(Ea/kb)*(1/(273.15 + tC)-1/(273.15 + Tr)));
  return Y;
}
