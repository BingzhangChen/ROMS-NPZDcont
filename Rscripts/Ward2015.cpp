#include <Rcpp.h>
#include <math.h>
#include <cmath>
using namespace Rcpp;

//Eq. 2 in Ward 2015:
//
// [[Rcpp::export]]
double Picofrac(double Ctot, double Temp) {
    double Csm = 0.2414;
    double Ds  = 1.0; 
    double As  = 0.668;
    double bs  = .302;
    double Cs;
    Cs = log10(Csm)+log10(1.-exp(-Ds/Csm*Ctot))-As*exp(-bs*Temp);
    Cs = pow(10, Cs);
    return(Cs/Ctot);
}

// [[Rcpp::export]]
double Nanofrac(double Ctot, double Temp) {
    double Csm = 0.921;
    double Ds  = 1.0; 
    double As  = 0.129;
    double bs  = 0.173;
    double Cs;
    double pico;
    Cs = log10(Csm)+log10(1.-exp(-Ds/Csm*Ctot))-As*exp(-bs*Temp);
    Cs = pow(10, Cs)/Ctot;
    pico = Picofrac(Ctot, Temp);
    return(Cs-pico);
}
