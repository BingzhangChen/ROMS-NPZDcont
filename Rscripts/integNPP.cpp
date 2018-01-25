#include <Rcpp.h>
#include <math.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP integC(int L, int M, int NMo, int NZ, SEXP Depth_, SEXP HZ_, SEXP NPPm_){
  //DESCRIPTION:
  //
  //INPUT PARAMETERS:
  //L: 
  //Transform NPPm and Depth into a vector:
  NumericVector  NPPm(clone(NPPm_));
  NumericVector Depth(clone(Depth_));
  NumericVector    Hz(clone(HZ_));

  int N = NPPm.size();
  int F = L * M * NMo;
  double npp;
  NumericVector NPPI(F);
  NumericVector dep(NZ);
  NumericVector  hz(NZ);

  if (N != (F * NZ)) {
      stop("Dimensions of NPPm do not match with input!");
  }
  //Integrated NPP to be calculated (unit: mgC m-2 d-1)
  //for (int i = 0; i < N; i++) {
   for (int j =0; j < M;  j++){
     for (int i =0; i < L; i++){
       for (int k =0; k < NZ; k++){
         dep[k] = Depth[k*L*M + j*L + i] ;  //Depth profile
          hz[k] =    Hz[k*L*M + j*L + i] ;  //Depth profile
         for (int t = 0; t < NMo; t++){
             //NPP depth profile
             if (dep[k] > -260.){
               npp = NPPm[t*L*M*NZ + k*L*M + j*L + i];  

               //Integrate through the water column
               NPPI[t*L*M + j*L + i]+= npp*hz[k]; //Unit: mgC m-2 d-1
             }
         }
       }
     }
   }
  return NPPI;
} 
//return Rcpp::List::create(Rcpp::Named("vec") = someVector,
//                          Rcpp::Named("lst") = someList,
//                          Rcpp::Named("vec2") = someOtherVector);
