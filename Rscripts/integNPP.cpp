#include <Rcpp.h>
#include <math.h>
#include <cmath>
using namespace Rcpp;

Function sum("sum");
// [[Rcpp::export]]
SEXP integC(integer L, integer M, integer NMo, integer NZ, SEXP mask_, SEXP Depth_, SEXP NPPm_){
  //DESCRIPTION:
  //
  //INPUT PARAMETERS:
  //L: 
  //Transform NPPm and Depth into a vector:
  NumericVector  NPPm(clone(NPPm_));
  NumericVector Depth(clone(Depth_));
  NumericMatrix  mask(clone(mask_));

  int N = NPPm.size();
  int M = L + M + NMo;
  NumericVector NPPI(M);

  //Integrated NPP to be calculated (unit: mgC m-2 d-1)
  //for (int i = 0; i < N; i++) {
   for (int i =0; i < L; i++){
       for (int j =0; j < M;  j++){
            dep = Depth[i,j,] ;  //Depth profile
            w   = dep > -260;    //Pick up indices above 260 m
           for (int k = 0; k < 12; k++){
               //NPP depth profile
               npp         = NPPm[i,j,w,k];  
               NPPI[i,j,k] = NPPI[i,j,k] +
                              sum(npp*Hz[i,j,w])*mask[i,j] //Unit: mgC m-2 d-1
           }
       }
   }
  return NPPI;
} 

//return Rcpp::List::create(Rcpp::Named("vec") = someVector,
//                          Rcpp::Named("lst") = someList,
//                          Rcpp::Named("vec2") = someOtherVector);
