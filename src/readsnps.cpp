// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// Routine to read SNPs from binary dosage files
// Currently does nothing but return a list with one value set to 1
// [[Rcpp::export]]
Rcpp::List readsnpsc() {
    int x = 1;
    
    return Rcpp::List::create(Rcpp::Named("x") = x);
}
