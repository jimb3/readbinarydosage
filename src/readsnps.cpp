// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

/*
// Temporary routines to write files that will produce errors in
// readbdformat routine
void writebderrorfiles() {
  std::ofstream bdfile;
  const char error1file[6] = {'b', 'o', 's', 'e', 0x0, 0x0};
  const char error2file[8] = {'a', 'b', 'c', 'd', 0x0, 0x0, 0x0, 0x0};
  const char error3file[8] = {'b', 'o', 's', 'e', 0x1, 0x0, 0x0, 0x0};
  const char error4file[8] = {'b', 'o', 's', 'e', 0x0, 0x0, 0x0, 0x0};
  const char error5file[8] = {'b', 'o', 's', 'e', 0x0, 0x1, 0x0, 0x0};
  const char error6file[8] = {'b', 'o', 's', 'e', 0x0, 0x3, 0x0, 0x0};
  
  bdfile.open("error1.bdose", std::ios_base::out | std::ios_base::binary);
  bdfile.write(error1file, sizeof(error1file));
  bdfile.close();

  bdfile.open("error2.bdose", std::ios_base::out | std::ios_base::binary);
  bdfile.write(error2file, sizeof(error2file));
  bdfile.close();
  
  bdfile.open("error3.bdose", std::ios_base::out | std::ios_base::binary);
  bdfile.write(error3file, sizeof(error3file));
  bdfile.close();
  
  bdfile.open("error4.bdose", std::ios_base::out | std::ios_base::binary);
  bdfile.write(error4file, sizeof(error4file));
  bdfile.close();
  
  bdfile.open("error5.bdose", std::ios_base::out | std::ios_base::binary);
  bdfile.write(error5file, sizeof(error5file));
  bdfile.close();
  
  bdfile.open("error6.bdose", std::ios_base::out | std::ios_base::binary);
  bdfile.write(error6file, sizeof(error6file));
  bdfile.close();
}
*/

// Routine to open a binary dosage file and read its format
void readbdformat(Rcpp::StringVector &status,
                  Rcpp::StringVector &filenames,
                  int &format,
                  int &subformat) {
  std::ifstream bdfile;
  const char bdoseheader[4] = { 'b', 'o', 's', 'e'}; 
  char header[8];

  bdfile.open(filenames[0], std::ios_base::in | std::ios_base::binary);
  if (!bdfile.is_open()) {
    status = "Failed to open file";
    bdfile.close();
    return;
  }
// Error 1
  bdfile.read(header, 8);
  if (bdfile.fail()) {
    status = "Error reading binary dosage file header";
    bdfile.close();
    return;
  }
  bdfile.close();
// Error 2  
  if (memcmp(bdoseheader, header, 4)) {
    status = "Not a binary dosage file";
    return;
  }
// Error 3  
  if (header[4] != 0 || header[6] != 0) {
    status = "Error reading format";
    return;
  }
// Error 4  
  if (header[5] < 1 || header[5] > 4) {
    status = "Unknown format";
    return;
  }
// Error 5  
  if (header[5] < 3 && (header[7] < 1 || header[7] > 2)) {
    status = "Unknown subformat for format 1 or 2";
    return;
  }
// Error 6  
  if (header[5] > 2 && (header[7] < 1 || header[7] > 4)) {
    status = "Unknown subformat for format 3 or 4";
    return;
  }

  format = header[5];
  subformat = header[7];
}

// Routine to read SNPs from binary dosage files
// Returns status of reading
// [[Rcpp::export]]
Rcpp::List readsnpsc(Rcpp::StringVector &filenames) {
    std::ifstream bdfile;
    int format = 0;
    int subformat = 0;
    Rcpp::StringVector status(1);

//    writebderrorfiles();    
    status[0] = "Good";
    
    readbdformat(status, filenames, format, subformat);
    return Rcpp::List::create(Rcpp::Named("status") = status);
}
