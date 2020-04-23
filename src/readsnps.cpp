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
extern const int NUMBEROFBASES = 3;
// 0x7ffe is 32,767 or 2^16 - 1
// 0xfffe is 65,534 or 2^32 - 1
// 0x2710 is 10,000
extern const unsigned short USBASE[NUMBEROFBASES] = {
  0x7ffe, // Used for format 1.1
  0xfffe, // Used for format 1.2
  0x2710  // Used for all other formats
};

// Values the short integers are multiplied by to get dosage and genetic
// probabilities
extern const double DBASE[NUMBEROFBASES] = {
  1. / USBASE[0],
  1. / USBASE[1],
  1. / USBASE[2]
};

// Find base used by given format
double findbase(int format, int subformat) {
  if (format == 1) {
    if (subformat == 1)
      return DBASE[0];
    return DBASE[1];
  }
  return DBASE[2];
}

// Routine to open a binary dosage file and read its format
void readbdformat(Rcpp::StringVector &status,
                  Rcpp::StringVector &filename,
                  int &format,
                  int &subformat) {
  std::ifstream bdfile;
  const char bdoseheader[4] = { 'b', 'o', 's', 'e'}; 
  char header[8];

  bdfile.open(filename[0], std::ios_base::in | std::ios_base::binary);
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

// Read dosages from formats 1.1, 2.1, 3.1, 3.3, 4.1, 4.3
int readdosages1(Rcpp::StringVector &filename,
                 Rcpp::NumericVector &indices,
                 Rcpp::StringVector &status,
                 Rcpp::NumericVector &dosage,
                 char *readbuffer,
                 int numsub,
                 double dbase) {
  std::ifstream bdfile;
  unsigned short *us = (unsigned short *)readbuffer;
  
  bdfile.open(filename[0], std::ios_base::in | std::ios_base::binary);
  bdfile.seekg(indices[0]);

  bdfile.read(readbuffer, 2 * numsub);
  for (int i = 0; i < 10; ++i) {
    dosage[i] = us[i] * dbase;
    Rcpp::Rcout << us[i] << '\t' << dosage[i] << std::endl;
  }

  bdfile.close();
  return 0;
}

// Read dosages only from formats 1.2, 2.2
int readdosages2(Rcpp::StringVector &filename,
                 Rcpp::NumericVector &indices,
                 Rcpp::StringVector &status,
                 Rcpp::NumericVector &dosage,
                 char *readbuffer,
                 int numsub,
                 double dbase) {
  std::ifstream bdfile;
  unsigned short *us = (unsigned short *)readbuffer;
  double p1, p2;
  
  bdfile.open(filename[0], std::ios_base::in | std::ios_base::binary);
  bdfile.seekg(indices[0]);
  
  bdfile.read(readbuffer, 2 * numsub);
  for (int i = 0; i < numsub; ++i) {
    p1 = us[i] * dbase;
    p2 = us[i + numsub] * dbase;
    dosage[i] = p1 + p2 + p2;
    dosage[i] = (dosage[i] > 2.) ? 2. : dosage[i];
  }
  for (int i = 0; i < 10; ++i)
    Rcpp::Rcout << dosage[i] << std::endl;
  
  bdfile.close();
  return 0;
}

// Read dosages and genetype probabilities from
// format 1.2, 2.2
int readbdall1(Rcpp::StringVector &filename,
               Rcpp::NumericVector &indices,
               Rcpp::StringVector &status,
               Rcpp::NumericVector &dosage,
               Rcpp::NumericVector &p0,
               Rcpp::NumericVector &p1,
               Rcpp::NumericVector &p2,
               char *readbuffer,
               int numsub,
               double dbase) {
  std::ifstream bdfile;
  unsigned short *us = (unsigned short *)readbuffer;
  
  bdfile.open(filename[0], std::ios_base::in | std::ios_base::binary);
  bdfile.seekg(indices[0]);
  
  bdfile.read(readbuffer, 4 * numsub);
  for (int i = 0; i < numsub; ++i) {
    p1[i] = us[i] * dbase;
    p2[i] = us[i + numsub] * dbase;
    p0[i] = 1. - p1[i] - p2[i];
    p0[i] = (p0[i] < 0.) ? 0. : p0[i];
    dosage[i] = p1[i] + p2[i] + p2[i];
    dosage[i] = (dosage[i] > 2.) ? 2. : dosage[i];
  }
  for (int i = 0; i < 10; ++i)
    Rcpp::Rcout << std::setw(15) << p0[i] << std::setw(15) << p1[i] << std::setw(15) << p2[i] << std::setw(15) << dosage[i] << std::endl;

  bdfile.close();
  return 0;
}

// Routine to read SNPs from binary dosage files
// Returns status of reading
// [[Rcpp::export]]
Rcpp::List readsnpsc(Rcpp::StringVector &filename,
                     Rcpp::NumericVector &indices,
                     int numsub,
                     Rcpp::LogicalVector dosageonly) {
    std::ifstream bdfile;
    int format = 0;
    int subformat = 0;
    char *readbuffer = NULL;
    Rcpp::StringVector status(1);
    Rcpp::NumericVector dosage(numsub);
    Rcpp::NumericVector p0(numsub);
    Rcpp::NumericVector p1(numsub);
    Rcpp::NumericVector p2(numsub);
    double dbase;
    
//    writebderrorfiles();    
    status[0] = "Good";
    
    readbdformat(status, filename, format, subformat);
    if (status[0] != "Good")
      return Rcpp::List::create(Rcpp::Named("status") = status);

    dbase = findbase(format, subformat);
    if (format == 1 && subformat == 1) {
      readbuffer = new char[2 * numsub];
      readdosages1(filename, indices, status,
                   dosage, readbuffer, numsub, dbase);
    } else if (format == 1 && subformat == 2) {
      readbuffer = new char[4 * numsub];
      readbdall1(filename, indices, status,
                 dosage, p0, p1, p2,
                 readbuffer, numsub, dbase);
    } else if (format == 2 && subformat == 1) {
      readbuffer = new char[2 * numsub];
      readdosages1(filename, indices, status,
                   dosage, readbuffer, numsub, dbase);
    } else if (format == 2 && subformat == 2) {
      readbuffer = new char[4 * numsub];
      readbdall1(filename, indices, status,
                 dosage, p0, p1, p2,
                 readbuffer, numsub, dbase);
    }    
    if (readbuffer)
      delete [] readbuffer;
    
    return Rcpp::List::create(Rcpp::Named("status") = status);
}
