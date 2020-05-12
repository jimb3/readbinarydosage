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
  0x2710  // Used for all other formats (0x2710 = 10000)
};

// Values the short integers are multiplied by to get dosage and genetic
// probabilities
extern const double DBASE[NUMBEROFBASES] = {
  1. / USBASE[0],
  1. / USBASE[1],
  1. / USBASE[2]
};

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

// **************************************************************************//
//                                                                           //
//               Convert data read to dosages and probabilities              //
//                                                                           //
// **************************************************************************//

void convertbd1(arma::uvec &subindex,
                arma::ivec &readbuffer,
                int start,
                int size,
                arma::ivec &subbuffer,
                arma::mat::iterator &dit,
                double dbase,
                unsigned short andbits) {
  char *p = (char *)readbuffer.memptr() + start;
  arma::Col<unsigned short> ureadbuffer((unsigned short *)p,
                                        size, false, true);
  arma::Col<unsigned short> usubbuffer((unsigned short *)subbuffer.memptr(),
                                       subindex.size(), false, true);

  usubbuffer = ureadbuffer.elem(subindex);
  arma::Col<unsigned short>::iterator u;
  for(u = usubbuffer.begin(); u != usubbuffer.end(); ++dit, ++u) {
    *dit = (*u & andbits) * dbase;
//    Rcpp::Rcout << std::setw(12) << *dit;
  }
//  Rcpp::Rcout << std::endl;
}

void convertbd2(arma::uvec &subindex,
                arma::ivec &readbuffer,
                int start,
                int size,
                arma::ivec &subbuffer,
                arma::mat::iterator &dit,
                double dbase,
                unsigned short andbits) {
}

// **************************************************************************//
//                                                                           //
//                          Reading in Blocks                                //
//                                                                           //
// **************************************************************************//


// Read the block
void readblock(Rcpp::StringVector &filename,
               double blkloc,
               double blksize,
               arma::ivec &readbuffer) {
  std::ifstream bdfile;
  
  bdfile.open(filename[0], std::ios_base::in | std::ios_base::binary);
  bdfile.seekg(blkloc);
  bdfile.read((char *)readbuffer.memptr(), blksize);
  bdfile.close();
}

// **************************************************************************//
//                                                                           //
//                          Reading in SNPs                                  //
//                                                                           //
// **************************************************************************//

// Routine to read SNPs from binary dosage files
// Returns status of reading
// [[Rcpp::export]]
Rcpp::List readsnpsc(Rcpp::StringVector &filename,
                     arma::uvec &subjects,
                     arma::uvec &snps,
                     int nsub,
                     arma::vec &snploc,
                     arma::vec &snpbytes,
                     arma::vec &blksnps,
                     arma::vec &firstsnp,
                     arma::vec &blkloc,
                     arma::vec &blkbytes,
                     Rcpp::LogicalVector &dosageonly) {
  std::ifstream bdfile;
  Rcpp::StringVector status(1);
  arma::uvec subindex;
  arma::uvec snpindex;
  arma::uvec block;
  arma::ivec readbuffer;
  arma::ivec subbuffer;
  arma::mat dosage;
  int format = 0;
  int subformat = 0;
  int maxblock;
  double dbase;
  unsigned short andbits;
  unsigned int currentblock;
  void (*convertbd)(arma::uvec &,
                    arma::ivec &,
                    int,
                    int,
                    arma::ivec &,
                    arma::mat::iterator &,
                    double,
                    unsigned short);

  status[0] = "Good";
  
  readbdformat(status, filename, format, subformat);
  if (status[0] != "Good")
    return Rcpp::List::create(Rcpp::Named("status") = status);

  maxblock = (int)blkbytes.max();
  readbuffer.set_size((maxblock + sizeof(int) - 1) / sizeof(int));
  
  subbuffer.zeros((2 * subjects.size() + sizeof(int) - 1) / sizeof(int));
  subindex = subjects - 1;
  
  snpindex = snps - 1;
  block = snpindex / blksnps[0];
  
  dosage.zeros(subindex.size(), snpindex.size());

  dbase = (format == 1) ? ((subformat == 1) ? DBASE[0] : DBASE[1]) : DBASE[2];
  andbits = (format == 1) ? 0xffff : 0x7fff;
  if (subformat == 1 || subformat == 3 || dosageonly[0] == TRUE)
    convertbd = convertbd1;
  else
    convertbd = convertbd2;

  currentblock = 0xffffffff;
  arma::uvec::iterator iblk = block.begin();
  arma::uvec::iterator isnp = snpindex.begin();
  arma::mat::iterator dit = dosage.begin();
  for (; iblk != block.end(); ++iblk, ++isnp) {
    if (currentblock != *iblk) {
      readblock(filename, blkloc[*iblk], blkbytes[*iblk], readbuffer);
      currentblock = *iblk;
    }
    convertbd(subindex, readbuffer,
              (int)(snploc[*isnp] - blkloc[*iblk]),
              nsub, subbuffer,
              dit, dbase, andbits);
  }
  
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("filename") = filename,
                            Rcpp::Named("subindex") = subindex,
                            Rcpp::Named("snpindex") = snpindex,
                            Rcpp::Named("format") = format,
                            Rcpp::Named("subformat") = subformat,
                            Rcpp::Named("snploc") = snploc,
                            Rcpp::Named("snpbytes") = snpbytes,
                            Rcpp::Named("block") = block,
                            Rcpp::Named("blksnps") = blksnps,
                            Rcpp::Named("firstsnp") = firstsnp,
                            Rcpp::Named("blkloc") = blkloc,
                            Rcpp::Named("blkbytes") = blkbytes,
                            Rcpp::Named("readbuffer") = readbuffer,
                            Rcpp::Named("subbuffer") = subbuffer,
                            Rcpp::Named("dosage") = dosage);
//  return readdosage(bdfile, subindices, dosageonly, snps, nsub,
//                    indices, dsize, blkbytes, readbuffer,
//                    subbuffer,
//                    status, format, subformat, dbase);
}
