// Utilities in this package
// Author: Hanyu WEI,   March 17, 2017

#ifndef UTIL_H
#define UTIL_H
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TROOT.h"

// interactive initilization of a matrix/vector
TMatrixD Matrix(Int_t row, Int_t column);
TVectorD Vector(Int_t row);

// print out two matrices M1, M2, and their product M1*M2
void MatrixMatrix(TMatrixD M1, TMatrixD M2);

// print out matrix M and vector V, and their product M*V
void MatrixVector(TMatrixD M, TVectorD V);

// test/validation on SVD decomposition implemented by ROOT
void SVD(TMatrixD M);

// converters for histrogram <---> matrix/vector
// former is input, latter is output

// rowcolumn is TRUE, histo(i+1, j+1) = matrix(i, j)
// rowcolumn is FALSE, histo(i+1, j+1) = matrix(j, i)
void H2M(const TH2D* histo, TMatrixD& mat, bool rowcolumn); 

void H2V(const TH1D* histo, TVectorD& vec);
void M2H(const TMatrixD mat, TH2D* histo);
void V2H(const TVectorD vec, TH1D* histo);

void CopyDir(TDirectory *source, bool blank_tree=false, bool verbose=true);
void CopyDir(TDirectory *source, TString TDirectory_name, bool blank_tree=false, bool verbose=true);
std::vector<TTree*>* CopyTrees(TDirectory *source, bool blank_tree=false, bool rename=false, TString TDirectory_name="", bool verbose=true);
std::vector<TTree*>* GetTrees(TDirectory *source,bool verbose=true);
#endif
