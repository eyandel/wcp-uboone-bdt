#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

template <typename T>
class HistReader
{
public:
  HistReader(std::string fin, std::string key, double limit = -1.)
    : fFile(new TFile(fin.c_str(), "read")), fLimit(limit)
  {
      fHist = (T*)fFile->Get(key.c_str());
  }
  template <typename... dtypes>
  double operator()(dtypes... values) const {
      int binidx = fHist->FindBin(static_cast<double>(values)...);
      if(fLimit < 0.)
        return (double)(fHist->GetBinContent(binidx));
      else
        return (fHist->GetBinContent(binidx) < fLimit) ? (double)fHist->GetBinContent(binidx) : fLimit;
  }
private:
  T* fHist;
  TFile* fFile;
  double fLimit;
};

template <typename T>
class HistReaderFile
{
public:
  HistReaderFile(TFile* fin, std::string key, double limit = -1.)
    : fFile(fin), fLimit(limit)
  {
      fHist = (T*)fFile->Get(key.c_str());
  }
  template <typename... dtypes>
  double operator()(dtypes... values) const {
      int binidx = fHist->FindBin(static_cast<double>(values)...);
      if(fLimit < 0.)
        return (double)(fHist->GetBinContent(binidx));
      else
        return (fHist->GetBinContent(binidx) < fLimit) ? (double)fHist->GetBinContent(binidx) : fLimit;
  }
private:
  T* fHist;
  TFile* fFile;
  double fLimit;
};

template <typename T>
class GraphReader
{
public:
  GraphReader(std::string fin, std::string key)
    : fFile(new TFile(fin.c_str(), "read"))
  {
      fGraph = (T*)fFile->Get(key.c_str());
  }
  template <typename... dtypes>
  double operator()(dtypes... values) const {
      return (double)(fGraph->Eval(static_cast<double>(values)...));
  }
private:
  T* fGraph;
  TFile* fFile;
};
