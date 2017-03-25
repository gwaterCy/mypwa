#ifndef DPF_PROPOGATOR_HH
#define DPF_PROPOGATOR_HH

#include <iostream>
#include <fstream>
#include <math.h>
#include "TNamed.h"

#if defined(USEROOT) || defined(__CINT__)
#include "RooComplex.h"
#else
#include "RooFitCore/RooComplex.hh"
#endif

#include "DPFPWAPoint.h"
#include "TComplex.h"

using namespace std;

class DPFPropogator {

public:
  DPFPropogator() {};

  virtual ~DPFPropogator() {}


  TComplex cro(Double_t sx, Double_t am1, Double_t am2)const;
  TComplex propogator980(Double_t mass, Double_t g11, Double_t g22,Double_t sx)const;
  TComplex pip(Double_t sx)const;
  TComplex propogator600(Double_t mass, Double_t b1, Double_t b2, Double_t b3, Double_t b4, Double_t b5, Double_t sx)const;
  TComplex propogator(Double_t mass,Double_t width,Double_t sx) const;
  TComplex propogator1270(Double_t mass,Double_t width,Double_t sx) const;


private:

  DPFPWAPoint* _dp;

//  ClassDef(DPFPropogator,1)
};

#endif
