#ifndef DPF_PWAPOINT
#define DPF_PWAPOINT

#include <iostream>
#include <math.h>
#include "assert.h"
#include "TString.h"
#include "DPFCoord.h"

using namespace std;
//using namespace DPFCoord;

class DPFPWAPoint {

public:

  DPFPWAPoint () {}
  inline ~DPFPWAPoint() {
  };
  DPFPWAPoint (Double_t, Double_t, Double_t, Double_t, Double_t, TString, TString);
  DPFPWAPoint (const DPFPWAPoint &other);

  Bool_t   valid (const Double_t dalX, const Double_t dalY) const;

  Double_t firstOfPair  (DPFCoord::Pair pair) const;
  Double_t secondOfPair (DPFCoord::Pair pair) const;
  Double_t notOfPair    (DPFCoord::Pair pair) const;

  Double_t firstOfPair  () const {return firstOfPair  (_dalitzX);}
  Double_t secondOfPair () const {return secondOfPair (_dalitzX);}
  Double_t notOfPair    () const {return notOfPair    (_dalitzX);}

  void print();

  DPFCoord::Index commonMember (const DPFCoord::Pair pair1, const DPFCoord::Pair pair2) const;
  DPFCoord::Index firstMember  (const DPFCoord::Pair pair)  const;
  DPFCoord::Index secondMember (const DPFCoord::Pair pair)  const;
  DPFCoord::Index notAMember   (const DPFCoord::Pair pair)  const;
  DPFCoord::Index notAMember   (const DPFCoord::Index i1,   const DPFCoord::Index i2) const;
  DPFCoord::Pair  otherPair    (const DPFCoord::Pair pair)  const;

  void qMaxMin(const Double_t q, Double_t &qmax, Double_t &qmin) const;
  void calBack(Double_t dalX, Double_t dalY,
       Double_t dalQ[], Double_t dalSQ[], Double_t pPi2[], Double_t pPi1[],
       Double_t a2[], Double_t a1[]) const;
  Double_t p2Cal(Double_t sx) const;
  Double_t p1Cal(Double_t sx) const;
  Double_t p0Cal(Double_t sx) const;
  Double_t calCos     (Double_t x, Double_t y) const;
  Double_t calY       (Double_t x, Double_t coshel) const;

  Double_t lowLimit (DPFCoord::Pair) const;
  Double_t highLimit(DPFCoord::Pair) const;
  void calcArea ();

  // define the Dalitz fitting plot: x vs y
  void definePlot (DPFCoord::Pair dalX, DPFCoord::Pair dalY) {
    _dalitzX = dalX;
    _dalitzY = dalY;
    if (dalX!=DPFCoord::AB && dalY!=DPFCoord::AB) _dalitzZ = DPFCoord::AB;
    if (dalX!=DPFCoord::BC && dalY!=DPFCoord::BC) _dalitzZ = DPFCoord::BC;
    if (dalX!=DPFCoord::CA && dalY!=DPFCoord::CA) _dalitzZ = DPFCoord::CA;
    calcArea();
  }

  DPFCoord::Pair _dalitzX;
  DPFCoord::Pair _dalitzY;
  DPFCoord::Pair _dalitzZ;

  Double_t _M, _m[4], _M2, _m2[4], _area, _sum2;
  TString _phspfile;
  TString _datafile;
  Double_t fm2Mmpi2,fm2Amds2,fm2Ampi2,mpi2Mmds2,D25Dfm2;
  Double_t xLower,xUpper,xDiff,yLower,yUpper,yDiff;
  Double_t fDel[4][4],fGel[4][4],E[4][4][4][4],G1[4][4][4][4],G3[4][4][4][4][4][4];

private:

//  ClassDef(DPFPWAPoint,0)
};

#endif
