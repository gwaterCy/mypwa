#ifndef COMMON_TOOLS_H
#define COMMON_TOOLS_H

#include <TString.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <algorithm>
#include <iomanip>

double CosTheta(TVector3 v1, TVector3 v2);

void V2LV(TLorentzVector &p4, double a[4]);

TString Int2Str(int i);

TString Double2Str(double d);

TString Float2Str(Float_t d);

bool InRange(double a, double min, double max);

bool OutRange(double a, double min, double max);

double valid_digits(double sp, int nDigits);

#endif
