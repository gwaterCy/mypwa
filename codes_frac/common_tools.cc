#ifndef COMMON_TOOLS_CC
#define COMMON_TOOLS_CC

#include <TString.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <sstream>
#include "common_tools.h"

using namespace std;

double CosTheta(TVector3 v1, TVector3 v2) {
    return v1.Dot(v2)/(v1.Mag()*v2.Mag());
}

void V2LV(TLorentzVector &p4, double a[4]) {
    p4 = TLorentzVector(a[1], a[2], a[3], a[0]);
}

TString Int2Str(int i) {
	stringstream ss;
	ss << i;
	return ss.str().c_str();
}
TString Double2Str(double d) {
	stringstream ss;
	ss << d;
	return ss.str().c_str();
}
TString Float2Str(Float_t d) {
	stringstream ss;
	ss << d;
	return ss.str().c_str();
}

bool InRange(double a, double min, double max) {
  return (a > min) && (a < max);
}

bool OutRange(double a, double min, double max) {
  return (a < min) || (a > max);
}

double valid_digits(double sp, int nDigits) {
    stringstream sst;
    sst << std::setprecision(nDigits) << sp;
    double rsp;
    sst >> rsp;
    cout << sp << " approx " << rsp << endl;
    return rsp;
}


#endif
