#ifndef DPF_ANGULAR_HH
#define DPF_ANGULAR_HH

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

#include "PWA_PARAS.h"

using namespace std;
const Double_t pion_mass(0.139570180);
const Double_t kaon_mass(0.493677);
const Double_t psi_mass(3.686);

class DPFAngular {

    public:
        DPFAngular() {};

        virtual ~DPFAngular() {}

        void setdp(DPFPWAPoint *dp) { _dp = dp; };

        Double_t scalar(Double_t *a1,Double_t *a2)const;

        //calculate0p中对每个事例计算PWA_PARAS结构中的参数
        Double_t calculate0p(
                Double_t, Double_t, Double_t, Double_t,
                Double_t, Double_t, Double_t, Double_t,
                Double_t, Double_t, Double_t, Double_t,
                Double_t, Double_t, Double_t, Double_t,
                Double_t, Double_t, Double_t, Double_t,
                PWA_PARAS &) const;
    private:

        DPFPWAPoint* _dp;

        //  ClassDef(DPFAngular,1)
};

#endif
