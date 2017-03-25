#ifndef DPF_PWAPDF
#define DPF_PWAPDF

#include <vector>
#include <iostream>
#include <fstream>
#include "TComplex.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TH1F.h"
#if defined(USEROOT) || defined(__CINT__)
#include "RooStringVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooComplex.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooRealConstant.h"
#else
#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooAbsReal.hh"
#include "RooFitCore/RooAbsData.hh"
#include "RooFitCore/RooFitResult.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooRealProxy.hh"
#include "RooFitCore/RooListProxy.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooComplex.hh"
#include "RooFitCore/RooRealConstant.hh"
#endif
#include "DPFCoord.h"
#include "DPFPWAPoint.h"
#include "DPFAngular.h"
#include "DPFPropogator.h"

#include "PWA_PARAS.h"
#include "PWA_CTRL.H"


class DPFPWAPdf : public RooAbsPdf {
    public:
        TString work_path;
        //RooAbsReal b;
        DPFPWAPdf(const char *name, const char *title,
                RooAbsReal&,
                DPFPWAPoint *dp);
        DPFPWAPdf(const DPFPWAPdf& other, const char* name=0) ;
        virtual TObject* clone(const char* newname) const {
            cout << "TObject clone" << endl;
            return new DPFPWAPdf(*this,newname); }
        inline virtual ~DPFPWAPdf() {
            for(Int_t i = 0; i < Nmc + Nmc_data; i++) {
                delete[] mcp1[i];
                delete[] mcp2[i];
                delete[] mcp3[i];
                delete[] mcp4[i];
                delete[] mcp5[i];
            }
            delete[] mcp1;
            delete[] mcp2;
            delete[] mcp3;
            delete[] mcp4;
            delete[] mcp5;
            cout << "delete mcp1~5!!!" << endl;
            fx.resize(0);
            pwa_paras.resize(0);
            _dp = NULL;

            for(int i = 0; i < Nmc + Nmc_data; i++) {
                delete[] mlk[i];
            }
            delete[] mlk;
            removeAll();
                delete _spinIter;
                delete _massIter;
                delete _mass2Iter;
                delete _widthIter;
                delete _g1Iter;
                delete _g2Iter;
                delete _b1Iter;
                delete _b2Iter;
                delete _b3Iter;
                delete _b4Iter;
                delete _b5Iter;
                delete _rhoIter;
                delete _fracIter;
                delete _phiIter;
                delete _propIter;
        };
        void removeAll();
        void setup_iter_vec();
        inline void showNumAll() {
            cout << _spinList.getSize() << endl;
            cout << _massList.getSize() << endl;
            cout << _mass2List.getSize() << endl;
            cout << _widthList.getSize() << endl;
            cout << _g1List.getSize() << endl;
            cout << _g2List.getSize() << endl;
            cout << _b1List.getSize() << endl;
            cout << _b2List.getSize() << endl;
            cout << _b3List.getSize() << endl;
            cout << _b4List.getSize() << endl;
            cout << _b5List.getSize() << endl;
            cout << _rhoList.getSize() << endl;
            cout << _fracList.getSize() << endl;
            cout << _phiList.getSize() << endl;
            cout << _propList.getSize() << endl;

        };
        void  project(const RooArgList& cPar, const RooArgList& fPar,const char* fname);
        void projectkk(const PWA_CTRL &);
        void projectpipi(const PWA_CTRL &);
        void addResonance600  (const Char_t* name, const Char_t* title, RooAbsReal& spin,
                RooAbsReal& mass, RooAbsReal& b1,  RooAbsReal& b2, RooAbsReal& b3, RooAbsReal& b4, RooAbsReal& b5, RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,RooAbsReal& propType);
        void addResonance980  (const Char_t* name, const Char_t* title, RooAbsReal& spin,
                RooAbsReal& mass, RooAbsReal& g1,  RooAbsReal& g2, RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,RooAbsReal& propType);
        void addResonance1680  (const Char_t* name, const Char_t* title, RooAbsReal& spin,
                RooAbsReal& mass1680, RooAbsReal& mass2,RooAbsReal& width,RooAbsReal& g1,  RooAbsReal& g2, RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,
                RooAbsReal& propType);
        void addResonance  (const TString name, const TString title,RooAbsReal& spin,
                RooAbsReal& mass, RooAbsReal& width,RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,RooAbsReal& propType);
        //  void addResonance  (const Char_t* name, const Char_t* title,RooAbsReal& spin,
        //                      RooAbsReal& mass, RooAbsReal& width,RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,RooAbsReal& propType);
        RooArgSet *fitFractions(const RooArgList& newPar, Bool_t print=kFALSE, ostream& os=std::cout);
        Double_t scalar(Double_t *a1, Double_t *a2)const;

        //Double_t getValV(const RooArgSet* nset) const;

        int count_lines(TString);
        void showAllParas();
        void writeToFile(TString fname = "ParaResults.root");
        void readFromFile(TString fname = "ParaResults.root");
        void change_value(TString rn, double xx);
        void paras_getval() const;
    protected:

        RooRealProxy idp ;

        RooRealProxy p11 ;
        RooRealProxy p12 ;
        RooRealProxy p13 ;
        RooRealProxy p14 ;
        RooRealProxy p21 ;
        RooRealProxy p22 ;
        RooRealProxy p23 ;
        RooRealProxy p24 ;
        RooRealProxy p31 ;
        RooRealProxy p32 ;
        RooRealProxy p33 ;
        RooRealProxy p34 ;
        RooRealProxy p41 ;
        RooRealProxy p42 ;
        RooRealProxy p43 ;
        RooRealProxy p44 ;
        RooRealProxy p51 ;
        RooRealProxy p52 ;
        RooRealProxy p53 ;
        RooRealProxy p54 ;
        Double_t evaluate() const;
        Double_t evaluate(int) const;
        //Double_t evaluate (
        //        Double_t _p11, Double_t _p12, Double_t _p13, Double_t _p14,
        //        Double_t _p21, Double_t _p22, Double_t _p23, Double_t _p24,
        //        Double_t _p31, Double_t _p32, Double_t _p33, Double_t _p34,
        //        Double_t _p41, Double_t _p42, Double_t _p43, Double_t _p44,
        //        Double_t _p51, Double_t _p52, Double_t _p53, Double_t _p54) const;
        Double_t calEva(const PWA_PARAS &pp, int) const;
        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
        Double_t analyticalIntegral (Int_t code, const char* rangeName) const;

        vector<PWA_PARAS> pwa_paras;
        mutable vector<Double_t> fx;
        //vector<PWA_PARAS> pwa_paras_data;
        mutable bool update_fx;
        void cu_inti_data(double * &h_float_pp,int * &h_parameter,double * &h_paraList,double *&h_fx,double* &h_mlk,int iEnd);
        void store_pwa_paras(); //用来进行内存换时间的操作，将所有需要的PWA_PARAS参数放到队列pwa_paras和pwa_paras_data中去
        void store_fx(int, int) const;
    private:
        void initialize();
        RooListProxy _spinList;
        RooListProxy _massList;
        RooListProxy _mass2List;
        RooListProxy _widthList;
        RooListProxy _g1List;
        RooListProxy _g2List;
        RooListProxy _b1List;
        RooListProxy _b2List;
        RooListProxy _b3List;
        RooListProxy _b4List;
        RooListProxy _b5List;
        RooListProxy _rhoList;
        RooListProxy _fracList;
        RooListProxy _phiList;
        RooListProxy _propList;

        TIterator* _spinIter;
        TIterator* _massIter;
        TIterator* _mass2Iter;
        TIterator* _widthIter;
        TIterator* _g1Iter;
        TIterator* _g2Iter;
        TIterator* _b1Iter;
        TIterator* _b2Iter;
        TIterator* _b3Iter;
        TIterator* _b4Iter;
        TIterator* _b5Iter;
        TIterator* _rhoIter;
        TIterator* _fracIter;
        TIterator* _phiIter;
        TIterator* _propIter;

        mutable vector<double> paraList;
        int _CN_spinList;
        int _CN_massList;
        int _CN_mass2List;
        int _CN_widthList;
        int _CN_g1List;
        int _CN_g2List;
        int _CN_b1List;
        int _CN_b2List;
        int _CN_b3List;
        int _CN_b4List;
        int _CN_b5List;
        int _CN_rhoList;
        int _CN_fracList;
        int _CN_phiList;
        int _CN_propList;
        int _CN_end;

        //  RooListProxy _nameList;
        TString *nameList;
        TString *titleList;
        TString *titleListT;
        Int_t Nmc;
        Int_t Nmc_data;
        Int_t nAmps;
        Int_t nStates;
        Int_t nStatesb;
        Int_t nStatesg1g2;
        Int_t nStateswidth;
        Int_t nStatesmass2;

        Double_t **mcp1;
        Double_t **mcp2;
        Double_t **mcp3;
        Double_t **mcp4;
        Double_t **mcp5;

        Double_t **mlk; //用来计算likelihood附加项的数组
        mutable Double_t penalty;
        mutable Double_t penalty_data;

        //Double_t **mcp1_data;
        //Double_t **mcp2_data;
        //Double_t **mcp3_data;
        //Double_t **mcp4_data;
        //Double_t **mcp5_data;

        DPFPWAPoint *_dp;
        DPFAngular _amp;
        DPFPropogator _prop;

        mutable Double_t anaIntegral;


        //  ClassDef(DPFPWAPdf,0) // Your description goes here...
};

#endif

