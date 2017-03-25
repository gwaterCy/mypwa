#ifndef DPF_PWAPDF
#define DPF_PWAPDF


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
class DPFPWAPdf : public RooAbsPdf {
public:
  DPFPWAPdf(const char *name, const char *title,
	      RooAbsReal& _p11,
	      RooAbsReal& _p12,
	      RooAbsReal& _p13,
	      RooAbsReal& _p14,
	      RooAbsReal& _p21,
	      RooAbsReal& _p22,
	      RooAbsReal& _p23,
	      RooAbsReal& _p24,
	      RooAbsReal& _p31,
	      RooAbsReal& _p32,
	      RooAbsReal& _p33,
	      RooAbsReal& _p34,
	      RooAbsReal& _p41,
	      RooAbsReal& _p42,
	      RooAbsReal& _p43,
	      RooAbsReal& _p44,
	      RooAbsReal& _p51,
	      RooAbsReal& _p52,
	      RooAbsReal& _p53,
	      RooAbsReal& _p54,
              DPFPWAPoint *dp);
  DPFPWAPdf(const DPFPWAPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new DPFPWAPdf(*this,newname); }
  inline virtual ~DPFPWAPdf() { };
  void  project(const RooArgList& cPar, const RooArgList& fPar,const char* fname); 
  void addResonance600  (const Char_t* name, const Char_t* title, RooAbsReal& spin,
                      RooAbsReal& mass, RooAbsReal& b1,  RooAbsReal& b2, RooAbsReal& b3, RooAbsReal& b4, RooAbsReal& b5, RooAbsReal& rho, RooAbsReal& phi,RooAbsReal& propType);
  void addResonance980  (const Char_t* name, const Char_t* title, RooAbsReal& spin,
                      RooAbsReal& mass, RooAbsReal& g1,  RooAbsReal& g2, RooAbsReal& rho, RooAbsReal& phi,RooAbsReal& propType);
  void addResonance1680  (const Char_t* name, const Char_t* title, RooAbsReal& spin,
                      RooAbsReal& mass1680, RooAbsReal& mass2,RooAbsReal& width,RooAbsReal& g1,  RooAbsReal& g2, RooAbsReal& rho, RooAbsReal& phi,
                      RooAbsReal& propType);
  void addResonance  (const Char_t* name, const Char_t* title,RooAbsReal& spin,
                      RooAbsReal& mass, RooAbsReal& width,RooAbsReal& rho, RooAbsReal& phi,RooAbsReal& propType);
  RooArgSet *fitFractions(const RooArgList& newPar, Bool_t print=kFALSE, ostream& os=std::cout);
  Double_t scalar(Double_t *a1, Double_t *a2)const;
  int count_lines();
protected:

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
  Double_t evaluate()const; 
  Double_t evaluate (Double_t _p11, Double_t _p12, Double_t _p13, Double_t _p14,
                     Double_t _p21, Double_t _p22, Double_t _p23, Double_t _p24,
                     Double_t _p31, Double_t _p32, Double_t _p33, Double_t _p34,
                     Double_t _p41, Double_t _p42, Double_t _p43, Double_t _p44,
                     Double_t _p51, Double_t _p52, Double_t _p53, Double_t _p54)const;
  Double_t calEva (Double_t _p11, Double_t _p12, Double_t _p13, Double_t _p14,
                     Double_t _p21, Double_t _p22, Double_t _p23, Double_t _p24,
                     Double_t _p31, Double_t _p32, Double_t _p33, Double_t _p34,
                     Double_t _p41, Double_t _p42, Double_t _p43, Double_t _p44,
                     Double_t _p51, Double_t _p52, Double_t _p53, Double_t _p54) const ;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral (Int_t code, const char* rangeName) const;

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
  RooListProxy _phiList;
  RooListProxy _propList;
//  RooListProxy _nameList;
  TString *nameList;
  TString *titleList;
  TString *titleListT;
  Int_t Nmc;
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
  DPFPWAPoint *_dp;
  DPFAngular _amp;
  DPFPropogator _prop;
//  ClassDef(DPFPWAPdf,0) // Your description goes here...
};
 
#endif

