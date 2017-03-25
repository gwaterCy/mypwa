/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * Copyright (c) 2000-2005, Regents of the University of California          * 
 *                          and Stanford University. All rights reserved.    * 
 *                                                                           * 
 * Redistribution and use in source and binary forms,                        * 
 * with or without modification, are permitted according to the terms        * 
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             * 
 *****************************************************************************/ 

// -- CLASS DESCRIPTION [PDF] -- 
// Your description goes here... 

#include "DPFPWAPdf.h" 

#include <stdio.h>  
#include <time.h>  
#include <stdlib.h>  
#include <omp.h>  

/*#include "RooAbsReal.h" 
#include "RooListProxy.h" 
#include "RooStringVar.h" 
#include "TFile.h" 
#include "TH1F.h" 
#include "TMath.h" 
#include<TLorentzVector.h>
#include<TStopwatch.h>
*/
Double_t rk=0.493677,rp=0.13957018 ;

//ClassImp(DPFPWAPdf) 

DPFPWAPdf::DPFPWAPdf(const char *name, const char *title, 
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
        DPFPWAPoint *dp) :
            RooAbsPdf(name,title), 
            p11("p11","p11",this,_p11),
            p12("p12","p12",this,_p12),
            p13("p13","p13",this,_p13),
            p14("p14","p14",this,_p14),
            p21("p21","p21",this,_p21),
            p22("p22","p22",this,_p22),
            p23("p23","p23",this,_p23),
            p24("p24","p24",this,_p24),
            p31("p31","p31",this,_p31),
            p32("p32","p32",this,_p32),
            p33("p33","p33",this,_p33),
            p34("p34","p34",this,_p34),
            p41("p41","p41",this,_p41),
            p42("p42","p42",this,_p42),
            p43("p43","p43",this,_p43),
            p44("p44","p44",this,_p44),
            p51("p51","p51",this,_p51),
            p52("p52","p52",this,_p52),
            p53("p53","p53",this,_p53),
            p54("p54","p54",this,_p54),
            _spinList ("spinList",  "list of spins", this),
            _massList ("massList",  "list of masses", this),
            _mass2List ("mass2List",  "list of mass2s", this),
            _widthList("widthList", "list of widths", this),
            _g1List ("g1List",  "list of g1s", this),
            _g2List ("g2List",  "list of g2s", this),
            _b1List ("b1List",  "list of b1s", this),
            _b2List ("b2List",  "list of b2s", this),
            _b3List ("b3List",  "list of b3s", this),
            _b4List ("b4List",  "list of b4s", this),
            _b5List ("b5List",  "list of b5s", this),
            _rhoList("rhoList", "list of rhos", this),
            _phiList("phiList", "list of phis", this),
            _propList("propList", "list of props", this)
{
    _dp         = dp;
    _amp.setdp(_dp);
    //  cout<<"haha: "<< __LINE__ << endl;
    initialize();
    //  cout<<"haha: "<< __LINE__ << endl;
} 

DPFPWAPdf::DPFPWAPdf(const DPFPWAPdf& other, const char* name) :  
    RooAbsPdf(other,name), 
    p11("p11",this,other.p11),
    p12("p12",this,other.p12),
    p13("p13",this,other.p13),
    p14("p14",this,other.p14),
    p21("p21",this,other.p21),
    p22("p22",this,other.p22),
    p23("p23",this,other.p23),
    p24("p24",this,other.p24),
    p31("p31",this,other.p31),
    p32("p32",this,other.p32),
    p33("p33",this,other.p33),
    p34("p34",this,other.p34),
    p41("p41",this,other.p41),
    p42("p42",this,other.p42),
    p43("p43",this,other.p43),
    p44("p44",this,other.p44),
    p51("p51",this,other.p51),
    p52("p52",this,other.p52),
    p53("p53",this,other.p53),
    p54("p54",this,other.p54),
    _spinList ("spinList",  this,other._spinList),
    _massList ("massList",  this,other._massList),
    _mass2List ("mass2List",  this,other._mass2List),
    _widthList("widthList", this,other._widthList),
    _g1List ("g1List",  this,other._g1List),
    _g2List ("g2List",  this,other._g2List),
    _b1List ("b1List",  this,other._b1List),
    _b2List ("b2List",  this,other._b2List),
    _b3List ("b3List",  this,other._b3List),
    _b4List ("b4List",  this,other._b4List),
    _b5List ("b5List",  this,other._b5List),
    _rhoList("rhoList", this,other._rhoList),
    _phiList("phiList", this,other._phiList),
    _propList("propList", this,other._propList)
{
    _dp         = other._dp;
    _amp.setdp(_dp);
    //  cout<<"haha: "<< __LINE__ << endl;
    Nmc=other.Nmc;
    nAmps=other.nAmps;
    nStates=other.nStates;
    nStatesb=other.nStatesb;
    nStatesg1g2=other.nStatesg1g2;
    nStateswidth=other.nStateswidth;
    nStatesmass2=other.nStatesmass2;
    nameList=new TString[60];
    titleList=new TString[60];
    titleListT=new TString[60];
    for(Int_t i=0;i<60;i++){
        nameList[i]=other.nameList[i];
        titleList[i]=other.titleList[i];
        titleListT[i]=other.titleListT[i];
    }
    mcp1=new Double_t*[Nmc];
    mcp2=new Double_t*[Nmc];
    mcp3=new Double_t*[Nmc];
    mcp4=new Double_t*[Nmc];
    mcp5=new Double_t*[Nmc];
    for(Int_t i=0;i<Nmc;i++){
        mcp1[i]=new Double_t[4];
        mcp2[i]=new Double_t[4];
        mcp3[i]=new Double_t[4];
        mcp4[i]=new Double_t[4];
        mcp5[i]=new Double_t[4];
    }
    for(Int_t i=0;i<Nmc;i++){
        mcp1[i][0]=other.mcp1[i][0];mcp1[i][1]=other.mcp1[i][1];mcp1[i][2]=other.mcp1[i][2];mcp1[i][3]=other.mcp1[i][3];
        mcp2[i][0]=other.mcp2[i][0];mcp2[i][1]=other.mcp2[i][1];mcp2[i][2]=other.mcp2[i][2];mcp2[i][3]=other.mcp2[i][3];
        mcp3[i][0]=other.mcp3[i][0];mcp3[i][1]=other.mcp3[i][1];mcp3[i][2]=other.mcp3[i][2];mcp3[i][3]=other.mcp3[i][3];
        mcp4[i][0]=other.mcp4[i][0];mcp4[i][1]=other.mcp4[i][1];mcp4[i][2]=other.mcp4[i][2];mcp4[i][3]=other.mcp4[i][3];
        mcp5[i][0]=other.mcp5[i][0];mcp5[i][1]=other.mcp5[i][1];mcp5[i][2]=other.mcp5[i][2];mcp5[i][3]=other.mcp5[i][3];
    }
}
int DPFPWAPdf::count_lines() {
    ifstream in(_dp->_phspfile);
    string line;
    int n = 0;
    while (getline(in, line)) {
        n++;
    }
    cout << "PHSP file have " << n << " lines." << endl;
    return n;
}

void DPFPWAPdf::initialize()
{
    //  cout<<"haha: "<< __LINE__ << endl;
//    Nmc=300000;
//    Nmc = (count_lines() / 5 / 80 + 1) * 80;
    Nmc = count_lines() / 5;
    nAmps=0;
    nStates=0;
    nStatesb=0;
    nStatesg1g2=0;
    nStateswidth=0;
    nStatesmass2=0;
    nameList=new TString[60];
    titleList=new TString[60];
    titleListT=new TString[60];
    mcp1=new Double_t*[Nmc];
    mcp2=new Double_t*[Nmc];
    mcp3=new Double_t*[Nmc];
    mcp4=new Double_t*[Nmc];
    mcp5=new Double_t*[Nmc];
    for(Int_t i=0;i<Nmc;i++){
        mcp1[i]=new Double_t[4];
        mcp2[i]=new Double_t[4];
        mcp3[i]=new Double_t[4];
        mcp4[i]=new Double_t[4];
        mcp5[i]=new Double_t[4];
    }
    Double_t fx1,fy1,fz1,ft1,fx2,fy2,fz2,ft2,fx3,fy3,fz3,ft3,fx4,fy4,fz4,ft4,fx5,fy5,fz5,ft5;
    FILE *fp;
    if((fp=fopen(_dp->_phspfile,"r"))==NULL)
    {printf("can't open input file");
        return;
    }
    Int_t i=0;
    while(fscanf(fp,"%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n",&fx1,&fy1,&fz1,&ft1,&fx2,&fy2,&fz2,&ft2,&fx3,&fy3,&fz3,&ft3,&fx4,&fy4,&fz4,&ft4,&fx5,&fy5,&fz5,&ft5)!=EOF)
    {
        //  cout<<"haha: "<< __LINE__ << endl;
        mcp1[i][0]=fx1;mcp1[i][1]=fy1;mcp1[i][2]=fz1;mcp1[i][3]=ft1;
        mcp2[i][0]=fx2;mcp2[i][1]=fy2;mcp2[i][2]=fz2;mcp2[i][3]=ft2;
        mcp3[i][0]=fx3;mcp3[i][1]=fy3;mcp3[i][2]=fz3;mcp3[i][3]=ft3;
        mcp4[i][0]=fx4;mcp4[i][1]=fy4;mcp4[i][2]=fz4;mcp4[i][3]=ft4;
        mcp5[i][0]=fx5;mcp5[i][1]=fy5;mcp5[i][2]=fz5;mcp5[i][3]=ft5;
        i++;
    }
    fclose(fp);
//    Nmc = count_lines() / 5;
    if (i != Nmc) {
        cout << "There is memory allocated error!" << endl;
        exit(1);
    } 
}

Double_t DPFPWAPdf::evaluate()const
{ 
    //  cout<<"haha: "<< __LINE__ << endl;
    return evaluate(p11,p12,p13,p14,p21,p22,p23,p24,p31,p32,p33,p34,p41,p42,p43,p44,p51,p52,p53,p54);
} 

Double_t DPFPWAPdf::evaluate (Double_t _p11, Double_t _p12, Double_t _p13, Double_t _p14,
        Double_t _p21, Double_t _p22, Double_t _p23, Double_t _p24,
        Double_t _p31, Double_t _p32, Double_t _p33, Double_t _p34,
        Double_t _p41, Double_t _p42, Double_t _p43, Double_t _p44,
        Double_t _p51, Double_t _p52, Double_t _p53, Double_t _p54)const
/////return value of object function

{
    //  cout<<"haha: "<< __LINE__ << endl;
    // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
    Double_t sum = 1e-20;
    // cout<<"p11="<<_p11<<endl;
    // cout<<"p21="<<_p21<<endl;
    Double_t dsigma=calEva(_p11,_p12,_p13,_p14,_p21,_p22,_p23,_p24,_p31,_p32,_p33,_p34,_p41,_p42,_p43,_p44,_p51,_p52,_p53,_p54);
    //  Double_t sigma=mcIntegral();
    //  sum=dsigma/sigma;
    sum=dsigma;
    return (sum <= 0) ? 1e-20 : sum;
} 


void DPFPWAPdf::addResonance600(const Char_t* name, const Char_t* title, RooAbsReal& spin, 
        RooAbsReal& mass, RooAbsReal& b1, RooAbsReal& b2, RooAbsReal& b3,RooAbsReal& b4,RooAbsReal& b5,RooAbsReal& rho, RooAbsReal& phi,RooAbsReal& propType)
{       _spinList.add(spin);
    //  cout<<"haha: "<< __LINE__ << endl;
    _massList.add(mass);
    _b1List.add(b1);
    _b2List.add(b2);
    _b3List.add(b3);
    _b4List.add(b4);
    _b5List.add(b5);
    _rhoList.add(rho);
    _phiList.add(phi);
    _propList.add(propType);
    nameList[nAmps]= name;
    titleListT[nAmps]= title;
    titleList[nStates]= title;
    nAmps++;
    nStates++;
    nStatesb++;
}   
void DPFPWAPdf::addResonance980(const Char_t* name, const Char_t* title, RooAbsReal& spin, 
        RooAbsReal& mass, RooAbsReal& g1,RooAbsReal& g2,RooAbsReal& rho, RooAbsReal& phi,RooAbsReal& propType)
{       _spinList.add(spin);
    //  cout<<"haha: "<< __LINE__ << endl;
    _massList.add(mass);
    _g1List.add(g1);
    _g2List.add(g2);
    _rhoList.add(rho);
    _phiList.add(phi);
    _propList.add(propType);
    nameList[nAmps]= name;
    titleListT[nAmps]= title;
    titleList[nStates]= title;
    nAmps++;
    nStates++;
    nStatesg1g2++;
}   
void DPFPWAPdf::addResonance1680(const Char_t* name, const Char_t* title, RooAbsReal& spin, 
        RooAbsReal& mass1680, RooAbsReal& mass2,RooAbsReal& width,RooAbsReal& g1,RooAbsReal& g2,RooAbsReal& rho, RooAbsReal& phi,
        RooAbsReal& propType)
{       _spinList.add(spin);
    //  cout<<"haha: "<< __LINE__ << endl;
    _massList.add(mass1680);
    _mass2List.add(mass2);
    _widthList.add(width);
    _g1List.add(g1);
    _g2List.add(g2);
    _rhoList.add(rho);
    _phiList.add(phi);
    _propList.add(propType);
    nameList[nAmps]= name;
    titleListT[nAmps]= title;
    titleList[nStates]= title;
    nAmps++;
    nStates++;
    nStatesg1g2++;
    nStatesmass2++;
    nStateswidth++;
}   
void DPFPWAPdf::addResonance(const Char_t* name, const Char_t* title, RooAbsReal& spin, 
        RooAbsReal& mass, RooAbsReal& width,RooAbsReal& rho, RooAbsReal& phi,RooAbsReal& propType)
{       _spinList.add(spin);
    //  cout<<"haha: "<< __LINE__ << endl;
    _massList.add(mass);
    _widthList.add(width);
    _rhoList.add(rho);
    _phiList.add(phi);
    _propList.add(propType);
    nameList[nAmps]= name;
    titleListT[nAmps]= title;
    titleList[nStates]= title;
    nAmps++;
    nStates++;
    nStateswidth++;
}   
Double_t DPFPWAPdf::calEva (Double_t _p11, Double_t _p12, Double_t _p13, Double_t _p14,
        Double_t _p21, Double_t _p22, Double_t _p23, Double_t _p24, 
        Double_t _p31, Double_t _p32, Double_t _p33, Double_t _p34, 
        Double_t _p41, Double_t _p42, Double_t _p43, Double_t _p44, 
        Double_t _p51, Double_t _p52, Double_t _p53, Double_t _p54) const 
////return square of complex amplitude
{
    //	static int A=0;
    //	A++;
    Double_t value = 0.;
    //  cout<<"haha: "<< __LINE__ << endl;
    Double_t ak1[4],ak2[4],ak3[4],ak4[4],ak5[4];
    ak1[0]=_p11; ak1[1]=_p12; ak1[2]=_p13; ak1[3]=_p14;
    ak2[0]=_p21; ak2[1]=_p22; ak2[2]=_p23; ak2[3]=_p24;
    ak3[0]=_p31; ak3[1]=_p32; ak3[2]=_p33; ak3[3]=_p34;
    ak4[0]=_p41; ak4[1]=_p42; ak4[2]=_p43; ak4[3]=_p44;
    ak5[0]=_p51; ak5[1]=_p52; ak5[2]=_p53; ak5[3]=_p54;
    //   Double_t q[4],akd[4],akx[4];
    Double_t wu[4],w0p22[4],ak23w[4],w2p2[4],w2p1[4],w2p3[2],w2p4[2],w2p5[2];
    Double_t b2qf2xx,b4qjvf2,b2qjv2,b2qjv3,b2qbv2,b2qbv3,b1qjv2,b1qjv3,b1qbv2,b1qbv3,q2r23,b1q2r23;
    Double_t sv,s23,sv2,sv3;
    Double_t ap23[4],apv2[4],apv3[4];
    Double_t wpf22[2];
    Double_t temp,b2qjvf2;
    TComplex fCF[nAmps][4];
    TComplex fCP[nAmps];
    TComplex pa[nAmps][nAmps];
    TComplex fu[nAmps][nAmps];
    Double_t w1p12_1[4],w1p13_1[4],w1p12_2[4],w1p13_2[4],w1p12_1u[4],w1p13_1u[4],w1p12_2u[4],w1p13_2u[4];
    Double_t w1p12_3[2],w1p13_3[2],w1p12_4[2],w1p13_4[2];
    Double_t w1m12[2],w1m13[2];
    TComplex crp1[nAmps],crp11[nAmps];
    TComplex cr0p11;
    TComplex ca2p1;
    TComplex cw2p11;
    TComplex cw2p12;
    TComplex cw2p15;
    TComplex cd0p1;
    TComplex cd0p2;
    TComplex cd2p1;
    TComplex cw;
    TComplex c1p12_12,c1p13_12,c1p12_13,c1p13_13,c1p12_14,c1p13_14;
    TComplex cr1m12_1,cr1m13_1;
    TComplex crpf1,crpf2;

    for(Int_t i=0;i<4;i++){
        ap23[i]=ak2[i]+ak3[i];
        apv2[i]=ak1[i]+ak2[i];
        apv3[i]=ak1[i]+ak3[i];
    }
    //	//  fSX=scalar(akx,akx);
    sv=scalar(ak1,ak1);
    //	cout<<"sv="<<sv<<endl;
    s23=scalar(ap23,ap23);
    //	cout<<"s23="<<s23<<endl;
    sv2=scalar(apv2,apv2);
    sv3=scalar(apv3,apv3);

    _amp.calculate0p(ak1, ak2, ak3, ak4, ak5, wu, w0p22, w2p1, w2p2, w2p3, w2p4, w2p5, ak23w, wpf22, b2qjvf2,b2qf2xx, b4qjvf2, b1q2r23);

    TIterator* _spinIter_ = _spinList.createIterator();
    TIterator* _massIter_ = _massList.createIterator();
    TIterator* _mass2Iter_ = _mass2List.createIterator();
    TIterator* _widthIter_ = _widthList.createIterator();
    TIterator* _g1Iter_ = _g1List.createIterator();
    TIterator* _g2Iter_ = _g2List.createIterator();
    TIterator* _b1Iter_ = _b1List.createIterator();
    TIterator* _b2Iter_ = _b2List.createIterator();
    TIterator* _b3Iter_ = _b3List.createIterator();
    TIterator* _b4Iter_ = _b4List.createIterator();
    TIterator* _b5Iter_ = _b5List.createIterator();
    TIterator* _rhoIter_ = _rhoList.createIterator();
    TIterator* _phiIter_ = _phiList.createIterator();
    TIterator* _propIter_ = _propList.createIterator();
    //  cout<<"haha: "<< __LINE__ << endl;
    //#pragma omp parallel for
    for(Int_t index=0; index<nAmps; index++) {
        //  cout<<"haha: "<< __LINE__ << endl;
        RooRealVar *mass  = (RooRealVar*)_massIter_->Next();
        //RooRealVar *width = (RooRealVar*)_widthIter->Next();
        RooRealVar *spin = (RooRealVar*)_spinIter_->Next();
        RooRealVar *rho = (RooRealVar*)_rhoIter_->Next();
        RooRealVar *phi = (RooRealVar*)_phiIter_->Next();
        RooRealVar *propType= (RooRealVar*)_propIter_->Next();
        Double_t rho0=rho->getVal();
        Double_t phi0=phi->getVal();
        Int_t spin_now=spin->getVal();
        Int_t propType_now=propType->getVal();

        fCP[index]=TComplex(rho0*TMath::Cos(phi0),rho0*TMath::Sin(phi0));
        //        cout<<"fCP[index]="<<fCP[index]<<endl;
        //        std::cout << __FILE__ << __LINE__ << " : " << propType_now << std::endl;
        switch(propType_now)
        {
            //  cout<<"haha: "<< __LINE__ << endl;
            //                     ordinary  Propagator  Contribution 
            case 1:
                {
                    RooRealVar *width = (RooRealVar*)_widthIter_->Next();
                    Double_t mass0=mass->getVal();
                    Double_t width0=width->getVal();
                    //					cout<<"mass0="<<mass0<<endl;
                    //					cout<<"width0="<<width0<<endl;
                    crp1[index]=_prop.propogator(mass0,width0,s23);
                }
                break;
                //	Flatte   Propagator Contribution
            case 2:
                {
                    RooRealVar *g1 = (RooRealVar*)_g1Iter_->Next();
                    RooRealVar *g2 = (RooRealVar*)_g2Iter_->Next();
                    Double_t mass980=mass->getVal();
                    Double_t g10=g1->getVal();
                    Double_t g20=g2->getVal();
                    //			cout<<"mass980="<<mass980<<endl;
                    //			cout<<"g10="<<g10<<endl;
                    //			cout<<"g20="<<g20<<endl;
                    crp1[index]=_prop.propogator980(mass980,g10,g20,s23);
                    //			cout<<"crp1[index]="<<crp1[index]<<endl;
                }
                break;
                // sigma  Propagator Contribution
            case 3:
                {
                    RooRealVar *b1 = (RooRealVar*)_b1Iter_->Next();
                    RooRealVar *b2 = (RooRealVar*)_b2Iter_->Next();
                    RooRealVar *b3 = (RooRealVar*)_b3Iter_->Next();
                    RooRealVar *b4 = (RooRealVar*)_b4Iter_->Next();
                    RooRealVar *b5 = (RooRealVar*)_b5Iter_->Next();
                    Double_t mass600=mass->getVal();
                    Double_t b10=b1->getVal();
                    Double_t b20=b2->getVal();
                    Double_t b30=b3->getVal();
                    Double_t b40=b4->getVal();
                    Double_t b50=b5->getVal();
                    crp1[index]=_prop.propogator600(mass600,b10,b20,b30,b40,b50,s23);
                    //			cout<<"crp1[index]3="<<crp1[index]<<endl;
                }
                break;
                // 1- or 1+  Contribution
            case 4:
                {
                    RooRealVar *width = (RooRealVar*)_widthIter_->Next();
                    Double_t mass0=mass->getVal();
                    Double_t width0=width->getVal();
                    crp1[index]=_prop.propogator(mass0,width0,sv2);
                    crp11[index]=_prop.propogator(mass0,width0,sv3);
                }
                break;
                //  phi(1650) f0(980) include flatte and ordinary Propagator joint Contribution
            case 5:
                {
                    RooRealVar *mass2  = (RooRealVar*)_mass2Iter_->Next();
                    RooRealVar *g1 = (RooRealVar*)_g1Iter_->Next();
                    RooRealVar *g2 = (RooRealVar*)_g2Iter_->Next();
                    Double_t mass980=mass2->getVal();
                    Double_t g10=g1->getVal();
                    Double_t g20=g2->getVal();
                    //					cout<<"mass980="<<mass980<<endl;
                    //					cout<<"g10="<<g10<<endl;
                    //					cout<<"g20="<<g20<<endl;
                    crp1[index]=_prop.propogator980(mass980,g10,g20,sv);
                    //					cout<<"crp1[index]="<<crp1[index]<<endl;
                    RooRealVar *width = (RooRealVar*)_widthIter_->Next();
                    Double_t mass1680=mass->getVal();
                    Double_t width1680=width->getVal();
                    //					cout<<"mass1680="<<mass1680<<endl;
                    //					cout<<"width1680="<<width1680<<endl;
                    crp11[index]=_prop.propogator(mass1680,width1680,s23);
                    //					cout<<"crp11[index]="<<crp11[index]<<endl;
                }
                break;
            case 6:
                {
                    RooRealVar *width = (RooRealVar*)_widthIter_->Next();
                    Double_t mass0=mass->getVal();
                    Double_t width0=width->getVal();
                    //					cout<<"mass0="<<mass0<<endl;
                    //					cout<<"width0="<<width0<<endl;
                    crp1[index]=_prop.propogator1270(mass0,width0,s23);
                    //			cout<<"crp1[index]6="<<crp1[index]<<endl;
                }
            default :
                ;
        }
        for(Int_t i=0;i<2;i++){
            //  cout<<"haha: "<< __LINE__ << endl;
            //		cout<<"spin_now="<<spin_now<<endl;
            switch(spin_now)
            { 
                case 11:
                    //1+_1 contribution
                    fCF[index][i]=w1p12_1[i]*crp1[index]+w1p13_1[i]*crp11[i];

                    break;
                case 12:
                    //1+_2 contribution
                    c1p12_12=crp1[index]/b2qbv2;
                    c1p13_12=crp11[index]/b2qbv3;
                    fCF[index][i]=w1p12_2[i]*c1p12_12+w1p13_2[i]*c1p13_12;

                    break;
                case 13:
                    //1+_3 contribution
                    c1p12_13=crp1[index]/b2qjv2;
                    c1p13_13=crp11[index]/b2qjv3;
                    fCF[index][i]=w1p12_3[i]*c1p12_13+w1p13_3[i]*c1p13_13;

                    break;
                case 14:
                    //1+_4 contribution
                    c1p12_12=crp1[index]/b2qbv2;
                    c1p13_12=crp11[index]/b2qbv3;
                    c1p12_14=c1p12_12/b2qjv2;
                    c1p13_14=c1p13_12/b2qjv3;
                    fCF[index][i]=w1p12_4[i]*c1p12_14+w1p13_4[i]*c1p13_14;

                    break;
                case 111:
                    //1-__1 contribution
                    cr1m12_1=crp1[index]/b1qjv2/b1qbv2;
                    cr1m13_1=crp11[index]/b1qjv3/b1qbv3;
                    fCF[index][i]=w1m12[i]*cr1m12_1+w1m13[i]*cr1m13_1;

                    break;
                case 191:
                    //phi(1650)f0(980)_1 contribution
                    //		cout<<"b1q2r23="<<b1q2r23<<endl;
                    crpf1=crp1[index]*crp11[index]/b1q2r23;
                    //		cout<<"crpf1="<<crpf1<<endl;
                    fCF[index][i]=ak23w[i]*crpf1;
                    //	cout<<"fCF[index][i]="<<fCF[index][i]<<endl;

                    break;
                case 192:
                    //phi(1650)f0(980)_2 contribution
                    crpf1=crp1[index]*crp11[index]/b1q2r23;
                    crpf2=crpf1/b2qjvf2;
                    fCF[index][i]=wpf22[i]*crpf2;

                    break;
                case 1:
                    //  cout<<"haha: "<< __LINE__ << endl;
                    //01 contribution
                    //	cout<<"wu[i]="<<wu[i]<<endl;
                    //	cout<<"crp1[index]="<<crp1[index]<<endl;
                    //	cout<<"index="<<index<<endl;
                    fCF[index][i]=wu[i]*crp1[index];
                    //	cout<<"fCF[index][i]="<<fCF[index][i]<<endl;
                    //	cout<<"i="<<i<<endl;

                    break;
                case 2:
                    //02 contribution
                    cr0p11=crp1[index]/b2qjvf2;
                    fCF[index][i]=w0p22[i]*cr0p11;
                    //	cout<<"fCF[index][i]02="<<fCF[index][i]<<endl;

                    break;
                case 21:	
                    //21 contribution
                    //	cout<<"b2qf2xx="<<b2qf2xx<<endl;
                    cw2p11=crp1[index]/b2qf2xx;
                    //	cout<<"cw2p11="<<cw2p11<<endl;
                    //	cout<<"w2p1[0]="<<w2p1[0]<<endl;
                    //	cout<<"w2p1[1]="<<w2p1[1]<<endl;
                    fCF[index][i]=w2p1[i]*cw2p11;
                    //	cout<<"fCF[index][i]21="<<fCF[index][i]<<endl;

                    break;
                case 22:
                    //22 contribution
                    cw2p11=crp1[index]/b2qf2xx;
                    cw2p12=cw2p11/b2qjvf2;
                    fCF[index][i]=w2p2[i]*cw2p12;

                    break;
                case 23:
                    //23 contribution
                    cw2p11=crp1[index]/b2qf2xx;
                    cw2p12=cw2p11/b2qjvf2;
                    fCF[index][i]=w2p3[i]*cw2p12;

                    break;
                case 24:
                    //24 contribution
                    cw2p11=crp1[index]/b2qf2xx;
                    cw2p12=cw2p11/b2qjvf2;
                    fCF[index][i]=w2p4[i]*cw2p12;

                    break;
                case 25:
                    //25 contribution
                    cw2p11=crp1[index]/b2qf2xx;
                    cw2p15=cw2p11/b4qjvf2;          
                    fCF[index][i]=w2p5[i]*cw2p15;

                default:		;
            }
        }

    }

    //#pragma omp parallel for reduction(+:value)
    for(Int_t i=0;i<nAmps;i++){
        //  cout<<"haha: "<< __LINE__ << endl;
        for(Int_t j=0;j<nAmps;j++){
            cw=fCP[i]*TComplex::Conjugate(fCP[j]);
            //    cout<<"cw="<<cw<<endl;
            if(i==j) pa[i][j]=cw.Re();
            else if(i<j) pa[i][j]=2*cw.Re();
            else pa[i][j]=2*cw.Im();

            cw=TComplex(0,0);
            for(Int_t k=0;k<2;k++){
                cw+=fCF[i][k]*TComplex::Conjugate(fCF[j][k])/2.0;
                //   cout<<"cwfu="<<cw<<endl;

            }
            if(i<=j) fu[i][j]=cw.Re();
            if(i>j) fu[i][j]=-cw.Im();
            //      cout<<"pa[i][j]="<<pa[i][j]<<endl;
            //      cout<<"fu[i][j]="<<fu[i][j]<<endl;

            temp=pa[i][j]*fu[i][j];
            value+=temp;
        }
    }

    delete _spinIter_;
    delete _massIter_;
    delete _mass2Iter_;
    delete _widthIter_;
    delete _g1Iter_;
    delete _g2Iter_;
    delete _b1Iter_;
    delete _b2Iter_;
    delete _b3Iter_;
    delete _b4Iter_;
    delete _b5Iter_;
    delete _rhoIter_;
    delete _phiIter_;
    delete _propIter_;
    return (value <= 0) ? 1e-20 : value;
}

Double_t DPFPWAPdf::scalar(Double_t *a1,Double_t *a2)const
{
    //  cout<<"haha: "<< __LINE__ << endl;
    Double_t scal=0;
    //	Double_t fDel[4][4];
    //	for(Int_t i=0;i<4;i++){
    //		for(Int_t j=0;j<4;j++){
    //			if(i==j){
    //				if(i<3) fDel[i][j]=-1;
    //				else fDel[i][j]=1;
    //			}
    //			else fDel[i][j]=0;
    //		}
    //	}

    for(Int_t i=0;i<4;i++){
        scal+=a1[i]*a2[i]*_dp->fDel[i][i];
    }
    return scal;
} 

Int_t DPFPWAPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName)const
{
    //  cout<<"haha: "<< __LINE__ << endl;
    if  (rangeName) cout << " " <<endl;
    RooArgSet theSet1, theSet2, theSet3;
    theSet1.add(RooArgSet(p11.arg(),p12.arg(),p13.arg(),p14.arg(),p21.arg(),p22.arg(),p23.arg(),p24.arg()));
    theSet2.add(RooArgSet(p31.arg(),p32.arg(),p33.arg(),p34.arg(),p41.arg(),p42.arg(),p43.arg(),p44.arg()));
    theSet3.add(RooArgSet(p51.arg(),p52.arg(),p53.arg(),p54.arg()));
    RooArgSet theSet4(theSet1,theSet2," "); 
    RooArgSet theSet(theSet4, theSet3," "); 
    if (matchArgs(allVars, analVars, theSet)) { cout<<"haha: " << __LINE__ <<endl;
        return 1;}
        return 0;
}

Double_t DPFPWAPdf::analyticalIntegral (Int_t code, const char* rangeName) const
{
    cout<<"haha: "<< __LINE__ << endl;
    assert (code==1);
    Double_t sum=0;
    cout << "Begin analyticalIntegral, Nmc = " << Nmc << endl;
    double startTime=omp_get_wtime();
#pragma omp parallel for reduction(+:sum)
    for(Int_t i=0;i<Nmc;i++)
    {
        //  cout<<"haha: "<< __LINE__ << endl;
        Double_t eva=calEva(mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3]);
        sum=sum+eva;
    }
    sum=sum/Nmc;
    double endTime=omp_get_wtime();   
    cout<<"Integral time : "<< endTime - startTime <<"\n";
    cout << "sum = " << sum << endl;
//    cout << "End analyticalIntegral" << endl;
    return (sum <= 0) ? 1e-20 : sum;
}
void DPFPWAPdf::project(const RooArgList& cPar, const RooArgList& fPar, const char* fname)
{
    //  cout<<"haha: "<< __LINE__ << endl;
    TFile *f=new TFile(fname,"recreate");
    TH1F *hm23=new TH1F("m23","m23",65,0.15,2.1);
    TH1F *hm45=new TH1F("m45","m45",100,0.9,1.1);
    //	TH1F *hmkk1=new TH1F("mkk1","mkk",100,1,3);
    //	TH1F *hmkk2=new TH1F("mkk2","mkk",100,1,3);
    //	TH1F *hmkk3=new TH1F("mkk3","mkk",100,1,3);
    TH1F *hmphi=new TH1F("mphi","mphi",100,0,3);
    TH1F *hmf0=new TH1F("mf0","mf0",100,0.9,1.1);
    TH1F *h23=new TH1F("costhetaX23","costhetaX23",50,-1,1);
    TH1F *h45=new TH1F("costhetaX45","costhetaX45",50,-1,1);
    //	Double_t e[60];
    //save the old result of rhos to e[] // 
    //	_rhoIter->Reset();
    //	for (Int_t i=0; i<nAmps; i++) {
    //		RooRealVar *rhos = (RooRealVar*)_rhoIter->Next();
    //		e[i] = rhos->getVal();
    //	}

//#pragma omp parallel for
    for(Int_t i=0;i<Nmc;i++)
    {
        TLorentzVector l1,l2,l3,l4,l5,l23,l45,ltot;
        l1.SetPx(mcp1[i][0]);l1.SetPy(mcp1[i][1]);l1.SetPz(mcp1[i][2]);l1.SetE(mcp1[i][3]);
        l2.SetPx(mcp2[i][0]);l2.SetPy(mcp2[i][1]);l2.SetPz(mcp2[i][2]);l2.SetE(mcp2[i][3]);
        l3.SetPx(mcp3[i][0]);l3.SetPy(mcp3[i][1]);l3.SetPz(mcp3[i][2]);l3.SetE(mcp3[i][3]);
        l4.SetPx(mcp4[i][0]);l4.SetPy(mcp4[i][1]);l4.SetPz(mcp4[i][2]);l4.SetE(mcp4[i][3]);
        l5.SetPx(mcp5[i][0]);l5.SetPy(mcp5[i][1]);l5.SetPz(mcp5[i][2]);l5.SetE(mcp5[i][3]);
        l23=l2+l3;
        l45=l4+l5;
        Double_t m23=l23.M();
        Double_t m45=l45.M();
        Double_t costhetax23=l23.CosTheta();
        Double_t costhetax45=l45.CosTheta();
        hmphi->Fill(l23.M(),1);
        hmf0->Fill(l45.M(),1);
        Double_t eva=calEva(mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3]);
        //std::cout << "eva: " << eva << std::endl;
        hm23->Fill(m23,eva);//total
        hm45->Fill(m45,eva);//total
        h23->Fill(costhetax23,eva);//total angle 
        h45->Fill(costhetax45,eva);//total angle 

        //		RooArgSet *Fb  = new RooArgSet();

        //		for (Int_t j=0; j<nStates; j++) {
        //			TString title  = titleList[j];
        //			TString nameF  = TString("F_")  + title;
        //			RooRealVar *rF = new RooRealVar(nameF, title,0,0,1e6);
        //			if (!Fb->add(*rF, kTRUE)) delete rF;
        //		}
        //calculate fraction one by one by setting other rho  to zero
        //		TIterator *iF = Fb->createIterator();
        //		while (RooRealVar *fb = (RooRealVar*)iF->Next()) {
        //			_rhoIter->Reset();
        //			TString fNameb = fb->GetTitle();
        //			for (Int_t j=0; j<nAmps; j++) {
        //				RooRealVar *rhob = (RooRealVar*)_rhoIter->Next();
        //				rhob->setVal(e[j]);
        //				TString nameb = titleListT[j];
        //				if (fNameb != nameb) {
        //					rhob->setVal(0.0);
        //				}
        //			}
        //			Double_t  numerator=calEva(mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3]);
        //			fb->setVal(numerator);
        //		}
        //		iF->Reset();
        //		RooRealVar *fb1 = (RooRealVar*)iF->Next();
        //		double n1=fb1->getVal();
        //		hmkk1->Fill(m23,n1); 
        //		RooRealVar *fb2 = (RooRealVar*)iF->Next();
        //		double n2=fb2->getVal();
        //		hmkk2->Fill(m23,n2); 
        //		RooRealVar *fb3 = (RooRealVar*)iF->Next();
        //		double n3=fb3->getVal();
        //		hmkk3->Fill(m23,n3); 
        ///recover
        //		_rhoIter->Reset();
        //		for (Int_t i=0; i<nAmps; i++) {
        //			RooRealVar *rhoc = (RooRealVar*)_rhoIter->Next();
        //			rhoc->setVal(e[i]);
        //		}
    }
    f->Write();
}

RooArgSet *DPFPWAPdf::fitFractions(const RooArgList& newPar, Bool_t print, ostream& os)
{
    cout << "Begin fitFraction" << endl;
    //  cout<<"haha: "<< __LINE__ << endl;
    Double_t a[100],b[100],c[100],d[100],e[100],f[100],g[100],h[100],j[100],k[100],l[100],m[100],n[100],o[100],p[100];
    //save the old result of parameters to a[],b[],c[],d[] // 
    TIterator* _spinIter = _spinList.createIterator();
    TIterator* _massIter = _massList.createIterator();
    TIterator* _mass2Iter = _mass2List.createIterator();
    TIterator* _widthIter = _widthList.createIterator();
    TIterator* _g1Iter = _g1List.createIterator();
    TIterator* _g2Iter = _g2List.createIterator();
    TIterator* _b1Iter = _b1List.createIterator();
    TIterator* _b2Iter = _b2List.createIterator();
    TIterator* _b3Iter = _b3List.createIterator();
    TIterator* _b4Iter = _b4List.createIterator();
    TIterator* _b5Iter = _b5List.createIterator();
    TIterator* _rhoIter = _rhoList.createIterator();
    TIterator* _phiIter = _phiList.createIterator();
    TIterator* _propIter = _propList.createIterator();
    //#pragma omp parallel for
    for (Int_t i=0; i<nAmps; i++) {
        RooRealVar *rhos = (RooRealVar*)_rhoIter->Next();
        a[i] = rhos->getVal();
        //cout << "rhos---------------- " << rhos->getVal() << endl;
    }
    for (Int_t i=0; i<nStates; i++) {
        RooRealVar *masss = (RooRealVar*)_massIter->Next();
        RooRealVar *phis = (RooRealVar*)_phiIter->Next();
        RooRealVar *spins = (RooRealVar*)_spinIter->Next();
        RooRealVar *propTypes= (RooRealVar*)_propIter->Next();
        b[i] = phis->getVal();
        c[i] = masss->getVal();
        f[i] = spins->getVal();
        g[i] = propTypes->getVal();
        //cout << "phis---------------- " << phis->getVal() << endl;
        //cout << "masss---------------- " << masss->getVal() << endl;
        //cout << "spins---------------- " << spins->getVal() << endl;
        //cout << "propTypes---------------- " << propTypes->getVal() << endl;
    }
    for (Int_t i=0; i<nStatesb; i++) {
        RooRealVar *b1s = (RooRealVar*)_b1Iter->Next();
        RooRealVar *b2s = (RooRealVar*)_b2Iter->Next();
        RooRealVar *b3s = (RooRealVar*)_b3Iter->Next();
        RooRealVar *b4s = (RooRealVar*)_b4Iter->Next();
        RooRealVar *b5s = (RooRealVar*)_b5Iter->Next();
        h[i] = b1s->getVal();
        j[i] = b2s->getVal();
        k[i] = b3s->getVal();
        l[i] = b4s->getVal();
        m[i] = b5s->getVal();
    }
    for (Int_t i=0; i<nStatesg1g2; i++) {
        RooRealVar *g1s = (RooRealVar*)_g1Iter->Next();
        RooRealVar *g2s = (RooRealVar*)_g2Iter->Next();
        o[i] = g1s->getVal();
        p[i] = g2s->getVal();
    }
    for (Int_t i=0; i<nStateswidth; i++) {
        RooRealVar *widths = (RooRealVar*)_widthIter->Next();
        d[i] = widths->getVal();
    }
    for (Int_t i=0; i<nStatesmass2; i++) {
        RooRealVar *mass2s = (RooRealVar*)_mass2Iter->Next();
        n[i] = mass2s->getVal();
    }

    //modify the value of rhos,phis,mass,width to disturbed value//
    TString in,ou;
    TIterator* parIter =newPar.createIterator() ;
    RooRealVar *par(0);
    while(0 != (par= (RooRealVar*)parIter->Next())) {
        TIterator* _spinIter = _spinList.createIterator();
        TIterator* _massIter = _massList.createIterator();
        TIterator* _mass2Iter = _mass2List.createIterator();
        TIterator* _widthIter = _widthList.createIterator();
        TIterator* _g1Iter = _g1List.createIterator();
        TIterator* _g2Iter = _g2List.createIterator();
        TIterator* _b1Iter = _b1List.createIterator();
        TIterator* _b2Iter = _b2List.createIterator();
        TIterator* _b3Iter = _b3List.createIterator();
        TIterator* _b4Iter = _b4List.createIterator();
        TIterator* _b5Iter = _b5List.createIterator();
        TIterator* _rhoIter = _rhoList.createIterator();
        TIterator* _phiIter = _phiList.createIterator();
        TIterator* _propIter = _propList.createIterator();
        RooRealVar *rhopar(0);
        in = par->GetName();
        while(0 != (rhopar= (RooRealVar*)_rhoIter->Next())) {
            ou = rhopar->GetName();
            if (in==ou) { rhopar->setVal(par->getVal()); }
        }
        RooRealVar *phipar(0);
        while(0 != (phipar= (RooRealVar*)_phiIter->Next())) {
            //cout << "phiparstart---------------- " << phipar->getVal() << endl;
            ou = phipar->GetName();
            if (in==ou) { phipar->setVal(par->getVal()); }
            //cout << "phiparlast---------------- " << phipar->getVal() << endl;
        }
        RooRealVar *spinpar(0);
        while(0 != (spinpar= (RooRealVar*)_spinIter->Next())) {
            ou = spinpar->GetName();
            if (in==ou) { spinpar->setVal(par->getVal()); }
        }
        RooRealVar *masspar(0);
        while(0 != (masspar= (RooRealVar*)_massIter->Next())) {
            ou = masspar->GetName();
            if (in==ou) { masspar->setVal(par->getVal()); }
        }
        RooRealVar *mass2par(0);
        while(0 != (mass2par= (RooRealVar*)_mass2Iter->Next())) {
            ou = mass2par->GetName();
            if (in==ou) { mass2par->setVal(par->getVal()); }
        }
        RooRealVar *widthpar(0);
        while(0 != (widthpar= (RooRealVar*)_widthIter->Next())) {
            ou = widthpar->GetName();
            if (in==ou) { widthpar->setVal(par->getVal()); }
        }
        RooRealVar *g1par(0);
        while(0 != (g1par= (RooRealVar*)_g1Iter->Next())) {
            ou = g1par->GetName();
            if (in==ou) { g1par->setVal(par->getVal()); }
        }
        RooRealVar *g2par(0);
        while(0 != (g2par= (RooRealVar*)_g2Iter->Next())) {
            ou = g2par->GetName();
            if (in==ou) { g2par->setVal(par->getVal()); }
        }
        RooRealVar *b1par(0);
        while(0 != (b1par= (RooRealVar*)_b1Iter->Next())) {
            ou = b1par->GetName();
            if (in==ou) { b1par->setVal(par->getVal()); }
        }
        RooRealVar *b2par(0);
        while(0 != (b2par= (RooRealVar*)_b2Iter->Next())) {
            ou = b2par->GetName();
            if (in==ou) { b2par->setVal(par->getVal()); }
        }
        RooRealVar *b3par(0);
        while(0 != (b3par= (RooRealVar*)_b3Iter->Next())) {
            ou = b3par->GetName();
            if (in==ou) { b3par->setVal(par->getVal()); }
        }
        RooRealVar *b4par(0);
        while(0 != (b4par= (RooRealVar*)_b4Iter->Next())) {
            ou = b4par->GetName();
            if (in==ou) { b4par->setVal(par->getVal()); }
        }
        RooRealVar *b5par(0);
        while(0 != (b5par= (RooRealVar*)_b5Iter->Next())) {
            ou = b5par->GetName();
            if (in==ou) { b5par->setVal(par->getVal()); }
        }
        RooRealVar *proppar(0);
        while(0 != (proppar= (RooRealVar*)_propIter->Next())) {
            ou = proppar->GetName();
            if (in==ou) { proppar->setVal(par->getVal()); }
        }
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
        delete _phiIter;
        delete _propIter; 
    }
    delete parIter;

    Double_t norm = analyticalIntegral(1,"norm");//total prob
    RooArgSet *Fb  = new RooArgSet();
    _rhoIter->Reset();
    //save modified rhos and phis to e[]
    //#pragma omp parallel for
    for (Int_t i=0; i<nAmps; i++) {
        RooRealVar *rhos = (RooRealVar*)_rhoIter->Next();
        e[i] = rhos->getVal();
    }

    for (Int_t i=0; i<nStates; i++) {
        TString title  = titleList[i];
        TString nameF  = TString("F_")  + title;
        RooRealVar *rF = new RooRealVar(nameF, title,0,0,1e6);
        if (!Fb->add(*rF, kTRUE)) delete rF;
    }
    //calculate fraction one by one by setting other rho  to zero
    TIterator *iF = Fb->createIterator();
    while (RooRealVar *fb = (RooRealVar*)iF->Next()) {
        _rhoIter->Reset();
        TString fNameb = fb->GetTitle();
        cout << "Fit Fraction ====> " << fNameb << endl;
        //#pragma omp parallel for
        for (Int_t i=0; i<nAmps; i++) {
            RooRealVar *rhob = (RooRealVar*)_rhoIter->Next();
            rhob->setVal(e[i]);
            TString nameb = titleListT[i];
            if (fNameb != nameb) {
                rhob->setVal(0.0);
            }
            //cout << "rho---------------- " << rhob->getVal() << endl;
        }
        Double_t  numerator= analyticalIntegral(1,"norm");
        fb->setVal(numerator/norm);
        //cout << "numer---------------- " << numerator << "  " << norm << endl;
    }
    //////////////////recover///////////////
    _spinIter->Reset();
    _massIter->Reset();
    _mass2Iter->Reset();
    _widthIter->Reset();
    _g1Iter->Reset();
    _g2Iter->Reset();
    _b1Iter->Reset();
    _b2Iter->Reset();
    _b3Iter->Reset();
    _b4Iter->Reset();
    _b5Iter->Reset();
    _rhoIter->Reset();
    _phiIter->Reset();
    _propIter->Reset();
    //#pragma omp parallel for
    for (Int_t i=0; i<nAmps; i++) {
        RooRealVar *rhoc = (RooRealVar*)_rhoIter->Next();
        rhoc->setVal(a[i]);
        //cout << "rhoc---------------- " << rhoc->getVal() << endl;
    }
    for (Int_t i=0; i<nStates; i++) {
        RooRealVar *massc = (RooRealVar*)_massIter->Next();
        RooRealVar *phic = (RooRealVar*)_phiIter->Next();
        RooRealVar *spinc = (RooRealVar*)_spinIter->Next();
        RooRealVar *propTypec= (RooRealVar*)_propIter->Next();
        phic->setVal(b[i]);
        massc->setVal(c[i]);
        spinc->setVal(f[i]);
        propTypec->setVal(g[i]);
        //cout << "phic---------------- " << phic->getVal() << endl;
        //cout << "massc---------------- " << massc->getVal() << endl;
        //cout << "spinc---------------- " << spinc->getVal() << endl;
        //cout << "propTypec---------------- " << propTypec->getVal() << endl;
    }
    for (Int_t i=0; i<nStatesb; i++) {
        RooRealVar *b1c = (RooRealVar*)_b1Iter->Next();
        RooRealVar *b2c = (RooRealVar*)_b2Iter->Next();
        RooRealVar *b3c = (RooRealVar*)_b3Iter->Next();
        RooRealVar *b4c = (RooRealVar*)_b4Iter->Next();
        RooRealVar *b5c = (RooRealVar*)_b5Iter->Next();
        b1c->setVal(h[i]);
        b2c->setVal(j[i]);
        b3c->setVal(k[i]);
        b4c->setVal(l[i]);
        b5c->setVal(m[i]);
        //cout << "b1c---------------- " << b1c->getVal() << endl;
        //cout << "b2c---------------- " << b2c->getVal() << endl;
        //cout << "b3c---------------- " << b3c->getVal() << endl;
        //cout << "b4c---------------- " << b4c->getVal() << endl;
        //cout << "b5c---------------- " << b5c->getVal() << endl;
    }
    for (Int_t i=0; i<nStatesg1g2; i++) {
        RooRealVar *g1c = (RooRealVar*)_g1Iter->Next();
        RooRealVar *g2c = (RooRealVar*)_g2Iter->Next();
        g1c->setVal(o[i]);
        g2c->setVal(p[i]);
        //cout << "g1c---------------- " << g1c->getVal() << endl;
        //cout << "g2c---------------- " << g2c->getVal() << endl;
    }
    for (Int_t i=0; i<nStateswidth; i++) {
        RooRealVar *widthc = (RooRealVar*)_widthIter->Next();
        widthc->setVal(d[i]);
        //cout << "widthc---------------- " << widthc->getVal() << endl;
    }
    for (Int_t i=0; i<nStatesmass2; i++) {
        RooRealVar *mass2c = (RooRealVar*)_mass2Iter->Next();
        mass2c->setVal(n[i]);
        //cout << "mass2c---------------- " << mass2c->getVal() << endl;
    }

    /////////////out put to file///////////////
    if (print) {
        TIterator *miter = Fb->createIterator();
        Double_t total=0;
        const RooRealVar *tval;
        //    cout<<"Fit Fractions :"<<endl;
        os<<endl;
        while ((tval = (RooRealVar*)miter->Next())) {
            cout<<"Fraction of "<<tval->GetTitle()<<"="<<tval->getVal()<<"\n"<<endl;
            os<<"Fraction of "<<tval->GetTitle()<<"="<<tval->getVal()<<"\n"<<endl;
            total += tval->getVal();
        }
        //   os<<endl;
        cout<<"total="<<total<<endl;
        os<<"Total = "<<total<<endl;
        delete miter;
    }
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
    delete _phiIter;
    delete _propIter; 
    delete parIter;
    return Fb;
}
