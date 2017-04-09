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

#include <TTree.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
#include <fstream>
#include "phikk_structure.h"
#include "phipipi_structure.h"
#include "kernel_calEva.h"
#include <sstream>
#include <iomanip>
#include "cu_PWA_PARAS.h"
#include "assert.h"
#include "MultDevice.h"
/*#include "RooAbsReal.h"
#include "RooListProxy.h"
#include "RooStringVar.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include<TLorentzVector.h>
#include<TStopwatch.h>
*/
//#define DEVICE_NUM 2
Double_t rk=0.493677,rp=0.13957018 ;

//ClassImp(DPFPWAPdf)
DPFPWAPdf::DPFPWAPdf(const char *name, const char *title,
        RooAbsReal& _idp,
        DPFPWAPoint *dp) :
            RooAbsPdf(name,title),
            idp("idp", "idp", this, _idp),
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
            _fracList("fracList", "list of fracs", this),
            _phiList("phiList", "list of phis", this),
            _propList("propList", "list of props", this)
{
    _dp         = dp;
    _amp.setdp(_dp);
    //  //cout<<"haha: "<< __LINE__ << endl;
    //cout << "(****************************)" << endl;
    idp.print(cout);
    //cout << "9999" << endl;
    initialize();
}

DPFPWAPdf::DPFPWAPdf(const DPFPWAPdf& other, const char* name) :
    RooAbsPdf(other,name),
    idp("idp", this, other.idp),
//    p11("p11",this,other.p11),
//    p12("p12",this,other.p12),
//    p13("p13",this,other.p13),
//    p14("p14",this,other.p14),
//    p21("p21",this,other.p21),
//    p22("p22",this,other.p22),
//    p23("p23",this,other.p23),
//    p24("p24",this,other.p24),
//    p31("p31",this,other.p31),
//    p32("p32",this,other.p32),
//    p33("p33",this,other.p33),
//    p34("p34",this,other.p34),
//    p41("p41",this,other.p41),
//    p42("p42",this,other.p42),
//    p43("p43",this,other.p43),
//    p44("p44",this,other.p44),
//    p51("p51",this,other.p51),
//    p52("p52",this,other.p52),
//    p53("p53",this,other.p53),
//    p54("p54",this,other.p54),
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
    _fracList("fracList", this,other._fracList),
    _phiList("phiList", this,other._phiList),
    _propList("propList", this,other._propList)
{
    //cout<<"haha: "<< __LINE__ << endl;
    ////cout << "idp=" << idp.getVal() << endl;
    idp.print(cout);
    //cout << endl;
    other.idp.print(cout);
    //cout << other.idp.max() << endl;
    _dp         = other._dp;
    _amp.setdp(_dp);
    //cout<<" DPFPWAPdf other haha: "<< __LINE__ << endl;
    Nmc=other.Nmc;
    Nmc_data = other.Nmc_data;
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
    mcp1=new Double_t*[Nmc + Nmc_data];
    mcp2=new Double_t*[Nmc + Nmc_data];
    mcp3=new Double_t*[Nmc + Nmc_data];
    mcp4=new Double_t*[Nmc + Nmc_data];
    mcp5=new Double_t*[Nmc + Nmc_data];
    for(Int_t i=0;i<Nmc + Nmc_data;i++){
        mcp1[i]=new Double_t[4];
        mcp2[i]=new Double_t[4];
        mcp3[i]=new Double_t[4];
        mcp4[i]=new Double_t[4];
        mcp5[i]=new Double_t[4];
    }
    for(Int_t i=0;i<Nmc + Nmc_data;i++){
        mcp1[i][0]=other.mcp1[i][0];mcp1[i][1]=other.mcp1[i][1];mcp1[i][2]=other.mcp1[i][2];mcp1[i][3]=other.mcp1[i][3];
        mcp2[i][0]=other.mcp2[i][0];mcp2[i][1]=other.mcp2[i][1];mcp2[i][2]=other.mcp2[i][2];mcp2[i][3]=other.mcp2[i][3];
        mcp3[i][0]=other.mcp3[i][0];mcp3[i][1]=other.mcp3[i][1];mcp3[i][2]=other.mcp3[i][2];mcp3[i][3]=other.mcp3[i][3];
        mcp4[i][0]=other.mcp4[i][0];mcp4[i][1]=other.mcp4[i][1];mcp4[i][2]=other.mcp4[i][2];mcp4[i][3]=other.mcp4[i][3];
        mcp5[i][0]=other.mcp5[i][0];mcp5[i][1]=other.mcp5[i][1];mcp5[i][2]=other.mcp5[i][2];mcp5[i][3]=other.mcp5[i][3];
    }
    mlk = new Double_t*[Nmc + Nmc_data];
    for(int i = 0; i < Nmc + Nmc_data; i++) {
        mlk[i] = new Double_t[nAmps];
    }
//    mcp1_data=new Double_t*[Nmc_data];
//    mcp2_data=new Double_t*[Nmc_data];
//    mcp3_data=new Double_t*[Nmc_data];
//    mcp4_data=new Double_t*[Nmc_data];
//    mcp5_data=new Double_t*[Nmc_data];
//    for(Int_t i=0;i<Nmc;i++){
//        mcp1_data[i]=new Double_t[4];
//        mcp2_data[i]=new Double_t[4];
//        mcp3_data[i]=new Double_t[4];
//        mcp4_data[i]=new Double_t[4];
//        mcp5_data[i]=new Double_t[4];
//    }
//    for(Int_t i=0;i<Nmc_data;i++){
//        mcp1_data[i][0]=other.mcp1_data[i][0];mcp1_data[i][1]=other.mcp1_data[i][1];mcp1_data[i][2]=other.mcp1_data[i][2];mcp1_data[i][3]=other.mcp1_data[i][3];
//        mcp2_data[i][0]=other.mcp2_data[i][0];mcp2_data[i][1]=other.mcp2_data[i][1];mcp2_data[i][2]=other.mcp2_data[i][2];mcp2_data[i][3]=other.mcp2_data[i][3];
//        mcp3_data[i][0]=other.mcp3_data[i][0];mcp3_data[i][1]=other.mcp3_data[i][1];mcp3_data[i][2]=other.mcp3_data[i][2];mcp3_data[i][3]=other.mcp3_data[i][3];
//        mcp4_data[i][0]=other.mcp4_data[i][0];mcp4_data[i][1]=other.mcp4_data[i][1];mcp4_data[i][2]=other.mcp4_data[i][2];mcp4_data[i][3]=other.mcp4_data[i][3];
//        mcp5_data[i][0]=other.mcp5_data[i][0];mcp5_data[i][1]=other.mcp5_data[i][1];mcp5_data[i][2]=other.mcp5_data[i][2];mcp5_data[i][3]=other.mcp5_data[i][3];
//    }
    //cout << "Begin === Reset!" << endl;
    //
    pwa_paras.resize(0);
    fx.resize(0);
    cout << "\npwa_paras workd!!!!!" << other.pwa_paras.size() << endl;
    //ofstream cout("data_pwa_paras");
    for(int i = 0; i <  other.pwa_paras.size(); i++) {
        pwa_paras.push_back(other.pwa_paras[i]);
      //  Double_t * temp_pwa_paras_pointer=(Double_t *)&other.pwa_paras[i];
       // int size_PWA_PARAS=sizeof(PWA_PARAS)/sizeof(Double_t);
        //for(int j=0;j<size_PWA_PARAS;j++)
          //  {
            //    cout << temp_pwa_paras_pointer[j] << " ";
        //    }
        //cout << endl;
        fx.push_back(other.fx[i]);
     }
    #ifdef GPU
    //init d_float_pp
    d_float_pp.resize(DEVICE_NUM);
    int i_End= other.pwa_paras.size();
    int array_num = sizeof(cu_PWA_PARAS) / sizeof(double);
    int array_size = array_num * i_End;
    h_float_pp = new double[array_size];
    for(int i=0;i<i_End;i++)
    {
        Double_t * k=(Double_t*)&pwa_paras[i];
        for(int j=0;j<array_num;j++)
            h_float_pp[i*array_num+j]=(double)k[j];
    }
    for(int i=0;i<DEVICE_NUM;i++)
    {
        double * temp_pp;
        cu_malloc_h_pp(h_float_pp,temp_pp,pwa_paras.size(),i);
        d_float_pp[i]=temp_pp;
    }
    free(h_float_pp);
    #endif
    //cout.close();
    setup_iter_vec();
    //cout << "LINE: " << __LINE__ << endl;
    paras_getval();

    //cout << "DPFPWAPdf other haha : = " << __LINE__ << endl;
}
void DPFPWAPdf::setup_iter_vec() {
    _spinIter = _spinList.createIterator();
    _massIter = _massList.createIterator();
    _mass2Iter = _mass2List.createIterator();
    _widthIter = _widthList.createIterator();
    _g1Iter = _g1List.createIterator();
    _g2Iter = _g2List.createIterator();
    _b1Iter = _b1List.createIterator();
    _b2Iter = _b2List.createIterator();
    _b3Iter = _b3List.createIterator();
    _b4Iter = _b4List.createIterator();
    _b5Iter = _b5List.createIterator();
    _rhoIter = _rhoList.createIterator();
    _fracIter = _fracList.createIterator();
    _phiIter = _phiList.createIterator();
    _propIter = _propList.createIterator();

    _CN_spinList = 0;
    _CN_massList = _CN_spinList + _spinList.getSize();
    _CN_mass2List = _CN_massList + _massList.getSize();
    _CN_widthList = _CN_mass2List + _mass2List.getSize();
    _CN_g1List = _CN_widthList + _widthList.getSize();
    _CN_g2List = _CN_g1List + _g1List.getSize();
    _CN_b1List = _CN_g2List + _g2List.getSize();
    _CN_b2List = _CN_b1List + _b1List.getSize();
    _CN_b3List = _CN_b2List + _b2List.getSize();
    _CN_b4List = _CN_b3List + _b3List.getSize();
    _CN_b5List = _CN_b4List + _b4List.getSize();
    _CN_rhoList = _CN_b5List + _b5List.getSize();
    _CN_fracList = _CN_rhoList + _rhoList.getSize();
    _CN_phiList = _CN_fracList + _fracList.getSize();
    _CN_propList = _CN_phiList + _phiList.getSize();
    _CN_end = _CN_propList + _propList.getSize();

    paraList.resize(_CN_end);
}

void DPFPWAPdf::paras_getval() const {
    _spinIter->Reset();
    _massIter->Reset();
    _mass2Iter->Reset();
    _widthIter->Reset();
    _g1Iter->Reset();
    _g2Iter->Reset();
    _b1Iter->Reset();
    _b3Iter->Reset();
    _b2Iter->Reset();
    _b4Iter->Reset();
    _b5Iter->Reset();
    _rhoIter->Reset();
    _fracIter->Reset();
    _phiIter->Reset();
    _propIter->Reset();
    for(int i = 0; i < _spinList.getSize(); i++) {
        paraList[i + _CN_spinList] = ((RooRealVar*)_spinIter->Next())->getVal();
    }
    for(int i = 0; i < _massList.getSize(); i++) {
        paraList[i + _CN_massList] = ((RooRealVar*)_massIter->Next())->getVal();
    }
    for(int i = 0; i < _mass2List.getSize(); i++) {
        paraList[i + _CN_mass2List] = ((RooRealVar*)_mass2Iter->Next())->getVal();
    }
    for(int i = 0; i < _widthList.getSize(); i++) {
        paraList[i + _CN_widthList] = ((RooRealVar*)_widthIter->Next())->getVal();
    }
    for(int i = 0; i < _g1List.getSize(); i++) {
        paraList[i + _CN_g1List] = ((RooRealVar*)_g1Iter->Next())->getVal();
    }
    for(int i = 0; i < _g2List.getSize(); i++) {
        paraList[i + _CN_g2List] = ((RooRealVar*)_g2Iter->Next())->getVal();
    }
    for(int i = 0; i < _b1List.getSize(); i++) {
        paraList[i + _CN_b1List] = ((RooRealVar*)_b1Iter->Next())->getVal();
    }
    for(int i = 0; i < _b2List.getSize(); i++) {
        paraList[i + _CN_b2List] = ((RooRealVar*)_b2Iter->Next())->getVal();
    }
    for(int i = 0; i < _b3List.getSize(); i++) {
        paraList[i + _CN_b3List] = ((RooRealVar*)_b3Iter->Next())->getVal();
    }
    for(int i = 0; i < _b4List.getSize(); i++) {
        paraList[i + _CN_b4List] = ((RooRealVar*)_b4Iter->Next())->getVal();
    }
    for(int i = 0; i < _b5List.getSize(); i++) {
        paraList[i + _CN_b5List] = ((RooRealVar*)_b5Iter->Next())->getVal();
    }
    for(int i = 0; i < _rhoList.getSize(); i++) {
        paraList[i + _CN_rhoList] = ((RooRealVar*)_rhoIter->Next())->getVal();
    }
    for(int i = 0; i < _fracList.getSize(); i++) {
        paraList[i + _CN_fracList] = ((RooRealVar*)_fracIter->Next())->getVal();
    }
    for(int i = 0; i < _phiList.getSize(); i++) {
        paraList[i + _CN_phiList] = ((RooRealVar*)_phiIter->Next())->getVal();
    }
    for(int i = 0; i < _propList.getSize(); i++) {
        paraList[i + _CN_propList] = ((RooRealVar*)_propIter->Next())->getVal();
    }
}


int DPFPWAPdf::count_lines(TString df) {
    ifstream in(df);
    string line;
    int n = 0;
    while (getline(in, line)) {
        n++;
    }
    //cout << df << " have " << n << " lines." << endl;
    return n;
}

void DPFPWAPdf::initialize()
{
    //  //cout<<"haha: "<< __LINE__ << endl;
    //    Nmc=300000;
    //    Nmc = (count_lines() / 5 / 80 + 1) * 80;
    //cout << "phsp file is " << _dp->_phspfile << endl;
    //cout << "data file is " << _dp->_datafile << endl;
    Nmc = count_lines(_dp->_phspfile) / 5;
    Nmc_data = count_lines(_dp->_datafile) / 5;
    //cout << "LINE" << __LINE__ << "   " << "Nmc = " << Nmc << endl;
    //cout << "LINE" << __LINE__ << "   " << "Nmc_data = " << Nmc_data << endl;
    nAmps=0;
    nStates=0;
    nStatesb=0;
    nStatesg1g2=0;
    nStateswidth=0;
    nStatesmass2=0;
    nameList=new TString[60];
    titleList=new TString[60];
    titleListT=new TString[60];

    //cout << "Start initialize mcp" << endl;
    mcp1=new Double_t*[Nmc + Nmc_data];
    mcp2=new Double_t*[Nmc + Nmc_data];
    mcp3=new Double_t*[Nmc + Nmc_data];
    mcp4=new Double_t*[Nmc + Nmc_data];
    mcp5=new Double_t*[Nmc + Nmc_data];
    for(Int_t i=0;i<Nmc + Nmc_data;i++){
        mcp1[i]=new Double_t[4];
        mcp2[i]=new Double_t[4];
        mcp3[i]=new Double_t[4];
        mcp4[i]=new Double_t[4];
        mcp5[i]=new Double_t[4];
    }
    Double_t fx1,fy1,fz1,ft1,fx2,fy2,fz2,ft2,fx3,fy3,fz3,ft3,fx4,fy4,fz4,ft4,fx5,fy5,fz5,ft5;
    Double_t dweight;
    FILE *fp;
    if((fp=fopen(_dp->_phspfile,"r"))==NULL)
    {printf("can't open input file");
        return;
    }
    //cout << "------->start input mcp(pshp)" << _dp->_phspfile << endl;
    Int_t i=0;
    while(fscanf(fp,"%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n",&fx1,&fy1,&fz1,&ft1,&fx2,&fy2,&fz2,&ft2,&fx3,&fy3,&fz3,&ft3,&fx4,&fy4,&fz4,&ft4,&fx5,&fy5,&fz5,&ft5)!=EOF)
    {
        //  //cout<<"haha: "<< __LINE__ << endl;
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
        //cout << "There is memory allocated error!" << endl;
        exit(1);
    }

//    mcp1_data=new Double_t*[Nmc_data];
//    mcp2_data=new Double_t*[Nmc_data];
//    mcp3_data=new Double_t*[Nmc_data];
//    mcp4_data=new Double_t*[Nmc_data];
//    mcp5_data=new Double_t*[Nmc_data];
//    for(Int_t i=0;i<Nmc;i++){
//        mcp1_data[i]=new Double_t[4];
//        mcp2_data[i]=new Double_t[4];
//        mcp3_data[i]=new Double_t[4];
//        mcp4_data[i]=new Double_t[4];
//        mcp5_data[i]=new Double_t[4];
//    }
    //Double_t fx1,fy1,fz1,ft1,fx2,fy2,fz2,ft2,fx3,fy3,fz3,ft3,fx4,fy4,fz4,ft4,fx5,fy5,fz5,ft5;
    //FILE *fp;

    mlk = new Double_t*[Nmc + Nmc_data];
    for(int i = 0; i < Nmc + Nmc_data; i++) {
        mlk[i] = new Double_t[nAmps];
    }

    if((fp=fopen(_dp->_datafile,"r"))==NULL)
    {printf("can't open input file");
        return;
    }
    //cout << "------->start input mcp(data)" << _dp->_datafile << endl;
    //Int_t i=0;
    i = Nmc;
    while(fscanf(fp,"%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n",&fx1,&fy1,&fz1,&ft1,&fx2,&fy2,&fz2,&ft2,&fx3,&fy3,&fz3,&ft3,&fx4,&fy4,&fz4,&ft4,&fx5,&fy5,&fz5,&ft5)!=EOF)
    {
        //  //cout<<"haha: "<< __LINE__ << endl;
        mcp1[i][0]=fx1;mcp1[i][1]=fy1;mcp1[i][2]=fz1;mcp1[i][3]=ft1;
        mcp2[i][0]=fx2;mcp2[i][1]=fy2;mcp2[i][2]=fz2;mcp2[i][3]=ft2;
        mcp3[i][0]=fx3;mcp3[i][1]=fy3;mcp3[i][2]=fz3;mcp3[i][3]=ft3;
        mcp4[i][0]=fx4;mcp4[i][1]=fy4;mcp4[i][2]=fz4;mcp4[i][3]=ft4;
        mcp5[i][0]=fx5;mcp5[i][1]=fy5;mcp5[i][2]=fz5;mcp5[i][3]=ft5;
        i++;
    }
    fclose(fp);

    if (i != Nmc + Nmc_data) {
        //cout << "There is memory allocated error!" << endl;
        exit(1);
    }
    //cout << "Begin Reset!" << endl;

    setup_iter_vec();
    paras_getval();

    fx.resize(Nmc + Nmc_data, 0);
    //for(int i = 0; i < Nmc + Nmc_data; i++) {
    //    fx.push_back(0);
    //}




    store_pwa_paras();
    //cout << "store pwa paras" << endl;
    //update_fx = true;
}

Double_t DPFPWAPdf::evaluate() const
{
    //  //cout<<"haha: "<< __LINE__ << endl;
    ////cout<<"haha: "<< __LINE__ << endl;
    //return evaluate(p11,p12,p13,p14,p21,p22,p23,p24,p31,p32,p33,p34,p41,p42,p43,p44,p51,p52,p53,p54);
    ////cout << idp << endl;
    //return DPFPWAPdf::evaluate((idp.absArg())->getVal());
    //RooRealProxy cc("cc", "cc", RooRealVar("b", "b", 1, 2));
    ////cout << idp.name() << endl;
    //ofstream ff("tt.dat");
    //idp.print(ff);
    //ff.close();
    return DPFPWAPdf::evaluate(idp);
//    //cout << "Begin evaluate" << endl;
//    Double_t sum = 0;
//    const Double_t ee = 1E-10;
//#pragma omp parallel for reduction(+:sum)
//    for(Int_t i = 0; i < Nmc_data; i++) {
//        ////cout << "i = " << i << endl;
//        ////cout << "length = " << pwa_paras_data.size() << endl;
//        ////cout << "haha: " << __LINE__ << "PWA_PARAS_data : " << pwa_paras_data[i].sv << endl;
//        Double_t eva = ((int)log(calEva(pwa_paras_data[i])) / ee ) * ee;
//    //    //cout << "eva = " << eva << endl;
//        sum = sum + eva;
//   //     //cout << "sum = " << sum << endl;
//    }
//    //cout << "final sum = " << sum << endl;
//    return - sum;
//    //return exp(- sum / analyticalIntegral(1, ""));
}
void DPFPWAPdf::cu_init_data(int * &h_parameter,double * &h_paraList,double *&h_fx,double * &h_mlk,int i_End) const
{
    //int array_num = sizeof(cu_PWA_PARAS) / sizeof(double);
    //int array_size = array_num * i_End;
    //int mem_size = array_size * sizeof(double);
    ////init h_float_pp
    //h_float_pp = (double *)malloc(mem_size);
    //for(int i=0;i<i_End;i++)
    //{
    //    Double_t * k=(Double_t*)&pwa_paras[i];
    //    for(int j=0;j<array_num;j++)
    ///        h_float_pp[i*array_num+j]=(double)k[j];
    //}
    //cout << h_float_pp << endl;
    //init h_parameter
    h_parameter=(int *)malloc(18*sizeof(int));

    h_parameter[0] =  _CN_spinList;
    h_parameter[1] =  _CN_massList;
    h_parameter[2] =  _CN_mass2List;
    h_parameter[3] =  _CN_widthList;
    h_parameter[4] =  _CN_g1List;
    h_parameter[5] =  _CN_g2List;
    h_parameter[6] =  _CN_b1List;
    h_parameter[7] =  _CN_b2List;
    h_parameter[8] =  _CN_b3List;
    h_parameter[9] =  _CN_b4List;
    h_parameter[10] =  _CN_b5List;
    h_parameter[11] =  _CN_rhoList;
    h_parameter[12] =  _CN_fracList;
    h_parameter[13] =  _CN_phiList;
    h_parameter[14] =  _CN_propList;
    h_parameter[15] =  nAmps;
    h_parameter[16] =  Nmc;
    h_parameter[17] = Nmc_data;
    //init h_paraList
    h_paraList=(double *)malloc(paraList.size()*sizeof(double));
    for(int i=0;i<paraList.size();i++)
    {
        h_paraList[i]=paraList[i];
    }
    //init h_fx
    //h_fx=(double *)malloc(i_End*sizeof(double));
    h_fx=new double[i_End];
    //init h_mlk
    h_mlk = new double[(Nmc + Nmc_data)*nAmps];
    //init h_paraList
    h_paraList=(double *)malloc(paraList.size()*sizeof(double));
    for(int i=0;i<paraList.size();i++)
    {
        h_paraList[i]=paraList[i];
    }
}
void DPFPWAPdf::store_fx(int iBegin, int iEnd) const {
    paras_getval();
//#pragma omp parallel
    //for(int i = 0; i < Nmc + Nmc_data; i++) {
#ifdef CPU
    clock_t start,end;
    start= clock();
    for(int i = iBegin; i < iEnd; i++) {
        double sum = calEva(pwa_paras[i], i);
        fx[i] = (sum <= 0) ? 1e-20 : sum;
        fx[i] = sum;
    }
    end=clock();
    cout << "cpu part  time :" <<(double)(end-start)/CLOCKS_PER_SEC << "S" << endl;
#endif
    //gpu part//
#ifdef GPU
    for(int i = iBegin; i < iEnd; i++) {
        double sum = calEva(pwa_paras[i], i);
        fx[i] = (sum <= 0) ? 1e-20 : sum;
        fx[i] = sum;
    }
    clock_t start,end;
    start= clock();
    int *h_parameter;
    double *h_paraList;
    double *h_fx;
    double *h_mlk;
    //cout << "\niEnd : " << iEnd << endl;
    cu_init_data(h_parameter,h_paraList,h_fx,h_mlk,iEnd);
    host_store_fx(d_float_pp,h_parameter,h_paraList,paraList.size(),h_fx,h_mlk,iEnd,iBegin);
    int error_num=0;
    double abs_error;
    double total_error=0.0;
    for(int i = 0; i < Nmc + Nmc_data; i++) {
        for(int j=0;j<nAmps;j++)
        {
            //if(abs(mlk[i][j]-h_mlk[i*nAmps+j])>0.0001) assert(0);
            //mlk[i][j]=h_mlk[i*nAmps+j];
            abs_error=abs(mlk[i][j]-h_mlk[i*nAmps+j]);
            if(abs_error>0.0000001)
            {
                error_num++;
                //if(i==413) printf("i : %d\n index : %d\n",i,j);
            }
            total_error+=abs_error;
        }
    }
    cout << "mlk error more than 0.000001 : " << error_num*100.0/(nAmps*(Nmc+Nmc_data)) << "\%  ave_error : "<< total_error/(nAmps*(Nmc+Nmc_data)) << endl;
    error_num=0;
    total_error=0.0;
    for(int i=iBegin;i<iEnd;i++)
    {
        h_fx[i] = (h_fx[i] <= 0) ? 1e-20 : h_fx[i];
        abs_error=abs(fx[i]-h_fx[i]);
        if(abs_error>0.000001)
        {
           //printf("%d\n",i);
            error_num++;
        }
        total_error+=abs_error;
        //if(abs(fx[i]-h_fx[i])>0.0001) assert(0);
        //fx[i]=(h_fx[i] <= 0)? 1e-20 : h_fx[i];
    }/*
    if(error_num>iEnd/2)
    {
        for(int i=iBegin;i<iEnd;i++)
        {
            cout << i <<": fx[i] : " << fx[i] << "   h_fx[i] : " << h_fx[i]<< endl;
        }
    }*/
    cout << "fx error more than 0.000001 : " << error_num*100.0/iEnd << "\%  ave_error : "<< total_error/iEnd << endl;

    //free memory
    free(h_parameter);
    free(h_paraList);
    free(h_fx);
    free(h_mlk);
    end=clock();
    cout << "gpu part  time :" <<(double)(end-start)/CLOCKS_PER_SEC << "S" << endl;
#endif
    //gpu part end!//
    Double_t sum = 0;
    Double_t carry = 0;

//#pragma omp parallel for private(carry) reduction(+:sum)
    for(int i = 0; i < Nmc; i++)
    {
        //  //cout<<"haha: "<< __LINE__ << endl;
        Double_t y = fx[i] - carry;
        Double_t t = sum + y;
        carry = (t - sum) - y;
        sum = t; // Kahan Summation
    }
    anaIntegral = sum;

    sum = 0;
    for(int i = 0; i < nAmps; i++)
    {
        double tt = 0;
//#pragma omp parallel for reduction(+:tt)
        for(int j = 0; j < Nmc; j++)
        {
            tt += mlk[j][i];
        }
        sum += sqrt(tt / Nmc);
    }
    penalty = sum;

    sum = 0;
    for(int i = 0; i < nAmps; i++)
    {
        double tt = 0;
//#pragma omp parallel for reduction(+:tt)
        for(int j = Nmc; j < Nmc + Nmc_data; j++)
        {
            tt += mlk[j][i];
        }
        sum += sqrt(tt);
    }
    penalty_data = sum;
}

Double_t DPFPWAPdf::evaluate(int _idp) const
{
    //if (_idp == 0 || update_fx) {
    if (_idp == 0) {
        ////cout << "_idp = " << _idp << "update_fx = " << update_fx << endl;
        store_fx(0, Nmc + Nmc_data);
        ////cout << "anaIntegral = " << anaIntegral << endl;
        ////cout << "anaIntegral / Nmc = " << anaIntegral / Nmc << endl;
        ////cout << "fx[0] = " << fx[0];
        //return fx[Nmc];
        double lambda = 1e0;
        ////cout << "penalty = " << penalty << endl;
        ////cout << "penalty_data = " << penalty_data << endl;
        ////cout << "penalty / Nmc = " << penalty / Nmc << endl;
        //cout << "exp(-lambda * penalty) = " << exp(-lambda * penalty) << endl;
        //cout << "exp(-lambda * penalty_data) = " << exp(-lambda * penalty_data) << endl;
        //cout << "exp(-lambda * penalty / Nmc_data) = " << exp(-lambda * penalty / Nmc_data) << endl;
        //cout << "exp(-lambda * penalty_data / Nmc_data) = " << exp(-lambda * penalty_data / Nmc_data) << endl;
        ////cout << "exp(-lambda * penalty_data / Nmc_data) = " << exp(-lambda * penalty_data / Nmc_data) << endl;
        //return exp(- 0.5 * lambda * penalty) * anaIntegral / Nmc;
        //return exp(- 0.5 * lambda * penalty * anaIntegral / Nmc) * fx[Nmc];
        //return fx[Nmc] * anaIntegral / Nmc * exp(- lambda * penalty_data / Nmc_data);

    }
//            double sum = calEva(pwa_paras[Nmc + _idp]);
//     return (sum <= 0) ? 1e-20 : sum;
//    //cout << "XXXXXX fx - calEva" << "idp = " << _idp << "--->>" << (fx[Nmc + _idp] - sum) << endl;

        double lambda = 1e0;
    return fx[Nmc + _idp] * exp(- lambda * penalty / Nmc_data);
}


void DPFPWAPdf::addResonance600(const Char_t* name, const Char_t* title, RooAbsReal& spin,
        RooAbsReal& mass, RooAbsReal& b1, RooAbsReal& b2, RooAbsReal& b3,RooAbsReal& b4,RooAbsReal& b5,RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,RooAbsReal& propType)
{       _spinList.add(spin);
    //  //cout<<"haha: "<< __LINE__ << endl;
    _massList.add(mass);
    _b1List.add(b1);
    _b2List.add(b2);
    _b3List.add(b3);
    _b4List.add(b4);
    _b5List.add(b5);
    _rhoList.add(rho);
    _fracList.add(frac);
    _phiList.add(phi);
    _propList.add(propType);
    nameList[nAmps]= name;
    titleListT[nAmps]= title;
    titleList[nStates]= title;
    nAmps++;
    nStates++;
    nStatesb++;
}
void DPFPWAPdf::removeAll() {
        //cout <<"Begin removeAll !!!" << endl;
        nAmps = 0;
        nStates = 0;
        nStatesb = 0;
        nStatesg1g2 = 0;
        nStateswidth = 0;
        nStatesmass2 = 0;
        _spinList.removeAll();
        _massList.removeAll();
        _mass2List.removeAll();
        _widthList.removeAll();
        _g1List.removeAll();
        _g2List.removeAll();
        _b1List.removeAll();
        _b2List.removeAll();
        _b3List.removeAll();
        _b4List.removeAll();
        _b5List.removeAll();
        _rhoList.removeAll();
        _fracList.removeAll();
        _phiList.removeAll();
        _propList.removeAll();
        //cout << "End removeAll!!!" << endl;
}
void DPFPWAPdf::addResonance980(const Char_t* name, const Char_t* title, RooAbsReal& spin,
        RooAbsReal& mass, RooAbsReal& g1,RooAbsReal& g2,RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,RooAbsReal& propType)
{       _spinList.add(spin);
    //  //cout<<"haha: "<< __LINE__ << endl;
    _massList.add(mass);
    _g1List.add(g1);
    _g2List.add(g2);
    _rhoList.add(rho);
    _fracList.add(frac);
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
        RooAbsReal& mass1680, RooAbsReal& mass2,RooAbsReal& width,RooAbsReal& g1,RooAbsReal& g2,RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,
        RooAbsReal& propType)
{       _spinList.add(spin);
    //  //cout<<"haha: "<< __LINE__ << endl;
    _massList.add(mass1680);
    _mass2List.add(mass2);
    _widthList.add(width);
    _g1List.add(g1);
    _g2List.add(g2);
    _rhoList.add(rho);
    _fracList.add(frac);
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
void DPFPWAPdf::addResonance(const TString name, const TString title, RooAbsReal& spin,
        RooAbsReal& mass, RooAbsReal& width,RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,RooAbsReal& propType)
{       _spinList.add(spin);
    //  //cout<<"haha: "<< __LINE__ << endl;
    _massList.add(mass);
    _widthList.add(width);
    _rhoList.add(rho);
    _fracList.add(frac);
    _phiList.add(phi);
    _propList.add(propType);
    nameList[nAmps]= name;
    titleListT[nAmps]= title;
    titleList[nStates]= title;
    nAmps++;
    nStates++;
    nStateswidth++;
}
Double_t DPFPWAPdf::calEva(const PWA_PARAS &pp, int idp) const
    ////return square of complex amplitude
{
    //	static int A=0;
    //	A++;
    int _N_spinList     = _CN_spinList;
    int _N_massList     = _CN_massList;
    int _N_mass2List    = _CN_mass2List;
    int _N_widthList    = _CN_widthList;
    int _N_g1List       = _CN_g1List;
    int _N_g2List       = _CN_g2List;
    int _N_b1List       = _CN_b1List;
    int _N_b2List       = _CN_b2List;
    int _N_b3List       = _CN_b3List;
    int _N_b4List       = _CN_b4List;
    int _N_b5List       = _CN_b5List;
    int _N_rhoList      = _CN_rhoList;
    int _N_fracList     = _CN_fracList;
    int _N_phiList      = _CN_phiList;
    int _N_propList     = _CN_propList;
    Double_t value = 0.;
    TComplex fCF[nAmps][4];
    TComplex fCP[nAmps];
    TComplex pa[nAmps][nAmps];
    TComplex fu[nAmps][nAmps];
    TComplex crp1[nAmps],crp11[nAmps];
    TComplex cr0p11;
    TComplex ca2p1;
    TComplex cw2p11;
    TComplex cw2p12;
    TComplex cw2p15;
    TComplex cw;
    TComplex c1p12_12,c1p13_12,c1p12_13,c1p13_13,c1p12_14,c1p13_14;
    TComplex cr1m12_1,cr1m13_1;
    TComplex crpf1,crpf2;

    for(Int_t index=0; index<nAmps; index++) {
        Double_t rho0 = paraList[_N_rhoList++];
        Double_t frac0 = paraList[_N_fracList++];
        Double_t phi0 = paraList[_N_phiList++];
        Int_t spin_now = paraList[_N_spinList++];
        Int_t propType_now = paraList[_N_propList++];
    //cout<<"haha: "<< __LINE__ << endl;

        rho0 *= TMath::Exp(frac0);
        fCP[index]=TComplex(rho0*TMath::Cos(phi0),rho0*TMath::Sin(phi0));
        //        //cout<<"fCP[index]="<<fCP[index]<<endl;
        //std::cout << __FILE__ << __LINE__ << " : " << propType_now << std::endl;
        switch(propType_now)
        {
            //  //cout<<"haha: "<< __LINE__ << endl;
            //                     ordinary  Propagator  Contribution
            case 1:
                {
                    //RooRealVar *width = (RooRealVar*)_widthIterV[omp_id]->Next();
                    Double_t mass0 = paraList[_N_massList++];
                    Double_t width0 = paraList[_N_widthList++];
                    //					//cout<<"mass0="<<mass0<<endl;
                    //					//cout<<"width0="<<width0<<endl;
                    crp1[index]=_prop.propogator(mass0,width0,pp.s23);
                }
                break;
                //	Flatte   Propagator Contribution
            case 2:
                {
                    //RooRealVar *g1 = (RooRealVar*)_g1IterV[omp_id]->Next();
                    //RooRealVar *g2 = (RooRealVar*)_g2IterV[omp_id]->Next();
                    Double_t mass980 = paraList[_N_massList++];
                    Double_t g10 = paraList[_N_g1List++];
                    Double_t g20 = paraList[_N_g2List++];
                    //Double_t g10=g1->getVal();
                    //Double_t g20=g2->getVal();
     //               			//cout<<"mass980="<<mass980<<endl;
     //               			//cout<<"g10="<<g10<<endl;
     //               			//cout<<"g20="<<g20<<endl;
     //                           //cout<<"pp.s23="<<pp.s23<< endl;
                    crp1[index]=_prop.propogator980(mass980,g10,g20,pp.s23);
     //               			//cout<<"crp1[index]="<<crp1[index]<<endl;
                }
                break;
                // sigma  Propagator Contribution
            case 3:
                {
                    //RooRealVar *b1 = (RooRealVar*)_b1IterV[omp_id]->Next();
                    //RooRealVar *b2 = (RooRealVar*)_b2IterV[omp_id]->Next();
                    //RooRealVar *b3 = (RooRealVar*)_b3IterV[omp_id]->Next();
                    //RooRealVar *b4 = (RooRealVar*)_b4IterV[omp_id]->Next();
                    //RooRealVar *b5 = (RooRealVar*)_b5IterV[omp_id]->Next();
                    //Double_t mass600=mass->getVal();
                    //Double_t b10=b1->getVal();
                    //Double_t b20=b2->getVal();
                    //Double_t b30=b3->getVal();
                    //Double_t b40=b4->getVal();
                    //Double_t b50=b5->getVal();
                    Double_t mass600 = paraList[_N_massList++];
                    Double_t b10 = paraList[_N_b1List++];
                    Double_t b20 = paraList[_N_b2List++];
                    Double_t b30 = paraList[_N_b3List++];
                    Double_t b40 = paraList[_N_b4List++];
                    Double_t b50 = paraList[_N_b5List++];
                    crp1[index]=_prop.propogator600(mass600,b10,b20,b30,b40,b50,pp.s23);
                    //			//cout<<"crp1[index]3="<<crp1[index]<<endl;
                }
                break;
                // 1- or 1+  Contribution
            case 4:
                {
                    //RooRealVar *width = (RooRealVar*)_widthIterV[omp_id]->Next();
                    //Double_t mass0=mass->getVal();
                    //Double_t width0=width->getVal();
                    Double_t mass0 = paraList[_N_massList++];
                    Double_t width0 = paraList[_N_widthList++];
                    crp1[index]=_prop.propogator(mass0,width0,pp.sv2);
                    crp11[index]=_prop.propogator(mass0,width0,pp.sv3);
                }
                break;
                //  phi(1650) f0(980) include flatte and ordinary Propagator joint Contribution
            case 5:
                {
                    //RooRealVar *mass2  = (RooRealVar*)_mass2IterV[omp_id]->Next();
                    //RooRealVar *g1 = (RooRealVar*)_g1IterV[omp_id]->Next();
                    //RooRealVar *g2 = (RooRealVar*)_g2IterV[omp_id]->Next();
                    //Double_t mass980=mass2->getVal();
                    //Double_t g10=g1->getVal();
                    //Double_t g20=g2->getVal();
                    Double_t mass980 = paraList[_N_mass2List++];
                    Double_t g10 = paraList[_N_g1List++];
                    Double_t g20 = paraList[_N_g2List++];
                    //					//cout<<"mass980="<<mass980<<endl;
                    //					//cout<<"g10="<<g10<<endl;
                    //					//cout<<"g20="<<g20<<endl;
                    crp1[index]=_prop.propogator980(mass980,g10,g20,pp.sv);
                    //					//cout<<"crp1[index]="<<crp1[index]<<endl;
                    //RooRealVar *width = (RooRealVar*)_widthIterV[omp_id]->Next();
                    //Double_t mass1680=mass->getVal();
                    //Double_t width1680=width->getVal();
                    Double_t mass1680 = paraList[_N_massList++];
                    Double_t width1680 = paraList[_N_widthList++];
                    //					//cout<<"mass1680="<<mass1680<<endl;
                    //					//cout<<"width1680="<<width1680<<endl;
                    crp11[index]=_prop.propogator(mass1680,width1680,pp.s23);
                    //					//cout<<"crp11[index]="<<crp11[index]<<endl;
                }
                break;
            case 6:
                {
                    //RooRealVar *width = (RooRealVar*)_widthIterV[omp_id]->Next();
                    //Double_t mass0=mass->getVal();
                    //Double_t width0=width->getVal();
                    Double_t mass0 = paraList[_N_massList++];
                    Double_t width0 = paraList[_N_widthList++];
                    //					//cout<<"mass0="<<mass0<<endl;
                    //					//cout<<"width0="<<width0<<endl;
                    crp1[index]=_prop.propogator1270(mass0,width0,pp.s23);
                    //			//cout<<"crp1[index]6="<<crp1[index]<<endl;
                }
            default :
                ;
        }
        //if(idp ==1) printf("crp1 : %f\n",crp1[index].Re());
    //cout << "LINE: " << __LINE__ << endl;
        for(Int_t i=0;i<2;i++){
            ////cout<<"haha: "<< __LINE__ << endl;
            //		//cout<<"spin_now="<<spin_now<<endl;
            //if(idp == 413) printf("spin_now : %d\n",spin_now);
            switch(spin_now)
            {
                case 11:
                    //1+_1 contribution
                    fCF[index][i]=pp.w1p12_1[i]*crp1[index]+pp.w1p13_1[i]*crp11[index];

                    break;
                case 12:
                    //1+_2 contribution
                    c1p12_12=crp1[index]/pp.b2qbv2;
                    c1p13_12=crp11[index]/pp.b2qbv3;
                    fCF[index][i]=pp.w1p12_2[i]*c1p12_12+pp.w1p13_2[i]*c1p13_12;

                    break;
                case 13:
                    //1+_3 contribution
                    c1p12_13=crp1[index]/pp.b2qjv2;
                    c1p13_13=crp11[index]/pp.b2qjv3;
                    fCF[index][i]=pp.w1p12_3[i]*c1p12_13+pp.w1p13_3[i]*c1p13_13;

                    break;
                case 14:
                    //1+_4 contribution
                    c1p12_12=crp1[index]/pp.b2qbv2;
                    c1p13_12=crp11[index]/pp.b2qbv3;
                    c1p12_14=c1p12_12/pp.b2qjv2;
                    c1p13_14=c1p13_12/pp.b2qjv3;
                    fCF[index][i]=pp.w1p12_4[i]*c1p12_14+pp.w1p13_4[i]*c1p13_14;

                    break;
                case 111:
                    //1-__1 contribution
                    cr1m12_1=crp1[index]/pp.b1qjv2/pp.b1qbv2;
                    cr1m13_1=crp11[index]/pp.b1qjv3/pp.b1qbv3;
                    fCF[index][i]=pp.w1m12[i]*cr1m12_1+pp.w1m13[i]*cr1m13_1;

                    break;
                case 191:
                    //phi(1650)f0(980)_1 contribution
                    //		//cout<<"b1q2r23="<<b1q2r23<<endl;
                    crpf1=crp1[index]*crp11[index]/pp.b1q2r23;
                    //		//cout<<"crpf1="<<crpf1<<endl;
                    fCF[index][i]=pp.ak23w[i]*crpf1;
                    //	//cout<<"fCF[index][i]="<<fCF[index][i]<<endl;

                    break;
                case 192:
                    //phi(1650)f0(980)_2 contribution
                    crpf1=crp1[index]*crp11[index]/pp.b1q2r23;
                    crpf2=crpf1/pp.b2qjvf2;
                    fCF[index][i]=pp.wpf22[i]*crpf2;

                    break;
                case 1:
                    //  //cout<<"haha: "<< __LINE__ << endl;
                    //01 contribution
                    //	//cout<<"wu[i]="<<wu[i]<<endl;
                    //	//cout<<"crp1[index]="<<crp1[index]<<endl;
                    //	//cout<<"index="<<index<<endl;
                    fCF[index][i]=pp.wu[i]*crp1[index];
                    //	//cout<<"fCF[index][i]="<<fCF[index][i]<<endl;
                    //	//cout<<"i="<<i<<endl;

                    break;
                case 2:
                    //02 contribution
                    cr0p11=crp1[index]/pp.b2qjvf2;
                    fCF[index][i]=pp.w0p22[i]*cr0p11;
                    //	//cout<<"fCF[index][i]02="<<fCF[index][i]<<endl;

                    break;
                case 21:
                    //21 contribution
                    //if(idp == 413) cout<<"crp1="<<crp1[index].Re()<<endl;
                    cw2p11=crp1[index]/pp.b2qf2xx;
                    //cout<<"cw2p11="<<cw2p11<<endl;
                    //	//cout<<"w2p1[0]="<<w2p1[0]<<endl;
                    //	//cout<<"w2p1[1]="<<w2p1[1]<<endl;
                    fCF[index][i]=pp.w2p1[i]*cw2p11;
                    //if(idp==413) printf("cw2p11 = %.10f fcf = %.10f\n",cw2p11.Im(),fCF[index][i].Im());
                    //cout<<"fCF[index][i]21="<<fCF[index][i]<<endl;

                    break;
                case 22:
                    //22 contribution
                    cw2p11=crp1[index]/pp.b2qf2xx;
                    cw2p12=cw2p11/pp.b2qjvf2;
                    fCF[index][i]=pp.w2p2[i]*cw2p12;

                    break;
                case 23:
                    //23 contribution
                    cw2p11=crp1[index]/pp.b2qf2xx;
                    cw2p12=cw2p11/pp.b2qjvf2;
                    fCF[index][i]=pp.w2p3[i]*cw2p12;

                    break;
                case 24:
                    //24 contribution
                    cw2p11=crp1[index]/pp.b2qf2xx;
                    cw2p12=cw2p11/pp.b2qjvf2;
                    fCF[index][i]=pp.w2p4[i]*cw2p12;

                    break;
                case 25:
                    //25 contribution
                    cw2p11=crp1[index]/pp.b2qf2xx;
                    cw2p15=cw2p11/pp.b4qjvf2;
                    fCF[index][i]=pp.w2p5[i]*cw2p15;

                default:		;
            }
        }

    }

    double carry(0);
    //#pragma omp parallel for reduction(+:value)
    for(Int_t i=0;i<nAmps;i++){
        //  //cout<<"haha: "<< __LINE__ << endl;
        for(Int_t j=0;j<nAmps;j++){
            cw=fCP[i]*TComplex::Conjugate(fCP[j]);
            //    //cout<<"cw="<<cw<<endl;
            if(i==j) pa[i][j]=cw.Re();
            else if(i<j) pa[i][j]=2*cw.Re();
            else pa[i][j]=2*cw.Im();

            cw=TComplex(0,0);
            for(Int_t k=0;k<2;k++){
                cw+=fCF[i][k]*TComplex::Conjugate(fCF[j][k])/2.0;
                //   //cout<<"cwfu="<<cw<<endl;

            }
            if(i<=j) fu[i][j]=cw.Re();
            if(i>j) fu[i][j]=-cw.Im();
            //      //cout<<"pa[i][j]="<<pa[i][j]<<endl;
            //      //cout<<"fu[i][j]="<<fu[i][j]<<endl;
            Double_t temp = pa[i][j]*fu[i][j];
            double y = temp - carry;
            double t = value + y;
            carry = (t - value) - y;

            value = t; // Kahan Summation
        }
    }

    for(Int_t i=0;i<nAmps;i++){
        TComplex cw=fCP[i]*TComplex::Conjugate(fCP[i]);
        double pa=cw.Re();

        cw=TComplex(0,0);
        for(Int_t k=0;k<2;k++){
            cw+=fCF[i][k]*TComplex::Conjugate(fCF[i][k])/2.0;
        }
        double fu=cw.Re();
        mlk[idp][i] = pa * fu;
        if(idp==413 && i==3) printf("pa: %.10f fu: %.10f mlk : %.10f\n",pa,fu,mlk[idp][i]);
    }
    //    delete _spinIter_;
    //    delete _massIter_;
    //    delete _mass2Iter_;
    //    delete _widthIter_;
    //    delete _g1Iter_;
    //    delete _g2Iter_;
    //    delete _b1Iter_;
    //    delete _b2Iter_;
    //    delete _b3Iter_;
    //    delete _b4Iter_;
    //    delete _b5Iter_;
    //    delete _rhoIter_;
    //    delete _phiIter_;
    //    delete _propIter_;
    ////cout << "value = " << value << endl;
    //if(idp==1) cout << pp.wu[0] <<" "<<_N_spinList<<" "<<paraList[0]<<" "<<endl;
    return (value <= 0) ? 1e-20 : value;
}

Double_t DPFPWAPdf::scalar(Double_t *a1,Double_t *a2)const
{
    //  //cout<<"haha: "<< __LINE__ << endl;
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
    //cout<<"getAnalyticalIntegral haha: "<< __LINE__ << endl;
    RooArgSet theSet1;
    theSet1.add(RooArgSet(idp.arg()));
         //ofstream cout("data_of_priveate_member");
         //cout << _CN_spinList << std::endl;
        //cout <<  _CN_massList << std::endl;
        //cout <<  _CN_mass2List << std::endl;
        //cout <<  _CN_widthList << std::endl;
        //cout <<  _CN_g1List << std::endl;
        //cout <<  _CN_g2List << std::endl;
        //cout <<  _CN_b1List << std::endl;
        //cout <<  _CN_b2List << std::endl;
        //cout <<  _CN_b3List << std::endl;
        //cout <<  _CN_b4List << std::endl;
        //cout <<  _CN_b5List << std::endl;
        //cout <<  _CN_rhoList << std::endl;
       // cout <<  _CN_fracList << std::endl;
        //cout <<  _CN_phiList << std::endl;
        ///cout <<  _CN_propList << std::endl;
        //cout <<  nAmps << std::endl;
        //cout <<  Nmc << std::endl;
        //cout << Nmc_data << std::endl;
        //cout << "paraList.size():" << paraList.size()<<endl;
        //for(int i=0;i<paraList.size();i++)
        //{
            //cout << paraList[i] << endl;
        //}
        //cout.close();
    //ofstream cout("data_fx_result");
    if (matchArgs(allVars, analVars, theSet1)) {
        store_fx(0, Nmc + Nmc_data);
        //cout << "haha: " << __LINE__ << endl;
       // for(int ik=0;ik<Nmc+Nmc_data;ik++)
        //{
          //  cout << fx[ik] << endl;
       // }
       // cout.close();
        return 1;
    }
    return 0;
    //if  (rangeName) //cout << " " <<endl;
    //RooArgSet theSet1, theSet2, theSet3;
    //theSet1.add(RooArgSet(p11.arg(),p12.arg(),p13.arg(),p14.arg(),p21.arg(),p22.arg(),p23.arg(),p24.arg()));
    //theSet2.add(RooArgSet(p31.arg(),p32.arg(),p33.arg(),p34.arg(),p41.arg(),p42.arg(),p43.arg(),p44.arg()));
    //theSet3.add(RooArgSet(p51.arg(),p52.arg(),p53.arg(),p54.arg()));
    //RooArgSet theSet4(theSet1,theSet2," ");
    //RooArgSet theSet(theSet4, theSet3," ");
    //if (matchArgs(allVars, analVars, theSet)) { //cout<<"haha: " << __LINE__ <<endl;
    //    return 1;}
    //    return 0;
}
// 下面这个部分如何加速？
void DPFPWAPdf::store_pwa_paras() {
    pwa_paras.resize(0);
    for(Int_t i = 0; i < Nmc + Nmc_data; i++) {
        PWA_PARAS the_pwa_paras;
        _amp.calculate0p(
                mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],
                mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],
                mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],
                mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],
                mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3],
                the_pwa_paras);
        pwa_paras.push_back(the_pwa_paras);
        //fx.push_back(0);
        ////cout << "PWA_PARAS has " << the_pwa_paras.sv << endl;
    }
//    store_fx();
//    //cout << "store_fx fx[0] = " << fx[0] << endl;
//    //cout << "store_fx fx[1] = " << fx[1] << endl;
//    for(Int_t i = 0; i < Nmc_data; i++) {
//        PWA_PARAS the_pwa_paras;
//        _amp.calculate0p(
//                mcp1_data[i][0],mcp1_data[i][1],mcp1_data[i][2],mcp1_data[i][3],
//                mcp2_data[i][0],mcp2_data[i][1],mcp2_data[i][2],mcp2_data[i][3],
//                mcp3_data[i][0],mcp3_data[i][1],mcp3_data[i][2],mcp3_data[i][3],
//                mcp4_data[i][0],mcp4_data[i][1],mcp4_data[i][2],mcp4_data[i][3],
//                mcp5_data[i][0],mcp5_data[i][1],mcp5_data[i][2],mcp5_data[i][3],
//                the_pwa_paras);
//        pwa_paras_data.push_back(the_pwa_paras);
        ////cout << "PWA_PARAS_data has " << the_pwa_paras.sv << endl;
        ////cout << "LINE:"<< __LINE__ << pwa_paras_data[i].sv << endl;
//    }
    //    //cout << "length = " << pwa_paras_data.size() << endl;
}

Double_t DPFPWAPdf::analyticalIntegral (Int_t code, const char* rangeName) const
{
    ////cout<<"haha: "<< __LINE__ << endl;
    //assert (code==1);
    //Double_t sum=0;
    const Double_t ee = 1e9;
    ////cout << "Begin analyticalIntegral, Nmc = " << Nmc << endl;
    //double startTime=omp_get_wtime();
    //store_fx(0, Nmc);
  //store_fx(0, Nmc + Nmc_data);
//#pragma omp parallel for reduction(+:sum)
//    for(int i = 0; i < Nmc; i++)
//    {
//        //  //cout<<"haha: "<< __LINE__ << endl;
//        sum=sum + fx[i];
//    }
//  double sum_ = 0;
//  //store_fx(0, Nmc + Nmc_data);
//#pragma omp parallel for reduction(+:sum_)
//    for(int i = 0; i < Nmc; i++)
//    {
//        sum_ += fx[i];
//    }
//    stringstream ss, ss_;
//    ss << std::setprecision(14) << sum;
//    ss >> sum;
//    ss_ << setprecision(14) << sum_;
//    ss_ >> sum_;


    //sum = floor(sum * ee) * ee;
    //sum_ = floor(sum_ * ee) * ee;
//    if (sum - sum_ != 0) {
//    //cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
//    //cout << sum << endl;
//    //cout << sum_ << endl;
//    //cout << sum - sum_ << endl;
//    //cout << (sum - sum_) / sum << endl;
//    //cout << 1.0 * (sum - sum_) / (1.0 * Nmc) << endl;
//    //cout << 1.0 * (sum - sum_) << endl;
//    }
//    sum_ = floor(sum_ * 1e7) * 1e-7;
//    //cout << "sum = " << sum << endl;
//    //cout << "uuuu" << sum / Nmc - sum_ / Nmc << endl;
//    //cout << "xxxx" << (sum - sum_) / Nmc << endl;
    //sum=sum/Nmc;
    anaIntegral = anaIntegral / Nmc;
//    //cout << "sum / Nmc = " << sum << endl;
//    //cout << Nmc << endl;
    //double endTime=omp_get_wtime();
    ////cout << "analyticalIntegral sum = " << sum << endl;
    //    //cout << "End analyticalIntegral" << endl;
    //update_fx = true;
    ////cout << "final sum = " << ((sum <= 0) ? 1e-20 : sum) << endl;
    return anaIntegral;
}
void DPFPWAPdf::projectkk(const PWA_CTRL & pwa_ctrl) {
    TString phsp_weight_file_name = "";
    if (pwa_ctrl.actResList.size() > 1) {
        phsp_weight_file_name = work_path + "phsp_pwa_kk_weight_" + "all.root";
    } else {
        phsp_weight_file_name = work_path + "phsp_pwa_kk_weight_" + pwa_ctrl.actResList[0] + ".root";
    }
    DATA_PWA_KK ss;
    DATA_ANA_KK aa;
    TFile *fout = new TFile(phsp_weight_file_name, "RECREATE");
    TTree *pwa_tr = new TTree("pwa_tr", "pwa information");
    TTree *ana_tr = new TTree("ana_tr", "ana information");
    pwa_tr->Branch("kk", &ss, MEM_PWA_KK);
    ana_tr->Branch("kk", &aa, MEM_ANA_KK);
    for(Int_t i = 0; i < Nmc; i++) {
        ss.Kp2X = mcp2[i][0]; ss.Kp2Y = mcp2[i][1]; ss.Kp2Z = mcp2[i][2]; ss.Kp2E = mcp2[i][3];
        ss.Km2X = mcp3[i][0]; ss.Km2Y = mcp3[i][1]; ss.Km2Z = mcp3[i][2]; ss.Km2E = mcp3[i][3];
        ss.Kp1X = mcp4[i][0]; ss.Kp1Y = mcp4[i][1]; ss.Kp1Z = mcp4[i][2]; ss.Kp1E = mcp4[i][3];
        ss.Km1X = mcp5[i][0]; ss.Km1Y = mcp5[i][1]; ss.Km1Z = mcp5[i][2]; ss.Km1E = mcp5[i][3];
        PWA_PARAS pp;
        _amp.calculate0p(
                mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],
                mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],
                mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],
                mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],
                mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3],
                pp);
        ss.weight = calEva(pp, i);
        pwa_tr->Fill();
        DATA_ORIG_KK tt;
        pwa_to_orig(ss, tt);
        orig_to_ana(tt, aa);
        ana_tr->Fill();
    }
    fout->Write();
    fout->Close();
    //cout << phsp_weight_file_name << " is created!!!!" << endl;
}
void DPFPWAPdf::projectpipi(const PWA_CTRL & pwa_ctrl) {
    TString phsp_weight_file_name = "";
    if (pwa_ctrl.actResList.size() > 1) {
        phsp_weight_file_name = work_path + "phsp_pwa_pipi_weight_" + "all.root";
    } else {
        phsp_weight_file_name = work_path + "phsp_pwa_pipi_weight_" + pwa_ctrl.actResList[0] + ".root";
    }
    DATA_PWA_PIPI ss;
    DATA_ANA_PIPI aa;
    TFile *fout = new TFile(phsp_weight_file_name, "RECREATE");
    TTree *pwa_tr = new TTree("pwa_tr", "pwa information");
    TTree *ana_tr = new TTree("ana_tr", "ana information");
    pwa_tr->Branch("pp", &ss, MEM_PWA_PIPI);
    ana_tr->Branch("pp", &aa, MEM_ANA_PIPI);
    for(Int_t i = 0; i < Nmc; i++) {
        ss.pipX = mcp2[i][0]; ss.pipY = mcp2[i][1]; ss.pipZ = mcp2[i][2]; ss.pipE = mcp2[i][3];
        ss.pimX = mcp3[i][0]; ss.pimY = mcp3[i][1]; ss.pimZ = mcp3[i][2]; ss.pimE = mcp3[i][3];
        ss.KpX = mcp4[i][0]; ss.KpY = mcp4[i][1]; ss.KpZ = mcp4[i][2]; ss.KpE = mcp4[i][3];
        ss.KmX = mcp5[i][0]; ss.KmY = mcp5[i][1]; ss.KmZ = mcp5[i][2]; ss.KmE = mcp5[i][3];
        PWA_PARAS pp;
        _amp.calculate0p(
                mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],
                mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],
                mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],
                mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],
                mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3],
                pp);
        ss.weight = calEva(pp, i);
        pwa_tr->Fill();
        DATA_ORIG_PIPI tt;
        pwa_to_orig(ss, tt);
        orig_to_ana(tt, aa);
        ana_tr->Fill();
    }
    fout->Write();
    fout->Close();
    //cout << phsp_weight_file_name << " is created!!!!" << endl;
}

void DPFPWAPdf::project(const RooArgList& cPar, const RooArgList& fPar, const char* fname)
{
    //  //cout<<"haha: "<< __LINE__ << endl;
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
        PWA_PARAS pp;
        _amp.calculate0p(
                mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],
                mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],
                mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],
                mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],
                mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3],
                pp);
        Double_t eva = calEva(pp, i);
        ////std::cout << "eva: " << eva << std::endl;
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
void DPFPWAPdf::writeToFile(TString fname) {
    TFile *ff = new TFile(fname, "RECREATE");
    RooListProxy bb;
    bb.Copy(_spinList);               bb.Write("_spinList_");
    bb.Copy(_massList);               bb.Write("_massList_");
    bb.Copy(_mass2List);              bb.Write("_mass2List_");
    bb.Copy(_widthList);              bb.Write("_widthList_");
    bb.Copy(_g1List);                 bb.Write("_g1List_");
    bb.Copy(_g2List);                 bb.Write("_g2List_");
    bb.Copy(_b1List);                 bb.Write("_b1List_");
    bb.Copy(_b2List);                 bb.Write("_b2List_");
    bb.Copy(_b3List);                 bb.Write("_b3List_");
    bb.Copy(_b4List);                 bb.Write("_b4List_");
    bb.Copy(_b5List);                 bb.Write("_b5List_");
    bb.Copy(_rhoList);                bb.Write("_rhoList_");
    bb.Copy(_fracList);               bb.Write("_fracList_");
    bb.Copy(_phiList);                bb.Write("_phiList_");
    bb.Copy(_propList);               bb.Write("_propList_");
    ff->Close();
}
void DPFPWAPdf::readFromFile(TString fname) {
    TFile *ff = new TFile(fname);
    RooListProxy bb;
    bb.Read("_spinList_");           _spinList.Copy(bb);
    bb.Read("_massList_");           _massList.Copy(bb);
    bb.Read("_mass2List_");          _mass2List.Copy(bb);
    bb.Read("_widthList_");          _widthList.Copy(bb);
    bb.Read("_g1List_");             _g1List.Copy(bb);
    bb.Read("_g2List_");             _g2List.Copy(bb);
    bb.Read("_b1List_");             _b1List.Copy(bb);
    bb.Read("_b2List_");             _b2List.Copy(bb);
    bb.Read("_b3List_");             _b3List.Copy(bb);
    bb.Read("_b4List_");             _b4List.Copy(bb);
    bb.Read("_b5List_");             _b5List.Copy(bb);
    bb.Read("_rhoList_");            _rhoList.Copy(bb);
    bb.Read("_fracList_");           _fracList.Copy(bb);
    bb.Read("_phiList_");            _phiList.Copy(bb);
    bb.Read("_propList_");           _propList.Copy(bb);
    ff->Close();

}

void DPFPWAPdf::showAllParas() {
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
    TIterator* _fracIter = _fracList.createIterator();
    TIterator* _phiIter = _phiList.createIterator();
    TIterator* _propIter = _propList.createIterator();
    //#pragma omp parallel for
    for (Int_t i=0; i<nAmps; i++) {
        RooRealVar *rhos = (RooRealVar*)_rhoIter->Next();
        RooRealVar *fracs = (RooRealVar*)_fracIter->Next();
        //cout << nameList[i] << " " << titleList[i] << " rho is " << rhos->getVal() << endl;
        //cout << nameList[i] << " " << titleList[i] << " frac is " << fracs->getVal() << endl;
        ////cout << "rhos---------------- " << rhos->getVal() << endl;
    }
    for (Int_t i=0; i<nStates; i++) {
        RooRealVar *masss = (RooRealVar*)_massIter->Next();
        RooRealVar *phis = (RooRealVar*)_phiIter->Next();
        RooRealVar *spins = (RooRealVar*)_spinIter->Next();
        RooRealVar *propTypes= (RooRealVar*)_propIter->Next();
        //cout << titleList[i] << " phi is " << phis->getVal() << endl;
        //cout << titleList[i] << " mass is " << masss->getVal() << endl;
        //cout << titleList[i] << " spin is " << spins->getVal() << endl;
        //cout << titleList[i] << " propType is " << propTypes->getVal() << endl;
        ////cout << "phis---------------- " << phis->getVal() << endl;
        ////cout << "masss---------------- " << masss->getVal() << endl;
        ////cout << "spins---------------- " << spins->getVal() << endl;
        ////cout << "propTypes---------------- " << propTypes->getVal() << endl;
    }
    //    for (Int_t i=0; i<nStatesb; i++) {
    //        RooRealVar *b1s = (RooRealVar*)_b1Iter->Next();
    //        RooRealVar *b2s = (RooRealVar*)_b2Iter->Next();
    //        RooRealVar *b3s = (RooRealVar*)_b3Iter->Next();
    //        RooRealVar *b4s = (RooRealVar*)_b4Iter->Next();
    //        RooRealVar *b5s = (RooRealVar*)_b5Iter->Next();
    //        h[i] = b1s->getVal();
    //        j[i] = b2s->getVal();
    //        k[i] = b3s->getVal();
    //        l[i] = b4s->getVal();
    //        m[i] = b5s->getVal();
    //    }
    //    for (Int_t i=0; i<nStatesg1g2; i++) {
    //        RooRealVar *g1s = (RooRealVar*)_g1Iter->Next();
    //        RooRealVar *g2s = (RooRealVar*)_g2Iter->Next();
    //        o[i] = g1s->getVal();
    //        p[i] = g2s->getVal();
    //    }
    //    for (Int_t i=0; i<nStateswidth; i++) {
    //        RooRealVar *widths = (RooRealVar*)_widthIter->Next();
    //        d[i] = widths->getVal();
    //    }
    //    for (Int_t i=0; i<nStatesmass2; i++) {
    //        RooRealVar *mass2s = (RooRealVar*)_mass2Iter->Next();
    //        n[i] = mass2s->getVal();
    //    }


}
RooArgSet *DPFPWAPdf::fitFractions(const RooArgList& newPar, Bool_t print, ostream& os)
{
    //cout << "Begin fitFraction" << endl;
    //  //cout<<"haha: "<< __LINE__ << endl;
    Double_t a[300],aa[300],b[300],c[300],d[300],e[300],ee[300],f[300],g[300],h[300],j[300],k[300],l[300],m[300],n[300],o[300],p[300];
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
    TIterator* _fracIter = _fracList.createIterator();
    TIterator* _phiIter = _phiList.createIterator();
    TIterator* _propIter = _propList.createIterator();
    //#pragma omp parallel for
    for (Int_t i=0; i<nAmps; i++) {
        RooRealVar *rhos = (RooRealVar*)_rhoIter->Next();
        RooRealVar *fracs = (RooRealVar*)_fracIter->Next();
        a[i] = rhos->getVal();
        aa[i] = fracs->getVal();
        ////cout << "rhos---------------- " << rhos->getVal() << endl;
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
        ////cout << "phis---------------- " << phis->getVal() << endl;
        ////cout << "masss---------------- " << masss->getVal() << endl;
        ////cout << "spins---------------- " << spins->getVal() << endl;
        ////cout << "propTypes---------------- " << propTypes->getVal() << endl;
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
        _fracIter->Reset();
        _phiIter->Reset();
        _propIter->Reset();
        RooRealVar *rhopar(0);
        in = par->GetName();
        while(0 != (rhopar= (RooRealVar*)_rhoIter->Next())) {
            ou = rhopar->GetName();
            if (in==ou) { rhopar->setVal(par->getVal()); }
        }
        RooRealVar *fracpar(0);
        while(0 != (fracpar= (RooRealVar*)_fracIter->Next())) {
            ou = fracpar->GetName();
            if (in==ou) { fracpar->setVal(par->getVal()); }
        }
        RooRealVar *phipar(0);
        while(0 != (phipar= (RooRealVar*)_phiIter->Next())) {
            ////cout << "phiparstart---------------- " << phipar->getVal() << endl;
            ou = phipar->GetName();
            if (in==ou) { phipar->setVal(par->getVal()); }
            ////cout << "phiparlast---------------- " << phipar->getVal() << endl;
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
    }

    Double_t norm = analyticalIntegral(1,"norm");//total prob
    RooArgSet *Fb  = new RooArgSet();
    _rhoIter->Reset();
    _fracIter->Reset();
    //save modified rhos and phis to e[]
    //#pragma omp parallel for
    for (Int_t i=0; i<nAmps; i++) {
        RooRealVar *rhos = (RooRealVar*)_rhoIter->Next();
        RooRealVar *fracs = (RooRealVar*)_fracIter->Next();
        e[i] = rhos->getVal();
        ee[i] = fracs->getVal();
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
        _fracIter->Reset();
        TString fNameb = fb->GetTitle();
        //cout << "Fit Fraction ====> " << fNameb << endl;
        //#pragma omp parallel for
        for (Int_t i=0; i<nAmps; i++) {
            RooRealVar *rhob = (RooRealVar*)_rhoIter->Next();
            RooRealVar *fracb = (RooRealVar*)_fracIter->Next();
            rhob->setVal(e[i]);
            fracb->setVal(ee[i]);
            TString nameb = titleListT[i];
            if (fNameb != nameb) {
                rhob->setVal(0.0);
                fracb->setVal(0.0);
            }
            ////cout << "rho---------------- " << rhob->getVal() << endl;
        }
        Double_t  numerator= analyticalIntegral(1,"norm");
        fb->setVal(numerator/norm);
        ////cout << "numer---------------- " << numerator << "  " << norm << endl;
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
    _fracIter->Reset();
    _phiIter->Reset();
    _propIter->Reset();
    //#pragma omp parallel for
    for (Int_t i=0; i<nAmps; i++) {
        RooRealVar *rhoc = (RooRealVar*)_rhoIter->Next();
        RooRealVar *fracc = (RooRealVar*)_fracIter->Next();
        rhoc->setVal(a[i]);
        fracc->setVal(aa[i]);
        ////cout << "rhoc---------------- " << rhoc->getVal() << endl;
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
        ////cout << "phic---------------- " << phic->getVal() << endl;
        ////cout << "massc---------------- " << massc->getVal() << endl;
        ////cout << "spinc---------------- " << spinc->getVal() << endl;
        ////cout << "propTypec---------------- " << propTypec->getVal() << endl;
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
        ////cout << "b1c---------------- " << b1c->getVal() << endl;
        ////cout << "b2c---------------- " << b2c->getVal() << endl;
        ////cout << "b3c---------------- " << b3c->getVal() << endl;
        ////cout << "b4c---------------- " << b4c->getVal() << endl;
        ////cout << "b5c---------------- " << b5c->getVal() << endl;
    }
    for (Int_t i=0; i<nStatesg1g2; i++) {
        RooRealVar *g1c = (RooRealVar*)_g1Iter->Next();
        RooRealVar *g2c = (RooRealVar*)_g2Iter->Next();
        g1c->setVal(o[i]);
        g2c->setVal(p[i]);
        ////cout << "g1c---------------- " << g1c->getVal() << endl;
        ////cout << "g2c---------------- " << g2c->getVal() << endl;
    }
    for (Int_t i=0; i<nStateswidth; i++) {
        RooRealVar *widthc = (RooRealVar*)_widthIter->Next();
        widthc->setVal(d[i]);
        ////cout << "widthc---------------- " << widthc->getVal() << endl;
    }
    for (Int_t i=0; i<nStatesmass2; i++) {
        RooRealVar *mass2c = (RooRealVar*)_mass2Iter->Next();
        mass2c->setVal(n[i]);
        ////cout << "mass2c---------------- " << mass2c->getVal() << endl;
    }

    /////////////out put to file///////////////
    if (print) {
        TIterator *miter = Fb->createIterator();
        Double_t total=0;
        const RooRealVar *tval;
        //    //cout<<"Fit Fractions :"<<endl;
        os<<endl;
        while ((tval = (RooRealVar*)miter->Next())) {
            //cout<<"Fraction of "<<tval->GetTitle()<<"="<<tval->getVal()<<"\n"<<endl;
            os<<"Fraction of "<<tval->GetTitle()<<"="<<tval->getVal()<<"\n"<<endl;
            total += tval->getVal();
        }
        //   os<<endl;
        //cout<<"total="<<total<<endl;
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
    delete _fracIter;
    delete _phiIter;
    delete _propIter;
    delete parIter;
    return Fb;
}
void DPFPWAPdf::change_value(TString rn, double xx) {
    ((RooRealVar*)_massList.find(rn))->setVal(xx);
}
//inline DPFPWAPdf::~DPFPWAPdf() {
//    //cout << "deconstruct DPFPWAPdf" << endl;
//}
