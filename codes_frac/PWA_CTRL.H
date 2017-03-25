#ifndef PWA_CTRL_H
#define PWA_CTRL_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TString.h"

typedef struct {
    std::vector<std::string> paraConstList;
    std::vector<std::string> actResList;
    std::vector<std::string> actResListKK;
    std::vector<std::string> actResListPP;
    std::string reloadFitParaFile;
    std::string storeFitParaFile;
    std::string workPath;
    std::string outf_phipp, outf_phikk, proj_phipp, proj_phikk;
    std::string phsp_phipp, phsp_phikk, data_phipp, data_phikk;
    std::string idx_pp, idx_kk;
} PWA_CTRL;

void string_to_vector(std::string, std::vector<std::string> &);

static bool readConfigFile(const std::string &, const std::string &, std::string &);


void refresh_pwa_ctrl(const std::string, PWA_CTRL &);

#endif
