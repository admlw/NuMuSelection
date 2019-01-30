// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraphErrors.h"

// local
#include "Configuration.h"
#include "EventCategoriser.h"
#include "DataTypes.h"
#include "AnalysisCuts.h"
#include "SelectionMaker.h"
#include "Efficiency.h"
#include "HistogramHandler.h"
#include "TreeHandler.h"
#include "TPaveText.h"

// cpp
#include <iostream>
#include <bitset>
#include <string>

// initialise classes
numusel::Configuration    _config;
numusel::TreeHandler      _treeHandler;

double proton_binVals_arr[8] = {0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 1.0};
double muon_binVals_arr[8] = {0.0, 0.2, 0.3, 0.4, 0.6, 0.80, 1.2};

std::vector<std::vector<double>> proton_binVals {
    {proton_binVals_arr[0], proton_binVals_arr[1]},
    {proton_binVals_arr[1], proton_binVals_arr[2]},
    {proton_binVals_arr[2], proton_binVals_arr[3]},
    {proton_binVals_arr[3], proton_binVals_arr[4]},
    {proton_binVals_arr[4], proton_binVals_arr[5]},
    {proton_binVals_arr[5], proton_binVals_arr[6]},
    {proton_binVals_arr[6], proton_binVals_arr[7]}
};

std::vector<std::vector<double>> muon_binVals {
    {muon_binVals_arr[0], muon_binVals_arr[1]},
    {muon_binVals_arr[1], muon_binVals_arr[2]},
    {muon_binVals_arr[2], muon_binVals_arr[3]},
    {muon_binVals_arr[3], muon_binVals_arr[4]},
    {muon_binVals_arr[4], muon_binVals_arr[5]},
    {muon_binVals_arr[5], muon_binVals_arr[6]}
};

