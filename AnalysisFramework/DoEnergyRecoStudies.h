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

std::vector<std::vector<double>> proton_binVals {
    {0.00, 0.05},
    {0.05, 0.10},
    {0.10, 0.15},
    {0.15, 0.20},
    {0.20, 0.30},
    {0.30, 0.40},
    {0.40, 1.00}
};

std::vector<std::vector<double>> muon_binVals {
    {0.00, 0.20},
    {0.20, 0.30},
    {0.30, 0.40},
    {0.40, 0.60},
    {0.60, 0.80},
    {0.80, 1.20}
};

