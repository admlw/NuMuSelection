// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPad.h"

// local
#include "EventCategoriser.h"
#include "DataTypes.h"
#include "AnalysisCuts.h"
#include "SelectionMaker.h"
#include "EfficiencyPurity.h"
#include "HistogramHandler.h"
#include "TreeHandler.h"
#include "Configuration.h"

// cpp
#include <iostream>
#include <bitset>
#include <string>

// initialise classes
numusel::Configuration    _config;
numusel::EventCategoriser _evcat;
numusel::SelectionMaker   _selmaker;
numusel::EfficiencyPurity _effpur;
numusel::HistogramHandler _histhandler;
numusel::TreeHandler      _treehandler;

std::vector<std::vector<hists_1d*>> plots_to_make;
std::vector<std::vector<eff_1d*>> eff_to_make;
