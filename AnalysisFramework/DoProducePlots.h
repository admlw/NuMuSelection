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
#include "TError.h"

// local
#include "Configuration.h"
#include "EventCategoriser.h"
#include "DataTypes.h"
#include "AnalysisCuts.h"
#include "SelectionMaker.h"
#include "Efficiency.h"
#include "HistogramHandler.h"
#include "TreeHandler.h"

// cpp
#include <iostream>
#include <bitset>
#include <string>

// initialise classes
numusel::Configuration    _config;
numusel::EventCategoriser _evcat;
numusel::SelectionMaker   _selmaker;
numusel::Efficiency       _eff;
numusel::HistogramHandler _histoHandler;
numusel::TreeHandler      _treeHandler;
numusel::AnalysisCuts     _anacuts;

std::vector<std::vector<hists_1d*>> plots_to_make;
std::vector<std::vector<trackhists_1d*>> trackplots_to_make;
std::vector<std::vector<eff_1d*>> eff_to_make;
std::vector<std::vector<hists_2d*>> plots_to_make_2D;

