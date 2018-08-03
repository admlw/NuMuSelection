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
#include "Configuration.h"
#include "EventCategoriser.h"
#include "DataTypes.h"
#include "AnalysisCuts.h"
#include "SelectionMaker.h"
#include "EfficiencyPurity.h"
#include "HistogramHandler.h"
#include "TreeHandler.h"

// cpp
#include <iostream>
#include <bitset>
#include <string>

std::vector<std::vector<hists_1d*>> plots_to_make;
std::vector<std::vector<eff_1d*>> eff_to_make;
std::vector<std::vector<pur_1d*>> pur_to_make;
