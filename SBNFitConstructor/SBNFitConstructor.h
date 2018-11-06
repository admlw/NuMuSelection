// cpp includes
#include <iostream>
#include <fstream>

// art includes
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"

// larsoft obj includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "uboone/EventWeight/MCEventWeight.h"

#include "TTree.h"
#include "TFile.h"

int sel_run;
int sel_subrun;
int sel_event;
double sel_resconstructed_neutrino_energy;
std::vector<bool>* eventCat = nullptr;

int ew_run;
int ew_subrun;
int ew_event;
double ew_nu_w;
double ew_nu_x;
double ew_nu_y;
double ew_nu_qsqr;

double out_reconstructed_neutrino_energy;

std::map<std::string, std::vector<double>> weights;

void initialiseTrees(TTree* sel, TTree* ew);

void setVariables(TTree* t);

