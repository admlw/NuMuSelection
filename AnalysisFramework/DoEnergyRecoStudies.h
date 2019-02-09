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
#include "TMath.h"

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

int nbins_p[6] = {50, 100, 150,150,200,150}; //250
double proton_binVals_arr[7] = {0.0, 0.05, 0.10, 0.20, 0.3, 0.4, 0.6};
std::vector<std::vector<double>> proton_binVals {
    {proton_binVals_arr[0], proton_binVals_arr[1]},
    {proton_binVals_arr[1], proton_binVals_arr[2]},
    {proton_binVals_arr[2], proton_binVals_arr[3]},
    {proton_binVals_arr[3], proton_binVals_arr[4]},
    {proton_binVals_arr[4], proton_binVals_arr[5]},
    {proton_binVals_arr[5], proton_binVals_arr[6]}
};

std::vector<std::vector<double>> proton_fitRange {
    {-0.4, 0.5},
    {-0.2, 0.15},    
    {-0.2, 0.1},
    {-0.1, 0.1},
    {-0.1, 0.03}, // -0.05, 0.025
    {-0.15, 0.05},
};

int nbins_mu_c = 120;
int nbins_mu_u[6] = {40, 80, 80, 150, 80, 80};
double muon_binVals_arr[8] = {0.0, 0.2, 0.3, 0.4, 0.6, 0.80, 1.2};
std::vector<std::vector<double>> muon_binVals {
    {muon_binVals_arr[0], muon_binVals_arr[1]},
        {muon_binVals_arr[1], muon_binVals_arr[2]},
        {muon_binVals_arr[2], muon_binVals_arr[3]},
        {muon_binVals_arr[3], muon_binVals_arr[4]},
        {muon_binVals_arr[4], muon_binVals_arr[5]},
        {muon_binVals_arr[5], muon_binVals_arr[6]}
};

Double_t CrystalBall(Double_t *x, Double_t *par){

    Double_t arg = 0;
    Double_t xv  = x[0];

    if ((xv-par[1])/par[2] >= -par[3]){
        arg = -0.5*(pow(((xv-par[1])/par[2]),2));
    }
    else{
        arg = (pow(par[3],2)/2.)+par[3]*((xv-par[1])/par[2]);
    }

    return par[0]*exp(arg);

}

Double_t LanGaus(Double_t *x, Double_t *par) {

    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation),
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

    // Numeric constants
    Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    Double_t mpshift  = -0.22278298;       // Landau maximum location

    // Control constants
    Double_t np = 100.0;      // number of convolution steps
    Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

    // Variables
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;

    // MP shift correction
    mpc = par[1] - mpshift * par[0];

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];

    step = (xupp-xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);

        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
};


Double_t ResolutionFunction(Double_t* x, Double_t *par){
 
    Double_t xv  = x[0];
    
    return par[0] + par[1]/sqrt(x[0]) + par[2]/x[0];

}
