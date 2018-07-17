#ifndef HISTOGRAMSTOPRODUCE_H
#define HISTOGRAMSTOPRODUCE_H

// ROOT
#include "TH1D.h"

std::vector<std::string> histoNames = {
  "nTracks",
  "nShowers",
  "trackLength",
  "log(lmipoverp)"
};

std::vector<std::string> histoLabels = {
  ";Number of tracks;",
  ";Number of showers;",
  ";Track length (cm);",
  ";Log(L_{MIP}/L_{p});"
};

std::vector<std::vector<double>> histoBins = {
  {10, 0, 10},
  {10, 0, 10},
  {50, 0, 700},
  {50, -10, 10}
};

#endif
