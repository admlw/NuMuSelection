#ifndef HISTOGRAMSTOPRODUCE_H
#define HISTOGRAMSTOPRODUCE_H

// ROOT
#include "TH1D.h"

std::vector<std::string> histoNames = {
  "nTracks",
  "nShowers"
};

std::vector<std::string> histoLabels = {
  ";Number of tracks;",
  ";Number of showers;"
};

std::vector<std::vector<double>> histoBins = {
  {10, 0, 10},
  {10, 0, 10}
};

#endif
