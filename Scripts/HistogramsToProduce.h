#ifndef HISTOGRAMSTOPRODUCE_H
#define HISTOGRAMSTOPRODUCE_H

// ROOT
#include "TH1D.h"

/**
 * These three vectors fully define all of the data/mc histograms to produce
 */
std::vector<std::string> histoNames = {
  "nTracks",
  "nShowers",
  "trackLength",
  "vertex_x",
  "vertex_y",
  "vertex_z",
  "track_startx",
  "track_endx",
  "track_starty",
  "track_endy",
  "track_startz",
  "track_endz",
  "track_theta",
  "track_costheta",
  "track_phi",
  "log(lmipoverp)",
  "track_mcs_fwd",
  "track_mcs_bwd",
  "muoncand_track_length",
  "muoncand_track_theta",
  "muoncand_track_costheta",
  "muoncand_track_phi",
  "muoncand_track_mcs_fwd",
  "muoncand_track_mcs_bwd"
};

std::vector<std::string> histoLabels = {
  ";Number of tracks;",
  ";Number of showers;",
  ";Track length (cm);",
  ";Vertex x position (cm);",
  ";Vertex y position (cm);",
  ";Vertex z position (cm);",
  ";Track start x (cm);",
  ";Track end x (cm);",
  ";Track start y (cm);",
  ";Track end y (cm);",
  ";Track start z (cm);",
  ";Track end z (cm);",
  ";Track theta (rad);",
  ";Track cos(theta);",
  ";Track phi;",
  ";Log(L_{MIP}/L_{p});",
  ";Track MCS (forward);",
  ";Track MCS (backward);",
  ";Muon candidate track length (cm);",
  ";Muon candidate track theta (rad);",
  ";Muon candidate track cos(theta);",
  ";Muon candidate track phi;",
  ";Muon candidate track MCS (forward);",
  ";Muon candidate track MCS (backward;",
};

std::vector<std::vector<double>> histoBins = {
  {10, 0, 10},
  {10, 0, 10},
  {50, 0, 700},
  {50, 0, 256},
  {50, -116.5, 116.5},
  {50, 0, 1040},
  {50, 0, 256},
  {50, 0, 256},
  {50, -116.5, 116.5},
  {50, -116.5, 116.5},
  {50, 0, 1040},
  {50, 0, 1040},
  {25, 0, 3.15},
  {25, -1, 1},
  {50, -3.15, 3.15},
  {50, -10, 10},
  {50, 0, 3},
  {50, 0, 3},
  {25, 0, 700},
  {25, 0, 3.15},
  {25, -1, 1},
  {25, -3.15, 3.15},
  {25, 0, 3},
  {25, 0, 3}
};

/*
 * These three vectors fully define all of the efficiency/purity histograms to produce
 */
std::vector<std::string> effpurNames = {
  "true_enu",
  "true_mu_p"
};

std::vector<std::string> effpurLabels = {
  ";E_{#nu}^{true} (Gev);",
  ";p_{#mu} (Gev);"
};

std::vector< std::vector<double> > effpurBins = {
  {25, 0, 3},
  {25, 0, 2.5}
};

#endif
