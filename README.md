# NuMu CC0PiNP Selection

This repository contains several pieces of code to aid in a CC0PiNP analysis, although in principle it could be extended to any analysis which inherits from the CC-inclusive selection.

The `Modules`, `Algorithms`, and `job` directories hold a LArSoft module, helper algorithms and fhicl files respectively. The module is designed to be bare-bones and essentially just pulls out important information from the art event and dumps it into a normal ROOT tree. 

The `AnalysisFramework` directory contains a C++ package for performing the analysis. As of right now there are 4 stages to the analysis. The CC-inclusive selected TPC objects are first required to pass a set of topology cuts (number of PFParticles, number of Tracks, number of Showers), and then there is a cut on the particle identification demanding that there be only one muon candidate in the TPC object. 
