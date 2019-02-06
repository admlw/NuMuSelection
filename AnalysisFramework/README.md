# Analysis README.

## Introduction

This framework is designed to run on files output from the NuMuSelection LArSoft module.

This module outputs two root files: one which contains selection information and one which contains eventweight information.

Running options are contained in `Configuration.h`, while cut values are contained within `AnalysisCuts.h`.

All scripts can be built with `make all`. This produces four executables:
- `DoProducePlots` Runs the selection and makes a large collection of plots depending on options defined in `Configuration.h`.
- `DoEnergyRecoStudies` takes in the output of `DoProducePlots` and finds the resolution and bias of the protons, contained muons, and uncontained muons separately (once you've done this, you should put the bias gradiant and intercept into the `Configuration.h` and rerun `DoProducePlots`.).
- `SortTree` MUST be run on the eventweight trees in order for the EventWeight matching to work.
- `SplitTree` can be sued to split the output files of `DoProducePlots` in order to run on the grid
