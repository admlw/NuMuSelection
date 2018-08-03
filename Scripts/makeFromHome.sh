# For reasons I don't understand, there are issues using the makefile with 
# my local install of ROOT, but this seems to work

g++ -o MakeDataSimulationComparisonPlots MakeDataSimulationComparisonPlots.cxx SelectionMaker.cxx EventCategoriser.cxx EfficiencyPurity.cxx AnalysisCuts.cxx `root-config --cflags --libs`

