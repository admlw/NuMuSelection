CPPFLAGS =-I $(ROOT_INC)

CXXFLAGS=-std=c++14
CXX=g++

LDFLAGS=$$(root-config --libs)

DoProducePlots: DoProducePlots.cxx EventCategoriser.cxx SelectionMaker.cxx AnalysisCuts.cxx Efficiency.cxx HistogramHandler.cxx TreeHandler.cxx
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

DoEnergyRecoStudies: DoEnergyRecoStudies.cxx EnergyReco.cxx TreeHandler.cxx
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

SplitTree: SplitTree.cxx TreeHandler.cxx
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

all: DoProducePlots DoEnergyRecoStudies SplitTree

clean:
	rm DoProducePlots
	rm DoEnergyRecoStudies
	rm SplitTree
	rm plots/*.png
	rm plots/*.pdf
