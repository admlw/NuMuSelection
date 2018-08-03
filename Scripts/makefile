CPPFLAGS =-I $(ROOT_INC)

CXXFLAGS=-std=c++14
CXX=g++

LDFLAGS=$$(root-config --libs)

ProducePlots: Main.cxx EventCategoriser.cxx SelectionMaker.cxx AnalysisCuts.cxx EfficiencyPurity.cxx HistogramHandler.cxx TreeHandler.cxx
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
	
all: ProducePlots

clean:
	rm *.o ProducePlots
