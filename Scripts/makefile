CPPFLAGS =-I $(ROOT_INC)

CXXFLAGS=-std=c++14
CXX=g++

LDFLAGS=$$(root-config --libs)

MakeDataSimulationComparisonPlots: MakeDataSimulationComparisonPlots.cxx
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
	
all: MakeDataSimulationComparisonPlots 

clean:
	rm *.o MakeDataSimulationComparisonPlots
