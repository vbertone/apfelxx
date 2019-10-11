
CXX = clang++

CXXFLAGS += -O3 -fPIC -std=c++11

# APFEL++
APFELPPINCS = $(shell apfelxx-config --cppflags)
APFELPPLIBS = $(shell apfelxx-config --ldflags)

# LHAPDF
LHAPDFINCS = $(shell lhapdf-config --cppflags)
LHAPDFLIBS = $(shell lhapdf-config --ldflags)

# NangaParbat
NANPARINCS = $(shell NangaParbat-config --cppflags)
NANPARLIBS = $(shell NangaParbat-config --ldflags)

# Now set up the compiler and link flags and libs
CXXFLAGS += $(LHAPDFINCS) $(APFELPPINCS) $(NANPARINCS)
LDFLAGS  += $(LHAPDFINCS) $(APFELPPINCS) $(NANPARINCS)


CLIBS += $(LHAPDFLIBS) $(APFELPPLIBS) $(NANPARLIBS)

install : all
all : logs FOvsAsy EvolutionFactor fNP LeptonicModulation

logs: logs.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

FOvsAsy: FOvsAsy.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

EvolutionFactor: EvolutionFactor.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

fNP: fNP.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

LeptonicModulation: LeptonicModulation.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cxx .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la logs FOvsAsy EvolutionFactor fNP *~

