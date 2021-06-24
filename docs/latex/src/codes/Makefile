
CXX = clang++

CXXFLAGS += -O3 -fPIC -std=c++11

# APFEL++
APFELPPINCS = $(shell apfelxx-config --cppflags)
APFELPPLIBS = $(shell apfelxx-config --ldflags)

# APFEL
APFELINCS = $(shell apfel-config --cppflags)
APFELLIBS = $(shell apfel-config --ldflags)

# LHAPDF
LHAPDFINCS = $(shell lhapdf-config --cppflags)
LHAPDFLIBS = $(shell lhapdf-config --ldflags)

# NangaParbat
#NANPARINCS = $(shell NangaParbat-config --cppflags)
#NANPARLIBS = $(shell NangaParbat-config --ldflags)

# Now set up the compiler and link flags and libs
CXXFLAGS += $(LHAPDFINCS) $(APFELPPINCS) $(NANPARINCS) $(APFELINCS)
LDFLAGS  += $(LHAPDFINCS) $(APFELPPINCS) $(NANPARINCS) $(APFELINCS)

CLIBS += $(LHAPDFLIBS) $(APFELPPLIBS) $(NANPARLIBS) $(APFELLIBS)

install : all
all : logs FOvsAsy FONLL EvolutionFactor LeptonicModulation PhaseSpaceRedFactor PhaseSpacePowerCorrs NewPrescription Alphas PDFs InterpolationFunctions Yterm #fNP

logs: logs.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

FOvsAsy: FOvsAsy.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

FONLL: FONLL.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

EvolutionFactor: EvolutionFactor.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

#fNP: fNP.o
#	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

LeptonicModulation: LeptonicModulation.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

PhaseSpaceRedFactor: PhaseSpaceRedFactor.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

PhaseSpacePowerCorrs: PhaseSpacePowerCorrs.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

NewPrescription: NewPrescription.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

Alphas: Alphas.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

PDFs: PDFs.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

InterpolationFunctions: InterpolationFunctions.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

Yterm: Yterm.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cxx .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la *~ logs FOvsAsy FONLL EvolutionFactor LeptonicModulation PhaseSpaceRedFactor PhaseSpacePowerCorrs NewPrescription Alphas PDFs InterpolationFunctions Yterm #fNP

