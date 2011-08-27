ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)

INCLUDES       = -I./include -I../TESCO/include -I${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/interface

CXX            = g++
CXXFLAGS       = -g -fPIC -Wno-deprecated -D_GNU_SOURCE -O2 $(INCLUDES) 
LD             = g++
LDFLAGS        = -g 
#SOFLAGS        = -O --no_exceptions  -dynamiclib -flat_namespace -undefined suppress
SOFLAGS        = -O --no_exceptions -shared   -flat_namespace 


CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer -L${CMSSW_RELEASE_BASE}/lib/slc5_amd64_gcc434 -lFWCoreFWLite -lCondFormatsJetMETObjects  
GLIBS          = $(filter-out -lNew, $(NGLIBS))

SRCS           = ../TESCO/src/base/TreeClassBase.C ../TESCO/src/base/TreeReader.cc ../TESCO/src/base/TreeAnalyzerBase.cc ../TESCO/src/base/UserAnalysisBase.cc \
                 src/MT2tree.cc src/MT2Analysis.cc \
                 ../TESCO/src/helper/TMctLib.cc ../TESCO/src/helper/mctlib.cc  \
                 ../TESCO/src/helper/Davismt2.cc ../TESCO/src/helper/Hemisphere.cc  src/MT2Analyzer.cc \
                 src/MassPlotter.cc ../TESCO/src/helper/PUWeight.C \

# We want dictionaries only for classes that have _linkdef.h files                                                               
DICTOBS =  $(patsubst %_linkdef.hh, %.o, \
                      $(patsubst dict/%, obj/dict_%, \
                          $(wildcard dict/*_linkdef.hh) ) )

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

OBJS += $(DICTOBS)

SHARED=shlib/libDiLeptonAnalysis.so


.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend PhysQC

# Rules ====================================
all:   RunMT2Analyzer MakeMassPlots shared

 shared: $(SHARED)
 $(SHARED): $(OBJS)
	@echo "Creating library $(SHARED)"
	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJS) -o $(SHARED)
	@echo "$(SHARED) successfully compiled!"

RunMT2Analyzer: src/exe/RunMT2Analyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^	

MakeMassPlots: src/exe/MakeMassPlots.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

obj/dict_%.o: include/%.hh ../TESCO/include/%.hh  dict/%_linkdef.hh 
	@echo "Generating dictionary for $<"
	$(ROOTSYS)/bin/rootcint -f $(patsubst %.o, %.C, $@) -c -Idict $(INCLUDES) $(notdir $^)
	$(CXX) -c $(CXXFLAGS) -o $@ $(patsubst %.o, %.C, $@)

clean:
	$(RM) $(OBJS)	
	$(RM) RunMT2Analyzer
	$(RM) MakeMassPlots

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

