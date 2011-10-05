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
		 src/QuickAnalysis.cc src/QuickAnalyzer.cc  src/MT2Shapes.cc

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
all:   RunMT2Analyzer RunQuickAnalyzer MakeMassPlots shared

 shared: $(SHARED)
 $(SHARED): $(OBJS)
	@echo "Creating library $(SHARED)"
	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJS) -o $(SHARED)
	@echo "$(SHARED) successfully compiled!"

RunMT2Analyzer: src/exe/RunMT2Analyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^	

RunQuickAnalyzer: src/exe/RunQuickAnalyzer.C $(OBJS)
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
	$(RM) RunQuickAnalyzer
	$(RM) MakeMassPlots

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

../TESCO/src/base/TreeClassBase.o: ../TESCO/include/base/TreeClassBase.h
../TESCO/src/base/TreeReader.o: ../TESCO/include/base/TreeReader.hh
../TESCO/src/base/TreeReader.o: ../TESCO/include/base/TreeClassBase.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/stdlib.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/features.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/sys/cdefs.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/wordsize.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/gnu/stubs.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/gnu/stubs-64.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/sys/types.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/types.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/typesizes.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/time.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/endian.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/endian.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/sys/select.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/select.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/sigset.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/time.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/sys/sysmacros.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/pthreadtypes.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/alloca.h
../TESCO/src/base/TreeAnalyzerBase.o: ../TESCO/include/base/TreeAnalyzerBase.hh
../TESCO/src/base/TreeAnalyzerBase.o: ../TESCO/include/base/TreeReader.hh
../TESCO/src/base/TreeAnalyzerBase.o: ../TESCO/include/base/TreeClassBase.h
../TESCO/src/base/TreeAnalyzerBase.o: ../TESCO/include/helper/Utilities.hh
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/stdio.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/libio.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/_G_config.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/wchar.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/wchar.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/gconv.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/stdio_lim.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/bits/sys_errlist.h
../TESCO/src/base/TreeAnalyzerBase.o: /usr/include/assert.h
../TESCO/src/base/UserAnalysisBase.o: ../TESCO/include/base/TreeReader.hh
../TESCO/src/base/UserAnalysisBase.o: ../TESCO/include/base/TreeClassBase.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/stdlib.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/features.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/sys/cdefs.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/wordsize.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/gnu/stubs.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/gnu/stubs-64.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/sys/types.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/types.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/typesizes.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/time.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/endian.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/endian.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/sys/select.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/select.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/sigset.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/time.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/sys/sysmacros.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/pthreadtypes.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/alloca.h
../TESCO/src/base/UserAnalysisBase.o: ../TESCO/include/helper/pdgparticle.hh
../TESCO/src/base/UserAnalysisBase.o: ../TESCO/include/helper/Monitor.hh
../TESCO/src/base/UserAnalysisBase.o: ../TESCO/include/base/UserAnalysisBase.hh
../TESCO/src/base/UserAnalysisBase.o: ../TESCO/include/base/TreeReader.hh
../TESCO/src/base/UserAnalysisBase.o: ../TESCO/include/helper/Utilities.hh
../TESCO/src/base/UserAnalysisBase.o: /usr/include/stdio.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/libio.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/_G_config.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/wchar.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/wchar.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/gconv.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/stdio_lim.h
../TESCO/src/base/UserAnalysisBase.o: /usr/include/bits/sys_errlist.h
../TESCO/src/base/UserAnalysisBase.o: ../TESCO/include/helper/PUWeight.h
src/MT2tree.o: ./include/MT2tree.hh ../TESCO/include/helper/Davismt2.h
src/MT2tree.o: /usr/include/math.h /usr/include/features.h
src/MT2tree.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
src/MT2tree.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
src/MT2tree.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/MT2tree.o: /usr/include/bits/mathcalls.h
src/MT2tree.o: ../TESCO/include/helper/TMctLib.h
src/MT2tree.o: ../TESCO/include/helper/mctlib.h
src/MT2tree.o: ../TESCO/include/helper/Hemisphere.hh
src/MT2tree.o: ../TESCO/include/helper/Utilities.hh /usr/include/stdio.h
src/MT2tree.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
src/MT2tree.o: /usr/include/libio.h /usr/include/_G_config.h
src/MT2tree.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/MT2tree.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
src/MT2tree.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/MT2tree.o: /usr/include/sys/types.h /usr/include/time.h
src/MT2tree.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MT2tree.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/MT2tree.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/MT2tree.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
src/MT2tree.o: /usr/include/alloca.h
src/MT2Analysis.o: ../TESCO/include/helper/Utilities.hh /usr/include/stdio.h
src/MT2Analysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/MT2Analysis.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/MT2Analysis.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/MT2Analysis.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/MT2Analysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/MT2Analysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/MT2Analysis.o: /usr/include/bits/stdio_lim.h
src/MT2Analysis.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/MT2Analysis.o: /usr/include/sys/types.h /usr/include/time.h
src/MT2Analysis.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MT2Analysis.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/MT2Analysis.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/MT2Analysis.o: /usr/include/sys/sysmacros.h
src/MT2Analysis.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/MT2Analysis.o: ./include/MT2Analysis.hh
src/MT2Analysis.o: ../TESCO/include/base/UserAnalysisBase.hh
src/MT2Analysis.o: ../TESCO/include/base/TreeReader.hh
src/MT2Analysis.o: ../TESCO/include/base/TreeClassBase.h
src/MT2Analysis.o: ../TESCO/include/helper/pdgparticle.hh
src/MT2Analysis.o: ../TESCO/include/helper/PUWeight.h
src/MT2Analysis.o: ../TESCO/include/base/TreeReader.hh
src/MT2Analysis.o: ../TESCO/include/helper/Davismt2.h /usr/include/math.h
src/MT2Analysis.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/MT2Analysis.o: /usr/include/bits/mathcalls.h
src/MT2Analysis.o: ../TESCO/include/helper/TMctLib.h
src/MT2Analysis.o: ../TESCO/include/helper/mctlib.h
src/MT2Analysis.o: ../TESCO/include/helper/Hemisphere.hh ./include/MT2tree.hh
src/MT2Analysis.o: /swshare/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_2_8/src/CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h
../TESCO/src/helper/TMctLib.o: ../TESCO/include/helper/TMctLib.h
../TESCO/src/helper/TMctLib.o: ../TESCO/include/helper/mctlib.h
../TESCO/src/helper/TMctLib.o: /usr/include/math.h /usr/include/features.h
../TESCO/src/helper/TMctLib.o: /usr/include/sys/cdefs.h
../TESCO/src/helper/TMctLib.o: /usr/include/bits/wordsize.h
../TESCO/src/helper/TMctLib.o: /usr/include/gnu/stubs.h
../TESCO/src/helper/TMctLib.o: /usr/include/gnu/stubs-64.h
../TESCO/src/helper/TMctLib.o: /usr/include/bits/huge_val.h
../TESCO/src/helper/TMctLib.o: /usr/include/bits/mathdef.h
../TESCO/src/helper/TMctLib.o: /usr/include/bits/mathcalls.h
../TESCO/src/helper/mctlib.o: ../TESCO/include/helper/mctlib.h
../TESCO/src/helper/mctlib.o: /usr/include/math.h /usr/include/features.h
../TESCO/src/helper/mctlib.o: /usr/include/sys/cdefs.h
../TESCO/src/helper/mctlib.o: /usr/include/bits/wordsize.h
../TESCO/src/helper/mctlib.o: /usr/include/gnu/stubs.h
../TESCO/src/helper/mctlib.o: /usr/include/gnu/stubs-64.h
../TESCO/src/helper/mctlib.o: /usr/include/bits/huge_val.h
../TESCO/src/helper/mctlib.o: /usr/include/bits/mathdef.h
../TESCO/src/helper/mctlib.o: /usr/include/bits/mathcalls.h
../TESCO/src/helper/Davismt2.o: ../TESCO/include/helper/Davismt2.h
../TESCO/src/helper/Davismt2.o: /usr/include/math.h /usr/include/features.h
../TESCO/src/helper/Davismt2.o: /usr/include/sys/cdefs.h
../TESCO/src/helper/Davismt2.o: /usr/include/bits/wordsize.h
../TESCO/src/helper/Davismt2.o: /usr/include/gnu/stubs.h
../TESCO/src/helper/Davismt2.o: /usr/include/gnu/stubs-64.h
../TESCO/src/helper/Davismt2.o: /usr/include/bits/huge_val.h
../TESCO/src/helper/Davismt2.o: /usr/include/bits/mathdef.h
../TESCO/src/helper/Davismt2.o: /usr/include/bits/mathcalls.h
../TESCO/src/helper/Hemisphere.o: ../TESCO/include/helper/Hemisphere.hh
../TESCO/src/helper/Hemisphere.o: ../TESCO/include/helper/Utilities.hh
../TESCO/src/helper/Hemisphere.o: /usr/include/stdio.h
../TESCO/src/helper/Hemisphere.o: /usr/include/features.h
../TESCO/src/helper/Hemisphere.o: /usr/include/sys/cdefs.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/wordsize.h
../TESCO/src/helper/Hemisphere.o: /usr/include/gnu/stubs.h
../TESCO/src/helper/Hemisphere.o: /usr/include/gnu/stubs-64.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/types.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/typesizes.h
../TESCO/src/helper/Hemisphere.o: /usr/include/libio.h
../TESCO/src/helper/Hemisphere.o: /usr/include/_G_config.h
../TESCO/src/helper/Hemisphere.o: /usr/include/wchar.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/wchar.h
../TESCO/src/helper/Hemisphere.o: /usr/include/gconv.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/stdio_lim.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/sys_errlist.h
../TESCO/src/helper/Hemisphere.o: /usr/include/stdlib.h
../TESCO/src/helper/Hemisphere.o: /usr/include/sys/types.h
../TESCO/src/helper/Hemisphere.o: /usr/include/time.h /usr/include/endian.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/endian.h
../TESCO/src/helper/Hemisphere.o: /usr/include/sys/select.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/select.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/sigset.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/time.h
../TESCO/src/helper/Hemisphere.o: /usr/include/sys/sysmacros.h
../TESCO/src/helper/Hemisphere.o: /usr/include/bits/pthreadtypes.h
../TESCO/src/helper/Hemisphere.o: /usr/include/alloca.h
src/MT2Analyzer.o: ./include/MT2Analyzer.hh
src/MT2Analyzer.o: ../TESCO/include/base/TreeAnalyzerBase.hh
src/MT2Analyzer.o: ../TESCO/include/base/TreeReader.hh
src/MT2Analyzer.o: ../TESCO/include/base/TreeClassBase.h
src/MT2Analyzer.o: ../TESCO/include/helper/Utilities.hh /usr/include/stdio.h
src/MT2Analyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/MT2Analyzer.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/MT2Analyzer.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/MT2Analyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/MT2Analyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/MT2Analyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/MT2Analyzer.o: /usr/include/bits/stdio_lim.h
src/MT2Analyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/MT2Analyzer.o: /usr/include/sys/types.h /usr/include/time.h
src/MT2Analyzer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MT2Analyzer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/MT2Analyzer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/MT2Analyzer.o: /usr/include/sys/sysmacros.h
src/MT2Analyzer.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/MT2Analyzer.o: ../TESCO/include/base/TreeReader.hh
src/MT2Analyzer.o: ./include/MT2Analysis.hh
src/MT2Analyzer.o: ../TESCO/include/base/UserAnalysisBase.hh
src/MT2Analyzer.o: ../TESCO/include/helper/pdgparticle.hh
src/MT2Analyzer.o: ../TESCO/include/helper/PUWeight.h
src/MT2Analyzer.o: ../TESCO/include/helper/Davismt2.h /usr/include/math.h
src/MT2Analyzer.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/MT2Analyzer.o: /usr/include/bits/mathcalls.h
src/MT2Analyzer.o: ../TESCO/include/helper/TMctLib.h
src/MT2Analyzer.o: ../TESCO/include/helper/mctlib.h
src/MT2Analyzer.o: ../TESCO/include/helper/Hemisphere.hh ./include/MT2tree.hh
src/MT2Analyzer.o: /swshare/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_2_8/src/CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h
src/MassPlotter.o: ./include/MassPlotter.hh ./include/MT2tree.hh
src/MassPlotter.o: ../TESCO/include/helper/Utilities.hh /usr/include/stdio.h
src/MassPlotter.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/MassPlotter.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/MassPlotter.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/MassPlotter.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/MassPlotter.o: /usr/include/_G_config.h /usr/include/wchar.h
src/MassPlotter.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/MassPlotter.o: /usr/include/bits/stdio_lim.h
src/MassPlotter.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/MassPlotter.o: /usr/include/sys/types.h /usr/include/time.h
src/MassPlotter.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MassPlotter.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/MassPlotter.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/MassPlotter.o: /usr/include/sys/sysmacros.h
src/MassPlotter.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/MassPlotter.o: ../TESCO/include/helper/Monitor.hh /usr/include/math.h
src/MassPlotter.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/MassPlotter.o: /usr/include/bits/mathcalls.h
../TESCO/src/helper/PUWeight.o: ../TESCO/include/helper/PUWeight.h
src/QuickAnalysis.o: ../TESCO/include/helper/Utilities.hh
src/QuickAnalysis.o: /usr/include/stdio.h /usr/include/features.h
src/QuickAnalysis.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
src/QuickAnalysis.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
src/QuickAnalysis.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
src/QuickAnalysis.o: /usr/include/libio.h /usr/include/_G_config.h
src/QuickAnalysis.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/QuickAnalysis.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
src/QuickAnalysis.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/QuickAnalysis.o: /usr/include/sys/types.h /usr/include/time.h
src/QuickAnalysis.o: /usr/include/endian.h /usr/include/bits/endian.h
src/QuickAnalysis.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/QuickAnalysis.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/QuickAnalysis.o: /usr/include/sys/sysmacros.h
src/QuickAnalysis.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/QuickAnalysis.o: ./include/QuickAnalysis.hh
src/QuickAnalysis.o: ../TESCO/include/base/TreeReader.hh
src/QuickAnalysis.o: ../TESCO/include/base/TreeClassBase.h
src/QuickAnalysis.o: ../TESCO/include/base/UserAnalysisBase.hh
src/QuickAnalysis.o: ../TESCO/include/base/TreeReader.hh
src/QuickAnalysis.o: ../TESCO/include/helper/pdgparticle.hh
src/QuickAnalysis.o: ../TESCO/include/helper/PUWeight.h
src/QuickAnalyzer.o: ./include/QuickAnalyzer.hh
src/QuickAnalyzer.o: ../TESCO/include/base/TreeAnalyzerBase.hh
src/QuickAnalyzer.o: ../TESCO/include/base/TreeReader.hh
src/QuickAnalyzer.o: ../TESCO/include/base/TreeClassBase.h
src/QuickAnalyzer.o: ../TESCO/include/helper/Utilities.hh
src/QuickAnalyzer.o: /usr/include/stdio.h /usr/include/features.h
src/QuickAnalyzer.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
src/QuickAnalyzer.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
src/QuickAnalyzer.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
src/QuickAnalyzer.o: /usr/include/libio.h /usr/include/_G_config.h
src/QuickAnalyzer.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/QuickAnalyzer.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
src/QuickAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/QuickAnalyzer.o: /usr/include/sys/types.h /usr/include/time.h
src/QuickAnalyzer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/QuickAnalyzer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/QuickAnalyzer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/QuickAnalyzer.o: /usr/include/sys/sysmacros.h
src/QuickAnalyzer.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/QuickAnalyzer.o: ../TESCO/include/base/TreeReader.hh
src/QuickAnalyzer.o: ./include/QuickAnalysis.hh
src/QuickAnalyzer.o: ../TESCO/include/base/UserAnalysisBase.hh
src/QuickAnalyzer.o: ../TESCO/include/helper/pdgparticle.hh
src/QuickAnalyzer.o: ../TESCO/include/helper/PUWeight.h
src/MT2Shapes.o: ./include/MT2Shapes.hh ./include/MT2tree.hh
src/MT2Shapes.o: ../TESCO/include/helper/Utilities.hh /usr/include/stdio.h
src/MT2Shapes.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/MT2Shapes.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/MT2Shapes.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/MT2Shapes.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/MT2Shapes.o: /usr/include/_G_config.h /usr/include/wchar.h
src/MT2Shapes.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/MT2Shapes.o: /usr/include/bits/stdio_lim.h
src/MT2Shapes.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/MT2Shapes.o: /usr/include/sys/types.h /usr/include/time.h
src/MT2Shapes.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MT2Shapes.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/MT2Shapes.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/MT2Shapes.o: /usr/include/sys/sysmacros.h
src/MT2Shapes.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/MT2Shapes.o: ../TESCO/include/helper/Monitor.hh /usr/include/math.h
src/MT2Shapes.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/MT2Shapes.o: /usr/include/bits/mathcalls.h
