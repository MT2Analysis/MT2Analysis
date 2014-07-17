CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src -I$(CMSSW_BASE)/src/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/CondFormats/JetMETObjects/interface -I$(CMSSW_BASE)/src/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/ -I$(CMSSW_BASE)/src/MT2Analysis/Code/MT2AnalysisCode/RootMacros/include/ -I$(CMSSW_BASE)/src/MT2Analysis/Code/MT2AnalysisCode/TESCO/include/

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -L$(CMSSW_BASE)/lib/slc5_amd64_gcc462/ -L/swshare/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/slc5_amd64_gcc462/ -L/swshare/cms/slc5_amd64_gcc462/cms/cmssw-patch/CMSSW_5_3_7_patch5/lib/slc5_amd64_gcc462 -L$(CMSSW_BASE)/src/MT2Analysis/Code/MT2AnalysisCode/MT2Code/shlib/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit -lFWCoreFWLite -lFWCoreUtilities -lDataFormatsCommon -lDataFormatsFWLite -lCondFormatsJetMETObjects -lDiLeptonAnalysis -lPhysics -L/swshare/cms/slc5_amd64_gcc462/external/zlib/1.2.3/lib -lz -lTreePlayer




all: computeLostLepton compareToHJ


prova: prova.cpp 
	$(CC) -Wall $(INCLUDES) -o prova prova.cpp $(CMSSW_BASE)/src/MT2Analysis/Code/MT2AnalysisCode/MT2Code/src/MT2tree.o $(ROOTFLAG) $(EXTRALIBS)

computeLostLepton: computeLostLepton.cpp MT2Region.o MT2LostLeptonUtilities.o MT2Common.o 
	$(CC) -Wall $(INCLUDES) -o computeLostLepton computeLostLepton.cpp MT2Region.o MT2LostLeptonUtilities.o MT2Common.o $(CMSSW_BASE)/src/MT2Analysis/Code/MT2AnalysisCode/MT2Code/src/MT2tree.o $(ROOTFLAG) $(EXTRALIBS)

compareToHJ: compareToHJ.cpp DrawBase.o fitTools.o
	$(CC) -Wall $(INCLUDES) -o compareToHJ compareToHJ.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)




MT2Common.o: src/MT2Common.cc
	$(CC) -Wall $(INCLUDES) -c src/MT2Common.cc $(ROOTFLAG) $(EXTRALIBS)


MT2Region.o: src/MT2Region.cc
	$(CC) -Wall $(INCLUDES) -c src/MT2Region.cc $(ROOTFLAG) $(EXTRALIBS)

MT2LostLeptonUtilities.o: src/MT2LostLeptonUtilities.cc MT2Region.o
	$(CC) -Wall $(INCLUDES) -c src/MT2LostLeptonUtilities.cc MT2Region.o $(ROOTFLAG) $(EXTRALIBS)



DrawBase.o: $(CMSSW_BASE)/src/CommonTools/DrawBase.C
	$(CC) -Wall $(INCLUDES) -c $(CMSSW_BASE)/src/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS)

fitTools.o: $(CMSSW_BASE)/src/CommonTools/fitTools.C
	$(CC) -Wall $(INCLUDES) -c $(CMSSW_BASE)/src/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)



clean:
	rm *.o
