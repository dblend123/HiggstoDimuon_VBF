CC=g++

IDIR=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.14.09-pafccj/
LDIR=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_0/external/slc7_amd64_gcc700/lib

#IDIR=/disk1/root/include/
#LDIR=/disk1/root/lib
CFLAGS=-I$(IDIR)


R_LDFLAGS=`root-config --ldflags`
R_LIBS=`root-config --libs`
R_CFLAGS=`root-config --cflags`
R_ALL = $(R_LADFLAGS) $(R_LIBS) $(R_CFLAGS)

#LIBS=-L$(LDIR)  -lCore  -lNew   \
#	-lGraf -lGraf3d -lHist -lHtml  -lMatrix -lMinuit -lNet \
#	-lPostscript -lTree -lPhysics -lRIO -lMathCore\
#	-lGpad -lGui -lGX11 -lX3d -lRint -lThread\
#	-L/usr/lib/   -lm -lncurses


myanalysis:
	$(CC)  $(CFLAGS)  MyAnaSkim.C MyMuMu.o -o MySkimming $(R_ALL) $(LIBS)
