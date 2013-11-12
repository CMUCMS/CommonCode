### MODIFY THE FOLLOWING LINE FOR YOUR SETUP ###
SUSYNTUPLIZER=$(CMSSW_BASE)/src/SUSYPhotonAnalysis/SusyNtuplizer/src
################################################

TARGET = libCMUCommon.so
STANDALONESRCFILES = Utilities.cc ObjectSelector.cc ObjectVars.cc ObjectTree.cc
RA3SRCFILES = SimpleEventProducer.cc PFParticleBugFix.cc
STANDALONEOBJECTS = $(patsubst %.cc,%.o,$(STANDALONESRCFILES))
OBJECTS = $(STANDALONEOBJECTS) $(patsubst %.cc,%.o,$(RA3SRCFILES))

CFLAGS = -c -O3 -Wall -fPIC
LFLAGS = -shared

INC = -I. -I$(shell root-config --incdir)
SUSYINC = -I$(SUSYNTUPLIZER)
LIBS = $(shell root-config --libs)

all: INC += $(SUSYINC)
all: $(OBJECTS)
	g++ $(LFLAGS) -o $(TARGET) $(LIBS) $^

standalone: CFLAGS += -DSTANDALONE
standalone: $(STANDALONEOBJECTS)
	g++ $(LFLAGS) -o $(TARGET) $(LIBS) $^
	touch STANDALONE_BUILD

clean:
	rm -f $(TARGET) *.o STANDALONE_BUILD > /dev/null 2>&1

%.o: %.cc %.h
	g++ $(CFLAGS) $(INC) -o $@ $< $(LIBS)

.PHONY: all standalone clean
