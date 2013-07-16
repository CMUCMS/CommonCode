### MODIFY THE FOLLOWING LINE FOR YOUR SETUP ###
SUSYNTUPLIZER=$(CMSSW_BASE)/src/SusyAnalysis/SusyNtuplizer/src
################################################

TARGET = libCMUCommon.so
SRCFILES = Utilities.cc ObjectSelector.cc ObjectVars.cc ObjectTree.cc SimpleEventProducer.cc
HEADERS = $(patsubst %.cc,%.h,$(SRCFILES))
OBJECTS = $(patsubst %.cc,%.o,$(SRCFILES))

CFLAGS = -c -O3 -Wall -fPIC
LFLAGS = -shared -Wl

INC = -I. -I$(shell root-config --incdir) -I$(SUSYNTUPLIZER)
LIBS = $(shell root-config --libs)

all: $(TARGET)

clean:
	rm -f $(TARGET) $(OLDTARGET) *.o > /dev/null 2>&1

$(TARGET): $(OBJECTS)
	g++ $(LFLAGS) -o $@ $(LIBS) $^

$(OLDTARGET): $(OLDOBJECTS) Dict.o
	g++ $(LFLAGS) -o $@ $(LIBS) $^

%.o: %.cc %.h
	g++ $(CFLAGS) $(INC) -o $@ $< $(LIBS)
