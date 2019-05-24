SOURCES=HEPTopTagger.cc example.cc
MYSOURCES=HEPTopTagger.cc EventAnalysis.cc   
MYSOURCES2=HEPTopTagger.cc JetAnalysis.cc 
OBJECTS=$(SOURCES:.cc=.o)
CC = g++ 

######### FASTJET ##################
FASTJETLIB  = $(shell fastjet-config --libs)
CXXFLAGS += $(shell root-config --cflags) 
CXXFLAGS += $(shell fastjet-config --cxxflags) 
CXXLIBS = $(shell root-config --glibs)
CXXLIBS += ${FASTJETLIB}

######### NSUBJETTINESS ############  
CXXLIBS += -L./Nsubjettiness  -lNsubjettiness
######### QJETS #################### 
CXXLIBS += -L./qjets/lib -lQjets

######### RULES ####################
all:
	cd qjets; make
	cd Nsubjettiness; make
	make example

example: $(OBJECTS)
	$(CC) $(CXXFLAGS) $(OBJECTS) $(CXXLIBS)  -o $@

EventAnalysis: $(MYSOURCES)
	$(CC) $(CXXFLAGS) $(MYSOURCES) $(CXXLIBS)  -o $@

JetAnalysis: $(MYSOURCES2)
	$(CC) $(CXXFLAGS) $(MYSOURCES2) $(CXXLIBS)  -o $@

.cc.o:
	$(CC) $(CXXFLAGS) $(INCLUDES) -c $<

clean:
	rm -f $(OBJECTS) example
	cd qjets; make clean
	cd Nsubjettiness; make clean
