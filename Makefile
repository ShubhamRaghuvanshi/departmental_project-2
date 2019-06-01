SOURCES=HEPTopTagger.cc example.cc
MYSOURCES=HEPTopTagger.cc EventAnalysis.cc   
MYSOURCES2=HEPTopTagger.cc Analyser.cc particleproperty.cc Func.cc zpttEventAnalysis.cc 
MYSOURCES3=HEPTopTagger.cc Analyser.cc particleproperty.cc Func.cc Analyser.h particleproperty.h Func.h /home/ehep/Downloads/products/fjcontrib-1.041/RecursiveTools/SoftDrop.cc /home/ehep/Downloads/products/fjcontrib-1.041/RecursiveTools/RecursiveSymmetryCutBase.cc pp2zp2tt2qqblv_EventAnalysis.cc 


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

zpttEventAnalysis: $(MYSOURCES2)
	$(CC) $(CXXFLAGS) $(MYSOURCES2) $(CXXLIBS)  -o $@

pp2zp2tt2qqblv_EventAnalysis: $(MYSOURCES3)
	$(CC) $(CXXFLAGS) $(MYSOURCES3) $(CXXLIBS)  -o $@

.cc.o:
	$(CC) $(CXXFLAGS) $(INCLUDES) -c $<

clean:
	rm -f $(OBJECTS) example
	cd qjets; make clean
	cd Nsubjettiness; make clean
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
