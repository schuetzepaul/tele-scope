
# user:
# set path to eudaq and GBL:
# export EUDAQ=/home/YOURID/eudaq
# export GBL=/home/YOURID/GBL/V01-17-00/cpp

ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)

ROOTLIBS = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)

# -g for gdb
# -pg for gprof
# -std=c++11

CXXFLAGS = -std=c++11 -O2 -Wall -Wextra $(ROOTCFLAGS) -I$(EUDAQ)/main/include -I$(LCIO)/include

scope: scope.cc
	g++ $(CXXFLAGS) -o scope scope.cc \
	$(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: scope'

tele: tele.cc
	g++ $(CXXFLAGS) -o tele tele.cc \
	$(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: tele'

quad: quad.cc
	g++ $(CXXFLAGS) -I$(GBL)/include -o quad quad.cc \
	-L$(GBL)/lib -lGBL $(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ -L$(LCIO)/lib -llcio
	@echo 'done: quad'

evd: evd.cc
	g++ $(CXXFLAGS) -o evd evd.cc \
	$(ROOTGLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: evd'
