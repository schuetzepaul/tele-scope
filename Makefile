
# user: set path to eudaq in EUDAQ, e.g.
# export EUDAQ=/home/YOU/eudaq

ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)

ROOTLIBS = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)

# -g for gdb
# -pg for gprof
# -std=c++11

CXXFLAGS = -std=c++11 -O2 -Wall -Wextra $(ROOTCFLAGS) -I$(EUDAQ)/main/include

scope: scope.cc
	g++ $(CXXFLAGS) scope.cc -o scope \
	$(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: scope'

tele: tele.cc
	g++ $(CXXFLAGS) tele.cc -o tele \
	$(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: tele'

evd: evd.cc
	g++ $(CXXFLAGS) evd.cc -o evd \
	$(ROOTGLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: evd'
