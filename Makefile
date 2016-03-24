
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
# root C flags =  -pthread -m64 -I/home/pitzl/ROOT/armin/root-cern/include

ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)

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
