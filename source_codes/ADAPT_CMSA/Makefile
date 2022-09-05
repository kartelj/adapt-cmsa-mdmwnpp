TARGET = mdmwnpp
CXXFLAGS = -ansi -O3
OBJS = Timer.o

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

CPLEXDIR      = /home/marko/Desktop/CPLEX_Studio127/cplex#adapt the path to "cplex" directory before use
CONCERTDIR    = /home/marko/Desktop/CPLEX_Studio127/concert#adapt the path to "concert" directory before use

CCC = g++ -std=c++17  #g++
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

all: ${TARGET}

mdmwnpp: mdmwnpp.o $(OBJS)
	$(CCC) $(CCFLAGS) mdmwnpp.o $(OBJS) -o mdmwnpp $(CCLNFLAGS)

mdmwnpp.o: mdmwnpp.cpp
	$(CCC) -c $(CCFLAGS) mdmwnpp.cpp -o mdmwnpp.o 

clean:
	@rm -f *~ *.o ${TARGET} core
