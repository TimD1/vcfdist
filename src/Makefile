TARGET = vcfdist
INSTALLDIR = /usr/local/bin
OBJS = globals.o print.o variant.o dist.o bed.o cluster.o phase.o timer.o

CXX ?= g++
CXXFLAGS ?= -O3
CXXFLAGS += -Wall -Wextra -std=c++17
LDFLAGS += -lz -lhts -lpthread -lm -lstdc++

$(TARGET): $(OBJS) main.cpp
	$(CXX) $(OBJS) main.cpp $(LDFLAGS) -o $(TARGET)

install:
	cp $(TARGET) $(INSTALLDIR)

clean:
	rm -f $(TARGET) *.o

bed.o: bed.cpp bed.h print.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) bed.cpp

cluster.o: cluster.cpp cluster.h variant.h fasta.h globals.h dist.h defs.h
	$(CXX) -c $(CXXFLAGS) cluster.cpp

dist.o: dist.cpp dist.h fasta.h variant.h cluster.h print.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) dist.cpp

globals.o: globals.cpp globals.h bed.h print.h defs.h timer.h
	$(CXX) -c $(CXXFLAGS) globals.cpp

phase.o: phase.cpp phase.h cluster.h print.h globals.h defs.h variant.h
	$(CXX) -c $(CXXFLAGS) phase.cpp

print.o: print.cpp print.h globals.h phase.h dist.h defs.h variant.h
	$(CXX) -c $(CXXFLAGS) print.cpp

timer.o: timer.cpp timer.h globals.h defs.h
	$(CXX) -c $(CXXFLAGS) timer.cpp

variant.o: variant.cpp variant.h print.h fasta.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) variant.cpp
