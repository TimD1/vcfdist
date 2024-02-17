TARGET = vcfdist
INSTALLDIR = /usr/local/bin
OBJS = globals.o print.o variant.o dist.o bed.o cluster.o phase.o edit.o timer.o

CXX ?= g++
CXXFLAGS ?= -O3
CXXFLAGS += -Wall -Wextra -std=c++17
LDFLAGS += -lz -lhts -lpthread -lm -lstdc++

all: $(TARGET)

install:
	cp $(TARGET) $(INSTALLDIR)

$(TARGET): $(OBJS) main.cpp
	$(CXX) $(OBJS) -o $(TARGET) main.cpp $(LDFLAGS)

globals.o: globals.cpp globals.h bed.h print.h defs.h timer.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) globals.cpp

print.o: print.cpp print.h globals.h phase.h dist.h edit.h defs.h variant.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) print.cpp

timer.o: timer.cpp timer.h globals.h defs.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) timer.cpp

variant.o: variant.cpp variant.h print.h fasta.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) variant.cpp

dist.o: dist.cpp dist.h fasta.h variant.h cluster.h print.h edit.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) dist.cpp

bed.o: bed.cpp bed.h print.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) bed.cpp

edit.o: edit.cpp edit.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) edit.cpp

cluster.o: cluster.cpp cluster.h variant.h fasta.h globals.h dist.h defs.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) cluster.cpp

phase.o: phase.cpp phase.h cluster.h print.h globals.h defs.h variant.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) phase.cpp

clean:
	rm -f $(TARGET) *.o
