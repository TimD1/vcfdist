CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O2
OBJS = globals.o print.o variant.o dist.o bed.o cluster.o phase.o edit.o
TARGET = vcfdist
LDLIBS = -lz -lhts -lstdc++fs

all: $(TARGET)

$(TARGET): $(OBJS) main.cpp
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) main.cpp $(LDLIBS)

globals.o: globals.cpp globals.h bed.h print.h defs.h
	$(CXX) -c $(CXXFLAGS) globals.cpp

print.o: print.cpp print.h globals.h phase.h dist.h edit.h defs.h
	$(CXX) -c $(CXXFLAGS) print.cpp

variant.o: variant.cpp variant.h print.h fasta.h defs.h
	$(CXX) -c $(CXXFLAGS) variant.cpp

dist.o: dist.cpp dist.h fasta.h variant.h cluster.h print.h edit.h defs.h
	$(CXX) -c $(CXXFLAGS) dist.cpp

bed.o: bed.cpp bed.h print.h defs.h
	$(CXX) -c $(CXXFLAGS) bed.cpp

edit.o: edit.cpp edit.h defs.h
	$(CXX) -c $(CXXFLAGS) edit.cpp

cluster.o: cluster.cpp cluster.h variant.h globals.h dist.h defs.h
	$(CXX) -c $(CXXFLAGS) cluster.cpp

phase.o: phase.cpp phase.h cluster.h print.h globals.h defs.h
	$(CXX) -c $(CXXFLAGS) phase.cpp

clean:
	rm -f $(TARGET) *.o
