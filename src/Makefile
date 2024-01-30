CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O3
OBJS = globals.o print.o variant.o dist.o bed.o cluster.o phase.o edit.o timer.o
TARGET = vcfdist
INSTALLDIR = /usr/local/bin

LDLIBS = -lz -lhts -lpthread
# g++ < 9.1 and llvm < 9.0 require linking to lstdc++fs
# if using g++ >= 9.1 or llvm >= 9.0, please comment the following line
LDLIBS += -lstdc++fs

all: $(TARGET)

install:
	cp $(TARGET) $(INSTALLDIR)

$(TARGET): $(OBJS) main.cpp
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) main.cpp $(LDLIBS)

globals.o: globals.cpp globals.h bed.h print.h defs.h timer.h
	$(CXX) -c $(CXXFLAGS) globals.cpp

print.o: print.cpp print.h globals.h phase.h dist.h edit.h defs.h variant.h
	$(CXX) -c $(CXXFLAGS) print.cpp

timer.o: timer.cpp timer.h globals.h defs.h
	$(CXX) -c $(CXXFLAGS) timer.cpp

variant.o: variant.cpp variant.h print.h fasta.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) variant.cpp

dist.o: dist.cpp dist.h fasta.h variant.h cluster.h print.h edit.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) dist.cpp

bed.o: bed.cpp bed.h print.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) bed.cpp

edit.o: edit.cpp edit.h defs.h globals.h
	$(CXX) -c $(CXXFLAGS) edit.cpp

cluster.o: cluster.cpp cluster.h variant.h fasta.h globals.h dist.h defs.h
	$(CXX) -c $(CXXFLAGS) cluster.cpp

phase.o: phase.cpp phase.h cluster.h print.h globals.h defs.h variant.h
	$(CXX) -c $(CXXFLAGS) phase.cpp

clean:
	rm -f $(TARGET) *.o
