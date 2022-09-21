CC = /opt/rh/devtoolset-8/root/bin/g++
CFLAGS = -g -Wall -std=c++2a -O2
IFLAGS = /usr/local/include/
OBJS = globals.o print.o variant.o dist.o bed.o cluster.o phase.o
TARGET = vcfdist
LDLIBS = -lz -lhts -lstdc++fs

all: $(TARGET)

$(TARGET): $(OBJS) main.cpp
	$(CC) $(CFLAGS) $(OBJS) -I$(IFLAGS) -o $(TARGET) main.cpp $(LDLIBS)

globals.o: globals.cpp globals.h bed.h print.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) globals.cpp

print.o: print.cpp print.h globals.h phase.h dist.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) print.cpp

variant.o: variant.cpp variant.h print.h fasta.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) variant.cpp

dist.o: dist.cpp dist.h fasta.h variant.h cluster.h print.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) dist.cpp

bed.o: bed.cpp bed.h print.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) bed.cpp

cluster.o: cluster.cpp cluster.h variant.h globals.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) cluster.cpp

phase.o: phase.cpp phase.h cluster.h print.h globals.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) phase.cpp

clean:
	rm -f $(TARGET) *.o
