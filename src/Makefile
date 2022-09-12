CC = /opt/rh/devtoolset-8/root/bin/g++
CFLAGS = -g -Wall -std=c++2a -O2 -fpermissive
IFLAGS = /usr/local/include/
OBJS = globals.o print.o vcf.o dist.o bed.o cluster.o phase.o
TARGET = vcfdist
LDLIBS = -lz -lhts -lstdc++fs

all: $(TARGET)

$(TARGET): $(OBJS) main.cpp
	$(CC) $(CFLAGS) $(OBJS) -I$(IFLAGS) -o $(TARGET) main.cpp $(LDLIBS)

globals.o: globals.cpp globals.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) globals.cpp

print.o: print.cpp print.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) print.cpp

vcf.o: vcf.cpp vcf.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) vcf.cpp

dist.o: dist.cpp dist.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) dist.cpp

bed.o: bed.cpp bed.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) bed.cpp

cluster.o: cluster.cpp cluster.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) cluster.cpp

phase.o: phase.cpp phase.h
	$(CC) -c $(CFLAGS) -I$(IFLAGS) phase.cpp

clean:
	rm -f $(TARGET) *.o
