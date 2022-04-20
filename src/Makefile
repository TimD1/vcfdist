CC = /opt/rh/devtoolset-8/root/bin/g++
CFLAGS = -g -Wall -std=c++2a -O2
IFLAGS = /usr/local/include/htslib/
TARGET = compare
LIB = /usr/local/lib
LDLIBS = -lz -lhts -lstdc++fs

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -I$(IFLAGS) -L $(LIB) -o $(TARGET) $(TARGET).cpp $(LDLIBS)

clean:
	rm $(TARGET)
