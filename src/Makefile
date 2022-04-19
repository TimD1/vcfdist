CC = g++
CFLAGS = -g -Wall -std=c++11 -O2
IFLAGS = /usr/local/include/htslib/
TARGET = compare
LIB = /usr/local/lib

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -I$(IFLAGS) -L $(LIB) -o $(TARGET) $(TARGET).cpp -lz -lhts

clean:
	rm $(TARGET)
