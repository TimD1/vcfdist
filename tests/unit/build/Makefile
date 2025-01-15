######################
# GOOGLE TEST CONFIG #
######################
GTEST_DIR = /usr/src/gtest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_FLAGS += -isystem $(GTEST_DIR)/include
GTEST_SRCS = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

###################
# MY TESTS CONFIG #
###################
SRC = ../../../src
TESTS = test_vcfdist
CXXFLAGS += -g -Wall -Wextra -std=c++17
LDFLAGS += -lz -lhts -lpthread -lm -lstdc++
OBJS = $(SRC)/globals.o $(SRC)/print.o $(SRC)/variant.o $(SRC)/dist.o $(SRC)/bed.o \
       $(SRC)/cluster.o $(SRC)/phase.o $(SRC)/timer.o test_dist.o

############
# MY TESTS #
############
$(TESTS): $(OBJS) gtest_main.a
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TESTS) $(LDFLAGS) gtest_main.a

clean:
	rm -f $(TESTS) gtest.a gtest_main.a *.o $(SRC)/*.o *.a

#####################
# BUILD GOOGLE TEST #
#####################
gtest-all.o: $(GTEST_SRCS)
	$(CXX) $(GTEST_FLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o: $(GTEST_SRCS)
	$(CXX) $(GTEST_FLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest_main.a: gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

##################
# BUILD MY TESTS #
##################
test_dist.o: ../dist/test_dist.cpp $(SRC)/dist.h $(GTEST_HEADERS)
	$(CXX) $(CXXFLAGS) -c ../dist/test_dist.cpp

#################
# BUILD VCFDIST #
#################
$(SRC)/bed.o: $(SRC)/bed.cpp $(SRC)/bed.h $(SRC)/print.h $(SRC)/defs.h $(SRC)/globals.h
	$(CXX) -c $(CXXFLAGS) $(SRC)/bed.cpp -o $(SRC)/bed.o

$(SRC)/cluster.o: $(SRC)/cluster.cpp $(SRC)/cluster.h $(SRC)/variant.h $(SRC)/fasta.h \
                  $(SRC)/globals.h $(SRC)/dist.h $(SRC)/defs.h
	$(CXX) -c $(CXXFLAGS) $(SRC)/cluster.cpp -o $(SRC)/cluster.o

$(SRC)/dist.o: $(SRC)/dist.cpp $(SRC)/dist.h $(SRC)/fasta.h $(SRC)/variant.h $(SRC)/cluster.h \
               $(SRC)/print.h $(SRC)/defs.h $(SRC)/globals.h
	$(CXX) -c $(CXXFLAGS) $(SRC)/dist.cpp -o $(SRC)/dist.o

$(SRC)/globals.o: $(SRC)/globals.cpp $(SRC)/globals.h $(SRC)/bed.h $(SRC)/print.h $(SRC)/defs.h \
                  $(SRC)/timer.h
	$(CXX) -c $(CXXFLAGS) $(SRC)/globals.cpp -o $(SRC)/globals.o

$(SRC)/phase.o: $(SRC)/phase.cpp $(SRC)/phase.h $(SRC)/cluster.h $(SRC)/print.h $(SRC)/globals.h \
                $(SRC)/defs.h $(SRC)/variant.h
	$(CXX) -c $(CXXFLAGS) $(SRC)/phase.cpp -o $(SRC)/phase.o

$(SRC)/print.o: $(SRC)/print.cpp $(SRC)/print.h $(SRC)/globals.h $(SRC)/phase.h $(SRC)/dist.h \
                $(SRC)/defs.h $(SRC)/variant.h
	$(CXX) -c $(CXXFLAGS) $(SRC)/print.cpp -o $(SRC)/print.o

$(SRC)/timer.o: $(SRC)/timer.cpp $(SRC)/timer.h $(SRC)/globals.h $(SRC)/defs.h
	$(CXX) -c $(CXXFLAGS) $(SRC)/timer.cpp -o $(SRC)/timer.o

$(SRC)/variant.o: $(SRC)/variant.cpp $(SRC)/variant.h $(SRC)/print.h $(SRC)/fasta.h $(SRC)/defs.h \
                  $(SRC)/globals.h
	$(CXX) -c $(CXXFLAGS) $(SRC)/variant.cpp -o $(SRC)/variant.o
