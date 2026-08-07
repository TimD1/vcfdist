######################
# GOOGLE TEST CONFIG #
######################
# Both platforms link a prebuilt Google Test: Ubuntu's libgtest-dev ships libgtest{,_main}.a on the
# default search path, Homebrew's googletest is prebuilt but needs its prefix passed explicitly.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    GTEST_PREFIX := $(shell brew --prefix googletest)
    GTEST_INC := -isystem $(GTEST_PREFIX)/include
    GTEST_LINK := -L$(GTEST_PREFIX)/lib -lgtest_main -lgtest
else
    GTEST_INC :=
    GTEST_LINK := -lgtest_main -lgtest
endif

###################
# MY TESTS CONFIG #
###################
SRC = ../../../src
TEST_SRC = ../src
TESTS = test_vcfdist
# Must stay identical to src/Makefile: the $(SRC) objects are shared between the two builds.
CXXFLAGS = -g -O3
CXXFLAGS += -Wall -Wextra -std=c++17
# var_fields members with no default are required; omitting one must not compile (GCC only)
CXXFLAGS += -Werror=missing-field-initializers
LDFLAGS += -lz -lhts -lpthread -lm -lstdc++
SRC_HEADERS = $(wildcard $(SRC)/*.h)
TEST_HEADERS = $(wildcard $(TEST_SRC)/*.h)
OBJS = $(SRC)/globals.o $(SRC)/print.o $(SRC)/variant.o $(SRC)/dist.o $(SRC)/bed.o \
       $(SRC)/cluster.o $(SRC)/phase.o $(SRC)/timer.o \
       test_helpers.o test_bed.o test_cluster.o test_defs.o test_dist.o test_globals.o \
       test_harness.o test_phase.o test_print.o test_timer.o test_variant.o

############
# MY TESTS #
############
$(TESTS): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TESTS) $(LDFLAGS) $(GTEST_LINK)

# $(SRC)/*.o belongs to src/Makefile; removing it here would recompile vcfdist on every test build.
clean:
	rm -f $(TESTS) *.o

##################
# BUILD MY TESTS #
##################
# One pattern rule covers every test file; adding a test file only requires adding it to OBJS.
# Every test object depends on all vcfdist and test headers, since each test may include any of
# them and the set is small enough that over-rebuilding costs less than a stale object.
%.o: $(TEST_SRC)/%.cpp $(SRC_HEADERS) $(TEST_HEADERS)
	$(CXX) $(CXXFLAGS) $(GTEST_INC) -c $< -o $@

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
