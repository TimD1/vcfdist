######################
# GOOGLE TEST CONFIG #
######################
# Portable across Linux CI (libgtest-dev sources at /usr/src/gtest) and macOS
# (Homebrew prebuilt Google Test, not on the default compiler search path).
UNAME_S := $(shell uname -s)

###################
# MY TESTS CONFIG #
###################
SRC = ../../../src
TEST_SRC = ../src
TESTS = test_vcfdist
CXXFLAGS += -g -Wall -Wextra -std=c++17
# var_fields members with no default are required; omitting one must not compile (GCC only)
CXXFLAGS += -Werror=missing-field-initializers
LDFLAGS += -lz -lhts -lpthread -lm -lstdc++
SRC_HEADERS = $(wildcard $(SRC)/*.h)
TEST_HEADERS = $(wildcard $(TEST_SRC)/*.h)
OBJS = $(SRC)/globals.o $(SRC)/print.o $(SRC)/variant.o $(SRC)/dist.o $(SRC)/bed.o \
       $(SRC)/cluster.o $(SRC)/phase.o $(SRC)/timer.o \
       test_helpers.o test_bed.o test_cluster.o test_defs.o test_dist.o test_globals.o \
       test_harness.o test_phase.o test_print.o test_timer.o test_variant.o

ifeq ($(UNAME_S),Darwin)
    # macOS: link Homebrew's prebuilt Google Test; no gtest sources to compile.
    GTEST_PREFIX := $(shell brew --prefix googletest)
    GTEST_INC := -isystem $(GTEST_PREFIX)/include
    GTEST_LINK := -L$(GTEST_PREFIX)/lib -lgtest -lgtest_main
    GTEST_LIB :=
    GTEST_HEADERS :=
else
    # Linux (CI): build Google Test from source at /usr/src/gtest.
    GTEST_DIR = /usr/src/gtest
    GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                    $(GTEST_DIR)/include/gtest/internal/*.h
    GTEST_INC := -isystem $(GTEST_DIR)/include
    GTEST_SRCS = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
    GTEST_LINK := gtest_main.a
    GTEST_LIB := gtest_main.a
endif

############
# MY TESTS #
############
$(TESTS): $(OBJS) $(GTEST_LIB)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TESTS) $(LDFLAGS) $(GTEST_LINK)

clean:
	rm -f $(TESTS) gtest.a gtest_main.a *.o $(SRC)/*.o *.a

#####################
# BUILD GOOGLE TEST #
#####################
# Linux only; on macOS GTEST_LIB is empty so these rules are never invoked.
gtest-all.o: $(GTEST_SRCS)
	$(CXX) $(GTEST_INC) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o: $(GTEST_SRCS)
	$(CXX) $(GTEST_INC) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest_main.a: gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

##################
# BUILD MY TESTS #
##################
# One pattern rule covers every test file; adding a test file only requires adding it to OBJS.
# Every test object depends on all vcfdist and test headers, since each test may include any of
# them and the set is small enough that over-rebuilding costs less than a stale object.
%.o: $(TEST_SRC)/%.cpp $(SRC_HEADERS) $(TEST_HEADERS) $(GTEST_HEADERS)
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
