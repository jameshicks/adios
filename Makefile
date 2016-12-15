EXEC = adios
# add -Wc++98-compat-pedantic to warn on c++11 features
CXXFLAGS = -std=c++11 -Wall -pedantic-errors -Wunreachable-code -O2 -g -Wno-unused-parameter -Wunused-function -Wextra
INCLUDES = -Iinclude
COMMON_SOURCES = ArgumentParser.cpp HiddenMarkov.cpp Linalg.cpp adios.cpp combinatorics.cpp common.cpp datamodel.cpp setops.cpp stringops.cpp vcf.cpp
COMMON_OBJECTS = $(COMMON_SOURCES:.cpp=.o)

UNITTEST_SOURCES = $(wildcard unittests/*.cpp)
UNITTEST_OBJECTS = $(UNITTEST_SOURCES:.cpp=.o)

CPPU_CC_FLAGS = $(shell pkg-config cpputest --cflags) -g -Wno-keyword-macro
CPPU_LD_FLAGS = $(shell pkg-config cpputest --libs)

LDFLAGS = -lm
export CXXFLAGS

adios: $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c main.cpp -o main.o
	$(CXX) main.o $(COMMON_OBJECTS) -o $(EXEC) $(LDFLAGS) 

$(COMMON_OBJECTS): %.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(UNITTEST_OBJECTS): %.o: %.cpp $(COMMON_OBJECTS) 
	$(CXX) $(CXXFLAGS) $(CPPU_CC_FLAGS) $(INCLUDES)  -c $< -o $@

clean:
	rm -rf $(EXEC) $(COMMON_OBJECTS)
	rm -rf unittests/adios_tester $(UNITTEST_OBJECTS)

unittest: $(COMMON_OBJECTS) $(UNITTEST_OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(CPPU_CC_FLAGS) -c unittests/AllTests.cpp -o AllTests.o
	$(CXX) $(INCLUDES) $(COMMON_OBJECTS) $(UNITTEST_OBJECTS) -o unittests/adios_tester $(LDFLAGS) $(CPPU_LD_FLAGS)
	./unittests/adios_tester -v

