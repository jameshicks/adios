EXEC = adios
# add -Wc++98-compat-pedantic to warn on c++11 features
CXXFLAGS = -std=c++11 -Wall -pedantic-errors -Wunreachable-code -O2 -g -Wno-unused-parameter -Wunused-function -Wextra
INCLUDES = -Iinclude
COMMON_SOURCES = ArgumentParser.cpp FileIOManager.cpp HiddenMarkov.cpp Linalg.cpp adios.cpp combinatorics.cpp common.cpp datamodel.cpp setops.cpp stringops.cpp vcf.cpp
COMMON_OBJECTS = $(COMMON_SOURCES:.cpp=.o)

LDFLAGS = -lm

UNITTEST_SOURCES = $(wildcard unittests/*.cpp)
UNITTEST_OBJECTS = $(UNITTEST_SOURCES:.cpp=.o)

CPPU_CXXFLAGS = $(CXXFLAGS) $(shell pkg-config cpputest --cflags) -Wno-keyword-macro
CPPU_LDFLAGS = $(LDFLAGS) $(shell pkg-config cpputest --libs)



.PHONY: $(EXEC) clean unittest all

all: $(EXEC)

$(EXEC): $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c main.cpp -o main.o
	$(CXX) main.o $(COMMON_OBJECTS) -o $(EXEC) $(LDFLAGS) 

$(COMMON_OBJECTS): %.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(UNITTEST_OBJECTS): %.o: %.cpp $(COMMON_OBJECTS) 
	$(CXX) $(CPPU_CXXFLAGS) $(INCLUDES)  -c $< -o $@

clean:
	rm -rf $(EXEC) $(COMMON_OBJECTS) main.o
	rm -rf unittests/unittester $(UNITTEST_OBJECTS) AllTests.o

unittest: $(COMMON_OBJECTS) $(UNITTEST_OBJECTS)
	$(CXX) $(CPPU_CXXFLAGS) $(INCLUDES) -c unittests/AllTests.cpp -o AllTests.o
	$(CXX) $(COMMON_OBJECTS) $(UNITTEST_OBJECTS) -o unittests/unittester $(CPPU_LDFLAGS)
	./unittests/unittester -v

