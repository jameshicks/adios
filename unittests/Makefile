EXEC = adios_tester
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
ADIOS_OBJECTS = $(filter-out ../main.o, $(wildcard ../*.o) )

CPPU_CC_FLAGS = $(shell pkg-config cpputest --cflags) -g -Wno-keyword-macro
CPPU_LD_FLAGS = $(shell pkg-config cpputest --libs)

$(EXEC): $(OBJECTS)
		$(CXX) -g $(OBJECTS) $(ADIOS_OBJECTS) -o $(EXEC) $(CPPU_LD_FLAGS)

%.o: %.cpp
	$(CXX) -g -c $(CXXFLAGS) $(CPPU_CC_FLAGS) -I../ $< -o $@

clean:
	rm -rf $(EXEC) $(OBJECTS)

