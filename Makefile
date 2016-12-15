EXEC = adios
# add -Wc++98-compat-pedantic to warn on c++11 features
CXXFLAGS = -std=c++11 -Wall -pedantic-errors -Wunreachable-code -O2 -g -Wno-unused-parameter -Wunused-function -Wextra
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
LDFLAGS = -lm
export CXXFLAGS

$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(EXEC) $(LDFLAGS) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(EXEC) $(OBJECTS)
	make -C unittests clean

test: $(OBJECTS)
	make -C unittests
	./unittests/adios_tester -v

