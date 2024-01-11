CFLAGS = -std=c++2a -O2 -I/usr/include/hdf5/serial/
LDFLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5 -lhdf5_cpp

# Define sources and objects
SOURCES = unitTest.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: unitTest


unitTest: unitTest.o
	g++ -o $@ $^ $(LDFLAGS)

%.o: %.cpp DECCOData.hpp
	g++ $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o
	rm -f *.hdf5
	rm -f *.h5
	rm -f unitTest
