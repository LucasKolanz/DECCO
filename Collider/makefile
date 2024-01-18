# # CC=g++-10
# # std=c++2a
# # CFLAGS=-I.
# # CPP=CC
# CPP=g++
# # CPPflags=-O3 -std=c++2a -fopenmp
# # CPPflags=-O3 -std=c++2a 
# CPPflags=-std=c++2a -O2

# ColliderSingleCore.x: ColliderSingleCore.cpp ../default_files/dust_const.hpp ball_group.hpp ../timing/timing.hpp
# 	${CPP} ${CPPflags} -o ColliderSingleCore.x ColliderSingleCore.cpp

# # ColliderSingleCore_COPY.x: ColliderSingleCore_COPY.cpp ../dust_const.hpp ../ball_group_COPY.hpp ../timing/timing.hpp
# # 	${CPP} ${CPPflags} -o ColliderSingleCore_COPY.x ColliderSingleCore_COPY.cpp

# clean:
# 	rm -f *.o
# 	rm -f *.x

CXX = g++
CXXFLAGS = -g -O2 -std=c++2a -fopenmp -I/usr/include/hdf5/serial/
LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5 -lhdf5_cpp

# Source files
SOURCES = ColliderSingleCore.cpp 
HPP = ../data/DECCOData.hpp ball_group.hpp



# Object files for both regular and MOC-generated sources
OBJECTS = $(SOURCES:.cpp=.o) 

# Executable name
EXEC = ColliderSingleCore.x

# Default target
all: $(HPP) $(EXEC) 

# Linking the executable
$(EXEC): $(OBJECTS) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJECTS) $(LIBS)

# Generic rule for compiling
%.o: %.cpp $(HPP)
	$(CXX) $(CXXFLAGS) -c $<



# Clean
clean:
	rm -f $(OBJECTS) $(EXEC) $(MOC_SOURCES) 