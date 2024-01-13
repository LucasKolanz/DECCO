CXX = CC
CXXFLAGS = -g -O2 -std=c++2a -fopenmp -DMPI_ENABLE -I${HDF5_DIR}/include
LIBS = -L${HDF5_DIR}/lib -lhdf5 -lhdf5_cpp

# Source files
SOURCES = Collider.cpp 
HPP = ../data/DECCOData.hpp ball_group.hpp



# Object files for both regular and MOC-generated sources
OBJECTS = $(SOURCES:.cpp=.o) 

# Executable name
EXEC = Collider.x

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