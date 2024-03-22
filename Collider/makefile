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

# CXX = g++
# CXXFLAGS = -g -O2 -std=c++2a -fopenmp -I${HDF5_DIR}/include
# LIBS = -L${HDF5_DIR}/lib -lhdf5 -lhdf5_cpp

# # Source files
# SOURCES = Collider.cpp 
# HPP = ../data/DECCOData.hpp ball_group.hpp



# # Object files for both regular and MOC-generated sources
# OBJECTS = $(SOURCES:.cpp=.o) 

# # Executable name
# EXEC = Collider.x

# # Default target
# all: $(HPP) $(EXEC) 

# # Linking the executable
# $(EXEC): $(OBJECTS) $(HPP)
# 	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJECTS) $(LIBS)

# # Generic rule for compiling
# %.o: %.cpp $(HPP)
# 	$(CXX) $(CXXFLAGS) -c $<



# # Clean
# clean:
# 	rm -f $(OBJECTS) $(EXEC) $(MOC_SOURCES) 





# Base CXXFLAGS, without HDF5-specific flags
CXXFLAGS = -g -O2 -std=c++2a -fopenmp

# Detect HDF5
# HDF5_TEST := $(shell ldconfig -p | grep hdf5 >/dev/null 2>&1 && echo 1 || echo 0)
HDF5_TEST = 1
MPI_TEST = 1
GPU_TEST = 1

CXX = CC
# Conditional HDF5 flags
ifeq ($(HDF5_TEST), 1)

    HDF5_CXXFLAGS = -I/usr/include/hdf5/serial #Perlmutter

#     HDF5_CXXFLAGS = -I/usr/include/hdf5/serial
#     HDF5_LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_cpp

    #For Perlmutter and Novus
    HDF5_CXXFLAGS = -I${HDF5_DIR}/include 
    HDF5_LIBS = -L${HDF5_DIR}/lib -lhdf5 -lhdf5_cpp

    CXXFLAGS += $(HDF5_CXXFLAGS) -DHDF5_ENABLE
    LIBS = $(HDF5_LIBS)
else
    LIBS =
endif

ifeq ($(MPI_TEST), 1)
    CXXFLAGS += -DMPI_ENABLE -I$(CRAY_MPICH_DIR)/include
    LIBS += -L$(CRAY_MPICH_DIR)/lib
endif

ifeq ($(GPU_TEST), 1)
    CXXFLAGS += -DGPU_ENABLE 
endif

# Source files
SOURCES = Collider.cpp 
HPP += ../data/DECCOData.hpp  ball_group.hpp ../utilities/Utils.hpp ../utilities/MPI_utilities.hpp

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