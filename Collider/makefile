



CXX = CC
CXX = nvcc
CXX = g++
CXX = mpic++
CXX = nvc++


# Base CXXFLAGS, without HDF5-specific flags
# CXXFLAGS = -g -O2 -std=c++17 -fopenmp -fopenacc -foffload=nvptx-none -#-fopenacc 
CXXFLAGS = -g -O2 -std=c++17 -fopenmp #-fopenacc 

# Detect HDF5
# HDF5_TEST := $(shell ldconfig -p | grep hdf5 >/dev/null 2>&1 && echo 1 || echo 0)
HDF5_TEST = 1
MPI_TEST = 1
GPU_TEST = 1
EXPERIMENTAL_FILESYSTEM = 1

LIBS += -L$(LD_LIBRARY_PATH)

# Conditional HDF5 flags
ifeq ($(HDF5_TEST), 1)

# 	HDF5_CXXFLAGS = -I/usr/include/hdf5/serial #Perlmutter

#     HDF5_CXXFLAGS = -I/usr/include/hdf5/serial
#     HDF5_LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_cpp

	#For Perlmutter and Novus

# 	HDF5_CXXFLAGS = -I${HDF5_DIR}/include 
	HDF5_LIBS = -lhdf5 -lhdf5_cpp
# 	HDF5_LIBS = -L${HDF5_DIR}/lib -lhdf5 -lhdf5_cpp

	CXXFLAGS += $(HDF5_CXXFLAGS) -DHDF5_ENABLE
	LIBS += $(HDF5_LIBS)
# else
# 	LIBS =
endif

ifeq ($(MPI_TEST), 1)
# 	#Works for Engr
# 	CXXFLAGS += -DMPI_ENABLE #-I$(CRAY_MPICH_DIR)/include
# 	LIBS += -lmpi
	#Works for Novus
	CXXFLAGS += -DMPI_ENABLE -I$(MPI_DIR)/include
	LIBS += -lmpi -L$(MPI_DIR)/lib


# 	LIBS += -L$(CRAY_MPICH_DIR)/lib
# 	CXXFLAGS += -DMPI_ENABLE -I/apps/intel/compilers_and_libraries_2017.4.196/linux/mpi/mic/lib/release/libmpi.so.12/include
# 	CXXFLAGS += -DMPI_ENABLE 
# 	LIBS += -lmpi
endif

ifeq ($(GPU_TEST), 1)
	#Novus
#     CXXFLAGS += -DGPU_ENABLE -Wl,-znoexecstack #-fopt-info-omp -gpu=cuda11.8
    CXXFLAGS += -DGPU_ENABLE -acc=gpu -Minfo=accel -Wl,-znoexecstack #-fopt-info-omp -gpu=cuda11.8
#     CXXFLAGS += -DGPU_ENABLE -acc=gpu -Minfo=accel -Wl,-znoexecstack -gpu=cuda11.8
#     CXXFLAGS += -DGPU_ENABLE -Wl,-znoexecstack #-fopt-info-omp -gpu=cuda11.8
# 	CXXFLAGS += -DGPU_ENABLE -Wl,-znoexecstack -fdump-tree-optimized
endif

ifeq ($(EXPERIMENTAL_FILESYSTEM),1)
	CXXFLAGS += -DEXPERIMENTAL_FILESYSTEM
	LIBS += -lstdc++fs
endif

# Source files
SOURCES = Collider.cpp ball_group.cpp ../utilities/MPI_utilities.cpp ../utilities/vec3.cpp ../data/DECCOData.cpp ../utilities/Utils.cpp ../utilities/simple_graph.cpp
TEST_SOURCES = testing.cpp ball_group.cpp ../utilities/MPI_utilities.cpp ../utilities/vec3.cpp ../data/DECCOData.cpp ../utilities/Utils.cpp ../utilities/simple_graph.cpp
HPP += ../data/DECCOData.hpp ball_group.hpp ../utilities/Utils.hpp ../utilities/MPI_utilities.hpp ../utilities/Utils.hpp ../utilities/simple_graph.hpp

# Object directory
OBJDIR = objs

# Extract the base names of the source files
SOURCE_NAMES = $(notdir $(SOURCES))
TEST_SOURCE_NAMES = $(notdir $(TEST_SOURCES))

# Object files in OBJDIR
OBJECTS = $(addprefix $(OBJDIR)/, $(SOURCE_NAMES:.cpp=.o))

# Object files in OBJDIR
TEST_OBJECTS = $(addprefix $(OBJDIR)/,  $(TEST_SOURCE_NAMES:.cpp=.o))

# Source directories
SRC_DIRS = $(sort $(dir $(SOURCES)))

# VPATH to find source files
VPATH = $(SRC_DIRS)

# Executable name
EXEC = Collider.x
TEST_EXEC = testing.x

# Default target
all: $(OBJDIR) $(HPP) $(EXEC)

# test target
# test: $(OBJDIR) $(HPP) $(TEST_EXEC)
test: $(TEST_EXEC)

# Create OBJDIR if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Linking the executable
$(EXEC): $(OBJECTS) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJECTS) $(LIBS)

# Linking test.x
$(TEST_EXEC): $(TEST_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TEST_EXEC) $(TEST_OBJECTS) $(LIBS)

# Compile rule
# $(OBJDIR)/%.o: %.cpp $(HPP)
$(OBJDIR)/%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@


# Clean
clean:
	rm -f $(OBJDIR)/*.o $(EXEC) $(TEST_EXEC)