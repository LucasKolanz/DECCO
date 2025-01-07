# Define compiler and flags
CXX := g++
CXXFLAGS := -std=c++11 -O2 -I/usr/local/include/open3d -I/usr/include/eigen3
LDFLAGS := -lopen3d
# LDFLAGS := -L/usr/local/lib -lopen3d

# Define target and source files
TARGET := fractDim
SRC := fractDim.cpp

# Define build rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Clean rule to remove binary
clean:
	rm -f $(TARGET)
