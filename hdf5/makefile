# CC=g++-10
# std=c++2a
# CFLAGS=-I.



# testhdf5.o: testhdf5.cpp hdf5.hpp
# 	g++ -std=c++2a -o testhdf5.o testhdf5.cpp -lhdf5 -lhdf5_cpp -L/home/kolanzl/Documents/hdf5-1.14.1-2-Std-ubuntu2204_64/hdf/HDF5-1.14.1-Linux/HDF_Group/HDF5/1.14.1/lib
h5cpp: testh5cpp.cpp h5CPP.hpp
	g++ -std=c++2a -o testh5CPP.o testh5cpp.cpp -lhdf5 -lhdf5_cpp -L/home/kolanzl/Documents/hdf5-1.14.1-2-Std-ubuntu2204_64/hdf/HDF5-1.14.1-Linux/HDF_Group/HDF5/1.14.1/lib -I/home/kolanzl/Documents/hdf5-1.14.1-2-Std-ubuntu2204_64/hdf/HDF5-1.14.1-Linux/HDF_Group/HDF5/1.14.1/include

test: unitTest.cpp h5CPP.hpp
	g++ -std=c++2a -o unitTest.o unitTest.cpp -lhdf5 -lhdf5_cpp -L/home/kolanzl/Documents/hdf5-1.14.1-2-Std-ubuntu2204_64/hdf/HDF5-1.14.1-Linux/HDF_Group/HDF5/1.14.1/lib -I/home/kolanzl/Documents/hdf5-1.14.1-2-Std-ubuntu2204_64/hdf/HDF5-1.14.1-Linux/HDF_Group/HDF5/1.14.1/include


clean:
	rm -f *.o
	rm -f *.hdf5
	rm -f *.h5


