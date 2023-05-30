#include <iostream>
#include "hdf5.hpp"
#include <sstream>

int main()
{
	std::stringstream energyBuffer;
	//make dummy data
	for (int i = 0; i < 1000; i++)
	{
		if (i < 999)
		{
			energyBuffer<<i;
		}
		else
		{
			energyBuffer<<i<<'\n';
		}
	}


	HDF5File f("testhdf5file.hdf5");
	f.write_sim_data(energyBuffer);
	// f.close();

	// HDF5File f1("testhdf5file.hdf5");
	// std::stringstream buffer;
	// f1.read_data("/simData","data",buffer);
	// std::cout<<buffer.str();
}