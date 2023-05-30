#include "h5CPP.hpp"
#include <iostream>
#include <vector>

int main() {
    std::cerr<<"START"<<std::endl;
    HDF5Handler handler("test.h5",10);
    std::cerr<<"CONSTRUCTED"<<std::endl;

    // Write data to file
    std::vector<double> data_to_write = {1.1, 2.2, 3.3, 4.4, 5.5};
    handler.createAppendFile(data_to_write,"simData");
    std::cerr<<"First Append"<<std::endl;

    // Append more data to file
    std::vector<double> more_data_to_write = {6.6, 7.7, 8.8, 9.9, 10.10};
    handler.createAppendFile(more_data_to_write,"constants");
    std::cerr<<"Second Append"<<std::endl;
    handler.createAppendFile(more_data_to_write,"constants",4);
    std::cerr<<"Third Append"<<std::endl;

    // Read data from file
    std::vector<double> data_read = handler.readFile("simData");

    // Print out the data read from the file
    for (double d : data_read) {
        std::cout << d << "  ";
    }
    std::cout << std::endl;

    data_read = handler.readFile("constants");
    // Print out the data read from the file
    for (double d : data_read) {
        std::cout << d << ", ";
    }
    std::cout << std::endl;

    return 0;
}
