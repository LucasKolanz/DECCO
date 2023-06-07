#include "h5CPP.hpp"
#include <iostream>
#include <vector>

bool uTest0()
{
    std::string ut_file = "unitTest0.h5";
    std::vector<double> data_to_write = {1.1, 2.2, 3.3, 4.4, 5.5};
    std::vector<double> more_data_to_write = {6.6, 7.7, 8.8, 9.9, 10.10};
    HDF5Handler handler(ut_file);
    // Write data to file
    handler.createAppendFile(data_to_write,"simData");
    // Append more data to file
    handler.createAppendFile(data_to_write,"constants");
    handler.createAppendFile(more_data_to_write,"constants");
    // Read data from file
    std::vector<double> data_read_simData = handler.readFile("simData");
    std::vector<double> data_read_constants = handler.readFile("constants");
    bool passed = true;
    for (int i = 0; i < data_to_write.size(); ++i)
    {
        if (data_to_write[i] != data_read_simData[i])
        {
            passed = false;
        }
    }
    for (int i = 0; i < more_data_to_write.size(); ++i)
    {
        if (data_to_write[i] != data_read_constants[i])
        {
            passed = false;
        }
        if (more_data_to_write[i] != data_read_constants[i+more_data_to_write.size()])
        {
            passed = false;
        }
    }

    return passed;
}

bool uTest1()
{
    std::string ut_file = "unitTest1.h5";
    std::vector<double> data_to_write = {1.1, 2.2, 3.3, 4.4, 5.5};
    std::vector<double> more_data_to_write = {6.6, 7.7, 8.8, 9.9, 10.10};
    HDF5Handler handler(ut_file,10);
    // Write data to file
    handler.createAppendFile(data_to_write,"simData");
    // Append more data to file
    handler.createAppendFile(data_to_write,"constants");
    handler.createAppendFile(more_data_to_write,"constants",data_to_write.size());
    // Read data from file
    std::vector<double> data_read_simData = handler.readFile("simData");
    std::vector<double> data_read_constants = handler.readFile("constants");
    bool passed = true;
    for (int i = 0; i < data_to_write.size(); ++i)
    {
        // std::cerr<<data_to_write[i]<<", "<<data_read_simData[i]<<std::endl;
        if (data_to_write[i] != data_read_simData[i])
        {
            passed = false;
        }
    }
    for (int i = 0; i < more_data_to_write.size(); ++i)
    {
        // std::cerr<<data_to_write[i]<<", "<<data_read_constants[i]<<std::endl;
        if (data_to_write[i] != data_read_constants[i])
        {
            passed = false;
        }
        std::cerr<<more_data_to_write[i]<<", "<<data_read_constants[i+more_data_to_write.size()]<<std::endl;
        if (more_data_to_write[i] != data_read_constants[i+more_data_to_write.size()])
        {
            passed = false;
        }
    }

    return passed;
}

int main() {

    int num_tests = 2;
    bool tests[num_tests];
    tests[0] = uTest0();
    tests[1] = uTest1();

    std::cerr<<"==========START UNIT TESTS=========="<<std::endl;
    std::cerr<<"Create (two different dataspaces for one file), append (to one of these), read (both dataspaces), unlimited size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[0]<<std::endl;
    std::cerr<<"Create (two different dataspaces for one file), append (to one of these), read (both dataspaces), fixed size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[1]<<std::endl;
    std::cerr<<"===========END UNIT TESTS==========="<<std::endl;
    


    return 0;
}
