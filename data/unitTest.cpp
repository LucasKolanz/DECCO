#include "DECCOData.hpp"
#include <iostream>
#include <vector>
#include <stdlib.h>

bool compareVecs(std::vector<double> v1, std::vector<double> v2,int v1i, int v2i, int len)
{
    for (int i = 0; i < len; ++i)
    {
        // std::cerr<<v1[v1i]<<", "<<v2[v2i]<<std::endl;
        if (v1[v1i] != v2[v2i])
            return false;
        v1i++;
        v2i++;
    }
    return true;
}

bool hdf5uTest0()
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

bool hdf5uTest1()
{
    std::string ut_file = "unitTest1.h5";
    std::vector<double> data_to_write = {1.1, 2.2, 3.3, 4.4, 5.5};
    std::vector<double> more_data_to_write = {6.6, 7.7, 8.8, 9.9, 10.10};
    HDF5Handler handler(ut_file,true);
    // Write data to file
    handler.createAppendFile(data_to_write,"simData",0,10);
    // Append more data to file
    handler.createAppendFile(data_to_write,"constants",0,10);
    handler.createAppendFile(more_data_to_write,"constants",data_to_write.size(),data_to_write.size()*2);
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
        // std::cerr<<more_data_to_write[i]<<", "<<data_read_constants[i+more_data_to_write.size()]<<std::endl;
        if (more_data_to_write[i] != data_read_constants[i+more_data_to_write.size()])
        {
            passed = false;
        }
    }


    return passed;
}

bool hdf5uTest2()
{
    std::string ut_file = "unitTest2.h5";
    
    HDF5Handler handler(ut_file);
    std::vector<double> data_to_write = {1.1, 2.2, 3.3, 4.4, 5.5};
    std::vector<double> more_data_to_write = {6.6, 7.7, 8.8, 9.9, 10.10};
    // Write data to file
    handler.createAppendFile(data_to_write,"simData");
    //Write metadata to file
    std::string metadata = "Hello, I am metadata. Who are you?";
    handler.attachMetadataToDataset(metadata,"simData");
    // Read metadata from file
    std::string readMetaData = handler.readMetadataFromDataset("simData");
    // Write some more data
    handler.createAppendFile(more_data_to_write,"simData");
    
    // std::cerr<<"read: "<<readMetaData<<"\tmeta: "<<metadata<<std::endl;
    std::vector<double> data_read_simData = handler.readFile("simData");
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
        // std::cerr<<data_to_write[i]<<", "<<data_read_simData[i]<<std::endl;
        // std::cerr<<more_data_to_write[i]<<", "<<data_read_simData[i+data_to_write.size()]<<std::endl;
        if (more_data_to_write[i] != data_read_simData[i+data_to_write.size()])
        {
            passed = false;
        }
    }
    return passed && (readMetaData == metadata);
}

bool hdf5uTest3()
{
    std::string ut_file = "unitTest3.h5";
    
    HDF5Handler handler(ut_file,true);
    std::vector<double> data_to_write = {1.1, 2.2, 3.3, 4.4, 5.5};
    std::vector<double> more_data_to_write = {6.6, 7.7, 8.8, 9.9, 10.10};
    // Write data to file
    handler.createAppendFile(data_to_write,"simData",0,10);
    //Write metadata to file
    std::string metadata = "Hello, I am metadata. Who are you?";
    handler.attachMetadataToDataset(metadata,"simData");
    // Write some more data
    handler.createAppendFile(more_data_to_write,"simData",data_to_write.size(),10);
    // Read metadata from file
    std::string readMetaData = handler.readMetadataFromDataset("simData");
    
    // std::cerr<<"read: "<<readMetaData<<"\tmeta: "<<metadata<<std::endl;
    std::vector<double> data_read_simData = handler.readFile("simData");
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
        // std::cerr<<more_data_to_write[i]<<", "<<data_read_simData[i+data_to_write.size()]<<std::endl;
        if (more_data_to_write[i] != data_read_simData[i+data_to_write.size()])
        {
            passed = false;
        }
    }
    return passed && (readMetaData == metadata);
}


//Write and read with hdf5 non-fixed size
bool DatauTest0()
{
    int num_particles = 10;
    int line_writes = 4;
    DECCOData data(std::string("datauTest0.h5"),num_particles);

    std::vector<double> simData_to_write(data.getWidth("simData")*line_writes);
    std::vector<double> constants_to_write(data.getWidth("constants")*line_writes);
    std::vector<double> energy_to_write(data.getWidth("energy")*line_writes);

    //make three lines of simData
    for (int i = 0; i < data.getWidth("simData")*line_writes; ++i)
    {
        simData_to_write[i] = i+i*0.1;
    }

    //make three lines of constants
    for (int i = 0; i < data.getWidth("constants")*line_writes; ++i)
    {
        constants_to_write[i] = (i)+100;
        // std::cerr<<constants_to_write[i]<<std::endl;
    }

    //make three lines of energy
    for (int i = 0; i < data.getWidth("energy")*line_writes; ++i)
    {
        energy_to_write[i] = (i+i*0.1)+1000;
    }

    bool b10 = data.Write(simData_to_write,"simData");
    bool b11 = data.Write(constants_to_write,"constants");
    bool b12 = data.Write(constants_to_write,"constants");
    bool b13 = data.Write(energy_to_write,"energy");

    //Verify metaData was written and is correct
    bool b50 = data.ReadMetaData("simData") == data.genMetaData(0);
    bool b51 = data.ReadMetaData("energy") == data.genMetaData(1);
    bool b52 = data.ReadMetaData("constants") == data.genMetaData(2);
    // bool b53 = data.ReadMetaData("timing") == data.genMetaData(3);

    //Read all of simData
    std::vector<double> data_read0 = data.Read("simData",true);
    std::vector<double> data_read4 = data.Read("energy");
    std::vector<double> data_read2 = data.Read("constants",false,1);
    std::vector<double> data_read3 = data.Read("constants",false,-1);
    std::vector<double> data_read1 = data.Read("constants",false,2);

    //verify simData read
    bool b1 = compareVecs(data_read0,simData_to_write,0,0,data_read0.size());
    //verify constants third line read
    bool b2 = compareVecs(data_read1,constants_to_write,0,data.getWidth("constants")*2,data_read1.size());
    //verify constants second line read
    bool b3 = compareVecs(data_read2,constants_to_write,0,data.getWidth("constants")*1,data_read2.size());
    //verify constants last line read
    bool b4 = compareVecs(data_read3,constants_to_write,0,constants_to_write.size()-data.getWidth("constants")*1,data_read3.size());
    //verify energy read
    bool b5 = compareVecs(data_read4,energy_to_write,0,0,data_read4.size());

    return b1 && b2 && b3 && b4 && b5 && b10 && b11 && b12 && b13 && b50 && b51 && b52;
}

//Write and read with hdf5 fixed size
bool DatauTest1()
{
    int num_particles = 11;
    int line_writes = 4;
    int writes = 4;
    int steps = 5;
    DECCOData data(std::string("datauTest1.h5"),num_particles,writes,steps);

    std::vector<double> simData_to_write(data.getWidth("simData")*writes*num_particles);
    std::vector<double> constants_to_write(data.getWidth("constants")*line_writes*num_particles/2);
    std::vector<double> energy_to_write(data.getWidth("energy")*line_writes*num_particles);

    //make three lines of simData
    for (int i = 0; i < data.getWidth("simData")*line_writes*num_particles; ++i)
    {
        simData_to_write[i] = i+i*0.1;
    }

    //make three lines of constants
    for (int i = 0; i < data.getWidth("constants")*line_writes*num_particles/2; ++i)
    {
        constants_to_write[i] = (i)+100;
    }

    //make three lines of energy
    for (int i = 0; i < data.getWidth("energy")*line_writes*num_particles; ++i)
    {
        energy_to_write[i] = (i+i*0.1)+1000;
    }

    // std::cerr<<data.max_size[2]<<std::endl;
    // std::cerr<<data.written_so_far[2]<<std::endl;
    // std::cerr<<"FIRSTWRITE"<<std::endl;
    bool b10 = data.Write(simData_to_write,"simData");
    // std::cerr<<"SECWRITE"<<std::endl;
    bool b11 = data.Write(constants_to_write,"constants");
    // std::cerr<<"THIRDWRITE"<<std::endl;
    bool b12 = data.Write(constants_to_write,"constants");
    // std::cerr<<"FOURWRITE"<<std::endl;
    bool b13 = data.Write(energy_to_write,"energy");
    // std::cerr<<"ENDWRITE"<<std::endl;

    // //Verify metaData was written and is correct
    bool b50 = data.ReadMetaData("simData") == data.genMetaData(0);
    bool b51 = data.ReadMetaData("energy") == data.genMetaData(1);
    bool b52 = data.ReadMetaData("constants") == data.genMetaData(2);
    // bool b53 = data.ReadMetaData("timing") == data.genMetaData(3);

    // //Read all of simData
    std::vector<double> data_read0 = data.Read("simData",true);
    std::vector<double> data_read4 = data.Read("energy");
    std::vector<double> data_read2 = data.Read("constants",false,1);
    std::vector<double> data_read3 = data.Read("constants",false,-1);
    std::vector<double> data_read1 = data.Read("constants",false,2);

    // //verify simData read
    bool b1 = compareVecs(data_read0,simData_to_write,0,0,data_read0.size());
    // //verify constants third line read
    bool b2 = compareVecs(data_read1,constants_to_write,0,data.getWidth("constants")*2,data_read1.size());
    // //verify constants second line read
    bool b3 = compareVecs(data_read2,constants_to_write,0,data.getWidth("constants")*1,data_read2.size());
    // //verify constants last line read
    bool b4 = compareVecs(data_read3,constants_to_write,0,constants_to_write.size()-data.getWidth("constants")*1,data_read3.size());
    // //verify energy read
    bool b5 = compareVecs(data_read4,energy_to_write,0,0,data_read4.size());

    return  b1 && b2 && b3 && b4 && b5 && b10 && b11 && b12 && b13 && b50 && b51 && b52;
}

bool DatauTest2()
{
    return false;
}

//Write and read with hdf5 non-fixed size
bool DatauTestCSV0()
{
    int num_particles = 10;
    int line_writes = 4;
    DECCOData data(std::string("datauTest0.csv"),num_particles);

    std::vector<double> simData_to_write(data.getWidth("simData")*line_writes);
    std::vector<double> constants_to_write(data.getWidth("constants")*line_writes);
    std::vector<double> energy_to_write(data.getWidth("energy")*line_writes);

    //make three lines of simData
    for (int i = 0; i < data.getWidth("simData")*line_writes; ++i)
    {
        simData_to_write[i] = i+i*0.1;
    }

    //make three lines of constants
    for (int i = 0; i < data.getWidth("constants")*line_writes; ++i)
    {
        constants_to_write[i] = (i)+100;
        // std::cerr<<constants_to_write[i]<<std::endl;
    }

    //make three lines of energy
    for (int i = 0; i < data.getWidth("energy")*line_writes; ++i)
    {
        energy_to_write[i] = (i+i*0.1)+1000;
    }

    bool b10 = data.Write(simData_to_write,"simData");
    bool b11 = data.Write(constants_to_write,"constants");
    bool b12 = data.Write(constants_to_write,"constants");
    bool b13 = data.Write(energy_to_write,"energy");

    //Verify metaData was written and is correct
    bool b50 = data.ReadMetaData("simData") == data.genMetaData(0);
    bool b51 = data.ReadMetaData("energy") == data.genMetaData(1);
    bool b52 = data.ReadMetaData("constants") == data.genMetaData(2);
    // bool b53 = data.ReadMetaData("timing") == data.genMetaData(3);

    //Read all of simData
    std::vector<double> data_read0 = data.Read("simData",true);
    std::vector<double> data_read4 = data.Read("energy");
    std::vector<double> data_read2 = data.Read("constants",false,1);
    std::vector<double> data_read3 = data.Read("constants",false,-1);
    std::vector<double> data_read1 = data.Read("constants",false,2);

    //verify simData read
    bool b1 = compareVecs(data_read0,simData_to_write,0,0,data_read0.size());
    //verify constants third line read
    bool b2 = compareVecs(data_read1,constants_to_write,0,data.getWidth("constants")*2,data_read1.size());
    //verify constants second line read
    bool b3 = compareVecs(data_read2,constants_to_write,0,data.getWidth("constants")*1,data_read2.size());
    //verify constants last line read
    bool b4 = compareVecs(data_read3,constants_to_write,0,constants_to_write.size()-data.getWidth("constants")*1,data_read3.size());
    //verify energy read
    bool b5 = compareVecs(data_read4,energy_to_write,0,0,data_read4.size());

    return b1 && b2 && b3 && b4 && b5 && b10 && b11 && b12 && b13 && b50 && b51 && b52;
}

//Write and read with hdf5 fixed size
bool DatauTestCSV1()
{
    int num_particles = 11;
    int line_writes = 4;
    int writes = 4;
    int steps = 5;
    DECCOData data(std::string("datauTest1.csv"),num_particles,writes,steps);

    std::vector<double> simData_to_write(data.getWidth("simData")*writes*num_particles);
    std::vector<double> constants_to_write(data.getWidth("constants")*line_writes*num_particles/2);
    std::vector<double> energy_to_write(data.getWidth("energy")*line_writes*num_particles);

    //make three lines of simData
    for (int i = 0; i < data.getWidth("simData")*line_writes*num_particles; ++i)
    {
        simData_to_write[i] = i+i*0.1;
    }

    //make three lines of constants
    for (int i = 0; i < data.getWidth("constants")*line_writes*num_particles/2; ++i)
    {
        constants_to_write[i] = (i)+100;
    }

    //make three lines of energy
    for (int i = 0; i < data.getWidth("energy")*line_writes*num_particles; ++i)
    {
        energy_to_write[i] = (i+i*0.1)+1000;
    }
    bool b10 = data.Write(simData_to_write,"simData");
    bool b11 = data.Write(constants_to_write,"constants");
    bool b12 = data.Write(constants_to_write,"constants");
    bool b13 = data.Write(energy_to_write,"energy");

    // //Verify metaData was written and is correct
    bool b50 = data.ReadMetaData("simData") == data.genMetaData(0);
    bool b51 = data.ReadMetaData("energy") == data.genMetaData(1);
    bool b52 = data.ReadMetaData("constants") == data.genMetaData(2);
    // bool b53 = data.ReadMetaData("timing") == data.genMetaData(3);

    // //Read all of simData
    std::vector<double> data_read0 = data.Read("simData",true);
    std::vector<double> data_read4 = data.Read("energy");
    std::vector<double> data_read2 = data.Read("constants",false,1);
    std::vector<double> data_read3 = data.Read("constants",false,-1);
    std::vector<double> data_read1 = data.Read("constants",false,2);

    // //verify simData read
    bool b1 = compareVecs(data_read0,simData_to_write,0,0,data_read0.size());
    // //verify constants third line read
    bool b2 = compareVecs(data_read1,constants_to_write,0,data.getWidth("constants")*2,data_read1.size());
    // //verify constants second line read
    bool b3 = compareVecs(data_read2,constants_to_write,0,data.getWidth("constants")*1,data_read2.size());
    // //verify constants last line read
    bool b4 = compareVecs(data_read3,constants_to_write,0,constants_to_write.size()-data.getWidth("constants")*1,data_read3.size());
    // //verify energy read
    bool b5 = compareVecs(data_read4,energy_to_write,0,0,data_read4.size());

    return  b1 && b2 && b3 && b4 && b5 && b10 && b11 && b12 && b13 && b50 && b51 && b52;
}

void hdf5Tests()
{
    int sys_call;
    int num_tests = 4;
    bool tests[num_tests];
    tests[0] = hdf5uTest0();
    sys_call = system("rm *.h5");
    tests[1] = hdf5uTest1();
    sys_call = system("rm *.h5");
    tests[2] = hdf5uTest2();
    sys_call = system("rm *.h5");
    tests[3] = hdf5uTest3();
    sys_call = system("rm *.h5");

    std::cerr<<"=======START hdf5 UNIT TESTS========"<<std::endl;
    std::cerr<<"Create (two different dataspaces for one file), append (to one of these), read (both dataspaces), unlimited size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[0]<<std::endl;
    std::cerr<<"Create (two different dataspaces for one file), append (to one of these), read (both dataspaces), fixed size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[1]<<std::endl;
    std::cerr<<"Create (metadata first and data second), read (metadata), unlimited size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[2]<<std::endl;
    std::cerr<<"Create (metadata first and data second), read (metadata), fixed size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[3]<<std::endl;
    std::cerr<<"===========END UNIT TESTS==========="<<std::endl;
}

void DataTests()
{
    int num_tests = 2;
    int sys_call;
    bool tests[num_tests];
    for (int i = 0; i < num_tests; ++i)
    {
        tests[i] = false;
    }
    tests[0] = DatauTest0();
    sys_call = system("rm *.h5");
    tests[1] = DatauTest1();
    sys_call = system("rm *.h5");
    // tests[2] = DatauTest2();
    // sys_call = system("rm *.h5");
    // tests[3] = DatauTest3();

    std::cerr<<"=======START data UNIT TESTS========"<<std::endl;
    std::cerr<<"Create (simData, constants, and energy), append (to constants), read (all), unlimited size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[0]<<std::endl;
    std::cerr<<"Create (simData, constants, and energy), append (to constants), read (all), fixed size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[1]<<std::endl;
    // std::cerr<<"Create (metadata first and data second), read (metadata), unlimited size"<<std::endl;
    // std::cerr<<"Unit test returned: "<<tests[2]<<std::endl;
    // std::cerr<<"Create (metadata first and data second), read (metadata), fixed size"<<std::endl;
    // std::cerr<<"Unit test returned: "<<tests[3]<<std::endl;
    std::cerr<<"===========END UNIT TESTS==========="<<std::endl;
}

void DataCSVTests()
{
    int num_tests = 2;
    int sys_call;
    bool tests[num_tests];
    for (int i = 0; i < num_tests; ++i)
    {
        tests[i] = false;
    }
    tests[0] = DatauTestCSV0();
    sys_call = system("rm *.h5");
    tests[1] = DatauTestCSV1();
    sys_call = system("rm *.h5");
    // tests[2] = DatauTest2();
    // tests[3] = DatauTest3();

    std::cerr<<"=======START data CSV UNIT TESTS========"<<std::endl;
    std::cerr<<"Create (simData, constants, and energy), append (to constants), read (all), unlimited size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[0]<<std::endl;
    std::cerr<<"Create (simData, constants, and energy), append (to constants), read (all), fixed size"<<std::endl;
    std::cerr<<"Unit test returned: "<<tests[1]<<std::endl;
    // std::cerr<<"Create (metadata first and data second), read (metadata), unlimited size"<<std::endl;
    // std::cerr<<"Unit test returned: "<<tests[2]<<std::endl;
    // std::cerr<<"Create (metadata first and data second), read (metadata), fixed size"<<std::endl;
    // std::cerr<<"Unit test returned: "<<tests[3]<<std::endl;
    std::cerr<<"===========END UNIT TESTS==========="<<std::endl;
}

int main() {

    hdf5Tests();
    DataTests();
    // DatauTestCSV0();
    


    return 0;
}
