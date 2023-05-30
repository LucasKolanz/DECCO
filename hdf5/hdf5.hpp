#include <iostream>
#include <sstream>
#include "/home/kolanzl/Documents/hdf5-1.14.1-2-Std-ubuntu2204_64/hdf/HDF5-1.14.1-Linux/HDF_Group/HDF5/1.14.1/include/hdf5.h"
#include "/home/kolanzl/Documents/hdf5-1.14.1-2-Std-ubuntu2204_64/hdf/HDF5-1.14.1-Linux/HDF_Group/HDF5/1.14.1/include/H5Cpp.h"

class HDF5File {
private:
    hid_t file_id;
    hid_t group_id;
    std::string filename;

public:
    HDF5File(const std::string& filename) : filename(filename) {
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        group_id = H5Gcreate2(file_id, "/simData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group_id);
        group_id = H5Gcreate2(file_id, "/constants", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group_id);
        group_id = H5Gcreate2(file_id, "/energy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group_id);
        H5Fclose(file_id);
    }

    // HDF5File(const std::string& filename) : filename(filename) {
    //     H5File file(filename, H5F_ACC_TRUNC);
    //     hsize_t     dimsf[2];              // dataset dimensions
    //     dimsf[0] = NX;
    //     DataSpace dataspace( RANK, dimsf );

    // }

    void write_sim_data(const std::stringstream& buffer) {
        file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        group_id = H5Gopen2(file_id, "/simData", H5P_DEFAULT);

        // Convert stringstream to buffer
        std::string str = buffer.str();
        const char* c_str = str.c_str();

        const hsize_t size = str.size();
        

        hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        hid_t dataset_id = H5Dcreate2(group_id, "simData", H5T_STD_U8LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, c_str);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
    }

    void write_constants(const std::stringstream& buffer) {
        file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        group_id = H5Gopen2(file_id, "/constants", H5P_DEFAULT);

        // Convert stringstream to buffer
        std::string str = buffer.str();
        const char* c_str = str.c_str();

        const hsize_t size = str.size();

        hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        hid_t dataset_id = H5Dcreate2(group_id, "constants", H5T_STD_U8LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, c_str);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
    }

    void write_energy(const std::stringstream& buffer) {
        file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        group_id = H5Gopen2(file_id, "/energy", H5P_DEFAULT);

        // Convert stringstream to buffer
        std::string str = buffer.str();
        const char* c_str = str.c_str();

        const hsize_t size = str.size();

        hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        hid_t dataset_id = H5Dcreate2(group_id, "energy", H5T_STD_U8LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, c_str);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
    }

    void read_data(const std::string& group_name, const std::string& dataset_name, std::stringstream& buffer) {
        file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        group_id = H5Gopen2(file_id, group_name.c_str(), H5P_DEFAULT);
        hid_t dataset_id = H5Dopen2(group_id, dataset_name.c_str(), H5P_DEFAULT);


        hid_t dataspace_id = H5Dget_space(dataset_id);
        hsize_t dims[1];
        H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
        size_t buffer_size = dims[0];

        // Resize the stringstream buffer
        buffer.str("");
        buffer.clear();
        buffer.str(std::string(buffer_size, '\0'));

        char* buffer_ptr = const_cast<char*>(buffer.str().c_str());
        H5Dread(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer_ptr);

        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
    }

    // void close()
    // {
    //     if (file_id >= 0) {
    //         H5Fclose(file_id);
    //     }
    // }

    // ~HDF5File() {
    //     // Close the HDF5 file if it's still open
    //     if (file_id >= 0) {
    //         H5Fclose(file_id);
    //     }
    // }
};