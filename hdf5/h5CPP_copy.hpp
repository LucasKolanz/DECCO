#include "H5Cpp.h"
#include <filesystem>
#include <iostream>
#include <vector>

class HDF5Handler {
    public:
        HDF5Handler(std::string filename) : filename(filename) {}

        void createAppendFile(std::vector<double>& data) {
		    H5::H5File file;
		    H5::DataSet dataset;
		    H5::DataSpace dataspace;

		    if(std::filesystem::exists(filename)) {
		        file = H5::H5File(filename, H5F_ACC_RDWR);
		        dataset = file.openDataSet("Dataset");
		        dataspace = dataset.getSpace();

		        // Get current size of the dataset
		        hsize_t dims_out[1];
		        dataspace.getSimpleExtentDims(dims_out, NULL);

		        // Extend the dataset to hold the new data
		        hsize_t size[1] = {dims_out[0] + data.size()};
		        dataset.extend(size);

		        // Select the extended part of the dataset
		        dataspace = dataset.getSpace();
		        hsize_t offset[1] = {dims_out[0]};
		        hsize_t dimsextend[1] = {data.size()};
		        dataspace.selectHyperslab(H5S_SELECT_SET, dimsextend, offset);

		        // Write the data to the extended part of the dataset
		        H5::DataSpace memspace(1, dimsextend);
		        dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
		    } else {
		        file = H5::H5File(filename, H5F_ACC_TRUNC);

		        hsize_t dims[1] = {data.size()};
		        hsize_t maxdims[1] = {H5S_UNLIMITED}; // Set maximum dimensions to unlimited
		        dataspace = H5::DataSpace(1, dims, maxdims);

		        H5::DSetCreatPropList plist;
		        hsize_t chunk_dims[1] = {std::min((hsize_t)1000, data.size())}; // Adjust chunk size as needed
		        plist.setChunk(1, chunk_dims);
		        dataset = file.createDataSet("Dataset", H5::PredType::NATIVE_DOUBLE, dataspace, plist);

		        dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
		    }
		}



        std::vector<double> readFile() {
            std::vector<double> data;
            if(std::filesystem::exists(filename)) {
                H5::H5File file(filename, H5F_ACC_RDONLY);
                H5::DataSet dataset = file.openDataSet("Dataset");
                H5::DataSpace dataspace = dataset.getSpace();
                
                hsize_t dims_out[1];
                dataspace.getSimpleExtentDims(dims_out, NULL);
                data.resize(dims_out[0]);
                
                dataset.read(&data[0], H5::PredType::NATIVE_DOUBLE);
            }
            return data;
        }

    private:
        std::string filename;
};