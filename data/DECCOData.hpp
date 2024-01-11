#include "H5Cpp.h"
#include "../default_files/dust_const.hpp"
#include "../utilities/vec3.hpp"
#include <filesystem>
#include <iostream>
#include <vector>
#include <sstream>

const int bufferlines = 10;
const int num_data_types = 4;
const std::string data_types[num_data_types] = {"simData","energy","constants","timing"};
const int single_ball_widths[num_data_types] = {11,6,3,2};

std::string getDataStringFromIndex(const int data_index)
{
	if (data_index < 0 || data_index > num_data_types-1)
	{
		std::cerr<<"DECCOData ERROR: data_index '"<<data_index<<"' is out of range."<<std::endl;
		return "DECCOData ERROR";
	}
	return data_types[data_index];
}


int getDataIndexFromString(std::string data_type)
{
	for (int i = 0; i < num_data_types; ++i)
	{
		if (data_type == data_types[i])
		{
			return i;
		}
	}
	std::cerr<<"DECCOData ERROR: dataType '"<<data_type<<"' not found in class."<<std::endl;
	return -1;
}

class CSVHandler
{
public:
	CSVHandler()=default;
	~CSVHandler(){};
	CSVHandler(std::string filename) : filename(filename) 
    {
    	// initialized = true;
    }

	bool writeSimData(std::vector<double> data, int width, std::string filename)
	{
		std::cerr<<"In writeSimData: "<<filename<<std::endl;
		filename += "simData.csv"; 
		try
		{
			std::string meta = "";
			if (!std::filesystem::exists(filename))
			{
				int num_particles = width/single_ball_widths[0];
				meta = genSimDataMetaData(num_particles);	
			}

			std::ofstream simWrite;
			simWrite.open(filename, std::ofstream::app);


			simWrite << meta;

			for (int i = 0; i < data.size(); ++i)
			{
				if (i%width == width-1)
				{
					simWrite << data[i] << '\n';
				}
				else
				{
					simWrite << data[i] << ',';
				}
			}
			simWrite.close();
		}
		catch(const std::exception& e)
		{
			return 0;
		}

		return 1;
	}

	bool writeEnergy(std::vector<double> data, int width, std::string filename)
	{
		filename += "energy.csv"; 
		try
		{
			std::string meta = "";

			if (!std::filesystem::exists(filename))
			{
				int num_particles = width/single_ball_widths[1];
				meta = genEnergyMetaData();	
			}
			
			std::ofstream energyWrite;
			energyWrite.open(filename, std::ofstream::app);

			energyWrite << meta;

			for (int i = 0; i < data.size(); ++i)
			{
				if (i%width == width-1)
				{
					energyWrite << data[i] << '\n';
				}
				else
				{
					energyWrite << data[i] << ',';
				}
			}
			energyWrite.close();
		}
		catch(const std::exception& e)
		{
			return 0;
		}

		return 1;
	}

	bool writeConstants(std::vector<double> data, int width, std::string filename)
	{
		filename += "constants.csv"; 
		try
		{
			//Consts has no meta data for csv
			// std::string meta = "";

			// if (!std::filesystem::exists(filename))
			// {
			// 	int num_particles = width/single_ball_widths[2];
			// 	meta = genConstantsMetaData();	
			// }
			
			std::ofstream constsWrite;
			constsWrite.open(filename, std::ofstream::app);

			for (int i = 0; i < data.size(); ++i)
			{
				if (i%width == width-1)
				{
					constsWrite << data[i] << '\n';
				}
				else
				{
					constsWrite << data[i] << ',';
				}
			}
			constsWrite.close();
		}
		catch(const std::exception& e)
		{
			return 0;
		}

		return 1;
	}

	// bool writeTiming(std::vector<double> data, int width, std::string filename)
	// {
	// 	filename += "timing.csv"; 
	// 	try
	// 	{
	// 		std::ofstream constsWrite;
	// 		constsWrite.open(filename, std::ofstream::app);

	// 		if (!std::filesystem::exists(filename))
	// 		{
	// 			int num_particles = width/single_ball_widths[1];
	// 			constsWrite << genEnergyMetaData();	
	// 		}

	// 		for (int i = 0; i < data.size(); ++i)
	// 		{
	// 			if (i%width == width-1)
	// 			{
	// 				constsWrite << data[i] << '\n';
	// 			}
	// 			else
	// 			{
	// 				constsWrite << data[i] << ',';
	// 			}
	// 		}
	// 	}
	// 	catch{
	// 		constsWrite.close();
	// 		return 0;
	// 	}

	// 	constsWrite.close();
	// 	return 1;
	// }

    std::string genSimDataMetaData(int num_particles)
	{
		std::ostringstream meta_data;
		meta_data << "x0,y0,z0,wx0,wy0,wz0,wmag0,vx0,vy0,vz0,bound0";

		for (int Ball = 1; Ball < num_particles; Ball++)  // Start at 2nd ball because first one was just written^.
	    {
	        // std::cout<<Ball<<','<<num_particles<<std::endl;
	        std::string thisBall = std::to_string(Ball);
	        meta_data << ",x" + thisBall << ",y" + thisBall << ",z" + thisBall << ",wx" + thisBall
	                  << ",wy" + thisBall << ",wz" + thisBall << ",wmag" + thisBall << ",vx" + thisBall
	                  << ",vy" + thisBall << ",vz" + thisBall << ",bound" + thisBall;
	        // std::cout << ",x" + thisBall << ",y" + thisBall << ",z" + thisBall << ",wx" + thisBall
	        //           << ",wy" + thisBall << ",wz" + thisBall << ",wmag" + thisBall << ",vx" + thisBall
	        //           << ",vy" + thisBall << ",vz" + thisBall << ",bound" + thisBall;

	    }
	    meta_data << '\n';

		return meta_data.str();
	}

	std::string genConstantsMetaData()
	{
		//There is no metadata for constant file with csv format. THis is just
		//here to be complete and stay compatable with the previous version.
		return "";
	}

	std::string genEnergyMetaData()
	{
		return "Time,PE,KE,E,p,L";
	}

	std::string genTimingMetaData()
	{
		//There is no metadata for timing file with csv format. THis is just
		//here to be complete and stay compatable with the previous version.
		return "";
	}

private:
	std::string filename;
};

void printVec(std::vector<double> v)
{
	int i;
	for (i = 0; i < v.size()-1; i++)
	{
		std::cout<<v[i]<<", ";
	}
	std::cout<<v[i]<<std::endl;
}

class HDF5Handler {
    public:
    	HDF5Handler()=default;
    	~HDF5Handler(){};
        HDF5Handler(std::string filename,bool fixed=false) : filename(filename),fixed(fixed) 
        {
        	initialized = true;
        }

        //Make start a class variable so you dont need to keep track
        void createAppendFile(std::vector<double>& data,const std::string datasetName,size_t start=0,hsize_t maxdim=0) {
		    // H5::H5File file;
		    // H5::DataSet dataset;
		    // H5::DataSpace dataspace;

		    if(std::filesystem::exists(filename)) {
		    	if (datasetExists(filename,datasetName))
		    	{
		    		// if (fixed)
		    		// 	appendDataSet(data,datasetName,start,);
		    		// else
		    			appendDataSet(data,datasetName,start);

		    	}
		    	else
		    	{
		    		createDataSet(data,datasetName,1,maxdim);
		    	}   
		    } else {
		        createDataSet(data,datasetName,0,maxdim);
		    }
		}

		//@params datasetName is the name of the dataset you want to read
		//@params offset specifies where to start reading in the dataset.
		//			If this value is negative, the offset is applied from the 
		//			end of the dataset.
		//@param len specifies the length of data to read starting from the offset.
		//@returns the data requested as a 1d vector of doubles. If the returned vector
		//			is empty, then the read failed, or you specified len=0.	
		// If you specify a length longer than the dataset, or an offset further away 
		// from zero then the length of the dataset, the whole dataset is returned.
		std::vector<double> readFile(const std::string datasetName, hsize_t start=0, hsize_t len=0) {
		    std::vector<double> data;
		    if (std::filesystem::exists(filename)) {
		        H5::H5File file(filename, H5F_ACC_RDONLY);
		        H5::DataSet dataset = file.openDataSet(datasetName);
		        H5::DataSpace dataspace = dataset.getSpace();

		        hsize_t dims_out[1];
		        dataspace.getSimpleExtentDims(dims_out, NULL);
		        hsize_t total_size = dims_out[0];


		        if (start > total_size)
		        {
		        	std::cerr<<"DECCOData ERROR: invalid start input"<<std::endl;
		        	return data;
		        }


		        if (len > total_size-start || len == 0)
		        {
		        	len = total_size-start;
		        }

		        data.resize(len);

		        // std::cerr<<"START: "<<start<<std::endl;
		        // std::cerr<<"LEN: "<<len<<std::endl;
		        // std::cerr<<"data.size(): "<<data.size()<<std::endl;
		        // std::cerr<<"total_size: "<<total_size<<std::endl;
		        // //start should be 
		        // //len should be 3

		        hsize_t offset[1] = {start};
		        hsize_t count[1] = {len};
		        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

		        hsize_t dimsm[1] = {len};              /* memory space dimensions */
 		        H5::DataSpace memspace(1, dimsm);

		        dataset.read(&data[0], H5::PredType::NATIVE_DOUBLE,memspace,dataspace);

		        dataspace.close();
		        memspace.close();
			    dataset.close();
			    file.close();
		    } else {
		        std::cerr << "File '" << filename << "' does not exist." << std::endl;
		    }

		    return data;
		}

		//@params datasetName is the name of the dataset you want to read
		//@params offset specifies where to start reading in the dataset.
		//			If this value is negative, the offset is applied from the 
		//			end of the dataset.
		//@param len specifies the length of data to read starting from the offset.
		//@param file is the .h5 file you want to read
		//@returns the data requested as a 1d vector of doubles. If the returned vector
		//			is empty, then the read failed, or you specified len=0.	
		// If you specify a length longer than the dataset, or an offset further away 
		// from zero then the length of the dataset, the whole dataset is returned.
		static std::vector<double> static_readFile(const std::string datasetName, hsize_t start=0, hsize_t len=0,bool neg_offset=false,std::string readfile="") {
		    std::vector<double> data;
		    if (std::filesystem::exists(readfile)) {
		        H5::H5File file(readfile, H5F_ACC_RDONLY);
		        H5::DataSet dataset = file.openDataSet(datasetName);
		        H5::DataSpace dataspace = dataset.getSpace();

		        hsize_t dims_out[1];
		        dataspace.getSimpleExtentDims(dims_out, NULL);
		        hsize_t total_size = dims_out[0];


		        if (start > total_size)
		        {
		        	std::cerr<<"DECCOData ERROR: invalid start input"<<std::endl;
		        	return data;
		        }

		        if (neg_offset)
		        {
		        	start = total_size - start;
		        }

		        if (len > total_size-start || len == 0)
		        {
		        	len = total_size-start;
		        }

		        data.resize(len);

		        // std::cerr<<"START: "<<start<<std::endl;
		        // std::cerr<<"LEN: "<<len<<std::endl;
		        // std::cerr<<"data.size(): "<<data.size()<<std::endl;
		        // std::cerr<<"total_size: "<<total_size<<std::endl;
		        // //start should be 
		        // //len should be 3

		        hsize_t offset[1] = {start};
		        hsize_t count[1] = {len};
		        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

		        hsize_t dimsm[1] = {len};              /* memory space dimensions */
 		        H5::DataSpace memspace(1, dimsm);

		        dataset.read(&data[0], H5::PredType::NATIVE_DOUBLE,memspace,dataspace);

		        dataspace.close();
		        memspace.close();
			    dataset.close();
			    file.close();
		    } else {
		        std::cerr << "File '" << readfile << "' does not exist." << std::endl;
		    }

		    return data;
		}

		

		static hsize_t get_data_length(std::string readfile,std::string datasetName)
		{
			hsize_t total_size;
			if (std::filesystem::exists(readfile)) {
		        H5::H5File file(readfile, H5F_ACC_RDONLY);
		        H5::DataSet dataset = file.openDataSet(datasetName);
		        H5::DataSpace dataspace = dataset.getSpace();

		        hsize_t dims_out[1];
		        dataspace.getSimpleExtentDims(dims_out, NULL);
		        total_size = dims_out[0];
		    } else {
		        std::cerr << "File '" << readfile << "' does not exist." << std::endl;
		    }
		    return total_size;
		}

		int readWrites(const std::string path, const std::string f)
		{
			std::string file_read = path + f;
			// std::cerr<<"readWrites file: "<<file_read<<std::endl;
			const std::string datasetName("writes");
	        H5::H5File file;
	        H5::DataSet dataset;
	        int value = 0; // uninitialized value

	        if(std::filesystem::exists(file_read)) {
	            file = H5::H5File(file_read, H5F_ACC_RDWR);
	        } else {
	            std::cerr << "File " << file_read << " doesn't exist." << std::endl;
	            exit(-1);
	        }

	        if (datasetExists(file_read,"writes"))
	        {

	            // Attempt to open the dataset
	            dataset = file.openDataSet(datasetName);

	            // Read the current value
	            dataset.read(&value, H5::PredType::NATIVE_INT);
	            // std::cerr<<"value in readWrites"
	        } 
	        else
	        {
	            // If the dataset doesn't exist, it hasn't been written to yet
	            value = 0;
	        }

	        // Write the (incremented or default) value back to the dataset
	        // dataset.write(&value, H5::PredType::NATIVE_INT);

	        // Close resources
	        dataset.close();
	        file.close();

	        return value;
		}

		void addWrites(int additional_writes) 
		{
			const std::string datasetName("writes");
	        H5::H5File file;
	        H5::DataSet dataset;
	        int value = additional_writes; // Default value

	        if(std::filesystem::exists(filename)) {
	            file = H5::H5File(filename, H5F_ACC_RDWR);
	        } else {
	            std::cerr << "File " << filename << " doesn't exist." << std::endl;
	            exit(-1);
	        }

	        if (datasetExists(filename,"writes"))
	        {
	            // Attempt to open the dataset
	            dataset = file.openDataSet(datasetName);

	            // Read the current value
	            dataset.read(&value, H5::PredType::NATIVE_INT);

	            // Increment the value
	            value+=additional_writes;
	            // std::cerr<<"WRITES IN INCREMENT WRITES: "<<value<<std::endl;
	        } 
	        else
	        {
	            // If the dataset doesn't exist, create it with a default value of 1
	            H5::DataSpace scalarSpace(H5S_SCALAR);
	            dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_INT, scalarSpace);
	            // std::cerr << "Dataset " << datasetName << " created with initial value: " << value << std::endl;
	        }

	        // Write the (incremented or default) value back to the dataset
	        dataset.write(&value, H5::PredType::NATIVE_INT);

	        // Close resources
	        dataset.close();
	        file.close();
	    }

		static int get_num_particles(std::string path, std::string file)
		{
			hsize_t total_size = get_data_length(path+file,"constants");

		    return (total_size*1.0)/(single_ball_widths[getDataIndexFromString("constants")]*1.0);
		}

		static void loadConsts(std::string path,std::string file,double *R,double *m,double *moi)
		{
			std::vector<double> constdata = static_readFile("constants",0,0,false,path+file);
		
			int j = -1;
			for (int i = 0; i < constdata.size(); i++)
			{
				if (i % 3 == 0)
				{
					j++;
					R[j] = constdata[i];
				}
				else if (i % 3 == 1)
				{
					m[j] = constdata[i];
				}
				else if (i % 3 == 2)
				{
					moi[j] = constdata[i];
				}
			}

		}

        // std::vector<double> readFile(const std::string datasetName) {
        //     std::vector<double> data;
        //     if(std::filesystem::exists(filename)) {
        //         H5::H5File file(filename, H5F_ACC_RDONLY);
        //         H5::DataSet dataset = file.openDataSet(datasetName);
        //         H5::DataSpace dataspace = dataset.getSpace();
                
        //         hsize_t dims_out[1];
        //         dataspace.getSimpleExtentDims(dims_out, NULL);
        //         data.resize(dims_out[0]);
                
        //         dataset.read(&data[0], H5::PredType::NATIVE_DOUBLE);
        //         file.close();
        //     }else{
        //     	std::cerr<<"File '"<<filename<<"' does not exist."<<std::endl;
        //     }
        //     return data;
        // }

        void attachMetadataToDataset(const std::string& metadata, const std::string& datasetName) 
        {
		    // Initialize the HDF5 library
        	if(std::filesystem::exists(filename)) {
			    H5::H5File file(filename, H5F_ACC_RDWR);
	    		H5::DataSet dataset;
			    // Open the specified dataset
			    if (datasetExists(filename,datasetName)){
			    	dataset = file.openDataSet(datasetName);
			    	// Check if the attribute's dataspace exists
				    if (attributeExists(dataset, "metadata")) {
				        // std::cerr << "DECCOHDF5 Warning: Attribute 'metadata' already exists for the dataset." << std::endl;
				        dataset.close();
				        file.close();
				        return;
				    }
			    }else{
			    	std::vector<double> dummy;
			    	dataset = createDataSet(dummy,datasetName);
			    }

			    // Create a string data type for the attribute
			    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);

			    

			    // Create a dataspace for the attribute
			    H5::DataSpace attrSpace = H5::DataSpace(H5S_SCALAR);

			    // Create the attribute and write the metadata
			    H5::Attribute metadataAttr = dataset.createAttribute("metadata", strType, attrSpace);
			    metadataAttr.write(strType, metadata);

			    // Close resources
			    metadataAttr.close();
			    attrSpace.close();
			    dataset.close();
			    file.close();
		    }else{
            	std::cerr<<"File '"<<filename<<"' does not exist."<<std::endl;
            }
            return;
		}

        void attachSimMetadataToDataset(const std::string& metadata, const std::string& metadataName, \
        							const std::string& datasetName) 
        {
		    // Initialize the HDF5 library
        	if(std::filesystem::exists(filename)) {
			    H5::H5File file(filename, H5F_ACC_RDWR);
	    		H5::DataSet dataset;
			    // Open the specified dataset
			    if (datasetExists(filename,datasetName)){
			    	dataset = file.openDataSet(datasetName);
			    	// Check if the attribute's dataspace exists
				    if (attributeExists(dataset, metadataName)) {
				        // std::cerr << "DECCOHDF5 Warning: Attribute 'metadata' already exists for the dataset." << std::endl;
				        dataset.close();
				        file.close();
				        return;
				    }
			    }else{
			    	std::vector<double> dummy;
			    	dataset = createDataSet(dummy,datasetName);
			    }

			    // Create a string data type for the attribute
			    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);

			    

			    // Create a dataspace for the attribute
			    H5::DataSpace attrSpace = H5::DataSpace(H5S_SCALAR);

			    // Create the attribute and write the metadata
			    H5::Attribute metadataAttr = dataset.createAttribute(metadataName, strType, attrSpace);
			    metadataAttr.write(strType, metadata);

			    // Close resources
			    metadataAttr.close();
			    attrSpace.close();
			    dataset.close();
			    file.close();
		    }else{
            	std::cerr<<"File '"<<filename<<"' does not exist."<<std::endl;
            }
            return;
		}

		//Non static version for reading from initialized class instance
		std::string readMetadataFromDataset(const std::string& datasetName,const std::string& metadataName="metadata") 
		{
		    std::string metadata;

		    // Initialize the HDF5 library
		    if(std::filesystem::exists(filename)) {
		    	H5::DataSet dataset;
                H5::H5File file(filename, H5F_ACC_RDONLY);
			    // Open the specified dataset
			    if (datasetExists(filename,datasetName)){
			    	dataset = file.openDataSet(datasetName);
			    }else{
			        std::cerr << "dataset '"<<datasetName<<"' does not exist for file '"<<filename<<"' ." << std::endl;
			    	return "ERROR RETRIEVING METADATA";
			    }

			    // Check if the attribute's dataspace exists
			    if (!attributeExists(dataset, metadataName)) {
			        std::cerr << "Attribute 'metadata' does not exist for dataset '"<<datasetName<<"' in file '"<<filename<<"' ." << std::endl;
			        dataset.close();
			        file.close();
			    	return "ERROR RETRIEVING METADATA";
			    }

			    // Open the metadata attribute
			    H5::Attribute metadataAttr = dataset.openAttribute(metadataName);

			    // Read the metadata attribute
			    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
			    metadataAttr.read(strType, metadata);
			    
		    	// Close resources
			    metadataAttr.close();
			    dataset.close();
			    file.close();
            }else{
            	std::cerr<<"File '"<<filename<<"' does not exist."<<std::endl;
            }



		    return metadata;
		}

		//static version for reading when no DECCOData class instance exists
		static std::string readMetadataFromDataset(const std::string& datasetName,const std::string& metafile,const std::string& metadataName="metadata") 
		{
		    std::string metadata;

		    // Initialize the HDF5 library
		    if(std::filesystem::exists(metafile)) {
		    	H5::DataSet dataset;
                H5::H5File file(metafile, H5F_ACC_RDONLY);
			    // Open the specified dataset
			    if (datasetExists(metafile,datasetName)){
			    	dataset = file.openDataSet(datasetName);
			    }else{
			        std::cerr << "dataset '"<<datasetName<<"' does not exist for file '"<<metafile<<"' ." << std::endl;
			    	return "ERROR RETRIEVING METADATA";
			    }

			    // Check if the attribute's dataspace exists
			    if (!attributeExists(dataset, metadataName)) {
			        std::cerr << "Attribute "<<metadataName<<" does not exist for dataset '"<<datasetName<<"' in file '"<<metafile<<"' ." << std::endl;
			        dataset.close();
			        file.close();
			    	return "ERROR RETRIEVING METADATA";
			    }

			    // Open the metadata attribute
			    H5::Attribute metadataAttr = dataset.openAttribute(metadataName);

			    // Read the metadata attribute
			    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
			    metadataAttr.read(strType, metadata);
			    
		    	// Close resources
			    metadataAttr.close();
			    dataset.close();
			    file.close();
            }else{
            	std::cerr<<"File '"<<metafile<<"' does not exist."<<std::endl;
            }



		    return metadata;
		}

		bool isInitialized()
		{
			return initialized;
		}

        static bool datasetExists(const std::string& f, const std::string& datasetName)
		{
		    // Open the HDF5 file
		    H5::H5File file(f, H5F_ACC_RDONLY);

		    // Check if the dataset exists
		    bool exists = H5Lexists(file.getId(), datasetName.c_str(), H5P_DEFAULT) > 0;
		    file.close();
		    return exists;
		}

        static bool attributeExists(const H5::H5Object& object, const std::string& attributeName) 
        {
    		return object.attrExists(attributeName);
    	}

    	std::string genSimDataMetaData(int width)
		{
			std::string meta_data = "";
			meta_data = "Columns: posx[ball],posy[ball],posz[ball],wx[ball],wy[ball],wz[ball],wtot[ball],velx[ball],vely[ball],velz[ball],bound[ball]\n";
			meta_data += "rows: writes\n";
			meta_data += "row width: " + std::to_string(width);

			return meta_data;
		}

		std::string genConstantsMetaData(int width)
		{
			std::string meta_data = "Columns: radius, mass, moment of inertia\n";
			meta_data += "rows: balls\n";	
			meta_data += "row width: " + std::to_string(width);

			return meta_data;
		}

		std::string genEnergyMetaData(int width)
		{
			std::string meta_data = "Columns: time, PE, KE, Etot, p, L\n";
			meta_data += "rows: writes\n";	
			meta_data += "row width: " + std::to_string(width);
			return meta_data;	
		}

		std::string genTimingMetaData(int width)
		{
			std::string meta_data = "Columns: number of balls, time spent in sim_one_step\n";
			meta_data += "rows: sims\n";	
			meta_data += "row width: " + std::to_string(width);
			return meta_data;	
		}

    private:
        std::string filename;
        bool fixed;
        bool initialized = false;
        




		// __attribute__((optimize("O0")))
		H5::DataSet createDataSet(std::vector<double>& data,const std::string datasetName,int fileExists=0,int max_size=-1)
		{
			H5::H5File file;
		    H5::DataSet dataset;
		    H5::DataSpace dataspace;
	        // std::cout<<"Filename in create: "<<filename<<std::endl;
	        // std::cout<<"datasetName: "<<datasetName<<std::endl;
	        // std::cout<<"data: "<<data[0]<<std::endl;
	        // int i;
			// for (i = 0; i < data.size()-1; i++)
			// {
			// 	std::cout<<data[i]<<", ";
			// }
			// std::cout<<data[i]<<std::endl;
	        // printVec(data);

		    if (fileExists == 1)//fileExists
		    {
				file = H5::H5File(filename, H5F_ACC_RDWR);
		    }
		    else
		    {
				file = H5::H5File(filename, H5F_ACC_TRUNC);
		    }
	        hsize_t maxdims[1];
	        hsize_t dims[1] = {data.size()};
	        if (fixed)
	        {
	        	maxdims[0] = max_size; // Set maximum dimensions to max_size
	        	// std::cout<<"maxsize: "<<max_size<<std::endl;
	        	dataspace = H5::DataSpace(1, maxdims);
		        dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE,dataspace);
	        	dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
	        }
	        else
	        {
	        	maxdims[0] = H5S_UNLIMITED; // Set maximum dimensions to unlimited
	        	dataspace = H5::DataSpace(1, dims, maxdims);
		        H5::DSetCreatPropList plist;
		        hsize_t chunk_dims[1] = {std::min((hsize_t)1000, (hsize_t)data.size())}; // Adjust chunk size as needed
		        plist.setChunk(1, chunk_dims);
		        dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace, plist);
	        	dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
	        }

	        dataspace.close();
	        dataset.close();
	        file.close();
	        return dataset;
		}

		void appendDataSet(std::vector<double>& data,const std::string datasetName,size_t start=0,int max_size=-1)
		{
			H5::H5File file;
		    H5::DataSet dataset;
		    H5::DataSpace dataspace;
		    if(std::filesystem::exists(filename)) 
		    {
				file = H5::H5File(filename, H5F_ACC_RDWR);
			}
			else
			{
				std::cerr<<"File "<<filename<<" doesn't exist."<<std::endl;
				exit(-1);
			}
	        // std::cout<<"datasetName: "<<datasetName<<std::endl;
			// int i;
			// for (i = 0; i < data.size()-1; i++)
			// {
			// 	std::cout<<data[i]<<", ";
			// }
			// std::cout<<data[i]<<std::endl;
	        // std::cout<<"Filename in append: "<<filename<<std::endl;
	        // dataset = file.openDataSet(datasetName);
	        // dataspace = dataset.getSpace();

	        try {
			    dataset = file.openDataSet(datasetName);
			} catch (const H5::Exception& error) {
			    std::cerr<<"H5 ERROR: "<<error.getDetailMsg()<<std::endl;
			}

			try {
			    dataspace = dataset.getSpace();
			} catch (const H5::Exception& error) {
			    std::cerr<<"H5 ERROR: "<<error.getDetailMsg()<<std::endl;
			}

	        // Get current size of the dataset
	        hsize_t dims_out[1];
	        if (fixed)
	        {
	        	dims_out[0] = max_size;
	        }
	        else
	        {
	        	dataspace.getSimpleExtentDims(dims_out, NULL);
	        }

	        // Extend the dataset to hold the new data
	        hsize_t size[1] = {dims_out[0] + data.size()};
	        if (!fixed)
	        {
	        	dataset.extend(size);
	        }

	        // Select the extended part of the dataset
	        dataspace = dataset.getSpace();
	        hsize_t offset[1];
	        if (fixed)
	        {

	        	offset[0] = {start};
	        }
	        else
	        {

	        	offset[0] = {dims_out[0]};
	        }
	        hsize_t dimsextend[1] = {data.size()};
	        dataspace.selectHyperslab(H5S_SELECT_SET, dimsextend, offset);

	        // Write the data to the extended part of the dataset
	        H5::DataSpace memspace(1, dimsextend);
	        dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
	        
	        dataspace.close();
	        dataset.close();
	        memspace.close();
	        file.close();
		}
};




class DECCOData
{
public:
	

	/*
	Mandatory class input
		storage_method: either "h5"/"hdf5" for an hdf5 file format or "csv" for a csv format. Will be whatever the file
						extension is in "fname" variable
	Optional class input
		num_particles : the number of particles in the DECCO simulation, only needed for a fixed hdf5 storage
		writes        : the number of writes that will happen in the DECCO simulation, only needed for a fixed hdf5 storage
		steps   	  : the number of sims for writing out timing
	*/
	DECCOData(std::string fname, int num_particles, int writes=-1, int steps=-1) : 
		filename(fname), 
		num_particles(num_particles),
		writes(writes),
		steps(steps)
	{
		//If writes is set we know we want a fixed hdf5 file storage
		if (writes > 0 || steps > 0)
		{
			fixed = true;
		}
		else
		{
			fixed = false;
		}

		int dot_index = filename.find('.');
		storage_type = filename.substr(dot_index+1,filename.length()-dot_index);
		//Transform storage_type to be not case sensitive
		transform(storage_type.begin(), storage_type.end(), storage_type.begin(), ::tolower);
		if (storage_type == "h5" || storage_type == "hdf5")
		{
			h5data = true;
			csvdata = false;
		}
		else if (storage_type == "csv")
		{
			csvdata = true;
			h5data = false;
			filename = filename.substr(0,filename.length()-4);
		}
		else
		{
			std::cerr<<"DECCOData ERROR: storage_type '"<<storage_type<<"' not available."<<std::endl;
			exit(-1);
		}


		//If user specified number of writes but a storage_type other than hdf5 then default to hdf5 and warn user
		// if (fixed && not h5data)
		// {
		// 	std::cerr<<"DECCOData warning: specified a non-negative number of writes (which indicates a fixed size hdf5 storage) and a storage_type other than hdf5. Defaulting to a fixed size hdf5 storage_type."<<std::endl;
		// }

	}

	//Constructor for restarting an existing job
	// DECCOData(std::string filename) : 
	// 	filename(filename)
	// {
		
	// }


	//copy constructor to handle the data handlers
	DECCOData& operator=(const DECCOData& rhs)
	{
		C=rhs.C;
		H=rhs.H;
		return *this;
	}

	//Write the data given with the method specified during class construction
	//This includes taking care of headers for csv and metadata for hdf5
	//@param data_type is one of "simData", "constants", "energy", "timing",
	//@param data is the data to be written. Can only be an std::stringstream 
	//		(for csv data_type) or std::vector<double> (for hdf5 data_type)
	//@returns true if write succeeded, false otherwise
	bool Write(std::vector<double> &data, std::string data_type,int add_writes=0)
	{
		bool retVal;
		if (h5data)
		{
			retVal = writeH5(data,data_type);
			if (add_writes > 0)
			{
				H.addWrites(add_writes);
			}
		}
		else if (csvdata)
		{
			std::cerr<<"IN WRITE: "<<filename<<std::endl;
			retVal = writeCSV(data,data_type,filename);
		}
		return retVal;
	}


	std::string ReadMetaData(int data_index) 
	{
		//if data_index is less than zero a bad data_type was input and write doesn't happen
		if (data_index < 0)
		{
			return "ERROR RETRIEVING METADATA";
		}
		std::string datasetName = getDataStringFromIndex(data_index);
		std::string readMetaData;
		if (h5data)
		{
			//Has the HDF5 handler been initiated yet?
			if (!H.isInitialized())
			{
				H = HDF5Handler(filename,fixed);
			}

			readMetaData = H.readMetadataFromDataset(datasetName);
		}
		else if (csvdata)
		{
			return "ERROR";
			// //Has the csv handler been initiated yet?
			// if (!C.isInitialized())
			// {
			// 	C = CSVHandler(filename,fixed);
			// }

			// readMetaData = C.readMetadataFromDataset(datasetName);
		}

		return readMetaData;
	}

	std::string ReadMetaData(std::string data_type) 
	{
		int data_index = getDataIndexFromString(data_type);
		std::string readMetaData = ReadMetaData(data_index);
		return readMetaData;
	}

	std::vector<double> Read(std::string data_type, bool all=true, int line=-1,std::string file="") 
	{
		if (file == "")
		{
			file = filename;
		}
		

		std::vector<double> data_read;
		if (file.substr(file.size()-3,file.size()) == ".h5")
		{

			//Has the HDF5 handler been initiated yet?
			if (!H.isInitialized())
			{
				H = HDF5Handler(file,fixed);
			}

			if (all)
			{
				data_read = H.readFile(data_type); 
			}
			else
			{
				int data_index = getDataIndexFromString(data_type);
				int start = 0;
				//if data_index is less than zero a bad data_type was input and write doesn't happen
				if (data_index < 0)
				{
					return data_read;
				}

				if (line < 0)
				{
					start = (writes+line)*widths[data_index];
				}
				else
				{
					start = (line)*widths[data_index];
					// start = line*widths[data_index];
				}

				data_read = H.readFile(data_type,start,widths[data_index]);
				// for (int i = 0; i < data_read.size(); i++)
				// {
				// 	std::cerr<<data_read[i]<<", ";
				// }
			}
		}
		else if (file.substr(file.size()-4,file.size()) == ".csv")
		{
			std::cerr<<"ERROR: csv file type not yet readable by DECCOData."<<std::endl;
			exit(-1);
		}
		return data_read;
	}

	void WriteMeta(const std::string& metadata, const std::string& metadataName, \
        			const std::string& datasetName) 
    {
    	H.attachSimMetadataToDataset(metadata,metadataName,datasetName);
    }


    //This function returns the status of how many writes have happened. It also sets the written_so_far values
    //	if there is more than zero and less then max_writes writes.
    //@return 0 if there is no writes so far (I don't think this should happen but if it does, more stuff needs to happen).
    //@return int >0 for how many writes there have been.
    //@return -1 if there are writes and the sim is already finished. 
    int setWrittenSoFar(const std::string path, const std::string file)
    {
    	hsize_t energy_size = H.get_data_length(path+file,"energy");
    	int energy_width = getSingleWidth("energy");
    	int writes_so_far=H.readWrites(path,file);
    	writes = writes_so_far;

    	if (writes_so_far*energy_width == energy_size)
    	{
    		return -1;
    	}
    	else if (writes_so_far == 0)
    	{
    		return 0;
    	}
    	
    	// std::cerr<<writes_so_far<<std::endl;
    	for (int i = 0; i < num_data_types; i++)
    	{
    		written_so_far[i] = writes_so_far*widths[i];
    		// std::cerr<<"widths: "<<widths[i]<<std::endl;
    		// std::cerr<<"written so far: "<<written_so_far[i]<<std::endl;
    	}

    	return writes_so_far;
    }


    void loadSimData(const std::string path, const std::string file,vec3 *pos,vec3 *w,vec3 *vel)
    {
    	std::vector<double> out = Read("simData",false,-1,path+file);
    	// printVec(out);
    	
    	for (int i = 0; i < num_particles; ++i)
    	{
    		int out_ind = i*getSingleWidth("simData");
    		pos[i].x = out[out_ind];
    		pos[i].y = out[out_ind+1];
    		pos[i].z = out[out_ind+2];

    		w[i].x = out[out_ind+3];
    		w[i].y = out[out_ind+4];
    		w[i].z = out[out_ind+5];

    		vel[i].x = out[out_ind+7];
    		vel[i].y = out[out_ind+8];
    		vel[i].z = out[out_ind+9];
    	}
    }

	int getNumTypes()
	{
		return num_data_types;
	}

	int getWidth(std::string data_type)
	{
		return widths[getDataIndexFromString(data_type)];
	}

	int getSingleWidth(std::string data_type)
	{
		return single_ball_widths[getDataIndexFromString(data_type)];
	}


	std::string genMetaData(int data_index)
	{
		if (data_index == 0)//simData
		{
			return genSimDataMetaData();
		}
		else if (data_index == 1)//energy
		{
			return genEnergyMetaData();
		}
		else if (data_index == 2)//constants
		{
			return genConstantsMetaData();
		}
		else if (data_index == 3)//timing
		{
			return genTimingMetaData();
		}
		std::cerr<<"DECCOData ERROR: data_index '"<<data_index<<"' is out of range."<<std::endl;
		return "DECCOData ERROR";
	}
	

	~DECCOData(){};
private:
	
	std::string storage_type;
	std::string filename;
	int num_particles;
	int writes;
	int steps;
	int fixed_width = -1;
	bool fixed;
	bool h5data;
	bool csvdata;
	HDF5Handler H; 
	CSVHandler C;


    //data types so everyone can see them
	
    int widths[num_data_types] = {  single_ball_widths[0]*num_particles,\
    								single_ball_widths[1],\
    								single_ball_widths[2],\
    								single_ball_widths[3]};
    int max_size[num_data_types] = {widths[0]*writes,\
    								widths[1]*writes,\
    								widths[2]*num_particles,\
    								widths[3]*steps};
    int written_so_far[num_data_types] = {0,0,0,0};

    //Write the data given with the csv method.
	//This includes taking care of headers for csv and metadata for hdf5
	//@param data_type is one of "simData", "constants", "energy", "timing",
	//@param data is the data to be written.  
	//@param filename is the absolute path to the base save file (everything except which data it is and the file extension)
	//@returns true if write succeeded, false otherwise
	bool writeCSV(std::vector<double> &data, std::string data_type, std::string filename)
	{
		std::cerr<<"FILE in writeCSV: "<<filename<<std::endl;
		if (getDataIndexFromString(data_type) == 0) //simData
		{
			return C.writeSimData(data,widths[0],filename);
		}
		else if (getDataIndexFromString(data_type) == 1) //energy
		{
			return C.writeEnergy(data,widths[1],filename);
		}
		else if (getDataIndexFromString(data_type) == 2) //consts
		{
			return C.writeConstants(data,widths[2],filename);
		}
		// else if (getDataIndexFromString(data_type) == 1) //timing
		// {
		// 	return C.writeTiming(data,filename);
		// }

		return 0;
	}

	//Write the data given with the hdf5 method.
	//This includes taking care of headers for csv and metadata for hdf5
	//@param data_type is one of "simData", "constants", "energy", "timing",
	//@param data is the data to be written.  
	//@returns true if write succeeded, false otherwise
	bool writeH5(std::vector<double> &data, std::string data_type)
	{
		int data_index = getDataIndexFromString(data_type);
		//if data_index is less than zero a bad data_type was input and write doesn't happen
		if (data_index < 0)
		{
			return 0;
		}

		//Has the HDF5 handler been initiated yet?
		if (!H.isInitialized())
		{
			H = HDF5Handler(filename,fixed);
		}		

		if (fixed)
		{
			H.createAppendFile(data,data_type,written_so_far[data_index],max_size[data_index]);
			written_so_far[data_index] += data.size();
			// std::cout<<data_type<<": "<<written_so_far[data_index]<<" / "<<max_size[data_index]<<std::endl;
		}
		else
		{
			H.createAppendFile(data,data_type);
		}

		H.attachMetadataToDataset(genMetaData(data_index),data_type);
		
		
		return 1;
	}

    


	std::string genSimDataMetaData()
	{
		std::string meta_data = "";
		if (h5data)
		{
			meta_data = H.genSimDataMetaData(widths[getDataIndexFromString("simData")]);
		}
		else if (csvdata)
		{
			meta_data = C.genSimDataMetaData(num_particles);
		}
		return meta_data;
	}

	std::string genConstantsMetaData()
	{
		std::string meta_data = "";
		if (h5data)
		{
			meta_data = H.genConstantsMetaData(widths[getDataIndexFromString("constants")]);
		}
		else if (csvdata)
		{
			meta_data = C.genConstantsMetaData();
		}
		return meta_data;
	}

	std::string genEnergyMetaData()
	{
		std::string meta_data = "";
		if (h5data)
		{
			meta_data = H.genEnergyMetaData(widths[getDataIndexFromString("energy")]);
		}
		else if (csvdata)
		{
			meta_data = C.genEnergyMetaData();
		}
		return meta_data;
	}

	std::string genTimingMetaData()
	{
		std::string meta_data = "";
		if (h5data)
		{
			meta_data = H.genTimingMetaData(widths[getDataIndexFromString("timing")]);
		}
		else if (csvdata)
		{
			meta_data = C.genTimingMetaData();
		}
		return meta_data;
	}
};