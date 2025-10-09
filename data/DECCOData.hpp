#pragma once
#ifndef DECCODATA_HPP
#define DECCODATA_HPP

#include <vector>

#include "../utilities/vec3.hpp"


#ifdef HDF5_ENABLE
	#include "H5Cpp.h"
#endif

// #ifdef MPI_ENABLE
// 	#include "mpi.h"
// #endif


#define ERR_RET_MET "ERROR RETRIEVING METADATA"

const int bufferlines = 10;
const int num_data_types = 4;
const std::string data_types[num_data_types] = {"simData","energy","constants","timing"};
const int single_ball_widths[num_data_types] = {11,6,3,2};

std::string getDataStringFromIndex(const int data_index);
int getDataIndexFromString(std::string data_type);
void printVec(std::vector<double> v);


class CSVHandler
{
public:
	CSVHandler();
	~CSVHandler();
	CSVHandler(std::string filename);

	bool writeSimData(std::vector<double> data, int width, std::string filename);

	bool writeEnergy(std::vector<double> data, int width, std::string filename);

	bool writeConstants(std::vector<double> data, int width, std::string filename);

	// bool writeTiming(std::vector<double> data, int width, std::string filename);

    std::string genSimDataMetaData(int num_particles);

	std::string genConstantsMetaData();

	std::string genEnergyMetaData();
	std::string genTimingMetaData();
	bool deleteData();
	
private:
	std::string filename;
};


#ifdef HDF5_ENABLE
class HDF5Handler {
    public:
    	HDF5Handler();
    	~HDF5Handler();
        HDF5Handler(std::string filename,bool fixed);

        //Make start a class variable so you dont need to keep track
        void createAppendFile(std::vector<double>& data,const std::string datasetName,
        						size_t start,hsize_t maxdim);
		//@params datasetName is the name of the dataset you want to read
		//@params offset specifies where to start reading in the dataset.
		//			If this value is negative, the offset is applied from the 
		//			end of the dataset.
		//@param len specifies the length of data to read starting from the offset.
		//@returns the data requested as a 1d vector of doubles. If the returned vector
		//			is empty, then the read failed, or you specified len=0.	
		// If you specify a length longer than the dataset, or an offset further away 
		// from zero then the length of the dataset, the whole dataset is returned.
		std::vector<double> readFile(const std::string datasetName, hsize_t start, hsize_t len);

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
		static std::vector<double> static_readFile(const std::string datasetName, hsize_t start=0, 
												hsize_t len=0,bool neg_offset=false,std::string readfile="");
		static void loadh5SimData(const std::string path, const std::string file,vec3 *pos,vec3 *w,vec3 *vel);

		static hsize_t get_data_length(std::string readfile,std::string datasetName);
		int readWrites(const std::string path, const std::string f);

		void addWrites(int additional_writes);

	    static bool sim_finished(std::string path, std::string file);

		static int get_num_particles(std::string path, std::string file);

		static void loadConsts(std::string path,std::string file,double *R,double *m,double *moi);

        // std::vector<double> readFile(const std::string datasetName);

        void attachMetadataToDataset(const std::string& metadata, const std::string& datasetName);

        void attachSimMetadataToDataset(const std::string& metadata, const std::string& metadataName, \
        							const std::string& datasetName);

		//Non static version for reading from initialized class instance
		std::string readMetadataFromDataset(const std::string& datasetName,const std::string& metadataName);

		//static version for reading when no DECCOData class instance exists
		static std::string readMetadataFromDataset(const std::string& datasetName,const std::string& metafile,
													const std::string& metadataName);
		void appendDataSet(std::vector<double>& data,const std::string datasetName,size_t start=0,int max_size=-1);

		bool isInitialized();

        static bool datasetExists(const std::string& f, const std::string& datasetName);

        static bool attributeExists(const H5::H5Object& object, const std::string& attributeName);
    	std::string genSimDataMetaData(int width);
		std::string genConstantsMetaData(int width);
		std::string genEnergyMetaData(int width);
		std::string genTimingMetaData(int width);
		bool deleteData();


		// __attribute__((optimize("O0")))
		H5::DataSet createDataSet(std::vector<double>& data,const std::string datasetName,
										int fileExists=0,int max_size=-1);

private:
	    std::string filename;
        bool fixed;
        bool initialized = false;
};
#endif




class DECCOData
{
public:
	
	~DECCOData();
	/*
	Mandatory class input
		storage_method: either "h5"/"hdf5" for an hdf5 file format or "csv" for a csv format. Will be whatever the file
						extension is in "fname" variable
	Optional class input
		num_particles : the number of particles in the DECCO simulation, only needed for a fixed hdf5 storage
		writes        : the number of writes that will happen in the DECCO simulation, only needed for a fixed hdf5 storage
		steps   	  : the number of sims for writing out timing
	*/
	DECCOData(std::string fname, int num_particles, int writes=-1, int steps=-1);

	//copy constructor to handle the data handlers
	DECCOData& operator=(const DECCOData& rhs);

	//Write the data given with the method specified during class construction
	//This includes taking care of headers for csv and metadata for hdf5
	//@param data_type is one of "simData", "constants", "energy", "timing",
	//@param data is the data to be written. Can only be an std::stringstream 
	//		(for csv data_type) or std::vector<double> (for hdf5 data_type)
	//@returns true if write succeeded, false otherwise
	bool Write(std::vector<double> &data, std::string data_type,int add_writes=0);
	std::string ReadMetaData(int data_index);

	std::string ReadMetaData(std::string data_type);
	std::vector<double> Read(std::string data_type, bool all=true, int line=-1,std::string file="");

	void WriteMeta(const std::string& metadata, const std::string& metadataName, \
        			const std::string& datasetName);
    //This function returns the status of how many writes have happened. It also sets the written_so_far values
    //	if there is more than zero and less then max_writes writes.
    //@return 0 if there is no writes so far (I don't think this should happen but if it does, more stuff needs to happen).
    //@return int >0 for how many writes there have been.
    //@return -1 if there are writes and the sim is already finished. 
    #ifdef HDF5_ENABLE
	    int setWrittenSoFar(const std::string path, const std::string file);
    #endif
  	void loadSimData(const std::string path, const std::string file,vec3 *pos,vec3 *w,vec3 *vel);
	int getNumTypes();
	int getWidth(std::string data_type);
	int getSingleWidth(std::string data_type);
	std::string genMetaData(int data_index);
	std::string getFileName();
	bool writeCheckpoint();
	bool deleteData();
	
	
	std::string filename;
private:
	
	std::string storage_type;
	int num_particles;
	int writes;
	int steps;
	int fixed_width = -1;
	bool fixed;
	bool h5data;
	bool csvdata;
	#ifdef HDF5_ENABLE
		HDF5Handler H; 
	#endif
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
	bool writeCSV(std::vector<double> &data, std::string data_type, std::string filename);

	//Write the data given with the hdf5 method.
	//This includes taking care of headers for csv and metadata for hdf5
	//@param data_type is one of "simData", "constants", "energy", "timing",
	//@param data is the data to be written.  
	//@returns true if write succeeded, false otherwise
	#ifdef HDF5_ENABLE
		bool writeH5(std::vector<double> &data, std::string data_type);
	#endif

	std::string genSimDataMetaData();
	std::string genConstantsMetaData();
	std::string genEnergyMetaData();
	std::string genTimingMetaData();

};

std::string find_restart_point(std::string path, const int index,bool relax=false);
std::string find_file_name(std::string path,int index,bool relax=false);
void delete_file(const std::string& delete_me);



#endif