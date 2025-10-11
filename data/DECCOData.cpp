
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <regex>
#include <cstring>


#include "DECCOData.hpp"
#include "../Collider/ball_group.hpp"
#include "../utilities/vec3.hpp"
#include "../utilities/Utils.hpp"
#include "../utilities/MPI_utilities.hpp"
#include "../utilities/simple_graph.hpp"


#ifdef EXPERIMENTAL_FILESYSTEM
	#include <experimental/filesystem>
	namespace fs = std::experimental::filesystem;
#else
	#include <filesystem>
	namespace fs = std::filesystem;
#endif

// #ifdef HDF5_ENABLE
// 	#include "H5Cpp.h"
// #endif

// #ifdef MPI_ENABLE
// 	#include "mpi.h"
// #endif

std::string trimFilename(const std::string& file)
{

    const auto pos = file.find_last_of('_');
    if (pos == std::string::npos) {
        MPIsafe_print(std::cerr, "ERROR: filename is not of correct format. Filename is '" + file + "'.\n");
        MPIsafe_exit(-1);
    }

    // If the chunk after '_' starts with "RELAX", keep "_RELAX" (exactly),
    // otherwise keep just up to and including the underscore.
    static constexpr char RELAX[] = "RELAX";
    if (file.compare(pos + 1, sizeof(RELAX) - 1, RELAX) == 0) {
        return file.substr(0, pos + 1 + (sizeof(RELAX) - 1));
    }

    return file.substr(0, pos + 1);
}
std::string getDataStringFromIndex(const int data_index)
{
	if (data_index < 0 || data_index > num_data_types-1)
	{
		MPIsafe_print(std::cerr,"DECCOData ERROR: data_index '"+std::to_string(data_index)+"' is out of range.\n");
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
	MPIsafe_print(std::cerr,"DECCOData ERROR: dataType '"+data_type+"' not found in class.\n");
	return -1;
}

void printVec(std::vector<double> v)
{
	if (getRank() == 0)
	{
		int i;
		for (i = 0; i < v.size()-1; i++)
		{
			std::cout<<v[i]<<", ";
		}
		std::cout<<v[i]<<std::endl;
	}
}


CSVHandler::CSVHandler()=default;
CSVHandler::~CSVHandler(){};
CSVHandler::CSVHandler(std::string filename) : filename(filename) 
{
	initialized = true;
}

bool CSVHandler::isInitialized()
{
	return initialized;
}

std::string CSVHandler::get_last_line(const std::string& path, const std::string& file)
{
	std::string filenameBase = trimFilename(file);
	std::string simDataFile = path + filenameBase + "simData.csv";

    if (auto simDataStream = std::ifstream(simDataFile, std::ifstream::in)) {
        MPIsafe_print(std::cerr,"\nParsing last line of data.\n");

        simDataStream.seekg(-1, std::ios_base::end);  // go to 
         // spot before the EOF

        bool keepLooping = true;
        bool first_run = true;
        while (keepLooping) {
            char ch = ' ';
            simDataStream.get(ch);  // Get current byte's data

            if (static_cast<int>(simDataStream.tellg()) <=
                1) {                     // If the data was at or before the 0th byte
                simDataStream.seekg(0);  // The first line is the last line
                keepLooping = false;     // So stop there
            } else if (ch == '\n' && not first_run) {     // If the data was a newline
                keepLooping = false;     // Stop at the current position (if we arent on the first character).
            } else {                     // If the data was neither a newline nor at the 0 byte
                simDataStream.seekg(-2, std::ios_base::cur);  // Move to the front of that data, then to
                                                              // the front of the data before it
            }
            first_run = false;
        }
        std::string line;
        std::getline(simDataStream, line);  // Read the current line
        return line;
    } else {

        std::string message("Could not open simData file: "+simDataFile+"... Exiting program.\n");
        MPIsafe_print(std::cerr,message);
        MPIsafe_exit(EXIT_FAILURE);
        return "ERROR"; //This shouldn't return but not returning anything is giving a warning
    }
}

void CSVHandler::loadCSVSimData(const std::string& path,const std::string& file,vec3 *pos,vec3 *w,vec3 *vel)
{
	std::string lineElement;
    // Get number of balls in file
    int num_particles = get_num_particles(path,file);

    std::string line = get_last_line(path,file);

    std::stringstream chosenLine(line);  // This is the last line of the read file, containing all data
                                         // for all balls at last time step
    
    int lineWidth = DECCOData::getSingleWidth("simData");
    // Get position and angular velocity data:
    for (int A = 0; A < num_particles; A++) 
    {
        for (int i = 0; i < 3; i++)  // Position
        {
            std::getline(chosenLine, lineElement, ',');
            pos[A][i] = std::stod(lineElement);
        }
        for (int i = 0; i < 3; i++)  // Angular Velocity
        {
            std::getline(chosenLine, lineElement, ',');
            w[A][i] = std::stod(lineElement);
        }
        std::getline(chosenLine, lineElement, ',');  // Angular velocity magnitude skipped
        for (int i = 0; i < 3; i++)                  // velocity
        {
            std::getline(chosenLine, lineElement, ',');
            vel[A][i] = std::stod(lineElement);
        }
        for (int i = 0; i < lineWidth - 10; i++)  // We used 10 elements. This skips the rest.
        {
            std::getline(chosenLine, lineElement, ',');
        }
    }
}

void CSVHandler::loadConsts(const std::string& path,const std::string& file,double *R,double *m,double *moi)
{
	std::string filenameBase = trimFilename(file);
	int num_particles = CSVHandler::get_num_particles(path,file);
	if (R != nullptr && m != nullptr && moi != nullptr)
	{
    	// Get radius, mass, moi:
	    std::string constantsFilename = path + filenameBase + "constants.csv";
	    if (auto ConstStream = std::ifstream(constantsFilename, std::ifstream::in)) {
	        std::string line, lineElement;
	        for (int A = 0; A < num_particles; A++) {
	            std::getline(ConstStream, line);  // Ball line.
	            std::stringstream chosenLine(line);
	            std::getline(chosenLine, lineElement, ',');  // Radius.
	            R[A] = std::stod(lineElement);
	            std::getline(chosenLine, lineElement, ',');  // Mass.
	            m[A] = std::stod(lineElement);
	            std::getline(chosenLine, lineElement, ',');  // Moment of inertia.
	            moi[A] = std::stod(lineElement);
	        }
	    } else {
	        MPIsafe_print(std::cerr,"Could not open constants file: " + constantsFilename + " ... Exiting program.\n");
	        MPIsafe_exit(EXIT_FAILURE);
	    }

	}
	else
	{
		MPIsafe_print(std::cerr,"ERROR: one or more of R, m, and/or moi are nullptr.");
	}
}

int CSVHandler::get_num_particles(std::string path, std::string filename)
{
	
	filename = trimFilename(filename);

	std::string constantsFilename = path + filename + "constants.csv";
    std::ifstream file(constantsFilename);
    if (!file.is_open()) return 0;

    int linecount = 0;
    std::string line;

    while (std::getline(file, line)) {
        // Check if the line has any non-whitespace characters
        bool nonempty = false;
        for (char c : line) {
            if (!std::isspace(static_cast<unsigned char>(c))) {
                nonempty = true;
                break;
            }
        }
        if (nonempty) ++linecount;
    }

    return linecount;
}

bool CSVHandler::writeSimData(std::vector<double> data, int width, std::string filename)
{
	filename += "simData.csv"; 
	try
	{
		std::string meta = "";
		if (!fs::exists(filename))
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

bool CSVHandler::writeEnergy(std::vector<double> data, int width, std::string filename)
{
	filename += "energy.csv"; 
	try
	{
		std::string meta = "";

		if (!fs::exists(filename))
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

bool CSVHandler::writeConstants(std::vector<double> data, int width, std::string filename)
{
	filename += "constants.csv"; 
	try
	{
		//Consts has no meta data for csv
		// std::string meta = "";

		// if (!fs::exists(filename))
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

// bool CSVHandler::writeTiming(std::vector<double> data, int width, std::string filename)
// {
// 	filename += "timing.csv"; 
// 	try
// 	{
// 		std::ofstream constsWrite;
// 		constsWrite.open(filename, std::ofstream::app);

// 		if (!fs::exists(filename))
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

bool CSVHandler::deleteData()
{
	bool simData_result = false;
	bool consts_result = false;
	bool energy_result = false;
	std::cerr<<"filename base: "<<filename<<std::endl;
	try {
		std::string simData_filename = filename+"simData.csv";
        if (std::filesystem::remove(simData_filename)) {
        	std::string message = "File '"+simData_filename+"' deleted successfully.\n";
        	MPIsafe_print(std::cerr,message);
        	simData_result = true;
        } else {
        	std::string message = "File '"+simData_filename+"' not found or could not be deleted.\n";
        	MPIsafe_print(std::cerr,message);
        	MPIsafe_exit(-1);
        }
    } catch (const std::filesystem::filesystem_error& ex) {
    	std::stringstream message;
    	message <<"Filesystem error: "<<ex.what()<<'\n';
        MPIsafe_print(std::cerr,message.str());
        MPIsafe_exit(-1);
    }

    try {
		std::string constants_filename = filename+"constants.csv";
        if (std::filesystem::remove(constants_filename)) {
        	std::string message = "File '"+constants_filename+"' deleted successfully.\n";
        	MPIsafe_print(std::cerr,message);
        	consts_result = true;
        } else {
        	std::string message = "File '"+constants_filename+"' not found or could not be deleted.\n";
        	MPIsafe_print(std::cerr,message);
        	MPIsafe_exit(-1);
        }
    } catch (const std::filesystem::filesystem_error& ex) {
    	std::stringstream message;
    	message <<"Filesystem error: "<<ex.what()<<'\n';
        MPIsafe_print(std::cerr,message.str());
        MPIsafe_exit(-1);
    }

    try {
		std::string energy_filename = filename+"energy.csv";
        if (std::filesystem::remove(energy_filename)) {
        	std::string message = "File '"+energy_filename+"' deleted successfully.\n";
        	MPIsafe_print(std::cerr,message);
        	energy_result = true;
        } else {
        	std::string message = "File '"+energy_filename+"' not found or could not be deleted.\n";
        	MPIsafe_print(std::cerr,message);
        	MPIsafe_exit(-1);
        }
    } catch (const std::filesystem::filesystem_error& ex) {
    	std::stringstream message;
    	message <<"Filesystem error: "<<ex.what()<<'\n';
        MPIsafe_print(std::cerr,message.str());
        MPIsafe_exit(-1);
    }

    return simData_result && consts_result && energy_result;
}

std::string CSVHandler::genSimDataMetaData(int num_particles)
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

std::string CSVHandler::genConstantsMetaData()
{
	//There is no metadata for constant file with csv format. THis is just
	//here to be complete and stay compatable with the previous version.
	return "";
}

std::string CSVHandler::genEnergyMetaData()
{
	return "Time,PE,KE,E,p,L\n";
}

std::string CSVHandler::genTimingMetaData()
{
	//There is no metadata for timing file with csv format. THis is just
	//here to be complete and stay compatable with the previous version.
	return "";
}

#ifdef HDF5_ENABLE
HDF5Handler::HDF5Handler()=default;
HDF5Handler::~HDF5Handler(){};
HDF5Handler::HDF5Handler(std::string filename,bool fixed=false) : filename(filename),fixed(fixed) 
{
	initialized = true;
}

//Make start a class variable so you dont need to keep track
void HDF5Handler::createAppendFile(std::vector<double>& data,const std::string datasetName,size_t start=0,hsize_t maxdim=0) {
    // H5::H5File file;
    // H5::DataSet dataset;
    // H5::DataSpace dataspace;

    if(fs::exists(filename)) {
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
std::vector<double> HDF5Handler::readFile(const std::string datasetName, hsize_t start=0, hsize_t len=0) {
    std::vector<double> data;
    if (fs::exists(filename)) {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        H5::DataSet dataset = file.openDataSet(datasetName);
        H5::DataSpace dataspace = dataset.getSpace();

        hsize_t dims_out[1];
        dataspace.getSimpleExtentDims(dims_out, NULL);
        hsize_t total_size = dims_out[0];

        if (start > total_size)
        {
        	MPIsafe_print(std::cerr, "DECCOData ERROR: invalid start input\n");
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
    	MPIsafe_print(std::cerr,"File '" + filename + "' does not exist.\n");
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
std::vector<double> HDF5Handler::static_readFile(const std::string datasetName, hsize_t start, hsize_t len,bool neg_offset,std::string readfile) 
{
    std::vector<double> data;
    if (fs::exists(readfile)) {
        H5::H5File file(readfile, H5F_ACC_RDONLY);
        H5::DataSet dataset = file.openDataSet(datasetName);
        H5::DataSpace dataspace = dataset.getSpace();

        hsize_t dims_out[1];
        dataspace.getSimpleExtentDims(dims_out, NULL);
        hsize_t total_size = dims_out[0];


        if (start > total_size)
        {
        	MPIsafe_print(std::cerr,"DECCOData ERROR: invalid start input\n");
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
    	MPIsafe_print(std::cerr,"File '" + readfile + "' does not exist.\n");
    }

    return data;
}

void HDF5Handler::loadh5SimData(const std::string path, const std::string file,vec3 *pos,vec3 *w,vec3 *vel)
{
	int n_particles = get_num_particles(path,file);
	int single_ball_width = single_ball_widths[getDataIndexFromString("simData")];
	int width = n_particles*single_ball_width;
	std::vector<double> out = static_readFile("simData",width,0,true,path+file); // get last line
	// printVec(out);
	
	for (int i = 0; i < n_particles; ++i)
	{
		int out_ind = i*single_ball_width;
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

hsize_t HDF5Handler::get_data_length(std::string readfile,std::string datasetName)
{
	hsize_t total_size;
	if (fs::exists(readfile)) {
        H5::H5File file(readfile, H5F_ACC_RDONLY);
        H5::DataSet dataset = file.openDataSet(datasetName);
        H5::DataSpace dataspace = dataset.getSpace();

        hsize_t dims_out[1];
        dataspace.getSimpleExtentDims(dims_out, NULL);
        total_size = dims_out[0];
    } else {
    	MPIsafe_print(std::cerr,"File '" + readfile + "' does not exist.\n");
    }
    return total_size;
}

int HDF5Handler::readWrites(const std::string path, const std::string f)
{
	std::string file_read = path + f;
	// std::cerr<<"readWrites file: "<<file_read<<std::endl;
	const std::string datasetName("writes");
    H5::H5File file;
    H5::DataSet dataset;
    int value = 0; // uninitialized value

    if(fs::exists(file_read)) {
        file = H5::H5File(file_read, H5F_ACC_RDWR);
    } else {
    	MPIsafe_print(std::cerr,"File " + file_read + " doesn't exist.\n");
        MPIsafe_exit(-1);
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

void HDF5Handler::addWrites(int additional_writes) 
{
	const std::string datasetName("writes");
    H5::H5File file;
    H5::DataSet dataset;
    int value = additional_writes; // Default value

    if(fs::exists(filename)) {
        file = H5::H5File(filename, H5F_ACC_RDWR);
    } else {
    	MPIsafe_print(std::cerr,"File " + filename + " doesn't exist.\n");
        MPIsafe_exit(-1);
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

bool HDF5Handler::sim_finished(std::string path, std::string file)
{
    std::vector<double> simdata = static_readFile("simData", 0, 0, false, path + file);
    if (!simdata.empty())
    {
        for (double value : simdata)
        {
            if (std::isnan(value))
            {
                return false;
            }
        }
        return true;
    }
    return false;
}


int HDF5Handler::get_num_particles(std::string path, std::string file)
{
	hsize_t total_size = get_data_length(path+file,"constants");

    return (total_size*1.0)/(single_ball_widths[getDataIndexFromString("constants")]*1.0);
}

void HDF5Handler::loadConsts(std::string path,std::string file,double *R,double *m,double *moi)
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
//     if(fs::exists(filename)) {
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

void HDF5Handler::attachMetadataToDataset(const std::string& metadata, const std::string& datasetName) 
{
    // Initialize the HDF5 library
	if(fs::exists(filename)) {
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
    	MPIsafe_print(std::cerr,"File '"+filename+"' does not exist.\n");
    }
    return;
}

void HDF5Handler::attachSimMetadataToDataset(const std::string& metadata, const std::string& metadataName, \
							const std::string& datasetName) 
{
    // Initialize the HDF5 library
	H5::H5File file;
	if(!fs::exists(filename)) 
	{
		file = H5::H5File(filename, H5F_ACC_TRUNC);
	}
	else
	{
    	file = H5::H5File(filename, H5F_ACC_RDWR);
	}

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
    // }else{
    // 	std::cout<<"File '"<<filename<<"' does not exist."<<std::endl;
    // }
    return;
}

//Non static version for reading from initialized class instance
std::string HDF5Handler::readMetadataFromDataset(const std::string& datasetName,const std::string& metadataName="metadata") 
{
    std::string metadata;

    // Initialize the HDF5 library
    if(fs::exists(filename)) {
    	H5::DataSet dataset;
        H5::H5File file(filename, H5F_ACC_RDONLY);
	    // Open the specified dataset
	    if (datasetExists(filename,datasetName)){
	    	dataset = file.openDataSet(datasetName);
	    }else{
	    	MPIsafe_print(std::cerr,"dataset '"+datasetName+"' does not exist for file '"+filename+"' .\n");
	    	return ERR_RET_MET;
	    }

	    // Check if the attribute's dataspace exists
	    if (!attributeExists(dataset, metadataName)) {
	    	MPIsafe_print(std::cerr,"Attribute 'metadata' does not exist for dataset '"+datasetName+"' in file '"+filename+"' .\n");
	        dataset.close();
	        file.close();
	    	return ERR_RET_MET;
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
    	MPIsafe_print(std::cerr,"File '"+filename+"' does not exist.\n");
    }

    return metadata;
}

//static version for reading when no DECCOData class instance exists
std::string HDF5Handler::readMetadataFromDataset(const std::string& datasetName,const std::string& metafile,const std::string& metadataName="metadata") 
{
    std::string metadata;

    // Initialize the HDF5 library
    if(fs::exists(metafile)) {
    	H5::DataSet dataset;
        H5::H5File file(metafile, H5F_ACC_RDONLY);
	    // Open the specified dataset
	    if (datasetExists(metafile,datasetName)){
	    	dataset = file.openDataSet(datasetName);
	    }else{
	    	MPIsafe_print(std::cerr,"dataset '"+datasetName+"' does not exist for file '"+metafile+"' .\n");
	    	return ERR_RET_MET;
	    }

	    // Check if the attribute's dataspace exists
	    if (!attributeExists(dataset, metadataName)) {
	    	MPIsafe_print(std::cerr,"Attribute "+metadataName+" does not exist for dataset '"+datasetName+"' in file '"+metafile+"' .\n");
	        dataset.close();
	        file.close();
	    	return ERR_RET_MET;
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
    	MPIsafe_print(std::cerr,"File '"+metafile+"' does not exist.\n");
    }



    return metadata;
}

bool HDF5Handler::isInitialized()
{
	return initialized;
}

bool HDF5Handler::datasetExists(const std::string& f, const std::string& datasetName)
{
    // Open the HDF5 file
    H5::H5File file(f, H5F_ACC_RDONLY);

    // Check if the dataset exists
    bool exists = H5Lexists(file.getId(), datasetName.c_str(), H5P_DEFAULT) > 0;
    file.close();
    return exists;
}

bool HDF5Handler::attributeExists(const H5::H5Object& object, const std::string& attributeName) 
{
	return object.attrExists(attributeName);
}

std::string HDF5Handler::genSimDataMetaData(int width)
{
	std::string meta_data = "";
	meta_data = "Columns: posx[ball],posy[ball],posz[ball],wx[ball],wy[ball],wz[ball],wtot[ball],velx[ball],vely[ball],velz[ball],bound[ball]\n";
	meta_data += "rows: writes\n";
	meta_data += "row width: " + std::to_string(width);

	return meta_data;
}

std::string HDF5Handler::genConstantsMetaData(int width)
{
	std::string meta_data = "Columns: radius, mass, moment of inertia\n";
	meta_data += "rows: balls\n";	
	meta_data += "row width: " + std::to_string(width);

	return meta_data;
}

std::string HDF5Handler::genEnergyMetaData(int width)
{
	std::string meta_data = "Columns: time, PE, KE, Etot, p, L\n";
	meta_data += "rows: writes\n";	
	meta_data += "row width: " + std::to_string(width);
	return meta_data;	
}

std::string HDF5Handler::genTimingMetaData(int width)
{
	std::string meta_data = "Columns: number of balls, time spent in sim_one_step\n";
	meta_data += "rows: sims\n";	
	meta_data += "row width: " + std::to_string(width);
	return meta_data;	
}

   



H5::DataSet HDF5Handler::createDataSet(std::vector<double>& data,const std::string datasetName,
										int fileExists,int max_size)
{
	H5::H5File file;
    H5::DataSet dataset;
    H5::DataSpace dataspace;

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
    	std::vector<double> data_write(max_size, std::nan("")); //initialize the data set with nan values (besides the initial conditions of course)
    	for (int i = 0; i < data.size(); ++i)
    	{
    		data_write[i] = data[i];
    	}

    	maxdims[0] = max_size; // Set maximum dimensions to max_size
    	// std::cout<<"maxsize: "<<max_size<<std::endl;
    	dataspace = H5::DataSpace(1, maxdims);
        dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE,dataspace);
    	dataset.write(&data_write[0], H5::PredType::NATIVE_DOUBLE);
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

void HDF5Handler::appendDataSet(std::vector<double>& data,const std::string datasetName,size_t start,int max_size)
{
	H5::H5File file;
    H5::DataSet dataset;
    H5::DataSpace dataspace;
    if(fs::exists(filename)) 
    {
		file = H5::H5File(filename, H5F_ACC_RDWR);
	}
	else
	{
		MPIsafe_print(std::cerr,"File "+filename+" doesn't exist.\n");
		MPIsafe_exit(-1);
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
		MPIsafe_print(std::cerr,"H5 ERROR: "+error.getDetailMsg()+'\n');
	}

	try {
	    dataspace = dataset.getSpace();
	} catch (const H5::Exception& error) {
		MPIsafe_print(std::cerr,"H5 ERROR: "+error.getDetailMsg()+'\n');
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

bool HDF5Handler::deleteData()
{
	bool result = false;
	try {
        if (std::filesystem::remove(filename)) {
        	std::string message = "File '"+filename+"' deleted successfully.\n";
        	MPIsafe_print(std::cerr,message);
        	result = true;
        } else {
        	std::string message = "File '"+filename+"' not found or could not be deleted.\n";
        	MPIsafe_print(std::cerr,message);
        	MPIsafe_exit(-1);
        }
    } catch (const std::filesystem::filesystem_error& ex) {
    	std::stringstream message;
    	message <<"Filesystem error: "<<ex.what()<<'\n';
        MPIsafe_print(std::cerr,message.str());
        MPIsafe_exit(-1);
    }

    return result;
}

#endif




DECCOData::~DECCOData(){};
/*
Mandatory class input
	storage_method: either "h5"/"hdf5" for an hdf5 file format or "csv" for a csv format. Will be whatever the file
					extension is in "fname" variable
Optional class input
	num_particles : the number of particles in the DECCO simulation, only needed for a fixed hdf5 storage
	writes        : the number of writes that will happen in the DECCO simulation, only needed for a fixed hdf5 storage
	steps   	  : the number of sims for writing out timing
*/
DECCOData::DECCOData(std::string fname, int num_particles, int writes, int steps) : 
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



	int dot_index = filename.find_last_of('.');
	storage_type = filename.substr(dot_index+1,filename.length()-dot_index);
	//Transform storage_type to be not case sensitive
	std::transform(storage_type.begin(), storage_type.end(), storage_type.begin(), ::tolower);
	if (storage_type == "h5" || storage_type == "hdf5")
	{
		h5data = true;
		csvdata = false;
		H = HDF5Handler(filename,fixed);
	}
	else if (storage_type == "csv")
	{
		filename = trimFilename(filename);
		csvdata = true;
		h5data = false;
		C = CSVHandler(filename);
		// filename = filename.substr(0,filename.length()-4);
	}
	else
	{
		MPIsafe_print(std::cerr,"DECCOData ERROR: storage_type '"+storage_type+"' not available.\n");
		MPIsafe_exit(-1);
	}
	#ifndef HDF5_ENABLE
		MPIsafe_print(std::cerr,"HDF5 not enabled. CSV format will be used for data.\n");
		csvdata = true;
		h5data = false;
		// if (storage_type != "csv")
		// {
		// 	// filename = filename.substr(0,filename.length()-4);
		// }
	#endif


	//If user specified number of writes but a storage_type other than hdf5 then default to hdf5 and warn user
	// if (fixed && not h5data)
	// {
	// 	std::cerr<<"DECCOData warning: specified a non-negative number of writes (which indicates a fixed size hdf5 storage) and a storage_type other than hdf5. Defaulting to a fixed size hdf5 storage_type."<<std::endl;
	// }

}



//copy constructor to handle the data handlers
DECCOData& DECCOData::operator=(const DECCOData& rhs)
{
	C=rhs.C;
	#ifdef HDF5_ENABLE
		H=rhs.H;
	#endif
	return *this;
}

//Write the data given with the method specified during class construction
//This includes taking care of headers for csv and metadata for hdf5
//@param data_type is one of "simData", "constants", "energy", "timing",
//@param data is the data to be written. Can only be an std::stringstream 
//		(for csv data_type) or std::vector<double> (for hdf5 data_type)
//@returns true if write succeeded, false otherwise
bool DECCOData::Write(std::vector<double> &data, std::string data_type,int add_writes)
{
	bool retVal;
	if (h5data)
	{
		#ifdef HDF5_ENABLE
			retVal = writeH5(data,data_type);
			if (add_writes > 0)
			{
				H.addWrites(add_writes);
			}
		#else
			return 0;
		#endif
	}
	else if (csvdata)
	{
		retVal = writeCSV(data,data_type,filename);
	}
	return retVal;
}


std::string DECCOData::ReadMetaData(int data_index) 
{
	//if data_index is less than zero a bad data_type was input and write doesn't happen
	if (data_index < 0)
	{
		return ERR_RET_MET;
	}
	std::string datasetName = getDataStringFromIndex(data_index);
	std::string readMetaData;
	if (h5data)
	{
		#ifdef HDF5_ENABLE
			//Has the HDF5 handler been initiated yet?
			if (!H.isInitialized())
			{
				H = HDF5Handler(filename,fixed);
			}

			readMetaData = H.readMetadataFromDataset(datasetName);
		#else
			return "ERROR";
		#endif
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

std::string DECCOData::ReadMetaData(std::string data_type) 
{
	int data_index = getDataIndexFromString(data_type);
	std::string readMetaData = ReadMetaData(data_index);
	return readMetaData;
}

std::vector<double> DECCOData::Read(std::string data_type, bool all, int line,std::string file) 
{
	if (file == "")
	{
		file = filename;
	}
	
	std::vector<double> data_read;
	if (file.substr(file.size()-3,file.size()) == ".h5")
	{
		#ifdef HDF5_ENABLE
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

			}
		#else // csv's are still being read by functions in ball_group
			MPIsafe_print(std::cerr,"ERROR: csv file type not yet readable by DECCOData.\n");
			MPIsafe_exit(-1);
		#endif	
	}
	else if (file.substr(file.size()-4,file.size()) == ".csv")
	{
		MPIsafe_print(std::cerr,"ERROR: csv file type not yet readable by DECCOData.\n");
		MPIsafe_exit(-1);
	}
	return data_read;
}

void DECCOData::WriteMeta(const std::string& metadata, const std::string& metadataName, \
    			const std::string& datasetName) 
{
	#ifdef HDF5_ENABLE
		if (h5data)
		{
    		H.attachSimMetadataToDataset(metadata,metadataName,datasetName);
		}
		else
		{
			MPIsafe_print(std::cerr,"Function WriteMeta only available for h5 ouput data format.\n");
			MPIsafe_exit(-1);
		}
	#else
		MPIsafe_print(std::cerr,"Function WriteMeta only available for h5 ouput data format.\n");
		MPIsafe_exit(-1);
    #endif
}


//This function returns the status of how many writes have happened. It also sets the written_so_far values
//	if there is more than zero and less then max_writes writes.
//@return 0 if there is no writes so far (I don't think this should happen but if it does, more stuff needs to happen).
//@return int >0 for how many writes there have been.
//@return -1 if there are writes and the sim is already finished. 
#ifdef HDF5_ENABLE
int DECCOData::setWrittenSoFar(const std::string path, const std::string file)
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
#endif

bool DECCOData::deleteData()
{
	bool result = false;
	if (h5data)
	{
		result = H.deleteData();
	}
	else
	{
		result = C.deleteData();
	}
	return result;
}


void DECCOData::loadSimData(const std::string path, const std::string file,vec3 *pos,vec3 *w,vec3 *vel)
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

int DECCOData::getNumTypes()
{
	return num_data_types;
}

int DECCOData::getWidth(std::string data_type)
{
	return widths[getDataIndexFromString(data_type)];
}

int DECCOData::getSingleWidth(std::string data_type)
{
	return single_ball_widths[getDataIndexFromString(data_type)];
}


std::string DECCOData::genMetaData(int data_index)
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
	MPIsafe_print(std::cerr,"DECCOData ERROR: data_index '"+std::to_string(data_index)+"' is out of range.\n");
	return "DECCOData ERROR";
}

std::string DECCOData::getFileName()
{
	return filename;
}

bool DECCOData::writeCheckpoint()
{
	std::string checkpt_file;
	if (csvdata)
	{
		checkpt_file = filename+"checkpoint.txt";
	}
	else
	{
		checkpt_file = filename.substr(0,filename.find_last_of('_')+1)+"checkpoint.txt";
	}

	// std::cerr<<"CREATING FILE: "<<checkpt_file<<std::endl;
	std::ofstream output(checkpt_file);
	if (not output.good()) 
	{
        auto const errNo = errno;
        std::string errMessage = "Failed to create file "+checkpt_file+": "+std::to_string(errNo)+": "+std::strerror(errNo);
        throw std::runtime_error(errMessage);
        return false;
	}
	return true;
}







//Write the data given with the csv method.
//This includes taking care of headers for csv and metadata for hdf5
//@param data_type is one of "simData", "constants", "energy", "timing",
//@param data is the data to be written.  
//@param filename is the absolute path to the base save file (everything except which data it is and the file extension)
//@returns true if write succeeded, false otherwise
bool DECCOData::writeCSV(std::vector<double> &data, std::string data_type, std::string filename)
{
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
#ifdef HDF5_ENABLE
bool DECCOData::writeH5(std::vector<double> &data, std::string data_type)
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
#endif



std::string DECCOData::genSimDataMetaData()
{
	std::string meta_data = "";
	if (h5data)
	{
		#ifdef HDF5_ENABLE
			meta_data = H.genSimDataMetaData(widths[getDataIndexFromString("simData")]);
		#else
			return "ERROR in genSimDataMetaData";
		#endif
	}
	else if (csvdata)
	{
		meta_data = C.genSimDataMetaData(num_particles);
	}
	return meta_data;
}

std::string DECCOData::genConstantsMetaData()
{
	std::string meta_data = "";
	if (h5data)
	{
		#ifdef HDF5_ENABLE
			meta_data = H.genConstantsMetaData(widths[getDataIndexFromString("constants")]);
		#else
			return "ERROR in genConstantsMetaData";
		#endif
	}
	else if (csvdata)
	{
		meta_data = C.genConstantsMetaData();
	}
	return meta_data;
}

std::string DECCOData::genEnergyMetaData()
{
	std::string meta_data = "";
	if (h5data)
	{
		#ifdef HDF5_ENABLE
			meta_data = H.genEnergyMetaData(widths[getDataIndexFromString("energy")]);
		#else
			return "ERROR in genEnergyMetaData";
		#endif
	}
	else if (csvdata)
	{
		meta_data = C.genEnergyMetaData();
	}
	return meta_data;
}

std::string DECCOData::genTimingMetaData()
{
	std::string meta_data = "";
	if (h5data)
	{
		#ifdef HDF5_ENABLE
			meta_data = H.genTimingMetaData(widths[getDataIndexFromString("timing")]);
		#else
			return "ERROR in genTimingMetaData";
		#endif
	}
	else if (csvdata)
	{
		meta_data = C.genTimingMetaData();
	}
	return meta_data;
}

//Returns the path + filename of the specified index
//If index < 0 (default is -1) then it will return the largest (completed) index for csv
//and returns largest index that passes allConnected for hdf5 
//If relax is true, it looks for the relax file instead of the regular file.
//WARNING: this function also deletes unfinished files that cannot be used
//		   as a restart point, making it a bit dangerous to edit.
std::string find_restart_point(std::string path, const int index,bool relax/*=false*/)
{
	std::string simDatacsv;
	std::string datah5;
	std::string relax_string;

	if (relax)
    {
    	simDatacsv = "_RELAXsimData.csv";
	    datah5 = "_RELAXdata.h5";
		relax_string = "RELAX";
    }
    else
    {
	    simDatacsv = "_simData.csv";
	    datah5 = "_data.h5";
		relax_string = "";
    }
	std::string largest_file_name;
	// std::string simDatacsv = "simData.csv";
    // std::string datah5 = "data.h5";
    int csv = -1;
    std::string dataType = data_type_from_input(path);
    if (dataType == "csv")
    {
        csv = 1;
    }
    else if (dataType == "h5" || dataType == "hdf5")
    {
        csv = 0;
    }
    else
    {
        std::string test_file = path+std::to_string(index);
        if (fs::exists(test_file+simDatacsv))
        {
            return test_file + simDatacsv;
        }
        else if (fs::exists(test_file+datah5))
        {
            return test_file + datah5;
        }
        // else if (fs::exists(test_file+"_data.hdf5"))
        // {
        //     return test_file + "_data.hdf5";
        // }
        else
        {
            MPIsafe_print(std::cerr,"ERROR in find_whole_file_name: dataType '"+dataType+"' not recognized.");
            MPIsafe_exit(-1);
        }
    }

    //if index is greater than zero we know if its csv or h5 so just return here
    if (index >= 0)
    {
        if (csv == 1)
        {
            return path + std::to_string(index) + simDatacsv;
        }
        else
        {
            return path + std::to_string(index) + datah5;
        }
    }

    bool empty = false;

    if (getRank() == 0)
    {

	    std::regex pattern;
	    if (csv)
	    {
	    	if (relax){pattern = std::regex(R"((\d+)_RELAXsimData\.csv)");}
	    	else {pattern = std::regex(R"((\d+)_simData\.csv)");}
	    }
	    else
	    {
	    	if (relax){pattern = std::regex(R"((\d+)_RELAXdata\.h5)");}
	    	else {pattern = std::regex(R"((\d+)_data\.h5)");}
	    	
	    }

	    // Store (index, filename) pairs
	    std::vector<std::pair<int, std::string>> files;

	    for (const auto& entry : fs::directory_iterator(path)) {
	        if (entry.is_regular_file()) {
	            std::string fname = entry.path().filename().string();
	            std::smatch match;
	            if (std::regex_match(fname, match, pattern)) {
	                int index = std::stoi(match[1].str());
	                files.emplace_back(index, fname);
	            }
	        }
	    }

	    if (files.empty()) {empty = true;}

	    // Sort by index descending
	    std::sort(files.begin(), files.end(),
	              [](auto& a, auto& b) { return a.first > b.first; });



	    // Print results
	    for (auto it = files.begin(); it != files.end(); /* no ++ here */) 
	    {	
	    	const auto idx  = it->first;
			const auto name = it->second;
			std::string checkpt_file = path+std::to_string(idx)+"_"+relax_string+"checkpoint.txt";

	    	bool checkpoint_exists = fs::exists(checkpt_file);

	    	if (csv == 1 && index < 0 && !checkpoint_exists) 
		    {

	            std::string file1 = path + name;
	            std::string file2 = path + name.substr(0,name.size()-simDatacsv.size()) + "_" + relax_string + "constants.csv";
	            std::string file3 = path + name.substr(0,name.size()-simDatacsv.size()) + "_" + relax_string + "energy.csv";

	            std::cerr<<"Removing the following files: \n\t"<<file1<<"\n\t"<<file2<<"\n\t"<<file3<<std::endl;

	            delete_file(file1);
	            delete_file(file2);
	            delete_file(file3);
		        
		        it = files.erase(it);
		    }
			else if (csv == 0 && index < 0 && !checkpoint_exists)
			{
		    	Ball_group temp(idx,/*JKR=*/false);
		    	temp.parse_input_file(path);
		    	temp.loadSim(path,name,/*verbose=*/false);
		        if (!isConnected(temp.pos,temp.R,temp.attrs.num_particles))
		        {
		            std::string file1 = path + name;
	            	std::cerr<<"Removing the following files: \n\t"<<file1<<std::endl;
		            delete_file(file1);
	            	// std::cerr<<"Removing the following files: \n\t"<<checkpt_file<<std::endl;
		            // delete_file(checkpt_file);

		            it = files.erase(it);
				}
				else
				{
					break;
				}
			}
			else
			{
				break;
			}
		}

		if (!files.empty()) {
            largest_file_name = files.front().second;  // pick the new top
        } else {
            largest_file_name = "";
        }
	}

	MPIsafe_bcast_string(largest_file_name, /*root=*/0);

	if (largest_file_name == "")
	{
		return "";
	}
    return path + largest_file_name;
}

//This function returns the name of the file in the directory path with the index index.
//If relax is true, it looks for the relax file instead of the regular file
//It differes from find_whole_file_name because it won't delete anything and you give an index of -1.
std::string find_file_name(std::string path,int index,bool relax/*=false*/)
{
    std::string simDatacsv;
	std::string datah5;

	if (relax)
    {
    	simDatacsv = "_RELAXsimData.csv";
	    datah5 = "_RELAXdata.h5";
    }
    else
    {
	    simDatacsv = "_simData.csv";
	    datah5 = "_data.h5";
    }

    std::string file;
    int file_index;


    for (const auto & entry : fs::directory_iterator(path))
    {

        file = entry.path();
        size_t slash_pos = file.find_last_of("/");
        file = file.erase(0,slash_pos+1);
        size_t _pos = file.find_first_of("_");

        if (_pos != std::string::npos) // only go in here if we found a data file
        {
            //Is the data in csv format? (but first verify the call to substr wont fail)
            if (file.size() >= simDatacsv.size() && file.substr(file.size()-simDatacsv.size(),file.size()) == simDatacsv)
            {
                file_index = stoi(file.substr(0,file.find_first_of("_")));
               
                if (file_index == index)
                {
                    return path+file;
                }

            }
            else if (file.size() >= datah5.size() && file.substr(file.size()-datah5.size(),file.size()) == datah5)
            {
                file_index = stoi(file.substr(0,file.find_first_of("_")));
                if (file_index == index)
                {
                    return path+file;
                }
            }
        }
    }
    
    MPIsafe_print(std::cerr,"ERROR: file at path '"+path+"' with index '"+std::to_string(index)+"' not found. Now exiting . . .\n");
    exit(-1);
}

//Deletes the given file. Should give this function a full file path.
void delete_file(const std::string& delete_me)
{
	if (fs::exists(delete_me))
	{
		int status = remove(delete_me.c_str());
		if (status != 0)
	    {
	        MPIsafe_print(std::cerr,"File: '"+delete_me+"' could not be removed, now exiting with failure.\n");
	        MPIsafe_exit(EXIT_FAILURE);
	    }
	}
	else
	{
	    MPIsafe_print(std::cerr,"File: '"+delete_me+"' does not exist. Skipping deletion.\n");
	}
}