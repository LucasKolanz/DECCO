#include <iostream>
#include <cmath>
#include <fstream>
#include "grid_group.hpp"
#include "vec3.hpp"

class wrapper
{
public:
	vec3 *pos = nullptr;
	int num_particles;
	grid *g = nullptr;
	
	wrapper() = default;

	wrapper(int num,std::string file,double size)
	{
		num_particles = num;
		pos = new vec3[num];
		parseSimData(getLastLine(file));

		

		g = new grid(200,1e-5,size,pos);
	}

	void parseSimData(std::string line)
    {
        std::string lineElement;
        // Get number of balls in file
        int count = 269;
        int properties = 11;
        
        // int count = std::count(line.begin(), line.end(), ',') / properties + 1;
        // allocate_group(count);
        std::stringstream chosenLine(line);  // This is the last line of the read file, containing all data
                                             // for all balls at last time step
        // Get position and angular velocity data:
        for (int A = 0; A < num_particles; A++) {
            for (int i = 0; i < 3; i++)  // Position
            {
                std::getline(chosenLine, lineElement, ',');
                pos[A][i] = std::stod(lineElement);
            }
            for (int i = 0; i < 3; i++)  // Angular Velocity
            {
                std::getline(chosenLine, lineElement, ',');
                // w[A][i] = std::stod(lineElement);
            }
            std::getline(chosenLine, lineElement, ',');  // Angular velocity magnitude skipped
            for (int i = 0; i < 3; i++)                  // velocity
            {
                std::getline(chosenLine, lineElement, ',');
                // vel[A][i] = std::stod(lineElement);
            }
            for (int i = 0; i < properties - 10; i++)  // We used 10 elements. This skips the rest.
            {
                std::getline(chosenLine, lineElement, ',');
            }
        }
    }

	/// Get last line of previous simData by filename.
    [[nodiscard]] static std::string getLastLine(const std::string& filename)
    {
        std::string simDataFilepath = filename;
        if (auto simDataStream = std::ifstream(simDataFilepath, std::ifstream::in)) {
            std::cerr << "\nParsing last line of data.\n";

            simDataStream.seekg(-1, std::ios_base::end);  // go to one spot before the EOF

            bool keepLooping = true;
            while (keepLooping) {
                char ch = ' ';
                simDataStream.get(ch);  // Get current byte's data

                if (static_cast<int>(simDataStream.tellg()) <=
                    1) {                     // If the data was at or before the 0th byte
                    simDataStream.seekg(0);  // The first line is the last line
                    keepLooping = false;     // So stop there
                } else if (ch == '\n') {     // If the data was a newline
                    keepLooping = false;     // Stop at the current position.
                } else {                     // If the data was neither a newline nor at the 0 byte
                    simDataStream.seekg(-2, std::ios_base::cur);  // Move to the front of that data, then to
                                                                  // the front of the data before it
                }
            }
            std::string line;
            std::getline(simDataStream, line);  // Read the current line
            return line;
        } else {
            std::cerr << "Could not open simData file: " << simDataFilepath << "... Existing program."
                      << '\n';
            exit(EXIT_FAILURE);
        }
        return "-1";
    }
};

void test_extrema()
{
	std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/198_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_simData.csv";
	wrapper w(200,file,4e-5);

	

	// std::cout<<"HEREERER"<<std::endl;
	// w.g.setMaxGridIndex();
	std::cout<<"maxGridIndex: "<<w.g -> maxGridIndex<<std::endl;
	//should be 2
}

int main(int argc, char const *argv[])
{
	test_extrema();

	return 0;
}

