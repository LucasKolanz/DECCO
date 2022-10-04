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
		

		g = new grid(num,1e-5,size,pos);
	}

	wrapper(int num,double size,vec3* positions)
	{
		num_particles = num;
		// pos = new vec3[num];
		pos = positions;


		g = new grid(num,1e-5,size,pos);
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

void test_groups1()
{
	// std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/198_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_simData.csv";
	vec3 *pos = nullptr;
	pos = new vec3[8];
	// pos[0] = {0.5,0.75,0.88}; // should be ID {0,0,0}
	// pos[1] = {-0.5,-0.75,-0.88}; // should be ID {-1,-1,-1}
	// pos[2] = {-0.5,0.75,0.88}; // should be ID {-1,0,0}
	// pos[3] = {0.5,-0.75,0.88}; // should be ID {0,-1,0}
	// pos[4] = {0.5,0.75,-0.88}; // should be ID {0,0,-1}
	// pos[5] = {-0.5,0.75,-0.88}; // should be ID {-1,0,-1}
	// pos[6] = {0.5,-0.75,-0.88}; // should be ID {0,-1,-1}
	// pos[7] = {-0.5,-0.75,0.88}; // should be ID {-1,-1,0}

	pos[0] = {1.5,1.75,1.88}; // should be ID {1,1,1}
	pos[1] = {-1.5,-1.75,-1.88}; // should be ID {-2,-2,-2}
	pos[2] = {-1.5,1.75,1.88}; // should be ID {-2,1,1}
	pos[3] = {1.5,-1.75,1.88}; // should be ID {1,-2,1}
	pos[4] = {1.5,1.75,-1.88}; // should be ID {1,1,-2}
	pos[5] = {-1.5,1.75,-1.88}; // should be ID {-2,1,-2}
	pos[6] = {1.5,-1.75,-1.88}; // should be ID {1,-2,-2}
	pos[7] = {-1.5,-1.75,1.88}; // should be ID {-2,-2,1}
	
	wrapper w(8,1,pos);

	w.g -> printGroups();

	// std::cout<<"HEREERER"<<std::endl;
	// w.g.setMaxGridIndex();
	// std::cout<<"maxGridIndex: "<<w.g -> maxGridIndex<<std::endl;
	//should be 2
}

void test_groups2()
{
	
	std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/198_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_simData.csv";
	wrapper w(200,file,2e-5);

	w.g -> printGroups();

	// std::cout<<"HEREERER"<<std::endl;
	// w.g.setMaxGridIndex();
	// std::cout<<"maxGridIndex: "<<w.g -> maxGridIndex<<std::endl;
	//should be 2
}

void test_map1()
{
	vec3 *pos = nullptr;
	pos = new vec3[16];
	pos[0] = {0.5,0.75,0.88}; // should be ID {0,0,0}
	pos[1] = {-0.5,-0.75,-0.88}; // should be ID {-1,-1,-1}
	pos[2] = {-0.5,0.75,0.88}; // should be ID {-1,0,0}
	pos[3] = {0.5,-0.75,0.88}; // should be ID {0,-1,0}
	pos[4] = {0.5,0.75,-0.88}; // should be ID {0,0,-1}
	pos[5] = {-0.5,0.75,-0.88}; // should be ID {-1,0,-1}
	pos[6] = {0.5,-0.75,-0.88}; // should be ID {0,-1,-1}
	pos[7] = {-0.5,-0.75,0.88}; // should be ID {-1,-1,0}

	pos[8] = {1.5,1.75,1.88}; // should be ID {1,1,1}
	pos[9] = {-1.5,-1.75,-1.88}; // should be ID {-2,-2,-2}
	pos[10] = {-1.5,1.75,1.88}; // should be ID {-2,1,1}
	pos[11] = {1.5,-1.75,1.88}; // should be ID {1,-2,1}
	pos[12] = {1.5,1.75,-1.88}; // should be ID {1,1,-2}
	pos[13] = {-1.5,1.75,-1.88}; // should be ID {-2,1,-2}
	pos[14] = {1.5,-1.75,-1.88}; // should be ID {1,-2,-2}
	pos[15] = {-1.5,-1.75,1.88}; // should be ID {-2,-2,1}
	
	wrapper w(16,1,pos);
	// w.g -> printMap();
	for (int i = 0; i < 16; i++)
	{
		w.g -> printVector(w.g -> getGroup(w.g -> gridIDs[i]));
	}
	// std::cout<<w.g -> IDToGrid["-2-21"][0]<<std::endl;
	return;
}

void test_map2()
{
	std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/198_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_simData.csv";
	wrapper w(200,file,2e-5);
	w.g -> printMap(w.g -> IDToGrid);
	return;
}

int main(int argc, char const *argv[])
{
	// test_groups2();
	test_map1();

	return 0;
}

