#include <iostream>
#include <cmath>
#include <fstream>
#include "grid_group.hpp"
#include "../vec3.hpp"



// class wrapper
// {
// public:
// 	vec3 *pos = nullptr;
// 	int num_particles;
// 	grid g();
	
// 	wrapper() = default;

// 	wrapper(int num,std::string file,double size)
// 	{
// 		num_particles = num;
// 		pos = new vec3[num];
// 		parseSimData(getLastLine(file));
		

// 		grid g(num,1e-5,size,pos);
// 	}

// 	wrapper(int num,double size,vec3* positions)
// 	{
// 		num_particles = num;
// 		// pos = new vec3[num];
// 		pos = positions;


// 		grid g(num,1e-5,size,pos);
// 	}

// 	wrapper(int num,double radius,double size,vec3* positions,double tol)
// 	{
// 		num_particles = num;
// 		// pos = new vec3[num];
// 		pos = positions;


// 		grid g(num,radius,size,pos,tol);
// 	}

	
// };


void parseSimData(std::string line,vec3 *pos,int num_particles)
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



void test_groups1()
{
	// std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/198_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_simData.csv";
	vec3 *pos = nullptr;
	pos = new vec3[8];
	double radius = 1;
	double grid_size = 1;
	double tol = 3*radius;
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
	
	// wrapper w(8,1,pos);
	grid g(8,radius,grid_size,pos,tol);

	g.printGroups();

	// std::cout<<"HEREERER"<<std::endl;
	// g.setMaxGridIndex();
	// std::cout<<"maxGridIndex: "<<g.maxGridIndex<<std::endl;
	//should be 2
}

void test_groups2()
{
	double radius = 1e-5;
	double grid_size = 1e-4;
	double tol = 3*radius;
	std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/198_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_simData.csv";
	vec3 *pos = nullptr;
	pos = new vec3[200];
	parseSimData(getLastLine(file),pos,200);
	// wrapper w(200,file,2e-5);
	grid g(200,radius,grid_size,pos,tol);

	g.printGroups();

	// std::cout<<"HEREERER"<<std::endl;
	// g.setMaxGridIndex();
	// std::cout<<"maxGridIndex: "<<g.maxGridIndex<<std::endl;
	//should be 2
}

void test_map1()
{
	vec3 *pos = nullptr;
	pos = new vec3[16];
	double radius = 1e-5;
	double grid_size = 1e-4;
	double tol = 3*radius;
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
	
	// wrapper w(16,1,pos);
	grid g(30,radius,grid_size,pos,tol);
	// g.printMap();
	for (int i = 0; i < 16; i++)
	{
		g.printVector(g.getBalls(i));
	}
	// std::cout<<g.IDToGrid["-2-21"][0]<<std::endl;
	return;
}

void test_map2()
{
	std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/198_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_simData.csv";
	// wrapper w(200,file,2e-5);
	vec3 *pos = nullptr;
	pos = new vec3[200];
	parseSimData(getLastLine(file),pos,200);
	double radius = 1e-5;
	double grid_size = 1e-4;
	double tol = 3*radius;
	grid g(200,radius,grid_size,pos,tol);
	// g.printMap(g.IDToGrid);
	return;
}

void test_balls()
{
	std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/198_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_simData.csv";
	// wrapper w(200,file,2e-5);
	vec3 *pos = nullptr;
	pos = new vec3[200];
	parseSimData(getLastLine(file),pos,200);
	double radius = 1e-5;
	double grid_size = 1e-4;
	double tol = 3*radius;
	grid g(200,radius,grid_size,pos,tol);
	g.printVector(g.getBalls(100));
}

void test_balls2()
{//put one ball in every box and move around pos[0] so it will give appropriate balls 
	//test v1
	int num = 27;
	vec3 *pos = nullptr;
	pos = new vec3[num];
	double radius = 1e-5;
	double grid_size = 4*radius;
	double mid = grid_size/2.0;
	double tol = radius*3;
	pos[0] = {mid,mid,mid};
	// pos[27] = pos[0];
	int l = 1;

	for (int i = -1; i < 2; i++)
	{
		for (int j = -1; j < 2; j++)
		{
			for (int k = -1; k < 2; k++)
			{
				if (!(i == 0 && j == 0 && k == 0))
				{
					pos[l] = pos[0];
					pos[l][0] += i*grid_size;
					pos[l][1] += j*grid_size;
					pos[l][2] += k*grid_size;
					// std::cout<<pos[l]<<std::endl;
					l++;
				}
			}
		}	
	}

	// pos[28] = {0,0,0};
	// pos[29] = {grid_size+mid,0,0};


	grid g(num,radius,grid_size,pos,tol);
	for (int i = 0; i < num; ++i)
	{
		std::cout<<"i="<<i<<" vec len="<<g.getBalls(i).size()<<'\n';
		// g.printVector(g.getBalls(i));
	}
	//should give [0,27]
	
	//all vertex tests should include balls 0 and 27 and should print 9 balls
	// pos[0] = {0,0,0};
	// g.printVector(g.getBalls(0));
	// pos[0] = {grid_size,0,0};
	// g.printVector(g.getBalls(0));
	// pos[0] = {0,grid_size,0};
	// g.printVector(g.getBalls(0));
	// pos[0] = {0,0,grid_size};
	// g.printVector(g.getBalls(0));
	// pos[0] = {grid_size,grid_size,0};
	// g.printVector(g.getBalls(0));
	// pos[0] = {0,grid_size,grid_size};
	// g.printVector(g.getBalls(0));
	// pos[0] = {grid_size,0,grid_size};
	// g.printVector(g.getBalls(0));
	// pos[0] = {grid_size,grid_size,grid_size};
	// g.printVector(g.getBalls(0));

	// // all side tests should include balls 0 and 27 and be 3 balls long
	// pos[0] = {grid_size,mid,mid};
	// g.printVector(g.getBalls(0));
	// pos[0] = {0,mid,mid};
	// g.printVector(g.getBalls(0));
	// pos[0] = {mid,grid_size,mid};
	// g.printVector(g.getBalls(0));
	// pos[0] = {mid,0,mid};
	// g.printVector(g.getBalls(0));
	// pos[0] = {mid,mid,grid_size};
	// g.printVector(g.getBalls(0));
	// pos[0] = {mid,mid,0};
	// g.printVector(g.getBalls(0));

}

void test_balls3()
{
	vec3 *pos = nullptr;
	pos = new vec3[3];
	double grid_size = 1e-4;
	// double mid = grid_size/2.0;
	double radius = 1e-5;
	double tolerance = radius*4;

	// pos[0] = {-7.83584e-06,2.52112e-05,6.49667e-06};
	// pos[1] = {-7.83584e-06,3.19117e-06,6.49667e-06};
	// pos[2] = {1.56717e-05,-2.84023e-05,-1.29933e-05};

	// pos[0] = {9.51319e-07,-2.90044e-06,-1.05297e-05};
	// pos[1] = {9.51319e-07,-2.49204e-05,-1.05297e-05};
	// pos[2] = {-1.90264e-06,2.78209e-05,2.10594e-05};

	pos[0] = {1.11215e-05,1.85625e-05,1.11596e-05};
	pos[1] = {1.11215e-05,-3.45745e-06,1.11596e-05};
	pos[2] = {-2.2243e-05,-1.51051e-05,-2.23191e-05};

	grid g(3,radius,grid_size,pos,tolerance);
	g.printVector(g.getBalls(0));
	g.printVector(g.getBalls(1));
	g.printVector(g.getBalls(2));
}

int main(int argc, char const *argv[])
{
	// test_groups2();
	// test_map1();
	// test_balls2();
	timey t;
	t.start_event("map init");
	test_map2();
	t.end_event("map init");
	t.print_events();

	return 0;
}

